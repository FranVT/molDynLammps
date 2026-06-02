"""
    File with all the functions needed for multiple.jl script
"""

###########################
# UPDATE DAT FILE
###########################

function getDat(path)
"""
    Creates a dataframe from the dat file of the experiment.
"""
    file_path=joinpath(path,"dataAssembly.dat");
    aux=split.(readlines(file_path),",");

    df_aux=DataFrame();
    for (col, val) in zip(aux[1], aux[2])
        # Convertir a Float64 si es posible, sino mantener como String
        parsed_val = tryparse(Float64, val)
        if parsed_val !== nothing
            df_aux[!, col] = [parsed_val]
        else
            df_aux[!, col] = [val]
        end
    end

    df_aux[!, :dir] = [last(split(path,"/"))];

    return df_aux
end


function createDatFiles()
"""
    Create a csv file with the information of the dat files and assign an id as a function of the system.
"""
# Get all directories in the data directory
MAIN_DIR=dirname(pwd());
DATA_DIR="/run/media/franvt/rogelio/DinMol/gel-batch-1/data";
STORE_DIR=joinpath(MAIN_DIR,"datFiles");

# Read the directory where data is stored
SIMS_DIR=filter(isdir,readdir(DATA_DIR,join=true));

# Get all data files in a dataframe
data_info=mapreduce(s->getDat(s),vcat,SIMS_DIR);

# Add related paths
data_info[!,:PARENT_DIR].=[DATA_DIR];

# Create an id for each experiment
id=[
data_info.phi,
data_info."CL-Con",
data_info.Npart,
data_info.Temperature,
data_info.damp,
data_info."time-step",
data_info."N_heat",
data_info."N_isot",
data_info."N_CL",
data_info."N_MO",
data_info."L"
]

# Selección de categorias
categories=[:phi,Symbol("CL-Con"),:Temperature,:Npart,:damp,Symbol("time-step")];

# Creación de los subdataframes
data_bySystem=groupby(data_info,categories);

# Obtención de los valores numéricos de las categorias
id_values=unique(data_info[:,categories]);

# Creación de los ids
ids=map(s->join(string.(collect(id_values[s,:]))),1:nrow(id_values));

# Asignar el id a los data frames
map(s->data_bySystem[s][!,:id].=[ids[s]],1:length(data_bySystem));

# Guardar el dataframe en un csv
CSV.write(joinpath(pwd(),"datFiles","systemDatfiles.csv"),reduce(vcat,data_bySystem));

end

function extractDatFiles()
"""
    Function that returns a dataframe with the dat file information of the systems
"""
# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","systemDatfiles.csv");
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

return CSV.read(DAT_PATH,DataFrame);

end



#######
#   STRUCTURE FACTOR RELATED FUNCTIONS
#######

function getPositions(time_instant,dump_paths)
    """
        Extract the positions of N dataframes for one time instant.
    """
        # Obtenemos los dumps de los N experimentos para un instante de tiempo
        dumps=[getDump(path,time_instant) for path in dump_paths];

        # Vector de N elementos. Cada elemento es la posición de N partículas de los N_exp.
        return [getPosition(df) for df in dumps];
end

function computeDensity(pp::AbstractMatrix)
    # Sumar cosenos y senos a lo largo de las filas (dimensión 1),
    # obteniendo un vector fila 1×D con los resultados por columna.
    sum_cos = sum(cos.(pp), dims=1)   # 1 × N_direcciones
    sum_sin = sum(sin.(pp), dims=1)
    # Elevar al cuadrado elemento a elemento y convertir a vector columna
    return vec(sum_cos.^2 + sum_sin.^2)
end

function createNvector(n_max)
"""
    Create wave vector considering the periodic Boundary conditions restriction
"""
    # Crear los vectores por componentes enteras evitando el vector nulo
    vectores = [[x,y,z] for x in -n_max:n_max, y in -n_max:n_max, z in -n_max:n_max if (x,y,z) != (0,0,0)]

    # Agrupar los vectores por magnitud
    return group(v -> v[1]^2 + v[2]^2 + v[3]^2, vectores)
end

function structureFactor(pos_exp, mag_n, n, qo, N_part; N_exp=length(pos_exp))
    n_mag = length(mag_n)
    sq = zeros(n_mag)          # salida prealocada

    for (i, mag) in enumerate(mag_n)
        n_vec = n[mag]
        n_dir = length(n_vec)
        # Matriz de proyección: columnas son qo .* dirección
        Q = reduce(hcat, [qo .* v for v in n_vec])   # tamaño (3, N_direcciones)

        rho = zeros(n_dir, N_exp)   # densidad por dirección y experimento
        for it_exp in 1:N_exp
            pos = pos_exp[it_exp]           # (N_part, 3)
            PP = pos * Q                    # (N_part, n_dir)
            rho[:, it_exp] = computeDensity(PP)
        end

        # Promedio de direcciones y experimentos
        # Sq = mean(rho, dims=2) ./ N_part;  mean(Sq) == mean(rho) / N_part
        sq[i] = mean(rho) / N_part
    end

    return sq
end

function dumpAnalysis(dat_DF)
"""
    Analysis that are done from the dump files information
"""

    # Parameters of the system relevant to the structure factor
    L=2*dat_DF.L[1];    # Longitud de la caja
    N_instants=2;       # Instantes temporales a analizar
    n_max=2*round(Int,first(dat_DF.L)/sqrt(3));
    #2^1;          # Magnitud máxima de cada componente
    N_exp=nrow(dat_DF); # Cantidad de experimentos por sistema
    N_part=dat_DF.Npart[1];
    qo=2*pi/L;          # Considera condiciones periódicas de frontera

    # Create time steps range
    aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));

    # Get the index for the time instants that we are interested to analyzed
    ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));

    # Get the time steps
    timeSteps=aux_timeStep[ind];

    # Create file names
    file_names=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in timeSteps];

    # Path to the dumps
    dump_paths=joinpath.(dat_DF.PARENT_DIR,dat_DF.dir,"traj");

    # Crear los vectores n
    n=createNvector(n_max);
    mag_n=sort(collect(keys(n)));

    #=
        Start the analysis of the structure factor
    =#

    # Alocar memoria para el factor de estructura
    Sq_expval=zeros(length(n),length(file_names));
    
    println("Start cycle of Structure factor")
    for (it_time,~) in enumerate(timeSteps)
        file=file_names[it_time];   # Select one time step
     
        r_exp = getPositions(file, dump_paths); # Get the position of all experiments of the same system at the same time step
        pos_exp=[reduce(hcat, r_exp[it_exp]) for it_exp in 1:N_exp];    # Get all positions of all experiments at the same time step
    
        # Calcular y almacenar S(q) para todas las magnitudes
        Sq_expval[:, it_time] = structureFactor(pos_exp, mag_n, n, qo, N_part; N_exp=N_exp);

    end

    println("One system done")
    df=DataFrame([mag_n Sq_expval],[:mag, Symbol.("Sq",timeSteps)... ])

    return df #[mag_n Sq_expval]

end


function storeAllSq(data_bySystem,Sq_all)
"""
    Function that stores all structure factores analyzed
"""
    file_path=joinpath(pwd(),"analyzedData");
    files_name=map(s->string("structureFactorPBC",first(data_bySystem[s].id),".csv"),1:length(data_bySystem))

    map(s->CSV.write(joinpath(file_path,files_name[s]),Sq_all[s]),1:length(data_bySystem))
    println("Files written")
end




#######
#   DUMP FUNCTIONS
#######

function getPosition(dump)
"""
    Get the position of the central particles of a given dump
"""

    # Filtrar
    mask=(dump.type .==1) .| (dump.type .== 2.0);
    dump_filtered=dump[mask,:];

   return [dump_filtered.x,dump_filtered.y,dump_filtered.z]

end

function getDump(dir,file_name)
"""
    Get the data from a single dump file that stores one timeste information
"""
    data = split.(readlines(joinpath(dir,file_name))," ")[9:end];
    HEADERS=data[1][3:end];
    INFO=parse.(Float64,reduce(hcat,data[2:end]))';

    return DataFrame(INFO,HEADERS)
end

#######
#   OTHER  FUNCTIONS
#######


function meanFixystem(dirs)
"""
    Function that returns a dataframe with the mean of N experiments of the observables stored in a fix file.
"""

    # Obtener la información de lo fix files
    data_fix=map(s->dataSystem=extractFixScalar(s,"system_assembly.fixf"),dirs);

    # Variables auxiliares
    n_row = nrow(data_fix[1]);
    n_col = ncol(data_fix[1]);
    cols = names(data_fix[1]);

    # Apilar matrices en un arreglo 3D
    arr = Array{Float64}(undef, n_row, n_col, length(data_fix))
    for (i, df) in enumerate(data_fix)
        arr[:, :, i] = Matrix(df)   # asume columnas numéricas
    end

    # Promedio a lo largo de la tercera dimensión
    promedio = dropdims(mean(arr, dims=3), dims=3);

    return DataFrame(promedio,cols);
end


function extractFixScalar(path_system,file_name)
"""
    Function that extracts the information of fix files that stores global scalar values
"""
    aux=split.(readlines(joinpath(path_system,file_name))," ");
    header=aux[2][2:end];
    info=reduce(hcat,map(s->parse.(Float64,s),aux[3:end]));

    return DataFrame(info',header)
end


