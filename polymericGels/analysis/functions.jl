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

function computeDensity(pp)
    """
        Compute the density for the structure factor
        pp is the dot product.
        pp is an array. Each row represent the dot product of onde direction with a distance
    """
        return mapreduce(s->cos(s),+,pp)^2 + mapreduce(s->sin(s),+,pp)^2
end

function computeStructureFactor_PBC(L, Rmax, N_lambda, lambda_o, lambda_f, time_instant, dump_paths, N_part)
    """
        Factor de estructura usando solo vectores de onda compatibles con PBC.
        L: lado de la caja periódica.
        Rmax: máximo entero para nx,ny,nz tal que q_max ≈ (2π/L)*Rmax.
    """
    # 1. Dominio de magnitudes q (exactamente como en tu código original)
    q_min = 2π / lambda_f
    q_max = 2π / lambda_o
    q_dom = range(q_min, q_max, length=N_lambda)

    # 2. Generar todos los vectores PBC-compatibles y su bin
    qx_all = Float64[]
    qy_all = Float64[]
    qz_all = Float64[]
    qbin   = Int[]            # índice del punto más cercano en q_dom

    for nx in -Rmax:Rmax, ny in -Rmax:Rmax, nz in -Rmax:Rmax
        n2 = nx^2 + ny^2 + nz^2
        n2 == 0 && continue
        n = sqrt(n2)
        n > Rmax && continue
        qvec = (2π/L) .* (nx, ny, nz)
        push!(qx_all, qvec[1])
        push!(qy_all, qvec[2])
        push!(qz_all, qvec[3])
        push!(qbin, argmin(abs.(norm(qvec) .- q_dom)))
    end

    # 3. Cargar posiciones de los experimentos (igual que antes)
    r_exp = getPositions(time_instant, dump_paths)

    # 4. Calcular S(q) instantáneo para cada experimento
    S_q = [zeros(N_lambda) for _ in eachindex(r_exp)]

    for (it_exp, dist0) in enumerate(r_exp)
        dist = reduce(hcat, dist0)          # matriz N_part x 3
        acum = zeros(Float64, N_lambda)
        cont = zeros(Int, N_lambda)

        # Iterar sobre los vectores precomputados
        for i in eachindex(qx_all)
            qv = [qx_all[i], qy_all[i], qz_all[i]]
            fase = dist * qv               # vector de productos punto (N_part)
            rho  = computeDensity(fase)    # |∑ exp(i·fase_j)|^2
            idx  = qbin[i]
            acum[idx] += rho
            cont[idx] += 1
        end

        # Normalizar por número de vectores en la cáscara y por N_part
        for k in 1:N_lambda
            if cont[k] > 0
                S_q[it_exp][k] = acum[k] / (N_part * cont[k])
            else
                S_q[it_exp][k] = 0.0
            end
        end
        println("Experiment ", it_exp, " done")
    end

    return mean(S_q)
end


function analyzeStructureFactorPBC(dat_DF)
"""
    Function that stores the structure factor of a system
"""

# Parametros para obtener el factor de estructura
lambda_o=1; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_part=Int(dat_DF.Npart[1]);
N_qu=2^9; # EXPONENTE DEBE SER IMPAR Cantidad de direcciones
N_lambda=2^11; # Cantidad de magnitudes
N_instants=2;

# Path to the dumps
dump_paths=joinpath.(dat_DF.PARENT_DIR,dat_DF.dir,"traj");

# Seleccion de time instants
aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));
aux_id=aux_timeStep[ind];
time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];

# Get the structure factor for different time intervals
S_q=[zeros(N_lambda) for _ in eachindex(time_instants)];

println("Start of Structure factor analysis")
for it_t in eachindex(time_instants)
    S_q[it_t]=computeStructureFactor_PBC(lambda_f, lambda_f, N_lambda, lambda_o, lambda_f, time_instants[it_t], dump_paths, N_part);
    #computeStructureFactor(N_qu,N_lambda,lambda_o,lambda_f,time_instants[it_t],dump_paths,N_part);
    println("One time step done")
end

Sq_df=DataFrame(;
    timeStep = aux_id,
    q_domain = fill(range(2*pi/lambda_f,2*pi/lambda_o,length=N_lambda), length(aux_id)),
    Sq       = S_q,
    lambda_o = fill(lambda_o, length(aux_id)),
    lambda_f = fill(lambda_f, length(aux_id)),
    id=unique(dat_DF.id)[1]
)

CSV.write(joinpath(joinpath(pwd(),"analyzedData"),string("structureFactorPBC",Sq_df.id[1],".csv")),Sq_df)
println("Structure factor analysis done. File written.")

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


