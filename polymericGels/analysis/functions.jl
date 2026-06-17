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


function createDat_ID(DATA_DIR)
"""
    Create a csv file with the information of the dat files and assign an id as a function of the system.
"""
# Get all directories in the data directory
MAIN_DIR=dirname(pwd());
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

function extract_FixScalar(path_file)
"""
    Function that extracts the information of fix files that stores global scalar values
"""
    aux=split.(readlines(path_file)," ");
    header=aux[2][2:end];
    info=reduce(hcat,map(s->parse.(Float64,s),aux[3:end]));

    return DataFrame(info',header)
end


function mean_FixSystem(dirs)
"""
    Function that returns a dataframe with the mean of N experiments of the observables stored in a fix file.
"""

    # Obtener la información de lo fix files
    data_fix=map(s->dataSystem=extract_FixScalar(s),dirs);

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


function storeAvg_fix(df)
"""
    Store a csv file with the assembly avergage of a system of the fix file observables 
"""

    # Extraer informacion del archivo de fix
    col=[:PARENT_DIR,:dir,:file0];

    # For each group of dataframes
    files=unique(df[:,col]);

    # Create the path to each fix file of each experiment of the same system
    files.path=joinpath.(files.PARENT_DIR, files.dir, files.file0);

    # Compute the assambly average of the fix file information
    fix_info=mean_FixSystem(files.path);

    # Get the id of the system
    id_system=unique(df[:,:id]);

    # Add the id 
    fix_info.id.=id_system;

    # Create file name to store the data
    filename=string("fix_avg_",first(id_system),".csv");

    # Directory to store the information
    file_dir=joinpath(pwd(),"analyzedData");

    # File path to store the data
    file_path=joinpath(file_dir,filename);

    # Save the average 
    CSV.write(file_path,fix_info);

end


