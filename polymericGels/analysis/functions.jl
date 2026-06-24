###########################
# UPDATE DAT FILE
###########################

"""
    Creates a dataframe from the dat file of the experiment.
"""
function getDat(path)
    file_path = joinpath(path, "dataAssembly.dat")
    aux = split.(readlines(file_path), ",")

    df_aux = DataFrame()
    for (col, val) in zip(aux[1], aux[2])
        # Convertir a Float64 si es posible, sino mantener como String
        parsed_val = tryparse(Float64, val)
        if parsed_val !== nothing
            df_aux[!, col] = [parsed_val]
        else
            df_aux[!, col] = [val]
        end
    end

    df_aux[!, :dir] = [last(split(path, "/"))]

    return df_aux
end

"""
    Create a csv file with the information of the dat files and assign an id as a function of the system.
"""
function createDat_ID(DATA_DIR)
    # Get all directories in the data directory
    MAIN_DIR   = dirname(pwd())
    STORE_DIR  = joinpath(MAIN_DIR, "datFiles")

    # Read the directory where data is stored
    SIMS_DIR = filter(isdir, readdir(DATA_DIR, join = true))

    # Get all data files in a dataframe
    data_info = mapreduce(s -> getDat(s), vcat, SIMS_DIR)

    # Add related paths
    data_info[!, :PARENT_DIR] .= [DATA_DIR]

    # Create an id for each experiment
    id = [
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
    categories = [:phi, Symbol("CL-Con"), :Temperature, :Npart, :damp, Symbol("time-step")]

    # Creación de los subdataframes
    data_bySystem = groupby(data_info, categories)

    # Obtención de los valores numéricos de las categorias
    id_values = unique(data_info[:, categories])

    # Creación de los ids
    ids = map(1:nrow(id_values)) do s
        join(string.(collect(id_values[s, :])))
    end

    # Asignar el id a los data frames
    for (i, df) in enumerate(data_bySystem)
        df[!, :id] .= ids[i]
    end

    # Guardar el dataframe en un csv
    CSV.write(joinpath(pwd(), "datFiles", "systemDatfiles.csv"), reduce(vcat, data_bySystem))
end

"""
    Function that returns a dataframe with the dat file information of the systems
"""
function extractDatFiles()
    # Get directories
    MAIN_DIR   = pwd()
    DAT_PATH   = joinpath(MAIN_DIR, "datFiles", "systemDatfiles.csv")
    SAVE_DIR   = joinpath(MAIN_DIR, "analyzedData")

    return CSV.read(DAT_PATH, DataFrame)
end

#######
#   STRUCTURE FACTOR RELATED FUNCTIONS
#######

"""
    Extract the positions of N dataframes for one time instant.
"""
function getPositions(time_instant, dump_paths)
    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps = [getDump(path, time_instant) for path in dump_paths]

    # Vector de N elementos. Cada elemento es la posición de N partículas de los N_exp.
    return [getPosition(df) for df in dumps]
end

"""
    Get the names files 
"""
function getpathfilesSq(dat_DF, N_instants)
    # Create time steps range
    aux_timeStep = Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])))

    # Get the index for the time instants that we are interested to analyzed
    ind = round.(Int, LinRange(1, length(aux_timeStep), N_instants))

    # Get the time steps
    timeSteps = aux_timeStep[ind]

    # Create file names
    file_names = [replace("traj_assembly.*.dumpf", "*" => string(it)) for it in timeSteps]

    # Path to the dumps
    dump_paths = joinpath.(dat_DF.PARENT_DIR, dat_DF.dir, "traj")
    return (timeSteps, dump_paths, file_names)
end

"""
    Create the reciprocal space of wave vectors
"""
function createqdom(qmax, rq0, dq0, bin0, numbin)
    # Crear el espacio del vector recíproco para no definirlo en cada experimento
    qx = (-qmax:qmax) .* dq0
    qy = (-qmax:qmax) .* dq0
    qz = (-qmax:qmax) .* dq0

    qhis   = [[] for _ in 1:numbin]
    qxhis  = [[] for _ in 1:numbin]
    qyhis  = [[] for _ in 1:numbin]
    qzhis  = [[] for _ in 1:numbin]

    # Calcular la magnitud cuadrada
    for x in qx
        for y in qy
            for z in qz
                # Excluir el origen del espacio recíproco
                if x == 0 && y == 0 && z == 0
                    continue
                end

                q = sqrt(x^2 + y^2 + z^2)

                # Si la magnitud es grande que el maximo se salta
                if q > rq0
                    continue
                end

                # Determinar el índice del bin (manejo especial del borde)
                sbin = floor(Int, q / bin0) + 1

                # Coso para transformar de un dominio continuo a uno discreto
                if q % bin0 == 0.0
                    sbin -= 1
                end

                # Si la magnitud del vector supera la magnitud de interés se lo salta
                if sbin > numbin
                    continue
                end

                append!(qhis[sbin], q)
                append!(qxhis[sbin], x)
                append!(qyhis[sbin], y)
                append!(qzhis[sbin], z)
            end
        end
    end

    # Compute the mean handeling the empty vectors
    qmean=[isempty(v) ? 0.0 : mean(v) for v in qhis];

    return (qxhis, qyhis, qzhis, qhis, qmean)
end

"""
    Compute the structure factor of a set of positions
"""
function computeSq(numbin, ntotav, qxhis, qyhis, qzhis, qhis, r)
    Sq = zeros(numbin, 2)
    rho = [[] for _ in 1:numbin]

    for it_bin in 1:numbin
        # Seleccion de las componentes del vector de onda
        qx = qxhis[it_bin]
        qy = qyhis[it_bin]
        qz = qzhis[it_bin]

        # Manage null vectors
        if isempty(qx) || isempty(qy) || isempty(qz)
            append!(rho[it_bin], 0.0)
        else
                    # Calcular la densidad para cada bin
            for it_q in eachindex(qx)
                vq = [qx[it_q], qy[it_q], qz[it_q]]
                dp = r * vq  # Producto punto del vector para cada partícula
                rho_re = sum(cos.(dp))
                rho_im = sum(sin.(dp))
                sq = (rho_re^2 + rho_im^2) / ntotav
                append!(rho[it_bin], sq)
            end
        end 

    end

    # Guardamos información
    # Valor esperado del factor de estructura
    Sq[:, 1] = sum.(rho) ./ numbin 
    smax = maximum(Sq[:, 1])
    Sq[:, 2] = Sq[:, 1] / smax

    return Sq
end

"""
    Compute the avg structure factor of a system at a given time step
"""
function computeSqmean(file_names, timeStep, file_path, dump_paths, dat_DF, numbin, ntotav,
                       qxhis, qyhis, qzhis, qhis, qmean, info)

    # Preparar informacion para el factor de estructura
    r_system = getPositions(file_names, dump_paths)

    # Calcular el factor de estructura para N experimentos del mismo sistema
    for it_exp in eachindex(r_system)
        r = reduce(hcat, r_system[it_exp])
        info[it_exp][:, 2:3] = computeSq(numbin, ntotav, qxhis, qyhis, qzhis, qhis, r)
        println("One experiment done")
    end

    # Se obtienen promedios
    info = mean(info)        # Compute the assembly mean of the Sq
    info[:, 1] = qmean       # Store the mean of the wave vectors

    # Paths and file names and stuff
    filename = string("structureFactorPBC", first(dat_DF.id), "time", timeStep, ".csv")

    # Create the data frame
    df = DataFrame([repeat([timeStep], numbin) info],
                   [:timeStep, :qmean, :Sqmean, :Sqmeannorm])

    # Save the data frame
    CSV.write(joinpath(file_path, filename), df)
end

"""
    Compute all structure factor of all system of N time instants
"""
function computeAllTimeSqmean(dat_DF, N_instants, qmax0)
    # Parameters of the system relevant to the structure factor
    L = 2 * dat_DF.L[1]            # Longitud de la caja
    N_exp = nrow(dat_DF)           # Cantidad de experimentos por sistema
    N_part = dat_DF.Npart[1]       # Numero de particulas centrales

    # Crea los dominios del vector de onda
    xc = L                 # longitud de la caja en x
    yc = L                 # longitud de la caja en y
    zc = L                 # longitud de la caja en z
    dq0 = 2 * pi / xc      # Δq fundamental
    qmax = Int(floor(qmax0 / dq0))   # número entero de pasos hasta qmax0
    rq0 = qmax * dq0       # valor máximo real de |q| usado
    bin0 = dq0             # nuevo ancho de bin (bin0 original * dq0)
    numbin = Int(floor(qmax * dq0 / bin0)) + 1  # número total de bines

    # Parametros para el factor de estructura
    ntotav = Int(N_part)   # Número total de partículas
    (qxhis, qyhis, qzhis, qhis, qmean) = createqdom(qmax, rq0, dq0, bin0, numbin)  # vector de onda
    file_path = joinpath(pwd(), "analyzedData")

    # Obtener los paths y los nombres de los archivos a analizar
    (timeSteps, dump_paths, file_names) = getpathfilesSq(dat_DF, N_instants)

    # Select one time step
    # Alocar memoria
    info = [zeros(numbin, 3) for _ in 1:N_exp]

    for it_time in eachindex(file_names)
        computeSqmean(file_names[it_time], timeSteps[it_time], file_path, dump_paths, dat_DF,
                      numbin, ntotav, qxhis, qyhis, qzhis, qhis, qmean, info)
        println("One time step done")
    end

    println("One system done")
end

#######
#   DUMP FUNCTIONS
#######

"""
    Get the position of the central particles of a given dump
"""
function getPosition(dump)
    # Filtrar
    mask = (dump.type .== 1) .| (dump.type .== 2.0)
    dump_filtered = dump[mask, :]

    return [dump_filtered.x, dump_filtered.y, dump_filtered.z]
end

"""
    Get the data from a single dump file that stores one timeste information
"""
function getDump(dir, file_name)
    data = split.(readlines(joinpath(dir, file_name)), " ")[9:end]
    HEADERS = data[1][3:end]
    INFO = parse.(Float64, reduce(hcat, data[2:end]))'

    return DataFrame(INFO, HEADERS)
end

#######
#   CLUSTER ANALYSIS
#######

"""
    Perform cluster analysis for a system over multiple time instants.
"""
function clusterAnalysis(dat_DF, N_instants)
    # Get the paths to the data
    (timeSteps, dump_paths, file_names) = getpathfilesSq(dat_DF, N_instants)

    # Get the position of all the experiments of the same system
    N_exp = length(dat_DF.Nexp)   # Amount of experiments of the same system

    # Save memory
    size_clusters = [[] for _ in 1:N_exp]
    N_clusters    = [zeros(N_instants) for _ in 1:N_exp]

    # Get the amount of clusters at each time step for all experiments
    # $ C_clusters is an id
    for it_exp in 1:N_exp
        # Get the position of all time steps for one experiment
        positions = map(s -> getDump(dump_paths[it_exp], s), file_names)

        # Get the amount of central particles in the clusters
        filtered = filter.(row -> row.type in (1, 2), positions)
        counts = map(s -> combine(groupby(s, :c_clusters), nrow => :count), filtered)

        # Store the information
        size_clusters[it_exp] = map(s -> s.count, counts)

        # Get the amount of clusters for each time step
        N_clusters[it_exp] = length.(size_clusters[it_exp])
    end

    # Compute the assemble average of amount of clusters in each timestep
    # row|col -> time|phi
    N_clusters_mean = mean(N_clusters)

    # Compute the assemble average of the amount of central particles in the biggest cluster in each timestep
    # row|col -> time|phi
    max_particles = mean(map(s -> maximum.(s), size_clusters))

    # For the histogram
    hist_size = map(l -> mapreduce(s -> size_clusters[s][l], vcat, 1:N_exp), eachindex(timeSteps))

    # Save the information per timeStep

    # Directory to store the information
    file_dir = joinpath(pwd(), "analyzedData")

    # Paths and file names and stuff
    filename = string("clusterAnalysis", first(dat_DF.id), ".csv")

    # File path to store the data
    file_path = joinpath(file_dir, filename)

    # Create the DataFrame
    df = DataFrame([timeSteps N_clusters_mean max_particles hist_size],
                   [:timeStep, :nClusters, :maxParticles, :hist])

    # Save the data frame
    CSV.write(file_path, df)

    println("One cluster analysis saved")
end


#######
#   FIX FUNCTIONS
#######

"""
    Function that extracts the information of fix files that stores global scalar values
"""
function extract_FixScalar(path_file)
    aux = split.(readlines(path_file), " ")
    header = aux[2][2:end]
    info = reduce(hcat, map(s -> parse.(Float64, s), aux[3:end]))

    return DataFrame(info', header)
end

"""
    Function that returns a dataframe with the mean of N experiments of the observables stored in a fix file.
"""
function mean_FixSystem(dirs)
    # Obtener la información de lo fix files
    data_fix = map(s -> extract_FixScalar(s), dirs)

    # Variables auxiliares
    n_row = nrow(data_fix[1])
    n_col = ncol(data_fix[1])
    cols  = names(data_fix[1])

    # Apilar matrices en un arreglo 3D
    arr = Array{Float64}(undef, n_row, n_col, length(data_fix))
    for (i, df) in enumerate(data_fix)
        arr[:, :, i] = Matrix(df)   # asume columnas numéricas
    end

    # Promedio a lo largo de la tercera dimensión
    promedio = dropdims(mean(arr, dims = 3), dims = 3)

    return DataFrame(promedio, cols)
end

"""
    Store a csv file with the assembly avergage of a system of the fix file observables
"""
function storeAvg_fix(df)
    # Extraer informacion del archivo de fix
    col = [:PARENT_DIR, :dir, :file0]

    # For each group of dataframes
    files = unique(df[:, col])

    # Create the path to each fix file of each experiment of the same system
    files.path = joinpath.(files.PARENT_DIR, files.dir, files.file0)

    # Compute the assambly average of the fix file information
    fix_info = mean_FixSystem(files.path)

    # Get the id of the system
    id_system = unique(df[:, :id])

    # Add the id
    fix_info.id .= id_system

    # Create file name to store the data
    filename = string("fix_avg_", first(id_system), ".csv")

    # Directory to store the information
    file_dir  = joinpath(pwd(), "analyzedData")

    # File path to store the data
    file_path = joinpath(file_dir, filename)

    # Save the average
    CSV.write(file_path, fix_info)
end


