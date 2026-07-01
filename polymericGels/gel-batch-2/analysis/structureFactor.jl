"""
    Script en el que se adapta el código de FORTRAN de Claudia a Julia
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll
using SplitApplyCombine

include("functions.jl")

#=
    Functions
=#

#=
    Script
=#

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Main Parameters of the analysis
N_systems=length(data_bySystem);
N_instants=2;               # Instantes temporales a analizar
qmax0 = 6;              # 3 es el min sin que cause problemas

# Select one sysem
it_system=1;
dat_DF=data_bySystem[it_system];

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
        # Preparar informacion para el factor de estructura
        r_system = getPositions(file_names[1], dump_paths)

        # Calcular el factor de estructura para N experimentos del mismo sistema
        for it_exp in eachindex(r_system)
            r = reduce(hcat, r_system[it_exp])
            info[it_exp][:, 2:3] = computeSq(numbin, ntotav, qxhis, qyhis, qzhis, qhis, r)
            println("One experiment done")
        end

        # Se obtienen promedios
        global info = mean(info)        # Compute the assembly mean of the Sq
        global info[:, 1] = qmean       # Store the mean of the wave vectors

        # Paths and file names and stuff
        filename = string("structureFactorPBC", first(dat_DF.id), "time", timeSteps[it_time], ".csv")

        # Create the data frame
        df = DataFrame([repeat([timeSteps[it_time]], numbin) info],
                   [:timeStep, :qmean, :Sqmean, :Sqmeannorm])

        #computeSqmean(file_names[it_time], timeSteps[it_time], file_path, dump_paths, dat_DF,
        #              numbin, ntotav, qxhis, qyhis, qzhis, qhis, qmean, info)
        println("One time step done")
    end

    println("One system done")

