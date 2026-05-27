"""
    Script to debug stuff
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll
using SplitApplyCombine

include("functions.jl")

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
    n_max=2^2;          # Magnitud máxima de cada componente
    N_exp=nrow(dat_DF); # Cantidad de experimentos por sistema
    N_part=dat_DF.Npart[1];
    qo=2*pi/L;

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

    for (it_time,~) in enumerate(timeSteps)
        file=file_names[it_time];   # Select one time step
     
        r_exp = getPositions(file, dump_paths); # Get the position of all experiments of the same system at the same time step
        pos_exp=[reduce(hcat, r_exp[it_exp]) for it_exp in 1:N_exp];    # Get all positions of all experiments at the same time step
    
        # Calcular y almacenar S(q) para todas las magnitudes
        Sq_expval[:, it_time] = structureFactor(pos_exp, mag_n, n, qo, N_part; N_exp=N_exp);

    end
    
    return [mag_n Sq_expval]

end


#function avgSq(dat_DF)
#    return mean([dumpAnalysis(dat_DF) for s in dat_DF])
#end


# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Compute the analysis from dump files for all the systems
Sq_all=[dumpAnalysis(dat_DF) for dat_DF in data_bySystem];

# Pasar los resultados a un DataFrame
df=map(s->DataFrame(s,[:nMag, :Sqo, :Sqf]),Sq_all);


nothing
