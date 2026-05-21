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


# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);


# Argumentos de la función para analizar el factor de estructura
id_system=5;

# Dat file of the system to be analyzed
dat_DF=data_bySystem[id_system];

# Parameters of the system relevant to the structure factor
L=2*dat_DF.L[1];    # Longitud de la caja
N_instants=2;       # Instantes temporales a analizar
n_max=2^3;          # Magnitud máxima de cada componente
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
    # Select one time step
    file=file_names[it_time];

    # Get the position of all experiments of the same system at the same time step 
    # [x,y,z for all experiments]
    r_exp = getPositions(file, dump_paths);

    # Get all positions of all experiments at the same time step
    pos_exp=[reduce(hcat, r_exp[it_exp]) for it_exp in 1:N_exp];
    
    # Iterar entre magnitudes
    for (it_mag,mag) in enumerate(mag_n)
        n_vec=n[mag];           # Get the vectors with same magnitude
        n_dir=length(n_vec);    # Number of directions per magnitude
        rho=zeros(n_dir,N_exp); # Alocar memoria para la densidad

    
        # Iterar entre experimentos para la misma magnitud
        for it_exp in 1:N_exp
            pos=pos_exp[it_exp];    # Position of particles of one experimen of one experimentt
            # Calcular la densidad para todos los vectores de misma magnitud para un experimento
            for it_vec in eachindex(n_vec)
                pp=pos*(qo.*n_vec[it_vec]);             # Calcular producto
                rho[it_vec,it_exp]=computeDensity(pp);  # Obtener la densidad
            end
        end

    # Sq = mean(rho, dims=2) ./ N_part;  mean(Sq) == mean(rho) / N_part
      
    Sq_expval[it_mag,it_time]=mean(rho)/N_part; # Store the structure factor at the same magnitude

    end

end

    # Valor esperado de misma magnitud, distintas direcciones 
    #rho_expval=mean(rho)


# Create the wave vector for periodic boundary conditions

# Factor escala del vector de onda

nothing
