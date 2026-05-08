"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra
using Distances, Random
using GLMakie, JLD

Random.seed!(1234)

# Include auxiliary files
include("functions.jl")

# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","experiments_dat.csv");
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

dat_files=CSV.read(DAT_PATH,DataFrame);

# Selection of the system by parameters
phi=0.01;
Temp=0.05;
N_part=5000.0;
CL_con=0.05;

# Se filtra el dataframe 
dat_DF = subset(dat_files,
    :phi => ByRow(==(phi)),
    :Temperature => ByRow(==(Temp)),
    :Npart => ByRow(==(N_part)),
    :"CL-Con" => ByRow(==(CL_con))
)

# Path to the dumps
dump_paths=joinpath.(dat_DF.PARENT_DIR,dat_DF.dir,"traj");

# Parametros para obtener el factor de estructura
N_qu=2^9; # EXPONENTE DEBE SER IMPAR Cantidad de direcciones
lambda_o=1; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=2^8; # Cantidad de magnitudes
N_instants=2;

# Seleccion de time instants
aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));
aux_id=aux_timeStep[ind];

time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];

time_instant=time_instants[end];

N_part=Int(dat_DF.Npart[1]);

S_q=computeStructureFactor(N_qu,N_lambda,lambda_o,lambda_f,time_instant,dump_paths);

q_min=2*pi/lambda_f;
q_max=2*pi/lambda_o;
q_dom=range(q_min,q_max,length=N_lambda);


fig=Figure()
ax=Axis(fig[1:1,1:1],
        limits=(nothing,nothing,0,10)
       )
#vlines!(ax,2*pi/(1.2),linestyle=:dash,color=:blue)
#vlines!(ax,2*pi/(2*1.2),linestyle=:dash,color=:blue)
#vlines!(ax,2*pi/(3*1.2),linestyle=:dash,color=:blue)

#vlines!(ax,2*pi/(lambda_f),linestyle=:solid,color=:black)
#vlines!(ax,2*pi/(0.5*lambda_f),linestyle=:solid,color=:black)
#vlines!(ax,2*pi/(0.25*lambda_f),linestyle=:solid,color=:black)
#vlines!(ax,2*pi/(0.125*lambda_f),linestyle=:solid,color=:black)
#vlines!(ax,2*pi/(0.0625*lambda_f),linestyle=:solid,color=:black)

scatterlines!(ax,(2pi)./q_dom,S_q)

