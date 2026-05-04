"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra

# Include auxiliary files
include("functions.jl")

# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","experiments_dat.csv");
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

dat_files=CSV.read(DAT_PATH,DataFrame);

# Selection of the system by parameters
phi=0.05;
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
N_qu=2^7; # EXPONENTE DEBE SER IMPAR Cantidad de direcciones
lambda_o=0.5; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=2^9; # Cantidad de magnitudes
N_instants=2;

# Seleccion de time instants
aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));
aux_id=aux_timeStep[ind];

time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];

Sq_t=[getTimeEvolSq(N_qu,dump_paths,time_instant) for time_instant in time_instants];

q_domain=[Sq_t[it][1] for it in eachindex(Sq_t)];
Sq=[Sq_t[it][2] for it in eachindex(Sq_t)];


Sq_df=DataFrame(;
    timeStep = aux_id,
    q_domain = q_domain,
    Sq       = Sq,
    lambda_o = fill(lambda_o, length(aux_id)),
    lambda_f = fill(lambda_f, length(aux_id)),
    id=unique(dat_DF.id)[1]
)

#CSV.write(joinpath(SAVE_DIR,string("structureFactor",Sq_df.id[1],".csv")),Sq_df)
