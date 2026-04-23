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
phi=0.02;
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
N_qu=2^2; # Cantidad de direcciones
lambda_o=1; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=2^3; # Cantidad de magnitudes

# Seleccion de time instants
aux_id=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
aux_id=aux_id[1:100:end];

time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];

Sq_t=[getTimeEvolSq(dump_paths,time_instant) for time_instant in time_instants];


Sq_df=DataFrame(;
    timeStep = aux_id,
    Sq       = reduce.(vcat, Sq_t),
    lambda_o = fill(lambda_o, length(aux_id)),
    lambda_f = fill(lambda_f, length(aux_id)),
    id=unique(dat_DF.id)[1]
)

# CSV.write(joinpath(SAVE_DIR,"structureFactor.csv"),Sq_df)


CSV.write(joinpath(SAVE_DIR,string("structureFactor",Sq_df.id[1],".csv")),Sq_df)
