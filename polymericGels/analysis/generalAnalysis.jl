"""
    Script that execute all the analysis from a given system
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra
using Random

Random.seed!(1234)

# Include auxiliary files
include("functions.jl")

# Get all directories in the data directory
MAIN_DIR=dirname(pwd());
DATA_DIR="/run/media/franvt/rogelio/DinMol/gel-batch-1/data";
#joinpath(MAIN_DIR,"data");

# Read the directory where data is stored
SIMS_DIR=filter(isdir,readdir(DATA_DIR,join=true));

# Get all data files in a dataframe
data_info=mapreduce(s->getDat(s),vcat,SIMS_DIR);

# Selección de categorias
categories=[:phi,Symbol("CL-Con"),:Temperature,:Npart,:damp,Symbol("time-step")];

# Creación de los subdataframes
data_bySystem=groupby(data_info,categories);

# Obtención de los valores numéricos de las categorias
id_values=unique(data_info[:,categories]);

# Creación de los ids


#=
# Selection of the system by parameters
phi=0.05;
Temp=0.05;
N_part=5000.0;
CL_con=0.05;

# Se filtra el dataframe 
data_filtrada = subset(data_info,
    :phi => ByRow(==(phi)),
    :Temperature => ByRow(==(Temp)),
    :Npart => ByRow(==(N_part)),
    :"CL-Con" => ByRow(==(CL_con))
)
=#

