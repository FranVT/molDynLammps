"""
    Multiple directories analysis
"""

using DataFrames, CSV
using Statistics, StatsBase
using Distances, Graphs
using GraphMakie, LinearAlgebra
using BenchmarkTools, UUIDs

# Include auxiliary files
include("functions.jl")

# Get all directories in the data directory
MAIN_DIR=dirname(pwd());
DATA_DIR=joinpath(MAIN_DIR,"data");

# Read the directory where data is stored
SIMS_DIR=filter(isdir,readdir(DATA_DIR,join=true));

# Get all data files in a dataframe
data_info=mapreduce(s->getDat(s),vcat,SIMS_DIR);

# Selection of the system by parameters
phi=0.02;
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

# Promedio de los N experimentos
data_system=meanFixystem(data_filtrada.dir);

# Save the data
id=string(uuid4());
file_name_dat=string("dat-$phi-$Temp-$N_part-$CL_con.csv");
file_name_info=string("system_info_$id.csv");

CSV.write(joinpath(pwd(),"data_mod",file_name_info),data_system);

data_filtrada[!,:file_system].=[file_name_info];
CSV.write(joinpath(pwd(),"data_mod",file_name_dat),data_filtrada)


