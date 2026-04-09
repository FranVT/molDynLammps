"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra

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
data_system=meanFixystem(joinpath.(DATA_DIR,data_filtrada.dir));

# Ids
id=string(uuid4());
file_name_dat=string("dat-$phi-$Temp-$N_part-$CL_con.csv");
file_name_info=string("system_info_$id.csv");

# Add related paths
data_filtrada[!,:file_system].=[file_name_info];
data_filtrada[!,:file_name].=[file_name_dat];
data_filtrada[!,:id_ref].=[id];

# Save the data
CSV.write(joinpath(pwd(),"data_mod",file_name_info),data_system);
CSV.write(joinpath(pwd(),"data_mod",file_name_dat),data_filtrada);
