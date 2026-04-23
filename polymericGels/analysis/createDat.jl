"""
    Script that search for simulation directories at a given path.
    Then create a dat.csv file.
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra

# Include auxiliary files
include("functions.jl")

# Get all directories in the data directory
MAIN_DIR=dirname(pwd());
DATA_DIR="/run/media/franvt/rogelio/DinMol/gel-batch-1";
STORE_DIR=joinpath(MAIN_DIR,"datFiles");



# Read the directory where data is stored
SIMS_DIR=filter(isdir,readdir(DATA_DIR,join=true));

# Get all data files in a dataframe
data_info=mapreduce(s->getDat(s),vcat,SIMS_DIR);

# Ids
id=string(uuid4());
file_name_dat=string("dat-$phi-$Temp-$N_part-$CL_con.csv");
file_name_info=string("system_info_$id.csv");

# Save the data
CSV.write(joinpath(pwd(),"data_mod",file_name_info),data_system);
CSV.write(joinpath(pwd(),"data_mod",file_name_dat),data_filtrada);
