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
MAIN_DIR=pwd();
DATA_DIR="/run/media/franvt/rogelio/DinMol/gel-batch-1/data/";
STORE_DIR=joinpath(MAIN_DIR,"datFiles");

# Read the directory where data is stored
SIMS_DIR=filter(isdir,readdir(DATA_DIR,join=true));

# Get all data files in a dataframe
data_info=mapreduce(s->getDat(s),vcat,SIMS_DIR);

# Add related paths
data_info[!,:PARENT_DIR].=[DATA_DIR];

# Create an id for each experiment
id=[
data_info.phi,
data_info."CL-Con",
data_info.Npart,
data_info.Temperature,
data_info.damp,
data_info."time-step",
data_info."N_heat",
data_info."N_isot",
data_info."N_CL",
data_info."N_MO",
data_info."L"
]

id=reduce(hcat,id);

id=[replace(join(string.(row)), "." => "") for row in eachrow(id)];

data_info[!,:id].=string.("id",id);

# Store the data
CSV.write(string(STORE_DIR,"/experiments_dat.csv"),data_info);

