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

# Store the data
CSV.write(string(STORE_DIR,"/experiments_dat.csv"),data_info);

