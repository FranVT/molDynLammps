#=
    Start of the Script
=#

using DataFrames, CSV

# Load the module
include("CreateParamsDF.jl")
using .CreateParamsDF

# Paths and directories
DIR_DATA = "/run/media/franvt/rogelio/DinMol/gel-batch-3/data/";
DIR_MAIN = pwd();
FILE_LOG = "params.log";
FILE_DAT = "dat.csv";

# Create a dataframe from all the experiments
df_dat=CreateParamsDF.collect_experiments_metadata(DIR_DATA)

# Save the dataframe of the metadata of the experiments 
CSV.write(joinpath(DIR_MAIN,FILE_DAT),df_dat)


