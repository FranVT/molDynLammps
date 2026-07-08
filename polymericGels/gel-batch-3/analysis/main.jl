#=
    Start of the Script
=#

include("CreateParamsDF.jl")

using DataFrames, CSV

# Load the module
using .CreateParamsDF

# Paths and directories
DIR_DATA = "/home/franvt/GitRepos/molDynLammps/polymericGels/gel-batch-3/data/";
#"/run/media/franvt/rogelio/DinMol/gel-batch-3/data/";
DIR_MAIN = pwd();
FILE_LOG = "params.log";
FILE_DAT = "dat.csv";

# Create a dataframe from all the experiments
df_dat=collect_experiments_metadata(DIR_DATA)

# Save the dataframe of the metadata of the experiments 
CSV.write(joinpath(DIR_MAIN,FILE_DAT),df_dat)

