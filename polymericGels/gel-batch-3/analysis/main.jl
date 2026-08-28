#=
    Start of the Script
=#

include("CreateParamsDF.jl")

using DataFrames, CSV

# Load the module
using .CreateParamsDF

# Paths and directories
#DIR_DATA = "/home/franvt/GitRepos/molDynLammps/polymericGels/gel-batch-3/data/";
#"/run/media/franvt/rogelio/DinMol/gel-batch-3/data/";
DIR_DATA = "/run/media/franvt/rogelio/DinMol/gel-batch-3-long/data/";
DIR_MAIN = pwd();
FILE_LOG = "params.log";
FILE_DAT = "dat.csv";

# Create a dataframe from all the experiments
df_dat=collect_experiments_metadata(DIR_DATA)

# Add the time heat interval 
df_dat[!,:time_heat].=df_dat.tstep.*df_dat.N_heat

# Add the time isothermal interval 
df_dat[!,:time_isothermal].=df_dat.tstep.*df_dat.N_isothermal

# Save the dataframe of the metadata of the experiments 
CSV.write(joinpath(DIR_MAIN,FILE_DAT),df_dat)

# Run the fix analysis
include("FixInfoAnalysis")

# Create the figures
include("FixInfoFigures")


