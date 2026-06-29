#=
    Script to analyze the infomration from the fix files
=#

using DataFrames, CSV

#=
    Functions
=#
"""
    extract_fix_scalar(path_file::String)
Function that extracts the information of fix files that stores global scalar values
"""
function extract_fix_scalar(path_file::String)
    aux = split.(readlines(path_file), " ")
    header = aux[2][2:end]
   
    aux_info = map(s -> parse.(Float64, s), aux[3:end])
    info = reduce(hcat, aux_info)

    return (header,info)
end

#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
FILE_DAT = "dat.csv";
FILE_FIX = "system_assembly.fixf";

# Read the dat file
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp];

# Group by system 
df_systems=groupby(df_dat,categories_system);

# Extract the data per system

# Select one system
df_system=df_systems[1];

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# Group by experiments
df_experiments=groupby(df_dat,categories_experiment);

# Select a set of experiments
df_set=df_experiments[1];

# Get the directories of all experiments of the system 
dir_set=String.(df_set.dir);

# Create the patchs to the files
path_fixf=joinpath.(dir_set,FILE_FIX);

# Extract

#extract_fix_scalar(path_file::String)
