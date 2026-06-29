#=
    Script to create dat files from the information stored 
=#

using DataFrames, CSV

# Directory where the data is
dir_data="/run/media/franvt/rogelio/DinMol/gel-batch-3/data/";

# This directory
dir_main=pwd();

# Read the directory and find subdirectories
dirs_info=filter(isdir, readdir(dir_data, join = true));
