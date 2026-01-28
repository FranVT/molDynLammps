"""
    Basic analysis
"""

using DataFrames, CSV
#using Plots, LaTeXStrings, Plots.PlotMeasures
#gr()
using Statistics
using GLMakie, LaTeXStrings, Typst_jll

# Load the functions
include("functions.jl")

# Selection of an specific simulation
date="2026-01-28-105503";

# Get the directory of the desire system
DIR=getDir(date);
DIR=DIR[1];


# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
data=extractFixScalar(DIR,FILE_NAME);

# Convert the array into a DataFrame
DATA_fix=DataFrame(data[2]',data[1]);

# Get the directory of the desire system
DIR=joinpath(DIR,"traj");

# Filename with the simulation data
FILE_NAME="traj_assembly.2000.dumpf";

# Get data
DATA_dump=getDump(DIR,FILE_NAME);

# Select parameters to filter dataframe (patches)
# type 3: PA
# type 4: PB


