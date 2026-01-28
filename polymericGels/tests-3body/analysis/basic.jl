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


# Parameter to select the system
T=0.05;
N_particles=500;
phi=0.3;
CL_con=0.05;

# Selection of an specific simulation
date="2026-01-23-150058";
#"2026-01-22-165623";
#"2026-01-22-163251";
#"2026-01-22-154929";
#"2026-01-21-161335";
#"2026-01-21-133721";
#"2026-01-20-151755";
#"2026-01-20-143651";
#"2026-01-20-135923";
#"2026-01-20-121005"; 
#"2026-01-20-111704";

#0.050.30.05500-2026-01-20-121005

# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=DIR[1];


# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
data=extractFixScalar(DIR,FILE_NAME);

# Convert the array into a DataFrame
DATA=DataFrame(data[2]',data[1]);

# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=joinpath(DIR[1],"traj");


# Filename with the simulation data
FILE_NAME="traj_assembly.4000000.dumpf";

# Get data
data=getDump(DIR,FILE_NAME);

# Select parameters to filter dataframe (patches)
# type 3: PA
# type 4: PB

patches=data[(data.type.==4.0).|(data.type.==3.0),:];

central=data[(data.type.==1.0).|(data.type.==2.0),:];

