"""
    Basic analysis script
"""

using DataFrames, CSV
#using Plots, LaTeXStrings, Plots.PlotMeasures
#gr()
using Statistics
using GLMakie, LaTeXStrings


# Get the directories
MAIN_DIR=dirname(pwd());
DATA_DIR=joinpath(MAIN_DIR,"data");

# Read the directory where data is stored
sims=filter(isdir,readdir(DATA_DIR,join=true));

# Parameter to select the system
T=0.05;
N_particles=500;
phi=0.3;
CL_con=0.05;

# Selection of an specific simulation
date="2026-01-20-111704";

# Creation of the id
id=string(T,phi,CL_con,N_particles)
id_c=string(T,phi,CL_con,N_particles,"-",date)

# filter!(s->s==id_c,readdir(DATA_DIR))

# Get the idex for the directory
indx=findall(!isempty,findall.(id_c,sims))

# Get the directory
DIR=sims[indx]

