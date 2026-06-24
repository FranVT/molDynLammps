"""
    Scripts to create figures
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll


include("functions.jl")

# Include the theme of the plots and usefull functions
include("functions_graphs.jl")

# Main directory
MAIN_DIR=pwd();

# Directory of the data already processed
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:id];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Select the graphs to create
figure_fix=0 
figure_Sq_t=0

# Fix graph
if figure_fix == 1
    # Get the information of all systems of the fix files 
    dfs_fix = map(s->get_fixInfo(data_bySystem[s],SAVE_DIR),eachindex(data_bySystem));

    # Store the figure of the energy 
    figure_fixEnergy(dfs_fix,dat_files)
end

# Time evolution of the structure factor of a system 
if figure_Sq_t == 1
    # Get the information of all systems of the structure factor files
    dfs_sf = map(s->get_sfInfo(data_bySystem[s],SAVE_DIR),eachindex(data_bySystem));

    # ids
    ids = map(s->unique(data_bySystem[s].id)[1],eachindex(data_bySystem));

    # Store the figure of the time evolution of the structure factor for each system
    map(s->figureSq_time_evol(dfs_sf[s],ids[s]),eachindex(ids))
end



