"""
    Script to debug stuff
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll
using SplitApplyCombine

include("functions.jl")

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Parameters
N_instants=25; # Select the amount of timesteps to analyze

# Cycle thru the different systems
map(s->clusterAnalysis(data_bySystem[s],N_instants),eachindex(data_bySystem))


nothing
