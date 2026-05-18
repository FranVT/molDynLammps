"""
    Script that execute all the analysis from a given system
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random

Random.seed!(1234)

# Include auxiliary files
include("functions.jl")

# Activate functions
up=0;
Sq=0;

# Update the data file with new systems
if up ==1 
    createDatFiles();
end

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Analyze the structure factor per system
if Sq == 1
    map(s->analyzeStructureFactor(data_bySystem[s]),1:length(data_bySystem));
end


