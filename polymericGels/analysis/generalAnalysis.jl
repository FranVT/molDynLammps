"""
    Script that execute all the analysis from a given system
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random

Random.seed!(1234)

# Include auxiliary files
include("functions.jl")


# Update the data file with new systems
up=0;
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
analyzeStructureFactor(data_bySystem[1])

