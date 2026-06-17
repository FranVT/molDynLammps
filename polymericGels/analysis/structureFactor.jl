"""
    Script en el que se adapta el código de FORTRAN de Claudia a Julia
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

# Main Parameters of the analysis
N_systems=length(data_bySystem);
N_instants=25;               # Instantes temporales a analizar
qmax0 = 6;              # 3 es el min sin que cause problemas

for it_system in (1,2,3,5) 
    computeAllTimeSqmean(data_bySystem[it_system],N_instants,qmax0)
end

