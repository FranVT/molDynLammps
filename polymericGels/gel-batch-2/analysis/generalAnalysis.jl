"""
    Script that execute all the analysis from a given system
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using SplitApplyCombine

Random.seed!(1234)

# Include auxiliary files
include("functions.jl")

# Activate functions
up=0;
selec_fixavg=0;
selec_Sq=0;
selec_Cluster=1;

# Assign a unique id to the dat files in the directory 
if up ==1 
    createDat_ID("/run/media/franvt/rogelio/DinMol/gel-batch-2/data");
end

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:id];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

if selec_fixavg == 1 
    # Store all the assembly avg of the systems
    #map(s->storeAvg_fix(data_bySystem[s]),eachindex(data_bySystem))
    storeAvg_fix(data_bySystem[2])
end

if selec_Sq == 1
    # Main Parameters of the analysis
    N_systems=length(data_bySystem);
    N_instants=25;               # Instantes temporales a analizar
    qmax0 = 6;              # 3 es el min sin que cause problemas

    for it_system in (2) 
        computeAllTimeSqmean(data_bySystem[it_system],N_instants,qmax0)
    end
end

if selec_Cluster == 1
    # Parameters
    N_instants=25; # Select the amount of timesteps to analyze

    # Cycle thru the different systems
    map(s->clusterAnalysis(data_bySystem[s],N_instants),eachindex(data_bySystem))
end
