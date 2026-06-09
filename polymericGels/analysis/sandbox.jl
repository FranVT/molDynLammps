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


dat_DF=data_bySystem[1];

# get the paths to the dump files foer each experiment
dump_paths=joinpath.(dat_DF.PARENT_DIR,dat_DF.dir,"traj");

paths=[readdir(it) for it in dump_paths];

# Time instants store in each directory
time_instants=[[split(it,".")[2] for it in p] for p in paths];












N_instants=1;


#function getpathfilesSq(dat_DF,N_instants)
"""
    Get the names files 
"""
     # Create time steps range
    #aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));

    # Get the index for the time instants that we are interested to analyzed
    #ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));

    # Get the time steps
    #timeSteps=aux_timeStep[ind];

    # Create file names
    #file_names=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in timeSteps];

    # Path to the dumps
#     return (timeSteps,dump_paths,file_names)
#end




nothing
