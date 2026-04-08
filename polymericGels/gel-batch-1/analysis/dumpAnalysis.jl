"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using Distances, Graphs
using GraphMakie, LinearAlgebra
using BenchmarkTools, UUIDs

# Include auxiliary files
include("functions.jl")

function getDump(dir,file_name)
"""
    Get the data from a single dump file that stores one timeste information
"""
    data = split.(readlines(joinpath(dir,file_name))," ")[9:end];
    HEADERS=data[1][3:end];
    INFO=parse.(Float64,reduce(hcat,data[2:end]))';

    return DataFrame(INFO,HEADERS)
end

function getDatInfo(DF_DIR)
"""
    This function get the dat.csv of experiments
"""

    archivos = filter(str -> occursin("dat-", str), DF_DIR)
    df_final = DataFrame()  # vacío inicial
    for (i, archivo) in enumerate(archivos)
        df_temp = CSV.read(archivo, DataFrame)
        append!(df_final, df_temp)
        df_temp = nothing  # liberar referencia
    end
    return df_final
end

# Get directories 
MAIN_DIR=dirname(pwd());
DATA_DIR=joinpath(MAIN_DIR,"data");
INFO_DIR=joinpath(pwd(),"data_mod");
DF_DIR=filter(isfile,readdir(INFO_DIR,join=true));

# Get the dat dataframes
dat_DF=getDatInfo(DF_DIR);

TRAJ_DIR=joinpath(dat_DF.dir[1],"traj");



