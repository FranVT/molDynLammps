"""
    Multiple directories analysis
"""

using DataFrames, CSV
using Statistics, StatsBase
using Distances, Graphs
using GraphMakie, LinearAlgebra
using BenchmarkTools

# Include auxiliary files
include("functions.jl")

# Get all directories in the data directory
MAIN_DIR=dirname(pwd());
DATA_DIR=joinpath(MAIN_DIR,"data");

# Read the directory where data is stored
SIMS_DIR=filter(isdir,readdir(DATA_DIR,join=true));

# Get all data files in a dataframe
df=mapreduce(s->getDat(s),vcat,SIMS_DIR);

