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

