"""
    Script en el que se adapta el código de FORTRAN de Claudia a Julia
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll
using SplitApplyCombine

include("functions.jl")

#=
    Functions
=#

#=
    Script
=#

# Paths and directories
DIR_MAIN = pwd();
DAT_PATH   = joinpath(DIR_MAIN, "datFiles")
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "systemDatfiles.csv";
FILE_FIX = "system_assembly.fixf";
FILE_DUMP = "traj_assembly.*.dumpf";

# Select the amount of time steps to analyze the structure factor
n_steps_Sq=8;

# Define the maximum wave vector
qmax_0=6;

# Select the categories that define a system
categories_system=[:phi,Symbol("CL-Con"),:Temperature,:damp];

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isot];

# Read the dat file
df_dat=CSV.read(joinpath(DAT_PATH,FILE_DAT), DataFrame);

# Group by system 
df_systems=groupby(df_dat,categories_system);



