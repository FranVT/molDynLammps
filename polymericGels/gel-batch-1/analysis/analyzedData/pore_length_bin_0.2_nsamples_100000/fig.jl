
using GLMakie
using CSV, DataFrames

# Paths and directories
DIR_MAIN = pwd();

# filename
filename_initial =  "pore_length_histogram_0.010.050.051.0500000.09.5e6_step_0_simulation_1.csv ";
filename_final =    "pore_length_histogram_0.010.050.051.0500000.09.5e6_step_10000000_simulation_1.csv";

# Get the data
df_initial=CSV.read(joinpath(DIR_MAIN,filename_initial), DataFrame);
df_final=CSV.read(joinpath(DIR_MAIN,filename_initial), DataFrame);

# Create categroies for the labels
categories_system=[:phi,Symbol("CL-Con"),:Temperature,:damp];    # Select the categories that define a system
categories_experiment=[:N_heat,:N_isot];  # Create categories to select different experiments (Just in case)



