
using GLMakie
using CSV, DataFrames


# Paths and directories
DIR_MAIN = pwd();



# filename
filename_initial =  "pl_0.csv";
filename_final =    "pl_10000000.csv";

# Get the data
df_initial=CSV.read(joinpath(DIR_MAIN,filename_initial), DataFrame);
df_final=CSV.read(joinpath(DIR_MAIN,filename_final), DataFrame);

# Create categroies for the labels
#categories_system=[:phi,Symbol("CL-Con"),:Temperature,:damp];    # Select the categories that define a system
#categories_experiment=[:N_heat,:N_isot];  # Create categories to select different experiments (Just in case)

fig=Figure()
ax=Axis(fig[1:1,1:1])

scatterlines!(ax,df_initial.bin,df_initial.proportion)
scatterlines!(ax,df_final."bin_center",df_final.fraction)






