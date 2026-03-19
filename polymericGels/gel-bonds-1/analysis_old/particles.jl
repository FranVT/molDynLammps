"""
    Script to analyse the cluster.
    Read the dump information
"""

using DataFrames, CSV
#using Plots, LaTeXStrings, Plots.PlotMeasures
#gr()
using Statistics
using GLMakie, LaTeXStrings, Typst_jll

# Load the functions
include("functions.jl")


# Parameter to select the system
T=0.05;
N_particles=500;
phi=0.3;
CL_con=0.05;

# Selection of and specific simualtion
date="2026-01-23-150058";

# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=joinpath(DIR[1],"traj");


# Filename with the simulation data
FILE_NAME="traj_assembly.4000000.dumpf";

# Get data
data=getDump(DIR,FILE_NAME);

# Select parameters to filter dataframe (patches)
# type 3: PA
# type 4: PB

patches=data[(data.type.==4.0).|(data.type.==3.0),:];

central=data[(data.type.==1.0).|(data.type.==2.0),:];


# Test of the potentials
id_ref=maximum(patches.id);

patch_ref=patches[patches.id.==id_ref,:];
patch_other=patches[patches.id.!=id_ref,:];

# compute the distances
delta_x=patch_ref.x.-patch_other.x;
delta_y=patch_ref.y.-patch_other.y;
delta_z=patch_ref.z.-patch_other.z;

delta_r = sqrt.(delta_x.^2 .+ delta_y.^2 .+ delta_z.^2);

# Filter which is near the reference patch
filt_ind=findall(x->x<0.6,delta_r);

# Get the indixes for the DataFrame management
ids=map(s->patch_other[s,:].id,filt_ind)

# Get the dataframes
df_select=map(s->patches[patches.id.==s,:],ids)



