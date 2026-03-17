"""
    Basic analysis script
"""

using DataFrames, CSV
using Statistics, StatsBase
using Distances, Graphs
using GraphMakie, LinearAlgebra

# https://juliagraphs.org/Graphs.jl/v1.5/
# https://graphsjl-docs.readthedocs.io/en/latest/

# Load the functions
include("functions.jl")
include("graphs-functions.jl") # Includes the graphical packages

# Parameter to select the system
T=0.05;
N_particles=500;
phi=0.01;
CL_con=0.05;

# Selection of an specific simulation
# 0.050.150.055000-2026-03-06-175535
date="2026-03-17-145727";

#"2026-03-06-175535";
#"2026-03-06-112337";
#"2026-03-05-124121";

#"2026-01-23-150058";
#"2026-01-22-165623";
#"2026-01-22-163251";
#"2026-01-22-154929";
#"2026-01-21-161335";
#"2026-01-21-133721";
#"2026-01-20-151755";
#"2026-01-20-143651";
#"2026-01-20-135923";
#"2026-01-20-121005"; 
#"2026-01-20-111704";

#0.050.30.05500-2026-01-20-121005


# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=DIR[1];

# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
data=extractFixScalar(DIR,FILE_NAME);

# Convert the array into a DataFrame
DATA=DataFrame(data[2]',data[1]);

"""
    P L O T S 
"""

damp=1;
dt=0.001;

# Create and save the graphics
fig_Temp(dt,DATA,T,DIR,id_c);
fig_Eng(dt,DATA,DIR,id_c);
fig_EngB(dt,DATA,DIR,id_c);
fig_EngSys(dt,DATA,DIR,id_c);
fig_EngPair(dt,DATA,DIR,id_c)


# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=joinpath(DIR[1],"traj");

# Filename with the simulation data
FILE_NAME="traj_assembly.9000000.dumpf";

# Get data
data=getDump(DIR,FILE_NAME);

L=2*14.984222; # Length of the box

# Get one dataframe per cluster in the system with its neighbors
clusters=getClusters(data,L); # Just central particles. Patches have been descarted

N_clusters=length(clusters);    # Amount of clusters in the system
cluster_Size=nrow.(clusters);   # Amount of central particles in each cluster in the system

# Get one graph for each cluster (strand length and loops)
graphs=map(s->createGraph(clusters[s]),eachindex(clusters))

# Structure factor


function structureFactor(vec_r,vec_K)
"""
    Get the factor of structure of a set of points
"""

    # Dot product between the wave vector and the position vector
    dot_kr=map(s->dot(vec_K,vec_r[s,:]),eachindex(vec_r[:,1]));
    I_k=sum(exp.(im.*dot_kr))

   return I_k
end

l=1.5; # "wave length"
vec_K=(2*pi).*[1/l,1/l,1/l]; # vector
# Structure factor per cluster
I_k=map(s->structureFactor([clusters[s].x clusters[s].y clusters[s].z],vec_K),eachindex(clusters));

#Structure factor of all the system
I_K=sum(I_k);




