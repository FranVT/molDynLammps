"""
    Basic analysis script
"""

using DataFrames, CSV
using Statistics, StatsBase

# Load the functions
include("functions.jl")
include("graphs-functions.jl") # Includes the graphical packages

# Parameter to select the system
T=0.05;
N_particles=5000;
phi=0.15;
CL_con=0.05;

# Selection of an specific simulation
# 0.050.150.055000-2026-03-06-175535
date="2026-03-06-175535";
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

clusters=map(collect(groupby(data,:c_clusters,sort=true))) do s
    s[(s.type.==1.0).|(s.type.==2.0),:]
end

# Start creating the list of neighbors to analuze distances and stuff
pos=[clusters[1].x clusters[1].y clusters[1].z]

sum(map(s->evaluate(Euclidean(),pos[1,:],pos[s,:]),2:length(pos[:,1])).<=1.2)



#=

# The amount of CL and MO
nrow.(typeCl1|>collect)





patches=data[(data.type.==4.0).|(data.type.==3.0),:];

# Count all same id of cl_clusters (Amount of clusters)
N_patchclusters=values(countmap(patches.c_clusters));

# Count the size of each cluster
cl_patches=countmap(N_patchclusters);
=#

