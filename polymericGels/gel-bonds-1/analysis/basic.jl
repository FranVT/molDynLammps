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

# Selection of an specific simulation
date="2026-03-19-130353";
path=getDir(date);
path_dump=joinpath(path,"traj");

# Read the data file
datConfig=getDatFile(path);

# Extract the data from the file
dataSystem=extractFixScalar(path,"system_assembly.fixf");


dt=datConfig."time-step"[1];



# Create and save the graphics
fig_Temp(dt,dataSystem,datConfig."Temperature"[1],path,date);
fig_Eng(dt,dataSystem,path,date);
fig_EngB(dt,dataSystem,path,date);
fig_EngSys(dt,dataSystem,path,date);
fig_EngPair(dt,dataSystem,path,date);



# Filename with the simulation data
FILE_NAME="traj_assembly.9000000.dumpf";

id_files=(0:Int(datConfig."save-dump"[1]):Int(datConfig."N_heat"[1]+datConfig."N_isot"[1]));
file_names=string.("traj_assembly.",id_files,".dumpf");

L=2*datConfig.L[1]; # Length of the box

# Cluster analysis
function clusterAnalysis(DIR,FILE_NAME,L)
"""
    Function that performs a quick cluster analysis.
    Return a dataframe with the position of the particles and neighbors and stuff
"""

    data=getDump(DIR,FILE_NAME);

    # Get one dataframe per cluster in the system with its neighbors
    clusters=getClusters(data,L); # Just central particles. Patches have been descarted

    N_clusters=length(clusters);    # Amount of clusters in the system
    cluster_Size=nrow.(clusters);   # Amount of central particles in each cluster in the system

    # Get one graph for each cluster (strand length and loops)
    graphs=map(s->createGraph(clusters[s]),eachindex(clusters))

    return (clusters,N_clusters,cluster_Size,graphs)

return 


(clusters,N_clusters,cluster_Size,graphs)=clusterAnalysis(path_dump,file_names[1],L)

#files=readdir(joinpath(path,"traj"),join=true);



#=
# Get data
# Structure factor
l=1.5; # "wave length"
vec_K=(2*pi).*[1/l,1/l,1/l]; # vector
# Structure factor per cluster
I_k=map(s->structureFactor([clusters[s].x clusters[s].y clusters[s].z],vec_K),eachindex(clusters));

#Structure factor of all the system
I_K=sum(I_k);





=#


