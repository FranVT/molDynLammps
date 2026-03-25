"""
    Basic analysis script
"""

using DataFrames, CSV
using Statistics, StatsBase
using Distances, Graphs
using GraphMakie, LinearAlgebra
using BenchmarkTools

# https://juliagraphs.org/Graphs.jl/v1.5/
# https://graphsjl-docs.readthedocs.io/en/latest/

# Load the functions
include("functions.jl")
include("graphs-functions.jl") # Includes the graphical packages

# Selection of an specific simulation
date="2026-03-25-093253";
#"2026-03-24-163630";
#"2026-03-24-102830";
#"2026-03-20-112940";
#"2026-03-20-105611";


#"2026-03-19-130353"; # COMPUTE CLUSTER
#"2026-03-20-105611";
#"2026-03-19-130353";
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

println("Imagenes del sistema listas")

# Cluster analysis
id_files=(0:Int(datConfig."save-dump"[1]):Int(datConfig."N_heat"[1]+datConfig."N_isot"[1]));
file_names=string.("traj_assembly.",id_files,".dumpf");

L=2*datConfig.L[1]; # Length of the box

# Reducción de la cantidad de time steps a analizar
N = 2^8  # cantidad deseada de puntos
log_min = log(1.0);
log_max = log(length(file_names));
pts = unique(round.(Int, exp.(range(log_min, log_max, length=N))))
pts = pts[pts .<= 100001]  # asegura límite superior

#=
println("Inicio de analisis de cluster")
info_cluster=map(pts) do s
    (N_clusters,cluster_Size)=clusterAnalysis(path_dump,file_names[s],L);   
end;
=#

println("Figuras de analisis de cluster")
fig_NumClusters(dt,id_files[pts],first.(info_cluster),path,date)
fig_MaxClusterPart(dt,id_files[pts],datConfig.Npart[1].-maximum.(last.(info_cluster)),path,date)


nothing



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


