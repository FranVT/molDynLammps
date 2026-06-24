"""
    Script to debug stuff
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll
using SplitApplyCombine

include("functions.jl")

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Parameters
N_instants=25; # Select the amount of timesteps to analyze

function clusterAnalysis(dat_DF)
# Get the paths to the data
(timeSteps, dump_paths, file_names) = getpathfilesSq(dat_DF, N_instants);

# Get the position of all the experiments of the same system 
N_exp=length(dat_DF.Nexp);  # Amount of experiments of the same system

# Save memory
size_clusters=[[] for _ in 1:N_exp];
N_clusters=[zeros(N_instants) for _ in 1:N_exp];

# Get the amount of clusters at each time step for all experiments
# $ C_clusters is an id
for it_exp in 1:N_exp
    # Get the position of all time steps for one experiment
    positions=map(s->getDump(dump_paths[it_exp], s),file_names);

    # Get the amount of central particles in the clusters 
    filtered = filter.(row -> row.type in (1, 2), positions);
    counts = map(s->combine(groupby(s, :c_clusters), nrow => :count),filtered);

    # Store the information
    size_clusters[it_exp]=map(s->s.count,counts);

    # Get the amount of clusters for each time step 
    N_clusters[it_exp]=length.(size_clusters[it_exp]);
end

# Compute the assemble average of amount of clusters in each timestep
# row|col -> time|phi
N_clusters_mean=mean(N_clusters);   

# Compute the assemble average of the amount of central particles in the biggest cluster in each timestep
# row|col -> time|phi
max_particles=mean(map(s-> maximum.(s),size_clusters));

# For the histogram 
hist_size=map(l->mapreduce(s->size_clusters[s][l],vcat,1:N_exp),eachindex(timeSteps))

# Save the information per timeStep

# Directory to store the information
file_dir  = joinpath(pwd(), "analyzedData")

# Paths and file names and stuff
filename = string("clusterAnalysis", first(dat_DF.id),".csv")

# File path to store the data
file_path = joinpath(file_dir, filename)

# Create the DataFrame
df = DataFrame([timeSteps N_clusters_mean max_particles hist_size],
                   [:timeStep, :nClusters, :maxParticles, :hist])
# Save the data frame
CSV.write(file_path, df)

println("One cluster analysisi saved")

end

# Cycle thru the different systems
map(s->clusterAnalysis(data_bySystem[s]),eachindex(data_bySystem))


nothing
