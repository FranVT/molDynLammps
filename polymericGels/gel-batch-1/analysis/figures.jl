"""
    Scripts to create figures
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll


include("functions.jl")

# Include the theme of the plots and usefull functions
include("functions_graphs.jl")

# Main directory
MAIN_DIR=pwd();

# Directory of the data already processed
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:id];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Select the graphs to create
figure_fix=0 
figure_Sq_t=0

# Fix graph
if figure_fix == 1
    # Get the information of all systems of the fix files 
    dfs_fix = map(s->get_fixInfo(data_bySystem[s],SAVE_DIR),eachindex(data_bySystem));

    # Store the figure of the energy 
    figure_fixEnergy(dfs_fix,dat_files)
end

# Time evolution of the structure factor of a system 
if figure_Sq_t == 1
    # Get the information of all systems of the structure factor files
    dfs_sf = map(s->get_sfInfo(data_bySystem[s],SAVE_DIR),eachindex(data_bySystem));

    # ids
    ids = map(s->unique(data_bySystem[s].id)[1],eachindex(data_bySystem));

    # Store the figure of the time evolution of the structure factor for each system
    map(s->figureSq_time_evol(dfs_sf[s],ids[s]),eachindex(ids))
end

# Cluster figure

# Get the cluster paths for each system
cluster_paths=map(s->get_cluster_paths(data_bySystem[s]),eachindex(data_bySystem));

# Read the information of the fix files
cluster_df=map(s->CSV.read(s,DataFrame),cluster_paths);

    # Prepare the color code/label
    phi_label=(100).*map(s->first(data_bySystem[s].phi),eachindex(data_bySystem))

    color_ref=maximum(phi_label);
    color_norm=phi_label./color_ref;
    color_min=minimum(color_norm);
    color_max=maximum(color_norm);

    # Prepare data for the figure
    time_domains=map(s->s.timeStep,cluster_df);
    dt=0.001;
    time_domains=(dt).*time_domains;



# Number of clusters 
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Time}~[\tau]",
                   ylabel = L"\langle\mathrm{Number~of~clusters}\rangle",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   #limits = (nothing,nothing,nothing,nothing),
                   limits = (10^2, 10^4.5, nothing, nothing),
                   xscale = log10,
                   yscale = log10
                  )

    # Prepare the y axis
    N_clusters=map(s->s.nClusters,cluster_df);

    for it_system in eachindex(phi_label)
        scatterlines!(ax_plot, time_domains[it_system], N_clusters[it_system],
                          color = color_norm[it_system],
                          colorrange = (color_min, color_max)
                         )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi\%", limits = (minimum(phi_label), maximum(phi_label)),
             colormap = cgrad(:viridis, length(phi_label), categorical = true), size = 25)

    save(string("fig_NumClusters_phiseries.png"), fig, px_per_unit = 300 / inch)

# Free particles
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Time}~[\tau]",
                   ylabel = L"\langle\mathrm{Free~particles}\rangle",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   #limits = (nothing,nothing,nothing,nothing),
                   limits = (10^2, 10^4.5, nothing, nothing),
                   xscale = log10,
                   yscale = log10
                  )

    # Prepare the y axis
    free_particles=map(s->s.maxParticles,cluster_df);

    for it_system in eachindex(phi_label)
        scatterlines!(ax_plot, time_domains[it_system], (5000).-free_particles[it_system],
                          color = color_norm[it_system],
                          colorrange = (color_min, color_max)
                         )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi\%", limits = (minimum(phi_label), maximum(phi_label)),
             colormap = cgrad(:viridis, length(phi_label), categorical = true), size = 25)

    save(string("fig_FreeParticles_phiseries.png"), fig, px_per_unit = 300 / inch)

# Histogram of number of particles in each clusters
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Size~of~cluster}",
                   ylabel = L"\mathrm{pdf}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (nothing,nothing,nothing,nothing),
                   #limits = (10^2, 10^4.5, nothing, nothing),
                   #xscale = log10,
                   #yscale = log10
                  )
    
    # Selection of the time step
    it_time=25; 
    time_label=dt*cluster_df[1].timeStep[it_time];

    # Prepare the range 
    hist_timestep=eval.(Meta.parse.(map(s->s.hist[it_time],cluster_df)));

    # Select the step for the historgrams
    bins_range=range(1,5000,length=120);

    for it_system in eachindex(phi_label)
        hist!(ax_plot, hist_timestep[it_system],
              bins=bins_range,
              color = color_norm[it_system],
              alpha=0.8,
              colorrange = (color_min, color_max),
              normalization=:pdf
              )
    end

    scatter!(ax_plot,0,0,color=(:white,0.1),label=latexstring("t=",time_label,"~[\\tau]"))

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi\%", limits = (minimum(phi_label), maximum(phi_label)),
             colormap = cgrad(:viridis, length(phi_label), categorical = true), size = 25)

    axislegend(position=:rt)

    save(string("fig_HistSize_phiseries.png"), fig, px_per_unit = 300 / inch)

# Time evolution of the histograms per concentration 

    # Selection of the concentration to plot
    it_phi=5;


    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Size~of~cluster}",
                   ylabel = L"\mathrm{pdf}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (nothing,nothing,nothing,nothing),
                   #limits = (10^2, 10^4.5, nothing, nothing),
                   #xscale = log10,
                   #yscale = log10
                  )

# Rename the colorbar
    time_domains=map(s->s.timeStep,cluster_df);
    dt=0.001;
    time_domains=(dt).*time_domains;

    color_ref=maximum(time_domains[1]);
    color_norm=time_domains[1]./color_ref;
    color_min=minimum(color_norm);
    color_max=maximum(color_norm);


    # Prepare the color code/label
    phi_labels=(100).*map(s->first(data_bySystem[s].phi),eachindex(data_bySystem));
    phi_label=phi_labels[it_phi];

    # Prepare the range 
    hist_timeevol=eval.(Meta.parse.( cluster_df[it_phi].hist )) #map(s->s.hist[it_time],cluster_df)));

    # Select the step for the historgrams
    bins_range=range(1,5000,length=120);

    for it_system in eachindex(color_norm)
        hist!(ax_plot, hist_timeevol[it_system],
              bins=bins_range,
              color = color_norm[it_system],
              alpha=0.8,
              colorrange = (color_min, color_max),
              normalization=:pdf
              )
    end

    scatter!(ax_plot,0,0,color=(:white,0.1),label=latexstring("\\phi=",phi_label,"\\%"))

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"t~[\tau]", limits = (minimum(time_domains[1]), maximum(time_domains[end])),
             colormap = cgrad(:viridis, length(color_norm), categorical = true), size = 25)

    axislegend(position=:rt)

    save(string("fig_HistSize_timeseries_phi",phi_label,".png"), fig, px_per_unit = 300 / inch)





