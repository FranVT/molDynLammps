#=
    Script to create the figures of the structure factor analysis
=#

using DataFrames, CSV
using GLMakie, LaTeXStrings 
#=
    Set figures 
=#
INCH = 96;
PT = 4/3;
CM = INCH / 2.54;

set_theme!(
    backgroundcolor = :white,
    fontsize = 16PT,                     # tamaño de letra base
    Axis = (
        titlesize = 20PT,                 # tamaño del título del eje
        xlabelsize = 14PT,                # tamaño etiqueta eje x
        ylabelsize = 14PT,                 # tamaño etiqueta eje y
        xticklabelsize = 12PT,             # tamaño números eje x
        yticklabelsize = 12PT,              # tamaño números eje y
        xgridstyle = :dash,               # estilo de la cuadrícula
        ygridstyle = :dash,
        spinewidth = 1.5PT,
    ),
    Legend = (
        labelsize = 10PT,                   # tamaño texto leyenda
        framewidth = 1.5PT,
    ),
    Colorbar = (
        labelsize = 14PT,
        ticklabelsize = 12PT,
    ),
    Figure = (
        size = (15CM, 12CM)                 # tamaño de la figura (ancho, alto)
    )
)

#=
    Functions 
=#
"""
    fig_sq_timeseries_loglog(df_systems::AbstractDataFrame, categories_figures::Vector{Symbol})

Function that creates a figure with the time evolution of a structure factor
"""
function fig_sq_timeseries_loglog(df_systems::GroupedDataFrame{DataFrame}, categories_figures::Vector{Symbol})

    # Create the figure for the time evolution for all systems 
for df_system in df_systems
    # Get the time domain
    system_time_domain=sort((DT).*unique(df_system.timeStep));

    # Get the data for the labels
    labels_plot=[df_system[1, col] for col in categories_figures];
        
    # Group the system by time domain
    df_system_evolution=groupby(df_system,:timeStep);

    # Prepare some labels
    time_min=minimum(system_time_domain);
    time_max=maximum(system_time_domain);

    # Prepare the color code
    color_norm  = system_time_domain ./ time_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

    # Starts the figure
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle ",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, 10^(-1.5), nothing),
                   xscale = log10,
                   yscale = log10
                  )

    for df_plot in df_system_evolution
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean[:];

        # Get the color
        color_label=first(DT*unique(df_plot.timeStep))/time_max;

        # Plot
        scatterlines!(ax_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end

        # Colobar to denote the time evolution 
        Colorbar(fig[1, 2], label = L"\tau", colormap = :viridis, limits = (time_min, time_max))

        

        # Plots to add systems labels
        plot!(ax_plot,[100],[100],label=latexstring("\\phi=",100*labels_plot[1],"\\%"),color=:black)
        plot!(ax_plot,[100],[100],label=latexstring("T=",labels_plot[3]),color=:black)
        axislegend(ax_plot,position=:lt)

        # Crate a file name with the labels
        file_name=string("fig_Sq_loglog",join(string.(labels_plot)),".png");

        save(file_name, fig, px_per_unit = 300 / INCH)
end

end






"""
    fig_sq_norm_timeseries_loglog(df_systems::AbstractDataFrame, categories_figures::Vector{Symbol})

Function that creates a figure with the time evolution of a structure factor normalize to a maximum
"""
function fig_sq_norm_timeseries_loglog(df_systems::GroupedDataFrame{DataFrame}, categories_figures::Vector{Symbol})
# Create the figure for the time evolution for all systems 
for df_system in df_systems

    # Get the time domain
    system_time_domain=sort((DT).*unique(df_system.timeStep));

    # Get the data for the labels
    labels_plot=[df_system[1, col] for col in categories_figures];
        
    # Group the system by time domain
    df_system_evolution=groupby(df_system,:timeStep);

    # Prepare some labels
    time_min=minimum(system_time_domain);
    time_max=maximum(system_time_domain);

    # Prepare the color code
    color_norm  = system_time_domain ./ time_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

    # Starts the figure
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle/S_{\max}(|\vec{q}|)",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, nothing, nothing),
                   xscale = log10,
                   yscale = log10
                  )

    for df_plot in df_system_evolution
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean_norm[:];

        # Get the color
        color_label=first(DT*unique(df_plot.timeStep))/time_max;

        # Plot
        scatterlines!(ax_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end

        # Colobar to denote the time evolution 
        Colorbar(fig[1, 2], label = L"\tau", colormap = :viridis, limits = (time_min, time_max))

        

        # Plots to add systems labels
        plot!(ax_plot,[100],[100],label=latexstring("\\phi=",100*labels_plot[1],"\\%"),color=:black)
        plot!(ax_plot,[100],[100],label=latexstring("T=",labels_plot[3]),color=:black)
        axislegend(ax_plot,position=:lt)

        # Crate a file name with the labels
        file_name=string("fig_Sq_norm_loglog",join(string.(labels_plot)),".png");

        save(file_name, fig, px_per_unit = 300 / INCH)
end

end




"""
    fig_sq_timeseries_semilog(df_systems::AbstractDataFrame, categories_figures::Vector{Symbol})

Function that creates a figure with the time evolution of a structure factor
"""
function fig_sq_timeseries_semilog(df_systems::GroupedDataFrame{DataFrame}, categories_figures::Vector{Symbol})
# Create the figure for the time evolution for all systems 
for df_system in df_systems

    # Get the time domain
    system_time_domain=sort((DT).*unique(df_system.timeStep));

    # Get the data for the labels
    labels_plot=[df_system[1, col] for col in categories_figures];
        
    # Group the system by time domain
    df_system_evolution=groupby(df_system,:timeStep);

    # Prepare some labels
    time_min=minimum(system_time_domain);
    time_max=maximum(system_time_domain);

    # Prepare the color code
    color_norm  = system_time_domain ./ time_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

    # Starts the figure
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle ",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, nothing, nothing),
                   xscale = log10
                  )

    for df_plot in df_system_evolution
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean[:];

        # Get the color
        color_label=first(DT*unique(df_plot.timeStep))/time_max;

        # Plot
        scatterlines!(ax_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end

        # Colobar to denote the time evolution 
        Colorbar(fig[1, 2], label = L"\tau", colormap = :viridis, limits = (time_min, time_max))

        

        # Plots to add systems labels
        plot!(ax_plot,[100],[100],label=latexstring("\\phi=",100*labels_plot[1],"\\%"),color=:black)
        plot!(ax_plot,[100],[100],label=latexstring("T=",labels_plot[3]),color=:black)
        axislegend(ax_plot,position=:lt)

        # Crate a file name with the labels
        file_name=string("fig_Sq_semilog",join(string.(labels_plot)),".png");

        save(file_name, fig, px_per_unit = 300 / INCH)
end

end




"""
    fig_sq_norm_timeseries_semilog(df_systems::AbstractDataFrame, categories_figures::Vector{Symbol})

Function that creates a figure with the time evolution of a structure factor normalize to a maximum
"""
function fig_sq_norm_timeseries_semilog(df_systems::GroupedDataFrame{DataFrame}, categories_figures::Vector{Symbol})
# Create the figure for the time evolution for all systems 
for df_system in df_systems

    # Get the time domain
    system_time_domain=sort((DT).*unique(df_system.timeStep));

    # Get the data for the labels
    labels_plot=[df_system[1, col] for col in categories_figures];
        
    # Group the system by time domain
    df_system_evolution=groupby(df_system,:timeStep);

    # Prepare some labels
    time_min=minimum(system_time_domain);
    time_max=maximum(system_time_domain);

    # Prepare the color code
    color_norm  = system_time_domain ./ time_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

    # Starts the figure
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle/S_{\max}(|\vec{q}|)",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, -0.5, 2),
                   xscale = log10
                  )

    for df_plot in df_system_evolution
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean_norm[:];

        # Get the color
        color_label=first(DT*unique(df_plot.timeStep))/time_max;

        # Plot
        scatterlines!(ax_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end

        # Colobar to denote the time evolution 
        Colorbar(fig[1, 2], label = L"\tau", colormap = :viridis, limits = (time_min, time_max))

        

        # Plots to add systems labels
        plot!(ax_plot,[100],[100],label=latexstring("\\phi=",100*labels_plot[1],"\\%"),color=:black)
        plot!(ax_plot,[100],[100],label=latexstring("T=",labels_plot[3]),color=:black)
        axislegend(ax_plot,position=:lt)

        # Crate a file name with the labels
        file_name=string("fig_Sq_norm_semilog",join(string.(labels_plot)),".png");

        save(file_name, fig, px_per_unit = 300 / INCH)
end

end

#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";

# Constants
DT=0.001;

# Read the directory 
files=readdir(DIR_SAVE);

# Get only those of the structure factor
files=filter(s -> occursin("structure_factor_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

# Group the Vector{DataFrame} into one DataFrame
df_combo=reduce(vcat,df_files);

# Select the categories that define a system
categories_figures=[:phi,:chi_4,:temp,:damp,:N_heat,:N_isothermal];

# Group by system and experiments
df_systems=groupby(df_combo,categories_figures);

#fig_sq_timeseries_loglog(df_systems,categories_figures)

#fig_sq_norm_timeseries_loglog(df_systems,categories_figures)

fig_sq_timeseries_semilog(df_systems,categories_figures)

fig_sq_norm_timeseries_semilog(df_systems,categories_figures)

