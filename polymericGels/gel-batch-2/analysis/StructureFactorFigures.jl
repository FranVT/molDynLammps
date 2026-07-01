#=
    Script for the figures of the new strufure factor analysis
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
        file_name=string("NEW_fig_Sq_loglog",join(string.(labels_plot)),".png");

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
        file_name=string("NEW_fig_Sq_norm_loglog",join(string.(labels_plot)),".png");

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
        file_name=string("NEW_fig_Sq_semilog",join(string.(labels_plot)),".png");

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
        file_name=string("NEW_fig_Sq_norm_semilog",join(string.(labels_plot)),".png");

        save(file_name, fig, px_per_unit = 300 / INCH)
end

end

"""
    fig_sq_phiseries_loglog(df_combo::DataFrame, categories_figures::Vector{Symbol})

Creates a figure that compares the final configuration structure factor for different packing fractions
"""
function fig_sq_phiseries_loglog(df_combo::DataFrame, categories_figures::Vector{Symbol})
# Get the biggest timestep
time_max= maximum(df_combo.timeStep);

# Get all rows at the max time
df_aux=df_combo[df_combo.timeStep .== time_max, :];

# Get the domain of the packing fraction
phi_domain=unique(df_aux.phi);

# Get the data for the labels
labels_plot=[df_aux[1, col] for col in categories_figures] # Waring: This only works for one case

    # Prepare some labels
    phi_min=minimum(phi_domain);
    phi_max=maximum(phi_domain);

    # Prepare the color code
    color_norm  = phi_domain ./ phi_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

# Group by system and experiments
df_systems=groupby(df_aux,categories_figures);

    # Starts the figure
    fig = Figure()
    ax_raw_plot = Axis(fig[1:3, 1:1],
                   #xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, 10^(-1.5), nothing),
                   xscale = log10,
                   yscale = log10
                  )

    for df_plot in df_systems
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean[:];

        # Get the color
        color_label=first(unique(df_plot.phi))/phi_max;

        # Plot
        scatterlines!(ax_raw_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end

    ax_norm_plot = Axis(fig[4:6, 1:1],
                   xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle/S_{\max}(|\vec{q}|)",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, nothing, nothing),
                   xscale = log10,
                   yscale = log10
                  )

    for df_plot in df_systems
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean_norm[:];

        # Get the color
        color_label=first(unique(df_plot.phi))/phi_max;

        # Plot
        scatterlines!(ax_norm_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end
        
        # Plots to add systems labels
        legend_aux_1=plot!(ax_raw_plot,[100],[100],color=:black);
        legend_aux_2=plot!(ax_raw_plot,[100],[100],color=:black);

        # Create and auxliary variable for the legend
        legend_aux=[legend_aux_1;legend_aux_2];

        # Create the labels array 
        legend_labels=[latexstring("t=",DT*time_max,"~\\tau");latexstring("T=",labels_plot[3])]

        # Create the label
        Legend(fig[1:1,2],legend_aux,legend_labels,
              halign=:left,
              valign=:top
             )

        # Colobar to denote the time evolution 
        Colorbar(fig[2:6, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

        # Crate a file name with the labels
        file_name=string("NEW_fig_Sq_phi_series_loglog",join(string.(labels_plot)),"_time_",time_max,".png");

        save(file_name, fig, px_per_unit = 300 / INCH)

end


"""
    fig_sq_phiseries_semilog(df_combo::DataFrame, categories_figures::Vector{Symbol})

Creates a figure that compares the final configuration structure factor for different packing fractions
"""
function fig_sq_phiseries_semilog(df_combo::DataFrame, categories_figures::Vector{Symbol})
# Get the biggest timestep
time_max= maximum(df_combo.timeStep);

# Get all rows at the max time
df_aux=df_combo[df_combo.timeStep .== time_max, :];

# Get the domain of the packing fraction
phi_domain=unique(df_aux.phi);

# Get the data for the labels
labels_plot=[df_aux[1, col] for col in categories_figures] # Waring: This only works for one case

    # Prepare some labels
    phi_min=minimum(phi_domain);
    phi_max=maximum(phi_domain);

    # Prepare the color code
    color_norm  = phi_domain ./ phi_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

# Group by system and experiments
df_systems=groupby(df_aux,categories_figures);

    # Starts the figure
    fig = Figure()
    ax_raw_plot = Axis(fig[1:3, 1:1],
                   #xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, nothing, nothing),
                   xscale = log10
                  )

    for df_plot in df_systems
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean[:];

        # Get the color
        color_label=first(unique(df_plot.phi))/phi_max;

        # Plot
        scatterlines!(ax_raw_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end

    ax_norm_plot = Axis(fig[4:6, 1:1],
                   xlabel = L"|\vec{q}|",
                   ylabel = L"\langle S(|\vec{q}|)\rangle/S_{\max}(|\vec{q}|)",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, nothing, nothing),
                   xscale = log10
                  )

    for df_plot in df_systems
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean_norm[:];

        # Get the color
        color_label=first(unique(df_plot.phi))/phi_max;

        # Plot
        scatterlines!(ax_norm_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end
        
        # Plots to add systems labels
        legend_aux_1=plot!(ax_raw_plot,[100],[100],color=:black);
        legend_aux_2=plot!(ax_raw_plot,[100],[100],color=:black);

        # Create and auxliary variable for the legend
        legend_aux=[legend_aux_1;legend_aux_2];

        # Create the labels array 
        legend_labels=[latexstring("t=",DT*time_max,"~\\tau");latexstring("T=",labels_plot[3])]

        # Create the label
        Legend(fig[1:1,2],legend_aux,legend_labels,
              halign=:left,
              valign=:top
             )

        # Colobar to denote the time evolution 
        Colorbar(fig[2:6, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

        # Crate a file name with the labels
        file_name=string("NEW_fig_Sq_phi_series_semilog",join(string.(labels_plot)),"_time_",time_max,".png");

        save(file_name, fig, px_per_unit = 300 / INCH)

end



#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzedData");

# Constants
DT=0.001;

# Read the directory 
files=readdir(DIR_SAVE);

# Get only those of the structure factor
files=filter(s -> occursin("NEW_structure_factor_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

# Group the Vector{DataFrame} into one DataFrame
df_combo=reduce(vcat,df_files);

# Select the categories that define a system
categories_figures=[:phi,Symbol("CL-Con"),:Temperature,:damp,:N_heat,:N_isot];

# Group by system and experiments
df_systems=groupby(df_combo,categories_figures);

#=
    Time series figures

=#
    fig_sq_timeseries_loglog(df_systems,categories_figures)

    fig_sq_norm_timeseries_loglog(df_systems,categories_figures)

    fig_sq_timeseries_semilog(df_systems,categories_figures)

    fig_sq_norm_timeseries_semilog(df_systems,categories_figures)


#=
    Compare different concentrations at the last time step

=#
    fig_sq_phiseries_loglog(df_combo,categories_figures)

    fig_sq_phiseries_semilog(df_combo,categories_figures)


