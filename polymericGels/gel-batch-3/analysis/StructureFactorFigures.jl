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

aspect_ratio = 1.6;
width = 18*CM;
height = width/aspect_ratio;

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
        size = (width, height)                 # tamaño de la figura (ancho, alto)
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

"""
    fig_sq_phiseries_loglog(df_combo::DataFrame, categories_figures::Vector{Symbol})

Creates a figure that compares the final configuration structure factor for different packing fractions
"""
function fig_sq_phiseries_loglog(df_combo::DataFrame, categories_figures::Vector{Symbol})
# Limits of the exponents
exp_min = -2.0; 
exp_max = 1.0; 

# Define length regimes
alpha_colors=0.3;
q_mins_regimes=[10^(exp_min), 10^(-0.5), 10^(0.5)];
q_maxs_regimes=[10^(-0.5), 10^(0.5), 10^(exp_max)];
q_colors_regimes=[(:black,alpha_colors),(:blue,alpha_colors),(:red,alpha_colors)];

# Exponents for the ticks domain
exp_ticks = range(exp_min,exp_max,step=1);

# Create ticks domain
x_ticks_domain = map(s->10^s,exp_ticks);

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
                   yscale = log10,
                   xticks = (x_ticks_domain,map(s->latexstring("10^{",s,"}"),exp_ticks))
                  )

    vspan!(ax_raw_plot,q_mins_regimes,q_maxs_regimes,color=q_colors_regimes)

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
                   yscale = log10,
                   xticks = (x_ticks_domain,map(s->latexstring("10^{",s,"}"),exp_ticks))
                  )

    vspan!(ax_norm_plot,q_mins_regimes,q_maxs_regimes,color=q_colors_regimes)

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
        file_name=string("fig_Sq_phi_series_loglog",join(string.(labels_plot)),"_time_",time_max,".png");

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
        file_name=string("fig_Sq_phi_series_semilog",join(string.(labels_plot)),"_time_",time_max,".png");

        save(file_name, fig, px_per_unit = 300 / INCH)

end


"""
    fig_sq_phiseries_loglog_initial_final(df_combo::DataFrame, categories_figures::Vector{Symbol})

Creates a figure that compares the initial and final configurations structure factor for different packing fractions
"""
function fig_sq_phiseries_semilog_loglog_final(df_combo::DataFrame, categories_figures::Vector{Symbol})

# Limits of the exponents
exp_min = -2.0; 
exp_max = 1.0; 

# Define length regimes
alpha_colors=0.3;
q_mins_regimes=[10^(exp_min), 10^(-0.5), 10^(0.5)];
q_maxs_regimes=[10^(-0.5), 10^(0.5), 10^(exp_max)];
q_colors_regimes=[(:black,alpha_colors),(:blue,alpha_colors),(:red,alpha_colors)];

# Exponents for the ticks domain
exp_ticks = range(exp_min,exp_max,step=1);

# Create ticks domain
x_ticks_domain = map(s->10^s,exp_ticks);

# Get the biggest timestep
time_max=maximum(df_combo.timeStep);

# Get all rows at the max time
df_aux_max=df_combo[df_combo.timeStep .== time_max, :];

# Group by system and experiments
df_systems=groupby(df_aux_max,categories_figures);

# Get the domain of the packing fraction
phi_domain_max=unique(df_aux_max.phi);

# Get the data for the labels
labels_plot=[df_aux_max[1, col] for col in categories_figures] # Waring: This only works for one case

    # Prepare some labels
    phi_min=minimum(phi_domain_max);
    phi_max=maximum(phi_domain_max);

    # Prepare the color code
    color_norm  = phi_domain_max ./ phi_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

    # Starts the figure
    fig = Figure()
    ax_semilog_plot = Axis(fig[1:3, 1:2],
                   xlabel = L"|\vec{q}|",
                   #ylabel = L"\langle S(|\vec{q}|)\rangle",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, 0, 250),
                   xticks = x_ticks_domain,
                   xscale = log10
                  )

    vspan!(ax_semilog_plot,q_mins_regimes,q_maxs_regimes,color=q_colors_regimes)

    for df_plot in df_systems
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean[:];

        # Get the color
        color_label=first(unique(df_plot.phi))/phi_max;

        # Plot
        scatterlines!(ax_semilog_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max))
    end

        plot!(ax_semilog_plot,[100],[100],color=:black,label=latexstring("t=",DT*time_max,"~\\tau"),markersize=0.0);
        plot!(ax_semilog_plot,[100],[100],color=:black,label=latexstring("T=",labels_plot[3]),markersize=0.0);

        axislegend(ax_semilog_plot,
                   position=:lt, 
                   framevisible = true,
                   framewidth = 0.0,
                   padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   patchsize=(0.0f0, 0.0f0)
                  )


    ax_loglog_plot = Axis(fig[1:3, 1:2],
                   #xlabel = L"|\vec{q}|",
                   #ylabel = L"\langle S(|\vec{q}|)\rangle",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.5), 10^1, 10^(-0.5), 10^(3)),
                   xscale = log10,
                   yscale = log10,
                   yaxisposition = :right
                  )
    hidespines!(ax_loglog_plot)
    hidexdecorations!(ax_loglog_plot)

    #vspan!(ax_loglog_plot,q_mins_regimes,q_maxs_regimes,color=q_colors_regimes)

    for df_plot in df_systems
        # Get the domain
        domain=df_plot.q_mean[:];
    
        # Get the range
        range=df_plot.Sq_mean[:];

        # Get the color
        color_label=first(unique(df_plot.phi))/phi_max;

        # Plot
        scatterlines!(ax_loglog_plot,domain,range,
                      color = color_label,
                      colorrange = (color_min, color_max),
                      marker = :utriangle,
                      markersize = 15
                     )
    end

        scatterlines!(ax_loglog_plot,[100],[100],color=:black,label=L"S(|\vec{q}|)");
        scatterlines!(ax_loglog_plot,[100],[100],color=:black,marker = :utriangle,label=L"\log(S(|\vec{q}|))");

        axislegend(ax_loglog_plot,
                   position=:rt, 
                   framevisible = true,
                   #framewidth = 0.0,
                   #padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   #patchsize=(0.0f0, 0.0f0)
                  )


       

        # Colobar to denote the time evolution 
        Colorbar(fig[1:3, 3], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

        # Crate a file name with the labels
        file_name=string("NEW_fig_Sq_phi_series__semilog_loglog_time_final.png");

        save(file_name, fig, px_per_unit = 300 / INCH)

end


#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_DATA = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
DIR_SAVE = joinpath(DIR_MAIN,"figures");

# Read the directory 
files=readdir(DIR_DATA);

# Get only those of the structure factor
files=filter(s -> occursin("structure_factor_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

# Group the Vector{DataFrame} into one DataFrame
df_group=reduce(vcat,df_files);

# NEED TO BE THE SAME AS THE FixInfoAnalysis.jl
# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp]; #,:tstep

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# Categories for time step
category_time = :timeStep;

# Group by time step
data_per_time_step = groupby(df_group,category_time);

# Select one time step
data_time_step = data_per_time_step[1];

for data_time_step in data_per_time_step

    # Get the time step for name and stuff
    time_step = first(data_time_step.timeStep);

    # Group by experiments
    data_per_experiment = groupby(data_time_step,categories_experiment);

    # Select one experiment
    data_experiment = data_per_experiment[1];

    for data_experiment in data_per_experiment

        # Figure same experiment, different systems
        fig_experiment = Figure();

        # For the legends
        legends = [];

        # Create an Axis
        ax_plot = Axis(fig_experiment[1:4,1],
                       xlabel = L"|\vec{q}|~[1/\sigma]",
                       ylabel = L"S(|\vec{q}|)",
                       xminorticksvisible = true,
                       xminorgridvisible = true,
                       limits = (nothing, nothing, nothing, nothing),
                       xscale = log10,
                       yscale = log10
                      )

        # Group by system 
        data_per_system = groupby(data_experiment,categories_system);

        # Select one system
        data_system = data_per_system[1];

        for (it_system,data_system) in enumerate(data_per_system)

            # Get the domain and range
            q_domain = data_system.q_mean;
            Sq_range = data_system.Sq_mean;

            # Add the plot
            lines!(ax_plot,q_domain,Sq_range,linewidth=3)

            # Create the legend
            aux=[];
            append!(aux,[latexstring("\\mathrm{System}~",it_system)])
            for cat in categories_system
                label=latexstring("\\mathrm{",string(cat),"}=~",first(data_system[!,cat]));
                append!(aux,[label])
            end
            append!(legends,[aux])
        end

    # Add legend in the figure
ax_legend = Axis(fig_experiment[1:1,1:1],
                limits = (0, 0, 0, 0)
                )

    # Hide everything
    hidespines!(ax_legend)
    hidedecorations!(ax_legend)

    # Create decoy plot per system to add the labels
    plots_decoy = mapreduce(s->[scatter!(ax_legend,1,1)],vcat,1:length(data_per_system))

    # Add the legend
    Legend(fig_experiment[5:6,1],
           plots_decoy,
           latexstring.(join.(legends,"~")),
           orientation = :vertical
          )

    # Display the figure
    display(fig_experiment)

    # Create the name by getting the experiment id
    aux_name = []
    for cat in categories_experiment
        label=string(cat,"_",first(data_experiment[!,cat]));
        append!(aux_name,[label])
    end

    ax_aux = Axis(fig_experiment[1:4,1],limits=(0,0,0,0))
    hidespines!(ax_aux)
    hidedecorations!(ax_aux)
    scatter!(ax_aux,1,1,markersize=0.1,color=:white,label=latexstring("N_{t}=",time_step))

    # add the time step
        axislegend(ax_aux,
                   position=:rt, 
                   framevisible = true,
                   framewidth = 0.0,
                   padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   patchsize=(0.0f0, 0.0f0)
                  )

    # Join the categories
    aux_name = join(aux_name,"_");

    # Create the file name
    figure_name = string("structure_factor_experiment_",aux_name,"_timestep_",time_step,".png")

    # save the figures
    save(joinpath(DIR_SAVE,figure_name), fig_experiment, px_per_unit = 300 / INCH)

    end
end

#=
# Select the categories that define a system
#categories_figures=[:phi,:chi_4,:temp,:damp,:N_heat,:N_isothermal];

# Group by system and experiments
#df_systems=groupby(df_combo,categories_figures);

#=
    Time series figures
    fig_sq_timeseries_loglog(df_systems,categories_figures)

    fig_sq_norm_timeseries_loglog(df_systems,categories_figures)

    fig_sq_timeseries_semilog(df_systems,categories_figures)

    fig_sq_norm_timeseries_semilog(df_systems,categories_figures)


=#

#=
    Compare different concentrations at the last time step
    fig_sq_phiseries_loglog(df_combo,categories_figures)

    fig_sq_phiseries_semilog(df_combo,categories_figures)

    fig_sq_phiseries_semilog_loglog_final(df_combo, categories_figures)
=#

# Limits of the exponents
exp_x_min = -2.0; 
exp_x_max = 1.0; 
exp_y_min = -1.0; 
exp_y_max = 3.0; 

# Define length regimes
alpha_colors=0.3;
q_mins_regimes=[10^(exp_x_min), 10^(-0.5), 10^(0.5)];
q_maxs_regimes=[10^(-0.5), 10^(0.5), 10^(exp_x_max)];
q_colors_regimes=[(:black,alpha_colors),(:blue,alpha_colors),(:red,alpha_colors)];

# Exponents for the ticks domain
exp_x_ticks = range(exp_x_min,exp_x_max,step=1);
exp_y_ticks = range(exp_y_min,exp_y_max,step=1);


# Create ticks domain
x_ticks_domain = map(s->10^s,exp_x_ticks);
y_ticks_domain = map(s->10^s,exp_y_ticks);


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
                   yscale = log10,
                   xticks = (x_ticks_domain,map(s->latexstring("10^{",s,"}"),exp_x_ticks)),
                   yticks = (y_ticks_domain,map(s->latexstring("10^{",s,"}"),exp_y_ticks))
                  )

    vspan!(ax_raw_plot,q_mins_regimes,q_maxs_regimes,color=q_colors_regimes)

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
      
        # Plots to add systems labels
        scatterlines!(ax_raw_plot,[10^2.0],[1],color=:black,linestyle=:solid,label=latexstring("t=",DT*time_max,"~\\tau"));
        #lines!(ax_plot,[-1],[-1],color=:black,linestyle=:dash,label=latexstring("t=",DT*time_step_final,"~\\tau"));

        axislegend(ax_raw_plot,
                   position=:rt, 
                   framevisible = true,
                   #framewidth = 0.0,
                   #padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   #patchsize=(0.0f0, 0.0f0)
                  )

        # Colobar to denote the time evolution 
        Colorbar(fig[1:3, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

        # Crate a file name with the labels
        file_name=string("fig_Sq_phi_series_loglog",join(string.(labels_plot)),"_time_",time_max,".png");

        #save(file_name, fig, px_per_unit = 300 / INCH)
=#



