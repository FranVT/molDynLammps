#=
    Script to create the figures of the structure factor analysis
=#

using DataFrames, CSV
using Statistics
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

# Define the color
color_map = :roma;

#=
    Functions 
=#

"""
    figures_per_time()

Create figures comparing different systems at the same time instant
"""
function figures_per_time(df_group, categories_experiment, categories_system)
# Categories for time step
category_time = :time;

# Group by time step
data_per_time_step = groupby(df_group,category_time);

# Select one time step
#data_time_step = data_per_time_step[1];

for data_time_step in data_per_time_step

    # Get the time step for name and stuff
    time_instant = first(data_time_step.time);

    # Group by experiments
    data_per_experiment = groupby(data_time_step,categories_experiment);

    # Select one experiment
    #data_experiment = data_per_experiment[1];

    for data_experiment in data_per_experiment

        # For the upper limit
        degree_limit_2 = ceil(log(10,maximum(data_experiment.Sq_mean)));

        # For the upper limit
        degree_limit_1 = ceil(log(10,minimum(data_experiment.Sq_mean)));

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
                       limits = (10^(-1.25), 10^(0.8), 10^(-0.5), 10^(degree_limit_2)),
                       xscale = log10,
                       yscale = log10
                      )

        # Group by system 
        data_per_system = groupby(data_experiment,categories_system);

        # Select one system
#        it_system = 1;
#        data_system = data_per_system[it_system];

        for (it_system,data_system) in enumerate(data_per_system)

        # Group by simulation
            # Group by simulation
            data_simulations = groupby(data_system,:Nsim);

            # Get the domain and range of all simulations
            q_domain = map(s->collect(data_simulations[s].q_mean),eachindex(data_simulations));
            Sq_range = map(s->collect(data_simulations[s].Sq_mean),eachindex(data_simulations));

            # Compute the average
            q_domain = reduce(vcat,mean(reduce(hcat,q_domain),dims=2));
            Sq_range = reduce(vcat,mean(reduce(hcat,Sq_range),dims=2));

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
    scatter!(ax_aux,1,1,markersize=0.1,color=:white,label=latexstring("t=",time_instant,"\\tau"))

    # add the time step
        axislegend(ax_aux,
                   position=:lb, 
                   framevisible = true,
                   framewidth = 0.0,
                   padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   patchsize=(0.0f0, 0.0f0)
                  )

    # Join the categories
    aux_name = join(aux_name,"_");

    # Create the file name
    figure_name = string("structure_factor_experiment_",aux_name,"_time_",time_instant,".png")

    # save the figures
    save(joinpath(DIR_SAVE,figure_name), fig_experiment, px_per_unit = 300 / INCH)

    end
end


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

# Add the category of time
df_group[!,:time] = df_group.timeStep .* df_group.tstep;

# For the upper limit
degree_limit = ceil(log(10,maximum(df_group.Sq_mean)));

# NEED TO BE THE SAME AS THE FixInfoAnalysis.jl
# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep]; #

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# List of linestyles for each system
label_systems = [
                 :solid,
                 :dot,
                 :dash,
                 :dashdot,
                 :dashdotdot
                ];

# Create figures comparing different systems at the same time instant
#figures_per_time(df_group, categories_experiment, categories_system)

# Figure of the time evolution of the same system

    # Group by experiments
    data_per_experiment = groupby(df_group,categories_experiment);

    # Select one experiment
    data_experiment = data_per_experiment[1];
    for data_experiment in data_per_experiment

        # Extract the time domain of all systems
        time_domain_experiment = unique(data_experiment.time);

        # Create the label for the time domain
        time_max = maximum(time_domain_experiment);
        time_min = minimum(time_domain_experiment);

        # Rime domain nromalize
        time_norm = time_domain_experiment./time_max;

        # get the min and max color
        color_min = time_min/time_max;
        color_max = time_max/time_max;

        # For the upper limit
        degree_limit_2 = ceil(log(10,maximum(data_experiment.Sq_mean)));

        # Figure same experiment, different systems
        fig_experiment = Figure();

        # For the legends
        legends = [];

        # Create an Axis
        ax_plot = Axis(fig_experiment[1:4,1:2],
                       xlabel = L"|\vec{q}|~[1/\sigma]",
                       ylabel = L"S(|\vec{q}|)",
                       xminorticksvisible = true,
                       xminorgridvisible = true,
                       limits = (10^(-1.25), 10^(0.8), 10^(-0.5), 10^(degree_limit_2)),
                       xscale = log10,
                       yscale = log10
                      )

        # Group by systems
        data_by_time = groupby(data_experiment,:time);

        # Select one time instant 
        # data_time = data_by_time[1]
        for data_time in data_by_time

            # Get the time and color label
            color_label = first(data_time.time);
            color_label = color_label/time_max;

            # Group by systems
            data_by_systems = groupby(data_time,categories_system);

            # Select one system
            #it_system = 1;
            #data_system = data_by_systems[it_system];

            for (it_system,data_system) in enumerate(data_by_systems)

                # Group by simulation
                data_simulations = groupby(data_system,:Nsim);

                # Get the domain and range of all simulations
                q_domain = map(s->collect(data_simulations[s].q_mean),eachindex(data_simulations));
                Sq_range = map(s->collect(data_simulations[s].Sq_mean),eachindex(data_simulations));

                # Compute the average
                q_domain = reduce(vcat,mean(reduce(hcat,q_domain),dims=2));
                Sq_range = reduce(vcat,mean(reduce(hcat,Sq_range),dims=2));

                # Add the plot
                lines!(ax_plot,q_domain,Sq_range,
                       linewidth=3,
                       colormap = color_map,
                       color=color_label,
                       colorrange = (color_min,color_max),
                       linestyle = label_systems[it_system]
                      )

                # Create the legend
                aux=[];
                append!(aux,[latexstring("\\mathrm{System:}~")])
                for cat in categories_system
                    label=latexstring("\\mathrm{",string(cat),"}=~",first(data_system[!,cat]));
                    append!(aux,[label])
                end
                append!(legends,[aux])
            end

        end

        # Make unique the stuff
        legends = unique(legends);

        # Add legend in the figure
        ax_legend = Axis(fig_experiment[1:1,1:1],
                    limits = (0, 0, 0, 0)
                    )

        # Hide everything
        hidespines!(ax_legend)
        hidedecorations!(ax_legend)

        # Create decoy plot per system to add the labels
        plots_decoy = mapreduce(s->[lines!(ax_legend,1,1,color=:black,linestyle=label_systems[s])],vcat,1:length(legends))

        # Add the legend
        Legend(fig_experiment[5:6,1:2],
                plots_decoy,
                latexstring.(join.(legends,"~")),
                orientation = :vertical
                ) 


        # Colobar to denote the time evolution 
        Colorbar(fig_experiment[1:4, 3], label = L"\tau", colormap = color_map, limits = (time_min, time_max))

        # Display the figure
        display(fig_experiment)

        # Create the name by getting the experiment id
        aux_name = []
        for cat in categories_experiment
            label=string(cat,"_",first(data_experiment[!,cat]));
            append!(aux_name,[label])
        end

        # Join the categories
        aux_name = join(aux_name,"_");

        # Create the file name
        figure_name = string("structure_factor_experiment_",aux_name,"_time_domain.png")

        # save the figures
        save(joinpath(DIR_SAVE,figure_name), fig_experiment, px_per_unit = 300 / INCH)
    
    end


