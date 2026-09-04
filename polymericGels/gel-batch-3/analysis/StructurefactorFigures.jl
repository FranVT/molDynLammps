#=
    Script that creates a figure of the structure factor
=#


using DataFrames, CSV
using GLMakie, LaTeXStrings 
using Statistics, LsqFit

#=
    Set figures 
=#
INCH = 96;
PT = 4/3;
CM = INCH / 2.54;

aspect_ratio = 1.6;
width = 20*CM;
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
color_map = :curl;


#=
    Functions 
=#

"""
    extract_Sq_analysis(DIR_DATA::String)
Get the averages
"""
function extract_Sq_analysis(DIR_DATA::String)

    # Read the directory 
    files=readdir(DIR_DATA);

    # Get only those of the structure factor
    files=filter(s -> occursin("structure_factor_", s), files);

    # Read the files
    df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

    # Create one dataframe
    df_files = reduce(vcat,df_files)
   
    return  df_files   

end



#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_DATA = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
DIR_SAVE = joinpath(DIR_MAIN,"figures");

# Combine the dataframes
df_group=extract_Sq_analysis(DIR_DATA);

# Add the time instant
df_group[!,:time] = df_group.timeStep.*df_group.tstep;


# NEED TO BE THE SAME AS THE FixInfoAnalysis.jl
# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:time_heat,:time_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# Compare the same experiment with different systems


# Group by experiments
data_per_experiment = groupby(df_group,categories_experiment);

# List of linestyles for each system
label_systems = [
                 :solid,
                 :dot,
                 :dash,
                 :dashdot,
                 :dashdotdot
                ];

# Select one experiment
#    data_experiment = data_per_experiment[1];
    
for (it_exp,data_experiment) in enumerate(data_per_experiment)

    # Start the figure
    fig = Figure()
    
    # Get the time domain
    time_domain = sort(unique(data_experiment.time));

    # Create the label for the time domain
    time_max = maximum(time_domain);
    time_min = minimum(time_domain);

    # Rime domain nromalize
    time_norm = time_domain./time_max;

    # get the min and max color
    color_min = time_min/time_max;
    color_max = time_max/time_max;

    # Group by system
    data_per_system = groupby(data_experiment,categories_system);

    # For the legends
    legends = [];

    # Add legend in the figure
    ax_legend = Axis(fig[1:4,1:5],
                    limits = (0, 0, 0, 0)
                    )

    # Hide everything
    hidespines!(ax_legend)
    hidedecorations!(ax_legend)

    # Prepare the ticks
    n_ticks = 10;
    q_aux_ticks = sort(unique(data_experiment.q_mean));
    l_domain = 2*pi./q_aux_ticks;
    ind_range = floor.(Int64,(10).^(range(log(10,1),log(10,length(q_aux_ticks)),length=n_ticks)));
    q_positions = round.(q_aux_ticks[ind_range],digits=2);
    q_ticks = latexstring.(q_positions);
    l_ticks = latexstring.(round.(l_domain[ind_range],digits=2));

    # --- Define tick positions (in q-space) and their top labels (λ = 2π/q) ---
    ax_bottom = Axis(fig[1:4, 1:5],
                         xlabel = L"|\vec{q}|",
                         ylabel = L"\mathrm{Intensity}",
                         xticks = (q_positions, q_ticks),
                             xscale = log10,
                             yscale = log10,
                             xticklabelrotation = pi/4
                            )

    # --- Top axis: wavelength λ ---
    ax_top = Axis(fig[1:4, 1:5],
                          xaxisposition = :top,
                          yaxisposition = :right,

        # Place ticks at the same data coordinates (q values),
        # but display the corresponding λ labels.
                          xticks = (q_positions, l_ticks),
                          xlabel = L"\mathrm{Wavelength}",

        # Spines: show only the top spine
                          topspinevisible = true,
                          bottomspinevisible = false,
                          leftspinevisible = false,
                          rightspinevisible = false,

                          xgridvisible = false,
                          #ygridvisible = false,

        # Hide all y‑axis decorations on the top axis
                          yticks = ([], []),
                          ylabelvisible = false,
                          ygridvisible = false,
                          yticklabelsvisible = false,
                             xscale = log10,
                             yscale = log10,
                             xticklabelrotation = pi/4
                         )

    # Synchronise limits and zoom/pan behaviour
    linkaxes!(ax_bottom, ax_top)

    # Select one system
    #it_system = 1
    #data_system = data_per_system[it_system];
    for (it_system,data_system) in enumerate(data_per_system)

        # Extract the q domain
        q_domain = unique(data_system.q_mean)[1:end-1];

        # Group by time instant 
        data_per_time = groupby(data_system,:time)

        # Select one time instant
        #data_time = data_per_time[1];
        for data_time in data_per_time

            # Get the time and color label
            color_label = first(data_time.time);
            color_label = color_label/time_max;

            # Group by simulation
            data_per_simulation = groupby(data_time,:Nsim);

            # Allocate for the mean
            Sq_mean = zeros(length(q_domain));

            # Compute the mean 
            for aux in data_per_simulation
                Sq_mean[:] += aux.Sq_mean[1:end-1]
            end
            Sq_mean = Sq_mean./length(data_per_simulation);

            # Take out the last value
            Sq_mean = Sq_mean;

            # Add the line to the plot
            lines!(ax_bottom, q_domain, Sq_mean,
                       colormap = color_map,
                       color=color_label,
                       colorrange = (color_min,color_max),
                       linestyle = label_systems[it_system],
                       linewidth = 4 
                      )
        end # time instants

        # Create the legend
        aux=[];
        append!(aux,[latexstring("\\mathrm{System}~",it_system)])
        for cat in categories_system
            label=latexstring("\\mathrm{",string(cat),"}=~",first(data_system[!,cat]));
            append!(aux,[label])
        end # labels 
        append!(legends,[aux])

    end #systems

    

    # Create decoy plot per system to add the labels
    plots_decoy = mapreduce(s->[lines!(ax_legend,1,1,
                                       linestyle = label_systems[s],
                                       color=:black,
                                      )],vcat,1:length(legends))

    # Add the legend
    Legend(fig[5:5,1:6],
            plots_decoy,
            latexstring.(join.(legends,"~")),
            orientation = :vertical
            )

    # Colobar to denote the time evolution 
    Colorbar(fig[1:4, 6:6], label = L"\tau", colormap = color_map, limits = (time_min, time_max))

    # Display the figure
    display(fig)

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
    save(joinpath(DIR_SAVE,figure_name), fig, px_per_unit = 300 / INCH)
 
end # experiments


#=
            # Model 
            model(q,p) = (p[1])./(q.^(p[2])) .+ p[3];

            # Set intial values for the fit
            p_initial = [mean(Sq_mean), 0.0, 0.0];

            # Add limits of the parameters
            p_lower = [0.0, 0.0, 0.0];
            p_upper = [Inf, Inf, Inf];

            # Fit the data
            fit = curve_fit(model, q_domain, Sq_mean, p_initial; lower=p_lower, upper=p_upper);

            # Get the parameters
            p_final = fit.param|>collect;

            # Get the standard error
            se = standard_errors(fit)|>collect;

            # Get the margin of error
            margin_of_error = margin_error(fit)|>collect;
=#



            # Adjust layout spacing for clarity
            #colgap!(fig, 0)
            #rowgap!(fig, 0)


