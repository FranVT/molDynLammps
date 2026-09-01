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

# Group by experiments
data_per_experiment = groupby(df_group,categories_experiment);

# Select one experiment
    data_experiment = data_per_experiment[1];
    
    # Group by system
    data_per_system = groupby(data_experiment,categories_system);

    # Select one system
    data_system = data_per_system[1];



        # Group by time instant 
        data_per_time = groupby(data_system,:time)

        # Select one time instant
        data_time = data_per_time[1];

            # Group by simulation
            data_per_simulation = groupby(data_time,:Nsim);

            # Extract the q domain
            q_domain = data_per_simulation[1].q_mean;

            # Allocate for the mean
            Sq_mean = zeros(length(q_domain));

            # Compute the mean 
            for aux in data_per_simulation
                Sq_mean[:] += aux.Sq_mean
            end
            Sq_mean = Sq_mean./length(data_per_simulation);

            # Take out the last value
            q_domain = q_domain[1:end-1];
            Sq_mean = Sq_mean[1:end-1];

            # Transform to length
            l_domain = 2*pi./q_domain;

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

            fig = Figure()
            
            # --- Define tick positions (in q-space) and their top labels (λ = 2π/q) ---
            ax_bottom = Axis(fig[1, 1],
                             xlabel = L"|\vec{q}|",
                             ylabel = L"\mathrm{Intensity}",
                             xticks = (q_domain[1:10:end], string.(q_domain[1:10:end])),
                            )
            lines!(ax_bottom, q_domain, Sq_mean)

            # --- Top axis: wavelength λ ---
            ax_top = Axis(fig[1, 1],
                          xaxisposition = :top,
                          yaxisposition = :right,

            # Place ticks at the same data coordinates (q values),
            # but display the corresponding λ labels.
                          xticks = (q_domain[1:10:end], string.(l_domain[1:10:end])),
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
                          yticklabelsvisible = false
                         )

            # Synchronise limits and zoom/pan behaviour
            linkaxes!(ax_bottom, ax_top)

            # Adjust layout spacing for clarity
            #colgap!(fig, 0)
            #rowgap!(fig, 0)


