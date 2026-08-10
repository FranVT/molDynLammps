#=
    To help me organize stuff

    Create one figure per system
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


# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
DIR_FIGURES = joinpath(DIR_MAIN,"figures");


# Read the directory 
files=readdir(DIR_SAVE);

# Select the analysis to graph 
files=filter(s -> occursin("fix_mean_", s), files);

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# Read the dat files
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Read the analysis files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

# Group the data into systems
df_systems=groupby(reduce(vcat,df_files),categories_system);

# Select one system
df_system = df_systems[1];

    # Get the categories information
    tuples_system = map(s-> (s , first(unique(df_system[!,s]))) ,categories_system);

    # Create directories. for labels
    cat_to_value=Dict(tuples_system);

    # Group by experiment
    df_experiments = groupby(df_system,categories_experiment);

    # Amount of experiments
    N_exp = length(df_experiments);

    # Create a list with line styles
    lines_styles = [:solid; :dot; :dash; :dashdot];

    # Prepare for the domain
    domains=[[] for _ in 1:N_exp];

    # Get the domain
    for (it,df) in enumerate(df_experiments)
        domains[it] =cat_to_value[:tstep].*df.TimeStep
    end

    # Prepare for the range 
    range=[[] for _ in 1:length(df_experiments)];

    # Get the range 
    for (it,df) in enumerate(df_experiments)
        range[it] = df.c_ep
    end

    # Start the figure
    fig = Figure();

    ax_plot = Axis(fig[1:1, 1:1],
                xlabel = L"t~[\tau]",
                ylabel = L"U(t)~[\epsilon]",
                xminorticksvisible = true,
                xminorgridvisible = true,
                #limits = (0, 10, 0, nothing)
                )

    # Add the plots
    for it in 1:N_exp
        lines!(ax_plot,domains[it],range[it],linestyle=lines_styles[it],color=:black)
    end
