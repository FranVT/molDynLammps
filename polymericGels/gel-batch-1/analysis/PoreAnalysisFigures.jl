#=
    Script to create figures of pore distribution
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

#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzedData");
FILE_DAT = "dat.csv";

# Constants
DT=0.001;

# Read the directory 
files=readdir(DIR_SAVE);

# Get only those of the structure factor
files=filter(s -> occursin("pore_analysis_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

# Group the Vector{DataFrame} into one DataFrame
df_combo=reduce(vcat,df_files);

# Get the time domain of the analysis
time_domain=unique(df_combo.timeStep);

# Group by the time domain 
df_time_domain=groupby(df_combo,:timeStep);

# Select one time step
time_step=time_domain[end];
df_time_step=df_time_domain[end];

# Get the domain of the packing fraction
phi_domain=unique(df_time_step.phi);

# Categories
categories_figures=[:phi,Symbol("CL-Con"),:Temperature,:damp,:N_heat,:N_isot];

# Get the data for the labels
labels_plot=[df_time_step[1, col] for col in categories_figures] # Waring: This only works for one case

    # Prepare some labels
    phi_min=minimum(phi_domain);
    phi_max=maximum(phi_domain);

    # Prepare the color code
    color_norm  = phi_domain ./ phi_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

# Group by systems
df_time_step_systems=groupby(df_time_step,categories_figures);

# Set the amount of bins
n_bins=64;

# Start the figure
    fig = Figure()
    ax_plot = Axis(fig[1:3, 1:1],
                   xlabel = L"r",
                   ylabel = L"\mathrm{pdf}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (0, nothing, 0, nothing)
                  )

        for df_system in df_time_step_systems
            # Get the mean value
            mean_pore=first(df_system.mean_pore_set);

            # Get the color
            color_label=first(unique(df_system.phi))/phi_max;

            # Plot the histogram
            hist!(ax_plot,df_system.histogram_pore[:], bins=n_bins, normalization = :pdf,
                color = color_label,
                colorrange = (color_min, color_max)
                )
        end 

    # Colobar to denote the time evolution 
    Colorbar(fig[2:3, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

        # Plots to add systems labels
        legend_aux_1=plot!(ax_plot,[-1],[-1],color=:black);
        legend_aux_2=plot!(ax_plot,[-1],[-1],color=:black);

        # Create and auxliary variable for the legend
        legend_aux=[legend_aux_1;legend_aux_2];

        # Create the labels array 
        legend_labels=[latexstring("t=",DT*time_step,"~\\tau");latexstring("T=",labels_plot[3])]

        # Create the label
        Legend(fig[1:1,2],legend_aux,legend_labels,
              halign=:left,
              valign=:top
             )

        # Crate a file name with the labels
        file_name=string("fig_pore_histogram_phiseries",join(string.(labels_plot)),"_time_",time_step,".png");

        save(file_name, fig, px_per_unit = 300 / INCH)
