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
files=filter(s -> occursin("pore_length_histogram_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

# Group the Vector{DataFrame} into one DataFrame
df_combo=reduce(vcat,df_files);

# Categories
categories_figures=[:phi,Symbol("CL-Con"),:Temperature,:damp,:N_heat,:N_isot];

# Get the time domain of the analysis
time_domain=unique(df_combo.timeStep);

# Group by the time domain 
df_time_domain=groupby(df_combo,:timeStep);

# Select final configuration
time_step_final=time_domain[end];
df_time_step_final=df_time_domain[end];

# Group by systems
df_time_step_systems_final=groupby(df_time_step_final,categories_figures);

# Select the initial configuration
time_step_initial=time_domain[1];
df_time_step_initial=df_time_domain[1];

# Group by systems
df_time_step_systems_initial=groupby(df_time_step_initial,categories_figures);

# Prepare the domains for the histogram
#max_total=maximum([maximum(df.histLength[:]) for df in df_time_domain]);

# Set the amount of bins
n_bins=45;

# Domain of bins
#bins_domain=range(0,max_total,length=n_bins);

# Get the domain of the packing fraction
phi_domain=unique(df_time_step_initial.phi);

    # prepare some labels
    phi_min=minimum(phi_domain);
    phi_max=maximum(phi_domain);

    # prepare the color code
    color_norm  = phi_domain ./ phi_max;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

    # create the color manually
    color_grad=cgrad(:viridis);

    # create the dictionary
    phi_to_color=Dict(phi_domain .=> trunc.(Int64,range(1,length(color_grad),length=length(phi_domain))))


# Get the data for the labels
labels_plot=[df_time_step_initial[1, col] for col in categories_figures] # Waring: This only works for one case


# Start the figure
    fig = Figure()


    ax_final_plot = Axis(fig[4:6, 1:1],
                   xlabel = L"r",
                   ylabel = L"\mathrm{pdf}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (0, nothing, 0, nothing)
                  )

        for df_system in df_time_step_systems_final
            # Get the mean value
            #mean_pore=first(df_system.mean_pore_set);

            # Get the color
            color_label=phi_to_color[first(unique(df_system.phi))];

            density!(ax_final_plot,df_system.histLength[:],
                     colormap = :viridis,
                     colorrange = (color_min, color_max),
                     #alpha=0.0,
                     color = (:white,0.0),
                     strokecolor = color_grad[color_label],
                     strokewidth = 1.5
                    )


        end 

        # Plots to add systems labels
        plot!(ax_final_plot,[-1],[-1],color=:black,markersize=0.0,label=latexstring("t=",DT*time_step_final,"~\\tau"));
        plot!(ax_final_plot,[-1],[-1],color=:black,markersize=0.0,label=latexstring("T=",labels_plot[3]));

        axislegend(ax_final_plot,
                   position=:rt, 
                   framevisible = true,
                   framewidth = 0.0,
                   padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   patchsize=(0.0f0, 0.0f0)
                  )


# Select final configuration
    ax_initial_plot = Axis(fig[1:3, 1:1],
                   xlabel = L"r",
                   ylabel = L"\mathrm{pdf}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (0, 7, 0, nothing)
                  )

        for df_system in df_time_step_systems_initial
            # Get the mean value
            #mean_pore=first(df_system.mean_pore_set);

            # Get the color
            color_label=phi_to_color[first(unique(df_system.phi))];

            density!(ax_initial_plot,df_system.histLength[:],
                     colormap = :viridis,
                     colorrange = (color_min, color_max),
                     #alpha=0.0,
                     color = (:white,0.0),
                     strokecolor = color_grad[color_label],
                     strokewidth = 1.5
                    )
        end 

    # Linkaxes
        linkyaxes!(ax_initial_plot,ax_final_plot)
        linkxaxes!(ax_initial_plot,ax_final_plot)

        # Plots to add systems labels
        plot!(ax_initial_plot,[-1],[-1],color=:black,markersize=0.0,label=latexstring("t=",DT*time_step_initial,"~\\tau"));
        plot!(ax_initial_plot,[-1],[-1],color=:black,markersize=0.0,label=latexstring("T=",labels_plot[3]));

        axislegend(ax_initial_plot,
                   position=:rt, 
                   framevisible = true,
                   framewidth = 0.0,
                   padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   patchsize=(0.0f0, 0.0f0)
                  )



    # Colobar to denote the time evolution 
    Colorbar(fig[1:6, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

        # Crate a file name with the labels
        #file_name=string("fig_pore_histogram_phiseries",join(string.(labels_plot)),"_time_initial_final.png");

        #save(file_name, fig, px_per_unit = 300 / INCH)


