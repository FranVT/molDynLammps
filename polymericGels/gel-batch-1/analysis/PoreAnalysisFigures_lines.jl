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
    collect_df(files::Vector{String})

Collect and combine the total histogram of N experiments
"""
function collect_df(files::Vector{String}, categories_figures::Vector{Symbol})

    # Read the files of one experiment
    df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

    # Get the poreDomain
    pore_domain=df_files[1].poreDomain;

    # Get the histogram
    hist=mapreduce(s->s.histLength,hcat,df_files);

    # Combine the results
    hist=sum(hist,dims=2)[:];

    # GEt ids
    ids=map(col->df_files[1][!, col],categories_figures)

    # Create a dataframe to plot
    df_to_store=DataFrame([pore_domain, hist],[:poreDomain, :histLength]);

    for (col, val) in zip(categories_figures, ids)
        df_to_store[!, col] = val 
    end
    
    # Add the time step to the dataframe
    df_to_store[!,:timeStep] .= df_files[1].timeStep

    return df_to_store
end

"""
    get_group_files(pattern::string,files::Vector{String})

Function that return and array of arrays of names of files of the same experiment at different time steps given an array of files names
"""
function get_group_files(pattern::string,files::Vector{String})

    # Read by experiments
    grupos = Dict{Tuple{String, Int}, Vector{String}}()
    for f in files 
        m = match(patron, f)
        if m !== nothing
            sistema = m.captures[1]   # cadena con los parámetros
            step = parse(Int, m.captures[2])
            clave = (sistema, step)
            # Agregar al grupo correspondiente
            if haskey(grupos, clave)
                push!(grupos[clave], f)
            else
                grupos[clave] = [f]
            end
        else
            @warn "Nombre no coincide con el patrón: $f"
        end
    end

    return grupos
end

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
files=filter(s -> occursin("pore_analysis_lines_NEW_", s), files);

# Prepare for reading the information
pattern = r"pore_analysis_lines_NEW_(.+)_step_(\d+)_simulation_\d+\.csv"

# Get the group the simulations by expriment and time domain
groups = get_group_files(pattern,files)

# Categories
categories_figures=[:phi,Symbol("CL-Con"),:Temperature,:damp,:N_heat,:N_isot];

# Collect the files names by experiment 
files = collect(values(groups));

# Collect and combine the analysis of N simulations
df_data=map(s->collect_df(s,categories_figures),files)

# Group the Vector{DataFrame} into one DataFrame
df_combo=reduce(vcat,df_data);

# Get the time domain of the analysis
time_domain=unique(df_combo.timeStep);

# Group by the time domain 
df_time_domain=groupby(df_combo,:timeStep);

# Select final configuration
    time_step_final=time_domain[1];
    df_time_step_final=df_time_domain[1];

    # Group by systems
    df_time_step_systems_final=groupby(df_time_step_final,categories_figures);

# Select the initial configuration
    time_step_initial=time_domain[end];
    df_time_step_initial=df_time_domain[end];

    # Group by systems
    df_time_step_systems_initial=groupby(df_time_step_initial,categories_figures);

# Get the domain of the packing fraction
phi_domain=sort(unique(df_time_step_initial.phi));

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

# Select final configuration
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Length~of~line}",
                   ylabel = L"\mathrm{Counts}/\mathrm{total}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (0, 55, 0, 0.01)
                  )

        for df_system in df_time_step_systems_initial

            # Get the color
            color_label=phi_to_color[first(unique(df_system.phi))];

            lines!(ax_plot,df_system.poreDomain,df_system.histLength[:]./sum(df_system.histLength[:]),
                      linewidth = 3,
                      linestyle = :solid,
                      colormap = :viridis,
                      colorrange = (color_min, color_max),
                      color = color_grad[color_label],
                      #normalization = :pdf
                     )
 
        end 

        for df_system in df_time_step_systems_final

            # Get the color
            color_label=phi_to_color[first(unique(df_system.phi))];

            lines!(ax_plot,df_system.poreDomain,df_system.histLength[:]./sum(df_system.histLength[:]),
                      linewidth = 3,
                      linestyle = :dash,
                      colormap = :viridis,
                      colorrange = (color_min, color_max),
                      color = color_grad[color_label],
                      #normalization = :pdf
                     )            
        end 

        # Plots to add systems labels
        lines!(ax_plot,[-1],[-1],color=:black,linestyle=:solid,label=latexstring("t=",DT*time_step_initial,"~\\tau"));
        lines!(ax_plot,[-1],[-1],color=:black,linestyle=:dash,label=latexstring("t=",DT*time_step_final,"~\\tau"));

        axislegend(ax_plot,
                   position=:rt, 
                   framevisible = true,
                   #framewidth = 0.0,
                   #padding=(0.0f0,4.0f0,1.0f0,1.0f0),
                   #patchsize=(0.0f0, 0.0f0)
                  )



    # Colobar to denote the time evolution 
    Colorbar(fig[1:1, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

        # Crate a file name with the labels
        file_name=string("fig_pore_histogram_phiseries_line_length",join(string.(labels_plot)),"_time_initial_final.png");

#        save(file_name, fig, px_per_unit = 300 / INCH)
