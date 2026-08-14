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
    pore_domain=df_files[1].bins;

    # Get the histogram
    hist=mapreduce(s->s.histogram_pore,hcat,df_files);

    # Combine the results
    hist=sum(hist,dims=2)[:];

    # GEt ids
    ids=map(col->df_files[1][!, col],categories_figures)

    # Create a dataframe to plot
    df_to_store=DataFrame([pore_domain, hist],[:poreDomain, :histRadius]);

    for (col, val) in zip(categories_figures, ids)
        df_to_store[!, col] = val 
    end
    
    # Add the time step to the dataframe
    df_to_store[!,:timeStep] .= df_files[1].timeStep

    return df_to_store
end

#=
function collect_df(DIR_SAVE::String, files::Vector{String})
    return [CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files]
end
=#

"""
    get_group_files(pattern::string,files::Vector{String})

Function that return and array of arrays of names of files of the same experiment at different time steps given an array of files names
"""
function get_group_files(patron::Regex,files::Vector{String})

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


"""
    create_histogram(low_bin::Float64,bin_size::Float64,high_bin::Float64,data::Vector{Float64})

Function that retuns the counts of each bins
"""
function create_histogram(low_bin::Float64,bin_size::Float64,high_bin::Float64,data::Vector{Float64})
    
    # Create the bounds of the bins
    low_bins = (low_bin:bin_size:high_bin);
    high_bins = (low_bin+bin_size:bin_size:high_bin+bin_size);


    return map(bin->length(findall(x-> low_bins[bin] .<=  x .<high_bins[bin],  sort(data))), eachindex(high_bins))
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
files=filter(s -> occursin("pore_analysis_spheres_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

# Combine the dataframes
df_group=reduce(vcat,df_files)

# NEED TO BE THE SAME AS THE FixInfoAnalysis.jl
# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

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

    # Group by experiments
    data_per_experiment = groupby(data_time_step,categories_experiment);

    data_experiment = data_per_experiment[1];

        # Figure same experiment, different systems
        fig_experiment = Figure();

        # For the legends
        legends = [];

        # Create an Axis
        ax_plot = Axis(fig_experiment[1:4,1],
                           xlabel = L"\mathrm{Time~}[\tau]",
                           ylabel = L"U(\tau)~[\epsilon]",
                           xminorticksvisible = true,
                           xminorgridvisible = true,
                           limits = (nothing, nothing, nothing, nothing),
                           xscale = log10
                          )

        # Group by system
        data_per_system = groupby(data_experiment,categories_system);

        data_system = data_per_system[1];


#=
#"pore_analysis_spheres_"

# Prepare for reading the information
pattern = r"pore_analysis_spheres_(.+)_step_(\d+)_simulation_\d+\.csv"

# Get the group the simulations by expriment and time domain
groups = get_group_files(pattern,files)

categories_system=[:phi,:chi_4,:temp,:damp,:tstep];    # Select the categories that define a system
categories_experiment=[:N_heat,:N_isothermal];  # Create categories to select different experiments (Just in case)

# Categories
categories_figures=[:phi,:chi_4,:temp,:damp,:N_heat,:N_isothermal,:timeStep,:tstep];

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

# Prepare tha labels
labels=unique(df_combo[!,categories_system]);

# Get the domain of the packing fraction
phi_domain=sort(unique(labels.phi));

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
                   xlabel = L"\mathrm{Radius}",
                   ylabel = L"\mathrm{Counts}/\mathrm{total}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (0, 10, 0, nothing)
                  )

        for df_system in df_time_step_systems_initial

            # Get the color
            color_label=phi_to_color[first(unique(df_system.phi))];

            lines!(ax_plot,df_system.poreDomain,df_system.histRadius[:]./sum(df_system.histRadius[:]),
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

            lines!(ax_plot,df_system.poreDomain,df_system.histRadius[:]./sum(df_system.histRadius[:]),
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
        file_name=string("fig_pore_histogram_phiseries_spheres",join(string.(labels_plot)),"_time_initial_final.png");
=#


#=
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"r",
                   ylabel = L"\mathrm{pdf}",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (0, high_bin, 0, nothing)
                  )

        for df_sets in df_experiments 
            for df_experiment in df_sets

                # Get the color
                color_label=first(unique(df_experiment.phi))/phi_max;

                # Get the time steps available
                time_domain=unique(df_experiment.timeStep);
            
                # Get the initial and final configurations
                time_initial=minimum(time_domain);
                time_final=maximum(time_domain);

                # Group by time step
                df_time=groupby(df_experiment,:timeStep);
         
                # Select the initial and final configuration
                df_time_intial = subset(df_time, :timeStep => (x -> x .== time_initial)) 
                df_time_final = subset(df_time, :timeStep => (x -> x .== time_final))
        
                hist_initial=create_histogram(low_bin,bin_size,high_bin,df_time_intial.histogram_pore)
                hist_final=create_histogram(low_bin,bin_size,high_bin,df_time_final.histogram_pore)

                lines!(ax_plot,bins_domain,hist_initial./sum(hist_initial),
                        linestyle=:solid,
                        linewidth=3,
                        color = color_label,
                        colorrange = (color_min, color_max)
                       )
                lines!(ax_plot,bins_domain,hist_final./sum(hist_final),
                        linestyle=:dash,
                        linewidth=3,
                        color = color_label,
                        colorrange = (color_min, color_max)
                       )
            end
        end 





#hist=create_histogram(low_bin,bin_size,high_bin,df_data[1].histogram_pore)

=#





#=
# Read the files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

# Group the Vector{DataFrame} into one DataFrame
df_combo=reduce(vcat,df_files);

# Get the time domain of the analysis
time_domain=unique(df_combo.timeStep);

# Group by the time domain 
df_time_domain=groupby(df_combo,:timeStep);

# Select one time step
time_step=time_domain[1];
df_time_step=df_time_domain[1];

# Get the domain of the packing fraction
phi_domain=unique(df_time_step.phi);

# Categories
categories_figures=[:phi,:chi_4,:temp,:damp,:N_heat,:N_isothermal];

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
        #save(file_name, fig, px_per_unit = 300 / INCH)
=#

