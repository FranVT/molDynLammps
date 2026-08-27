#=
    Script to create figures from the fix files info
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

"""
function eval_model(t,params)
    p = first(unique(params));

    return (p[1])./(t.^(p[2])) .+ p[3]
end

"""
    figure_potential_energy(data_experiment,categories_system,categories_experiment)
"""
function figure_potential_energy(data_experiment,categories_system,categories_experiment,DIR_SAVE)

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

    # For label
#    plot_fig=[Any for _ in 1:length(data_per_system)];

    for (it_system,data_system) in enumerate(data_per_system)

        # Get time step
        dt = first(data_system.tstep);

        # Start of the isothermal process
        start_isothermal = dt*first(data_system.N_heat);

        # Prepare the data for the figure
        time_domain = dt.*collect(data_system.TimeStep);

        # Extract the potential energy
        energy_range = collect(data_system.c_ep);

        # Create a mask to consider only isothermal process
        mask_isothermal = time_domain .> start_isothermal;

        # Apply the mask to domain and range
        #time_domain = time_domain[mask_isothermal];
        #energy_range = energy_range[mask_isothermal];

        # Add the energy to the figure
        lines!(ax_plot,time_domain,energy_range,
               #label=latexstring("\\mathrm{System}~",it_system),
               linewidth = 3
              )

        # Add the theoretical lines
        lines!(ax_plot,time_domain[mask_isothermal],eval_model(time_domain[mask_isothermal],data_system.fit_ep_parameters),
               color = :black,
               linestyle = :dash
              )

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

    # Join the categories
    aux_name = join(aux_name,"_");

    # Create the file name
    figure_name = string("energy_experiment_",aux_name,".png")

    save(joinpath(DIR_SAVE,figure_name), fig_experiment, px_per_unit = 300 / INCH)
end

"""
    parse_vector(str::String)

Functino to parse a string like "[1.0, 2.0]" to Vector{Float64}
"""
function parse_vector(str::String)
    # Quitar corchetes y split por comas
    cleaned = replace(str, r"\[|\]" => "")  # elimina [ y ]
    if isempty(cleaned)
        return Float64[]
    else
        return parse.(Float64, split(cleaned, ","))
    end
end

"""
    extract_fix_avg(DIR_DATA::String)
Get the averages
"""
function extract_fix_avg(DIR_DATA::String)

    # Read the directory 
    files=readdir(DIR_DATA);

    # Get only those of the structure factor
    files=filter(s -> occursin("fix_mean_", s), files);

    # Read the files
    df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

    # Create one dataframe
    df_files = reduce(vcat,df_files)

    # PArse the vectors of strings into vectors of floats
    df_files[!, :fit_ep_parameters] = parse_vector.(df_files[!, :fit_ep_parameters])
    df_files[!, :fit_ep_margin_error] = parse_vector.(df_files[!, :fit_ep_margin_error])
    df_files[!, :fit_ep_standard_error] = parse_vector.(df_files[!, :fit_ep_standard_error])
    
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
df_group=extract_fix_avg(DIR_DATA);

# NEED TO BE THE SAME AS THE FixInfoAnalysis.jl
# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# Group by experiments
data_per_experiment = groupby(df_group,categories_experiment);

for data_experiment in data_per_experiment
    figure_potential_energy(data_experiment,categories_system,categories_experiment,DIR_SAVE)
end



#=

# Prepare for a separation by system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Separate by systems
df_systems=groupby(df_group,categories_system);

# Prepare labels for the systems

    # Packing fraction domain
    phi_domain=unique(df_group.phi);

    # Time step domain
    tstep_domain=unique(df_group.tstep);

# Create labels

    # Packing fraction 
    phi_min   = 100*minimum(phi_domain);
    phi_max   = 100*maximum(phi_domain);
    color_phi_norm  = phi_domain ./ length(phi_domain);
    color_phi_min   = first(color_phi_norm);
    color_phi_max   = last(color_phi_norm);
    color_phi_label = Dict(phi_domain .=> color_phi_norm)

    # Time step
    line_styles=[:solid, :dash];
    tstep_label=Dict(tstep_domain .=> line_styles)

# Prepare the labels
    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Time~}[\tau]",
                   ylabel = L"K(\tau)~[\epsilon]",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (nothing, nothing, nothing, nothing),
                   xscale = log10
                  )

    # Plot the kinetick energy for each system
    for df_system in df_systems 
        # Prepate time domain
        time_domain=(df_system.tstep).*df_system.TimeStep

        aux_phi = first(df_system.phi);
        aux_tstep = first(df_system.tstep);

        lines!(ax_plot, time_domain, df_system.c_ek,
                      color = color_phi_label[aux_phi],
                      colorrange = (color_phi_min, color_phi_max),
                      linestyle = tstep_label[aux_tstep],
                      linewidth = 3 
                     )
    end

# System legends

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi\%", colormap = :viridis, limits = (phi_min, phi_max))

    # Time step lebel 
     for df_system in df_systems
        aux_tstep = first(df_system.tstep);

        lines!(ax_plot,[0],[0],
                color = :black,
                linestyle = tstep_label[aux_tstep],
                linewidth = 3,
                label=latexstring("\\Delta t=",aux_tstep,"~\\tau")
              );
     end
    
axislegend(ax_plot,
                   merge = true,
                   position=:rt, 
                   framevisible = false,
                  )

    save(string("fig_KineticEnergy_phiseries.png"), fig, px_per_unit = 300 / INCH)

    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Time~}[\tau]",
                   ylabel = L"U(\tau)~[\epsilon]",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (nothing, nothing, nothing, nothing),
                   xscale = log10
                  )

    # Plot the potential energy for each system
    for df_system in df_systems 
        # Prepate time domain
        time_domain=(df_system.tstep).*df_system.TimeStep

        aux_phi = first(df_system.phi);
        aux_tstep = first(df_system.tstep);

        lines!(ax_plot, time_domain, df_system.c_ep,
                      color = color_phi_label[aux_phi],
                      colorrange = (color_phi_min, color_phi_max),
                      linestyle = tstep_label[aux_tstep],
                      linewidth = 3 
                     )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi\%", colormap = :viridis, limits = (phi_min, phi_max))

    # Time step lebel 
     for df_system in df_systems
        aux_tstep = first(df_system.tstep);

        lines!(ax_plot,[0],[0],
                color = :black,
                linestyle = tstep_label[aux_tstep],
                linewidth = 3,
                label=latexstring("\\Delta t=",aux_tstep,"~\\tau")
              );
     end
    
        axislegend(ax_plot,
                   merge = true,
                   position=:rt, 
                   framevisible = false,
                  )

    save(string("fig_PotentialEnergy_phiseries.png"), fig, px_per_unit = 300 / INCH)

    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Time~}[\tau]",
                   ylabel = L"T(\tau)",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (nothing, nothing, nothing, nothing),
                   xscale = log10
                  )

    # Plot the temperature for each system
    for df_system in df_systems 
        # Prepate time domain
        time_domain=(df_system.tstep).*df_system.TimeStep

        aux_phi = first(df_system.phi);
        aux_tstep = first(df_system.tstep);

        lines!(ax_plot, time_domain, df_system.c_t,
                      color = color_phi_label[aux_phi],
                      colorrange = (color_phi_min, color_phi_max),
                      linestyle = tstep_label[aux_tstep],
                      linewidth = 3 
                     )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi\%", colormap = :viridis, limits = (phi_min, phi_max))

    # Time step lebel 
     for df_system in df_systems
        aux_tstep = first(df_system.tstep);

        lines!(ax_plot,[0],[0],
                color = :black,
                linestyle = tstep_label[aux_tstep],
                linewidth = 3,
                label=latexstring("\\Delta t=",aux_tstep,"~\\tau")
              );
     end
    
        axislegend(ax_plot,
                   merge = true,
                   position=:rt, 
                   framevisible = false,
                  )

    save(string("fig_Temperature_phiseries.png"), fig, px_per_unit = 300 / INCH)
=#
