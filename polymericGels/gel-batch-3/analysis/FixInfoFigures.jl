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
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";

# Constants
DT=0.001;

# Read the directory 
files=readdir(DIR_SAVE);

# Get only those of the structure factor
files=filter(s -> occursin("fix_mean_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_SAVE,file), DataFrame) for file in files];

# Combine the dataframes
df_group=reduce(vcat,df_files)

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

