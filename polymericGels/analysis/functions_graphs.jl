# these are relative to 1 CSS px
inch = 96;
pt = 4/3;
cm = inch / 2.54;

# Condigure the theme for ALL figures
set_theme!(
    backgroundcolor = :white,
    fontsize = 16pt,                     # tamaño de letra base
    Axis = (
        titlesize = 20pt,                 # tamaño del título del eje
        xlabelsize = 14pt,                # tamaño etiqueta eje x
        ylabelsize = 14pt,                 # tamaño etiqueta eje y
        xticklabelsize = 12pt,             # tamaño números eje x
        yticklabelsize = 12pt,              # tamaño números eje y
        xgridstyle = :dash,               # estilo de la cuadrícula
        ygridstyle = :dash,
        spinewidth = 1.5pt,
    ),
    Legend = (
        labelsize = 10pt,                   # tamaño texto leyenda
        framewidth = 1.5pt,
    ),
    Colorbar = (
        labelsize = 14pt,
        ticklabelsize = 12pt,
    ),
    Figure = (
        size = (15cm, 12cm)                 # tamaño de la figura (ancho, alto)
    )
)


####
#   INFORMATION ORGANIZATION
####

function classify_files(SAVE_DIR, id_str, pattern)
    # Filtrar archivos CSV que contengan el id (escapando puntos para seguridad)
    all_files = readdir(SAVE_DIR)
    
    # Get all the files with the id of the system
    matched = filter(f -> occursin(id_str, f) && endswith(f, ".csv"), all_files)
    
    # Separar por prefijo
    files = filter(f -> startswith(f, pattern), matched)
    
    return files 
end

"""
    Get the information of the fix files
"""
function get_fixInfo(df,SAVE_DIR)
    # Get the id of the system
    id_str = string(first(unique(df[:, :id])));

    # Get the list of all files
    fix_list = classify_files(SAVE_DIR, id_str, "fix_avg_");

    # Create the paths to the files
    fix_paths = joinpath.(SAVE_DIR, fix_list);

    # Read the information of the fix files
    df_fix=CSV.read(fix_paths[1],DataFrame);

    # Select the categories needed to do the plots
    categories_plot=[:id,:TimeStep,:c_t,:c_ep,:c_ek];

    # Select the information to plot
    return df_fix[:,categories_plot]
end

"""
    Get the information of the structure factor analysis
"""
function get_sfInfo(df,SAVE_DIR)
    # Get the id of the system
    id_str = string(first(unique(df[:, :id])));

    # Get the list of all files
    sf_list = classify_files(SAVE_DIR, id_str, "structureFactorPBC");

    sf_paths = joinpath.(SAVE_DIR, sf_list);

    # Read the information of the fix files
    df_sf=map(s->CSV.read(s,DataFrame),sf_paths);

    return df_sf
end

###
#   CLUSTER ANALYSIS
###

function get_cluster_paths(df)
    # Get the id of the system
    id_str = string(first(unique(df[:, :id])));

    # Get the file paths
    cluster_list=classify_files(SAVE_DIR, id_str, "clusterAnalysis")

    if length(cluster_list) == 1
        # Get the paths for the files
        return joinpath(SAVE_DIR, cluster_list[1])
    else
        println("Error, no correct dimensions")
    end

end


"""
    Plot potential and kinetic energy for a series of systems.
"""
function figure_fixEnergy(dfs_fix, dat_files)
    # Get the ids to create labels
    ids_plot = map(s -> unique(s.id)[1], dfs_fix);

    # Dictionaries to map id to info.
    dict_phi = Dict(dat_files.id .=> dat_files.phi);
    dict_CL  = Dict(dat_files.id .=> dat_files."CL-Con");
    dict_T   = Dict(dat_files.id .=> dat_files.Temperature);

    # Labels
    phis      = 100 * map(s -> dict_phi[s], ids_plot);
    phi_min   = minimum(phis);
    phi_max   = maximum(phis);

    # Prepare the color code
    sys_domain  = eachindex(dfs_fix);
    N_systems   = length(dfs_fix);
    color_norm  = phis ./ N_systems;
    color_min   = first(color_norm);
    color_max   = last(color_norm);

    # Prepare the time domain
    dt = 0.001;
    time_domains = map(s -> dt * (s.TimeStep), dfs_fix);

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
    for it_sys in sys_domain
        scatterlines!(ax_plot, time_domains[it_sys], dfs_fix[it_sys].c_ep,
                      color = color_norm[it_sys],
                      colorrange = (color_min, color_max)
                     )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

    save(string("fig_PotEnergy_phiseries.png"), fig, px_per_unit = 300 / inch)


    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"\mathrm{Time~}[\tau]",
                   ylabel = L"K(\tau)~[\epsilon]",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (nothing, nothing, nothing, nothing),
                   xscale = log10
                  )

    # Plot the potential energy for each system
    for it_sys in sys_domain
        scatterlines!(ax_plot, time_domains[it_sys], dfs_fix[it_sys].c_ek,
                      color = color_norm[it_sys],
                      colorrange = (color_min, color_max)
                     )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\phi", colormap = :viridis, limits = (phi_min, phi_max))

    save(string("fig_KineticEnergy_phiseries.png"), fig, px_per_unit = 300 / inch)

end

###
# STRUCTURE FACTOR
###

function figureSq_time_evol(system,id)

# Get the time steps of the S(q)
time_aux=map(s->unique(s.timeStep)[1],system);

# Sort the time steps 
time_domain=sort(time_aux);

# Get the index from the sorted order to the original file
ind_sort=map(s->findall(x->x==s,time_aux)[1],time_domain);

# Create the graph of the time evolution of the structure factor

# Time domain label
dt=0.001;
time_domain=(dt).*time_domain;

# Labels
to=first(time_domain);
tf=last(time_domain);

# Prepare the color code
color_ref=maximum(time_domain);
color_norm=time_domain./color_ref;
color_min=minimum(color_norm);
color_max=maximum(color_norm);

    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"|\vec{q}|",
                   ylabel = L"S(|\vec{q}|)",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (10^(-1.0), 10^1, 10^(-1.0), nothing),
                   xscale = log10,
                   yscale = log10
                  )

    # Plot the potential energy for each system
    for (it_sort,it_og) in enumerate(ind_sort)
        x = system[it_og].qmean
        y = system[it_og].Sqmean
        
        #idx = findall(x .!= 0.0)
        #idy = findall(y .!= 0.0)

        #println(idx)

        scatterlines!(ax_plot, x, y,
                      color = color_norm[it_sort],
                      colorrange = (color_min, color_max)
                     )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\tau", colormap = :viridis, limits = (to, tf))

    save(string("fig_Sq_",id,"_timeDomain.png"), fig, px_per_unit = 300 / inch)

    println("One time evolution of structure factor figure saved")

end


