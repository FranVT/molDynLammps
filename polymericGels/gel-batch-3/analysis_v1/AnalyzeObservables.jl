#=
    Script to analyze observables and stuff
=# 

using DataFrames, CSV
using Statistics, LsqFit
using GLMakie, LaTeXStrings

#=
 Functions
=#
"""
    directories_to_analyze(DIR_MAIN::String, FILE_DAT::String)

Return a String vector with the directories to compute the observables
"""
function directories_to_analyze(DIR_MAIN::String, FILE_DAT::String)
    # Read the dat file
    df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

    # Get the directories
    simulations_dir = df_dat.dir;

    # Get the ids
    simulations_ids = df_dat.id;

    # Create a directroy to store analyzed data
    # Do nothing id already existed
    mkpath.(joinpath.(simulations_dir,"observables"));

    # Get the time steps
    time_steps = df_dat.tstep;

    # Test to see if the time step is smaller or equal than 0.001
    test_time_steps = time_steps.<=0.001;

    # Compute the total amount of time steps for each simulation
    N_total = df_dat.N_heat .+ df_dat.N_isothermal;

    # Get the interval at which each time step a dump is stored
    N_save_dump = df_dat.N_dump;

    # Create the range of time steps for each space diagram
    time_steps_range = [range(0,N_total[s],step=N_save_dump[s]) for s in eachindex(N_total)];

    # Get the amount of time steps theoretically stored
    N_time_steps_to_store = length.(time_steps_range);

    # Get the amount of actually stored in the directory
    N_time_steps_stored = count.(isfile,readdir.(joinpath.(simulations_dir,"traj");join=true));

    # Test if the simulation is finished
    test_finish_simulation = N_time_steps_to_store.==N_time_steps_stored

    # Combine the tests
    test_filter = test_time_steps .& test_finish_simulation

    # Filter the directories of the simulations that al ready finish 
    simulations_dir = simulations_dir[test_filter];

    return String.(simulations_dir), time_steps_range[test_filter], df_dat[test_filter,:]
end

"""
Recibe un vector de paths y devuelve un vector con los **nombres** de
los archivos que sí existen en el sistema.

# Argumentos
- `paths` : vector de strings con rutas a archivos.

# Retorna
- `Vector{String}` con los nombres (sin directorio) de los archivos existentes.
"""
function files_exist(paths::Vector{<:AbstractString})
    return [p for p in paths if isfile(p)]
end

"""
Lee un archivo con bloques "Timestep N" / "q Sq" / datos
y devuelve un DataFrame con columnas: timestep, q, Sq
"""
function read_Sqtimesteps(ruta::String)
    dfs = DataFrame[]           # acumulador de bloques
    timestep = -1
    q_vals  = Float64[]
    Sq_vals = Float64[]

    function cerrar_bloque!()
        if timestep >= 0 && !isempty(q_vals)
            push!(dfs, DataFrame(timestep = fill(timestep, length(q_vals)),
                                 q        = q_vals,
                                 Sq       = Sq_vals))
        end
        empty!(q_vals)
        empty!(Sq_vals)
    end

    for linea in eachline(ruta)
        linea = strip(linea)

        if isempty(linea)
            continue
        elseif startswith(linea, "Timestep")
            cerrar_bloque!()                        # cierra el bloque anterior
            timestep = parse(Int, split(linea)[2])
        elseif startswith(linea, "q")               # encabezado "q Sq"
            continue
        else
            partes = split(linea)                   # ["1.24562", "3.7684764"]
            push!(q_vals,  parse(Float64, partes[1]))
            push!(Sq_vals, parse(Float64, partes[2]))
        end
    end
    cerrar_bloque!()                                # último bloque

    return vcat(dfs...)                             # un solo DataFrame
end

"""
    compute_Sq_mean_time_series(paths::Vector{String})

Return time domain, mean of the the wave vector and mean of the Sq given N simulations paths
"""
function compute_Sq_mean_time_series(paths::Vector{String})
    # Amount of simulations per system in one experiment
    N_sim = length(paths);

    # This only works if all files have the same timesteps stored.
    
    # Extract all dataframes
    df_Sq_all = read_Sqtimesteps.(paths);

    # Prepare to get time domain
    aux = [unique(df.timestep) for df in df_Sq_all];

    # Time domain
    time_domain = unique(reduce(vcat,aux));

    q_t_domain = [];
    Sq_t_domain = [];

    for s in eachindex(time_domain)

        q_mean = [];
        Sq_mean = [];

        for it_sim in 1:N_sim
            # Select one simulation
            df_aux = df_Sq_all[it_sim];
    
            # Create the time mask
            time_mask = df_aux.timestep .== time_domain[s];

            # Extract the data
            q_sim = df_aux.q[time_mask];
            Sq_sim = df_aux.Sq[time_mask];

            # Prepare for the mean 
            append!(q_mean,[q_sim])
            append!(Sq_mean,[Sq_sim])
        end

        # Compute the mean
        q_mean = reduce(+,q_mean)/N_sim;
        Sq_mean = reduce(+,Sq_mean)/N_sim;

        append!(q_t_domain,[q_mean]);
        append!(Sq_t_domain,[Sq_mean]);
    end

    return time_domain, q_t_domain, Sq_t_domain
end

"""
    moving_mean(y::Vector{Float64}; w::Int=3)

Return a "smother" set of points
"""
function moving_mean(y::Vector{Float64}; w::Int=5)
    n = length(y)
    @assert w ≥ 1 "w debe ser ≥ 1"
    @assert isodd(w) "w debe ser impar"
    @assert w ≤ n   "w no puede exceder la longitud de y"

    y_s = similar(y, Float64)
    r = w ÷ 2  # radio

    @inbounds for i in 1:n
        lo = max(1, i - r)
        hi = min(n, i + r)
        s = 0.0
        for j in lo:hi
            s += y[j]
        end
        y_s[i] = s / (hi - lo + 1)
    end

    return y_s
end

"""
    derivate(x::AbstractVector, y::AbstractVector)

Por DeepSeek
"""
function derivate(x::AbstractVector, y::AbstractVector)
    n = length(x)
    @assert n == length(y) "x e y deben tener la misma longitud"
    @assert n ≥ 3 "Se necesitan al menos 3 puntos"

    dy = zeros(eltype(y), n)

    # --- Primer punto: hacia adelante (3 puntos) ---
    h1 = x[2] - x[1]
    h2 = x[3] - x[2]
    dy[1] = y[1] * (-(2h1 + h2) / (h1 * (h1 + h2))) +
            y[2] * ((h1 + h2) / (h1 * h2)) +
            y[3] * (-h1 / (h2 * (h1 + h2)))

    # --- Puntos interiores: centrada (3 puntos) ---
    @inbounds for i in 2:n-1
        ha = x[i]   - x[i-1]   # espaciado hacia atrás
        hb = x[i+1] - x[i]     # espaciado hacia adelante
        dy[i] = y[i-1] * (-hb / (ha * (ha + hb))) +
                y[i]   * ((hb - ha) / (ha * hb)) +
                y[i+1] * ( ha / (hb * (ha + hb)))
    end

    # --- Último punto: hacia atrás (3 puntos) ---
    h1 = x[n-1] - x[n-2]
    h2 = x[n]   - x[n-1]
    dy[n] = y[n-2] * ( h2 / (h1 * (h1 + h2))) +
            y[n-1] * (-(h1 + h2) / (h1 * h2)) +
            y[n]   * ((h1 + 2h2) / (h2 * (h1 + h2)))

    return dy
end

"""
    norm_robust()

Normalize a set of point to identify easily the maximums and minimus
"""
function norm_robust(y)
    m = median(y)
    s = median(abs.(y .- m))  # MAD
    s == 0 && (s = std(y))
    return (y .- m) ./ s, m, s
end


function eval_model_log(t,p)
    return p[1].*t.+p[2]
end



#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_DATA = "/run/media/franvt/rogelio/DinMol/gel-batch-3-long/data/";
FILE_DAT = "dat.csv";
FILE_FIX = "system_assembly.fixf";

# Get the directories to analyze
simulations_dir, time_steps_range, df_dat = directories_to_analyze(DIR_MAIN,FILE_DAT);

# Organize the dat file into experiments and that stuff

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments
categories_experiment=[:time_heat,:time_isothermal];

# Group by experiment 
df_dat_experiments = groupby(df_dat,categories_experiment);

# Select one experiment
df_dat_experiment = df_dat_experiments[2];

    # Group by system
    df_dat_systems = groupby(df_dat_experiment,categories_system);

    # Select one system
    #df_dat_system = df_dat_systems[1];

    # To store the data from each system at the time series
    q_systems = [];
    Sq_systems = [];



    for df_dat_system in df_dat_systems

        # Get the paths to the files 
        paths = files_exist(joinpath.(df_dat_system.dir,"structure_factor.txt"));

        # Compute the mean of a time series of the structure factor
        time_domain, q_t_series, Sq_t_series = compute_Sq_mean_time_series(paths);

        # Create the fit for fractal dimension
        it_time = length(time_domain);

            # Get the data
            q_mean = q_t_series[it_time];
            Sq_mean = Sq_t_series[it_time];

            append!(q_systems,[q_mean]);
            append!(Sq_systems,[Sq_mean]);

    end
            
fig = Figure()

    # Prepare the ticks
    n_ticks = 10;
    q_aux_ticks = q_mean;
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

#=
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

ax2 = Axis(fig[1:4, 1:5], yticklabelcolor = :red, yaxisposition = :right)
hidespines!(ax2)
hidexdecorations!(ax2)


#scatterlines!(ax2, log.(10,q_mean), dSq_mean_smooth, color = (:red,0.5))
#scatterlines!(ax2, log.(10,q_mean), ddSq_mean_smooth, color = :red)


    lines!(ax_bottom,q_mean,fit_eval,linestyle=:dash, color =:black,
       linewidth=2.5)
=#
    
    # Add the line to the plot
    foreach(s->scatterlines!(ax_bottom,q_systems[s], Sq_systems[s],linewidth = 2),1:length(df_dat_systems)  )





            #=

            # Smooth the the information
            Sq_mean_smooth = moving_mean(Sq_mean);

            # First derivate
            dSq_mean_smooth = derivate(q_mean,Sq_mean_smooth);

            # Second derivative
            ddSq_mean_smooth = derivate(q_mean,dSq_mean_smooth);

            # Find the region to do the 1/q fit
            d2, m_aux, s_aux = norm_robust(ddSq_mean_smooth);

            #=
            # Threshold
            umbral = 10; 

            # picos locales en |d2|
            ind_peaks = Int[]
            for i in 2:length(d2)-1
                if abs(d2[i]) > abs(d2[i-1]) && abs(d2[i]) > abs(d2[i+1]) && abs(d2[i]) > umbral
                    push!(ind_peaks, i)  # índice en q
                end
            end
            =#

            # Modify the peaks to get the second derivative
            
            # Get the cut near the particle size
            q_fractal = q_mean[argmax(Sq_mean_smooth)];
            q_particle = 2*pi*0.25; # Bond distance between central particles 

            # Get the index at the middle
            ind_network = q_fractal .< q_mean .< q_particle

            # Select the region for the linear fit
            q_network = deepcopy(q_mean[ind_network])
            Sq_network = deepcopy(Sq_mean[ind_network])

            # Create the fit
            model(t,p) = (p[2])./t.^(p[1]) 

            # Set intial values for the fit
            p_initial = [1.0, 1.0];

            p_lower = [0.0, 1.0];
            p_upper = [Inf, Inf];

            # Fit the data
            fit = curve_fit(model, q_network, Sq_network, p_initial; lower=p_lower, upper=p_upper);

            # Get the parameters
            params_final = fit.param|>collect;

            # Evaluate the fit at the domain
            fit_eval = model(q_mean,params_final);


            

            #q_mean
            #Sq_mean
            #fit_eval


            #params_final



fig = Figure()

    # Prepare the ticks
    n_ticks = 10;
    q_aux_ticks = q_mean;
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

ax2 = Axis(fig[1:4, 1:5], yticklabelcolor = :red, yaxisposition = :right)
hidespines!(ax2)
hidexdecorations!(ax2)


#scatterlines!(ax2, log.(10,q_mean), dSq_mean_smooth, color = (:red,0.5))
#scatterlines!(ax2, log.(10,q_mean), ddSq_mean_smooth, color = :red)


    lines!(ax_bottom,q_mean,fit_eval,linestyle=:dash, color =:black,
       linewidth=2.5)

    # Add the line to the plot
    scatterlines!(ax_bottom, q_mean, Sq_mean,
        #linestyle = label_systems[it_system],
        linewidth = 2 
       )
    
#    scatterlines!(ax_bottom, q_mean, Sq_mean_smooth,
        #linestyle = label_systems[it_system],
#        linewidth = 2,
#        color=:orange
#       )

=#

#df_dat_categories = groupby(df_dat,[categories_experiment; categories_system]);

# Group each experiment by systems
#df_dat_experiment_system = [groupby(s,categories_system) for s in df_dat_experiment];


# Paths
#
