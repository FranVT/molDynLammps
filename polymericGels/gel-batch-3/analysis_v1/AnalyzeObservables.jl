#=
    Script to analyze observables and stuff
=# 

using DataFrames, CSV
using Statistics, LsqFit

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
df_dat_experiment = df_dat_categories[1];

    # Group by system
    df_dat_systems = groupby(df_dat_experiment,categories_system);

    # Select one system
    df_dat_system = df_dat_systems[1];

        # Get the paths to the files 
        paths = files_exist(joinpath.(df_dat_system.dir,"structure_factor.txt"));
      
        # Amount of simulations per system in one experiment
        N_sim = length(paths);

        # This only works if all files have the same timesteps stored.

        # To compute the mean
        #timestep = Array{Float64}[];
        q_mean = [];
        Sq_mean = [];

        # Exctract the Structure factor information
        for path in paths
            # Get the information
            df_Sq = read_Sqtimesteps(path)

            # Get the timestep
            #timestep = df_Sq.timestep;
            q_sim = df_Sq.q;
            Sq_sim = df_Sq.Sq;

            # Prepare for the mean 
            append!(q_mean,[q_sim])
            append!(Sq_mean,[Sq_sim])

        end
   
        # Compute the mean
        q_mean = reduce(+,q_mean)/N_sim;
        Sq_mean = reduce(+,Sq_mean)/N_sim

#=
        # Smooth the the information
        Sq_mean_smooth = moving_mean(Sq_mean);

        # First derivate
        dSq_mean_smooth = derivate(q_mean,Sq_mean_smooth);

        # Second derivative
        ddSq_mean_smooth = derivate(q_mean,dSq_mean_smooth);

        # Find the region to do the 1/q fit
        d2 = ddSq_mean_smooth;

        # Threshold
        umbral = 10; 

        # picos locales en |d2|
        ind_peaks = Int[]
        for i in 2:length(d2)-1
            if abs(d2[i]) > abs(d2[i-1]) && abs(d2[i]) > abs(d2[i+1]) && abs(d2[i]) > umbral
                push!(ind_peaks, i)  # índice en q
            end
        end

        # Modify the peaks to get the second derivative
        
        # Get the cut near the particle size
        q_fractal = q_sim[last(ind_peaks)];
        q_particle = 2*pi/1.6; # Bond distance between central particles 

        # Get the index at the middle
        ind_network = q_fractal .< q_sim .< q_particle

        # Select the region for the linear fit
        q_network = deepcopy(q_mean[ind_network])
        Sq_network = deepcopy(Sq_mean[ind_network])

        # Create the fit
        model(t,p) = (p[2])./t.^(p[1])

        # Set intial values for the fit
        p_initial = [1.0, 1.0];

        p_lower = [0.0, 0.0];
        p_upper = [Inf, Inf];

        # Fit the data
        fit = curve_fit(model, q_network, Sq_network, p_initial; lower=p_lower, upper=p_upper);

        # Get the parameters
        params_final = fit.param|>collect;

        # Evaluate the fit at the domain
        fit_eval = model(q_mean[ind_network],params_final);
=#

#df_dat_categories = groupby(df_dat,[categories_experiment; categories_system]);

# Group each experiment by systems
#df_dat_experiment_system = [groupby(s,categories_system) for s in df_dat_experiment];


# Paths
#
