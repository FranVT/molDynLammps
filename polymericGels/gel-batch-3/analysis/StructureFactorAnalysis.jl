#=
    Script for the structure factor analysis
=#

using DataFrames, CSV
using Statistics

#=
    Functions
=#

"""
    getDump(path::String)

Get the data from a single dump file that stores one timeste information
"""
function get_dump(path::String)
    data = split.(readlines(path), " ")[9:end]
    HEADERS = data[1][3:end]
    INFO = parse.(Float64, reduce(hcat, data[2:end]))'

    return DataFrame(INFO, HEADERS)
end

"""
    get_position_simulation(dump::DataFrame)

Get the position of the central particles of a given dump
"""
function get_position_simulation(dump::DataFrame)
    mask = (dump.type .== 1) .| (dump.type .== 2.0)
    dump_filtered = dump[mask, :]

    return [dump_filtered.x, dump_filtered.y, dump_filtered.z]
end

"""
    get_paths_simulation(path_dumpf_simulation::String,N_steps_Sq::Integer)

Get N_steps_Sq time steps positions paths
"""
function get_paths_simulation(path_dumpf_simulation::String,N_steps_Sq::Integer)

    # Read the directory to get the time steps stored 
    files_traj_simulation=readdir(path_dumpf_simulation);
    
    # Stored timesteps
    time_step=[parse(Int, split(s, ".")[2]) for s in files_traj_simulation]
   
    # Sort the time steps from low to high
    time_step=sort(time_step);

    # Create an array with the ids of the elements to analyzed
    ids_files_simulation=range(1,length(files_traj_simulation),length=N_steps_Sq);

    # Ensure integer numbers
    ids_files_simulation=Int64.(floor.(ids_files_simulation));

    # Get the time steps to analyse
    time_step_analyze=time_step[ids_files_simulation];

    # Re-create the filenane
    files_time_step = [replace(FILE_DUMP, "*" => string(it)) for it in time_step_analyze];

    # Create the paths
    files_traj_simulation=joinpath.(path_dumpf_simulation,files_time_step);

    #return get_position_simulation.(df_trajectories)
    return files_traj_simulation
end

"""
    createqdom(qmax::Real, rq0::Real, dq0::Real, bin0::Real, numbin::Integer)
    
Create the reciprocal space of wave vectors
"""
function createqdom(qmax::Real, rq0::Real, dq0::Real, bin0::Real, numbin::Integer)
    # Crear el espacio del vector recíproco para no definirlo en cada experimento
    qx = (-qmax:qmax) .* dq0
    qy = (-qmax:qmax) .* dq0
    qz = (-qmax:qmax) .* dq0

    qhis   = [[] for _ in 1:numbin]
    qxhis  = [[] for _ in 1:numbin]
    qyhis  = [[] for _ in 1:numbin]
    qzhis  = [[] for _ in 1:numbin]

    # Calcular la magnitud cuadrada
    for x in qx
        for y in qy
            for z in qz
                # Excluir el origen del espacio recíproco
                if x == 0 && y == 0 && z == 0
                    continue
                end

                q = sqrt(x^2 + y^2 + z^2)

                # Si la magnitud es grande que el maximo se salta
                if q > rq0
                    continue
                end

                # Determinar el índice del bin (manejo especial del borde)
                sbin = floor(Int, q / bin0) + 1

                # Coso para transformar de un dominio continuo a uno discreto
                if q % bin0 == 0.0
                    sbin -= 1
                end

                # Si la magnitud del vector supera la magnitud de interés se lo salta
                if sbin > numbin
                    continue
                end

                append!(qhis[sbin], q)
                append!(qxhis[sbin], x)
                append!(qyhis[sbin], y)
                append!(qzhis[sbin], z)
            end
        end
    end

    # Compute the mean handeling the empty vectors
    qmean=[isempty(v) ? 0.0 : mean(v) for v in qhis];

    return (qxhis, qyhis, qzhis, qhis, qmean)
end

"""
    computeSq(numbin::Integer, ntotav::Integer, qxhis, qyhis, qzhis, qhis, r)

Compute the structure factor of a set of positions
"""
function computeSq(numbin::Integer, ntotav::Integer, qxhis, qyhis, qzhis, qhis, r)
    Sq = zeros(numbin, 2)
    rho = [[] for _ in 1:numbin]

    for it_bin in 1:numbin
        # Seleccion de las componentes del vector de onda
        qx = qxhis[it_bin]
        qy = qyhis[it_bin]
        qz = qzhis[it_bin]

        # Manage null vectors
        if isempty(qx) || isempty(qy) || isempty(qz)
            append!(rho[it_bin], 0.0)
        else
                    # Calcular la densidad para cada bin
            for it_q in eachindex(qx)
                vq = [qx[it_q], qy[it_q], qz[it_q]]
                dp = r * vq  # Producto punto del vector para cada partícula
                rho_re = sum(cos.(dp))
                rho_im = sum(sin.(dp))
                sq = (rho_re^2 + rho_im^2) / ntotav
                append!(rho[it_bin], sq)
            end
        end 

    end

    # Guardamos información
    # Valor esperado del factor de estructura
    Sq[:, 1] = sum.(rho) ./ numbin 
    smax = maximum(Sq[:, 1])
    Sq[:, 2] = Sq[:, 1] / smax

    return Sq
end

"""
    save_Sq_timestep_mean()

Store a dataframe with the collective mean of S(q) of a set of simulations at a time step
"""
function save_Sq_timestep_mean(df_set::AbstractDataFrame, paths_traj_timestep::Vector{String}, n_bin::Integer, n_exp::Integer, n_tot_av::Integer, qx_his::Vector{Vector{Any}}, qy_his::Vector{Vector{Any}}, qz_his::Vector{Vector{Any}}, q_his::Vector{Vector{Any}}, q_mean::Vector{Float64}, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol}, DIR_SAVE::String)
    println("Start one save\n")

    # Get the dump
    df_timestep=get_dump.(paths_traj_timestep);

    # Get the positions
    position_timestep_simulations=get_position_simulation.(df_timestep);

    # Reduce
    position_timestep_simulations=reduce.(hcat,position_timestep_simulations);

    # Save memory space to save the structure factor and compute the mean 
    info = [zeros(n_bin, 3) for _ in 1:n_exp]

    # Compute the structure factor for all the simulations at a given time step
    println("Start the cycle")
    for it_sim in eachindex(paths_traj_timestep)
        r = position_timestep_simulations[it_sim]
        info[it_sim][:, 2:3] = computeSq(n_bin, n_tot_av, qx_his, qy_his, qz_his, q_his, r)
        println("One experiment done")
    end

    # Compute the mean over the simulations
    info = mean(info);
    
    # Store the mean of the wave vectors
    info[:, 1] = q_mean;

# Preparation to store the data

    # Get the time step analyzed from the files
    id_time_step=[parse(Int, match(r"traj_assembly\.(\d+)\.dumpf", s).captures[1]) for s in paths_traj_timestep];
    
    # Get the number
    id_time_step=first(unique(id_time_step));

    # Create the dataframe to be stored
    df_Sq=DataFrame([repeat([id_time_step], n_bin) info],
                   [:timeStep, :q_mean, :Sq_mean, :Sq_mean_norm]);

    # Define a set of categories
    categories_total=[categories_system; categories_experiment];

    # Create a file name from the categories_total
    ids_set_info=[df_set[1, col] for col in categories_total];

    # Add the values of the categories to the dataframe 
    for (col, val) in zip(categories_total, ids_set_info)
        df_Sq[!, col] .= val 
    end

    # Create a file name from the ids 
    file_name=string("structure_factor_",join(string.(ids_set_info)),"_step_",id_time_step,".csv");

    # Save the information
    CSV.write(joinpath(DIR_SAVE, file_name), df_Sq)

    println("Experiment ",file_name," saved\n")
end

"""
    save_Sq_set()

Save the collective mean structure factor of a set of simulations at a time domain 
"""
function save_Sq_set(df_set::AbstractDataFrame,n_steps_Sq::Integer, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol}, DIR_SAVE::String)
    println("Start one time series")

    # Get the directories of all experiments of the system 
    dir_set=String.(df_set.dir);

    # Create the paths to the files
    path_dumpf=joinpath.(dir_set,"traj");

    # Get all central particles position of all simulations 
    paths_dumpf_simulations=get_paths_simulation.(path_dumpf,n_steps_Sq);

    # Change the group: Each row is a time step
    paths_dumpf_simulations=reduce(hcat,paths_dumpf_simulations);

# Start the analysis of the Structure factor

    # Parameters
    l=2*first(unique(df_set.L));    # Compute the length of the simulation box of the experiments
    n_exp=nrow(df_set);    # Compute the amount of simulations
    n_cp=first(unique(df_set.N_PP));    # Amount of central particles

    # Crea los dominios del vector de onda
    x_c = l;                 # longitud de la caja en x
    y_c = l;                 # longitud de la caja en y
    z_c = l;                 # longitud de la caja en z
    dq_0 = 2 * pi / x_c;      # Δq fundamental
    qmax = Int(floor(qmax_0 / dq_0));   # número entero de pasos hasta qmax0
    rq_0 = qmax * dq_0;       # valor máximo real de |q| usado
    bin_0 = dq_0;             # nuevo ancho de bin (bin0 original * dq0)
    n_bin = Int(floor(qmax * dq_0 / bin_0)) + 1;  # número total de bines

    # Parametros para el factor de estructura
    n_tot_av = Int(n_cp);   # Número total de partículas

    # Compute the wave vectors for the experiment
    (qx_his, qy_his, qz_his, q_his, q_mean) = createqdom(qmax, rq_0, dq_0, bin_0, n_bin);

    # Store the assembly mean of a structure factor for the time steps selected
    foreach(s->save_Sq_timestep_mean(df_set,paths_dumpf_simulations[s,:],n_bin,n_exp,n_tot_av,qx_his,qy_his,qz_his,q_his,q_mean,categories_system,categories_experiment,DIR_SAVE),1:n_steps_Sq)

    println("End one time series")
end
#=
    Start of the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
FILE_FIX = "system_assembly.fixf";
FILE_DUMP = "traj_assembly.*.dumpf";

# Select the amount of time steps to analyze the structure factor
n_steps_Sq=2;

# Define the maximum wave vector
qmax_0=2;

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp];

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# Read the dat file
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Group by system 
df_systems=groupby(df_dat,categories_system);

# Save the mean of S(q) of a system given a set of experiments and a time domain
for df_aux in df_systems
    # Group by experiments
    df_experiments=groupby(df_aux,categories_experiment);

    # Save the mean of S(q) of an experiment given a set of simulation and a time domain 
    foreach(df_set->save_Sq_set(df_set,n_steps_Sq,categories_system,categories_experiment,DIR_SAVE) ,df_experiments)
end


