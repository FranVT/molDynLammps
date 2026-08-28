#=
    Script that uses the dump information to compute observables
=#


using DataFrames, CSV
using Statistics, LsqFit

#=
    Functions
=#

"""
    get_steps_compute()

Return the time instants to analyze the structure factor for each assembly process.
"""
function get_steps_compute(df_set::AbstractDataFrame, n_compute::Int64)

    # Ge the total amount of time steps
    total_steps = first(df_set.N_heat .+ df_set.N_isothermal);

    # Create the time steps domain stored
    steps_domain = range(0,total_steps,step=first(df_set.N_dump));

    # time domain
    time_domain = first(df_set.tstep).*steps_domain;

    # Mask heat
    mask_heat  = time_domain.<=first(df_set.time_heat);

    # Time domains per process
    time_heat = time_domain[mask_heat];
    time_isothermal = time_domain[.!mask_heat];

    # Steps domain
    steps_heat = steps_domain[mask_heat];
    steps_isothermal = steps_domain[.!mask_heat];

    # Index of the time domain
    idx_heat = round.(Int, LinRange(1, length(steps_heat), n_compute));
    idx_isothermal = round.(Int, LinRange(1, length(steps_isothermal), n_compute));

    # Apply the filter
    steps_heat_compute = steps_heat[idx_heat];
    steps_isothermal_compute = steps_isothermal[idx_isothermal];

    return [steps_heat_compute; steps_isothermal_compute]
end

"""
    correct_steps(files_traj_simulation::Vector{String}, steps_heat_compute::Vector{Int64}, steps_isothermal_compute::Vector{Int64})

If simulation hasn't finish. Extract available time steps
"""
function correct_steps(files_traj_simulation::Vector{String}, steps_heat_compute::Vector{Int64}, steps_isothermal_compute::Vector{Int64})
    # Get the time step analyzed from the files
    id_time_steps=[parse(Int, match(r"traj_assembly\.(\d+)\.dumpf", s).captures[1]) for s in files_traj_simulation];
   
    # sort them
    id_max_time_step = maximum(id_time_steps);

    # Select the time that are available
    mask_heat = steps_heat_compute .<= id_max_time_steps;
    mask_isothermal = steps_isothermal_compute .<= id_max_time_steps;

    steps_heat_compute = steps_heat_compute[mask_heat];
    steps_isothermal_compute = steps_isothermal_compute[mask_isothermal];

    return [steps_heat_compute; steps_isothermal_compute]
end
 
"""
    create_paths_compute(path_traj_timestep::String, df_set::AbstractDataFrame, steps_heat_compute::Vector{Int64}, steps_isothermal_compute::Vector{Int64})

Create the paths to the dump files to analyze
"""
function create_paths_compute(path_traj_timestep::String, df_set::AbstractDataFrame, steps_analyze::Vector{Int64})
    # read the directory
    files_traj_simulation=readdir(path_traj_timestep);

    # Test if the simulation is complete
    aux = Int64(first(df_set.N_heat .+ df_set.N_isothermal)/first(df_set.N_dump)) + 1

    # for mid simulation analysis
    if length(files_traj_simulation) != aux
        steps_analyze = correct_steps(files_traj_simulation,steps_heat_compute,steps_isothermal_compute)
    end

    # Get the names of the files
    files_compute = string.("traj_assembly.",steps_analyze,".dumpf");

    # Create the paths
    paths_compute = joinpath.(path_traj_timestep,files_compute);

    return paths_compute
end

"""
    createqdom(qmax::Real, rq0::Real, dq0::Real, bin0::Real, numbin::Integer)
    
Create the reciprocal space of wave vectors
"""
function createqdom(df_set::AbstractDataFrame, qmax_0::Real)

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

    # Crear el espacio del vector recíproco para no definirlo en cada experimento
    qx = (-qmax:qmax) .* dq_0
    qy = (-qmax:qmax) .* dq_0
    qz = (-qmax:qmax) .* dq_0

    qhis   = [[] for _ in 1:n_bin]
    qxhis  = [[] for _ in 1:n_bin]
    qyhis  = [[] for _ in 1:n_bin]
    qzhis  = [[] for _ in 1:n_bin]

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
                if q > rq_0
                    continue
                end

                # Determinar el índice del bin (manejo especial del borde)
                sbin = floor(Int, q / bin_0) + 1

                # Coso para transformar de un dominio continuo a uno discreto
                if q % bin_0 == 0.0
                    sbin -= 1
                end

                # Si la magnitud del vector supera la magnitud de interés se lo salta
                if sbin > n_bin
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

    return (qxhis, qyhis, qzhis, qhis, qmean, n_bin, n_tot_av)
end

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
function get_position_simulation(path::String)
    
    # Extract the dump
    dump = get_dump(path);

    # Create a mask to only considera central particles
    mask = (dump.type .== 1) .| (dump.type .== 2.0)
    dump_filtered = dump[mask, :]

    positions = [dump_filtered.x, dump_filtered.y, dump_filtered.z]

    return reduce(hcat,positions) 
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
    Sq[:, 1] = sum.(rho) ./ length.(rho) 
    smax = maximum(Sq[:, 1])
    Sq[:, 2] = Sq[:, 1] / smax

    return Sq
end

"""
    save_Sq_analysis(info::Matrix{Float64}, steps_analyze::Vector{Int64}, n_bin::Inte64, categories_id::Vector{Symbol}, ids_set_info::Vector{Float64}, it_sim::Int64, simulation_id::String)

Function that stores the structure factor analysis
"""
function save_Sq_analysis(info::Matrix{Float64}, steps_analyze::Int64, n_bin::Int64, categories_id::Vector{Symbol}, ids_set_info::Vector{Float64}, it_sim::Int64, simulation_id::String, DIR_SAVE::String)
        
    # Create the dataframe to be stored
    df_Sq=DataFrame([repeat([steps_analyze], n_bin) info],[:timeStep, :q_mean, :Sq_mean, :Sq_mean_norm]);

    # Add the values of the categories to the dataframe 
    for (col, val) in zip(categories_id, ids_set_info)
        df_Sq[!, col] .= val 
    end

    # Add a simulation identification
    df_Sq[!,:Nsim] .= it_sim;

    # Create a file name from the ids 
    file_name=string("structure_factor_",simulation_id,"_step_",steps_analyze[1],".csv");

    # Save the information
    CSV.write(joinpath(DIR_SAVE, file_name), df_Sq)
end


#=
    Start script
=#

# Paths and directories
DIR_DATA = "/run/media/franvt/rogelio/DinMol/gel-batch-3-long/data/";
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
PATTERN_DUMP = r"traj_assembly\.(\d+)\.dumpf";

# Read the dat file
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:time_heat,:time_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# Group the metadata into categories 
df_groups = groupby(df_dat,categories_id);

# Define the amoun of steps to compute the structure factor
n_compute = 5;
q_max = 6;

# Select one set
#df_set = df_groups[1];
for df_set in df_groups
    ids_set_info=[df_set[1, col] for col in categories_id];
    simulation_id = df_set.id;

    # Get the ids
    simulation_id = df_set.id;

    # Compute the tiome steps to analyze
    steps_analyze = get_steps_compute(df_set,n_compute)

    # Create the reciprocal domain
    (qx_his,qy_his,qz_his,q_his,q_mean,n_bin,n_tot_av) = createqdom(df_set,q_max);

    # Create the paths to the dump files
    paths_traj_timestep = joinpath.(df_set.dir,"traj")

    # Select one path
    #it_sim = 1
    #path_traj_timestep = paths_traj_timestep[it_sim];
    for (it_sim,path_traj_timestep) in enumerate(paths_traj_timestep)

        # Create the paths to the files
        paths_compute = create_paths_compute(path_traj_timestep,df_set,steps_analyze);

        # Select one path 
        for (it_step,path_to_compute) in enumerate(paths_compute)

            # Exctract the position
            r = get_position_simulation(path_to_compute);

            # Save memory space to save the structure factor and compute the mean 
            info = zeros(n_bin, 3)

            # Store the structure factor
            info[:, 2:3] = computeSq(n_bin, n_tot_av, qx_his, qy_his, qz_his, q_his, r)

            # Store the analysis 
            save_Sq_analysis(info,steps_analyze[1],n_bin,categories_id,ids_set_info,it_sim,simulation_id[it_sim],DIR_SAVE)
        end # time step
    end # simulation
end # set



