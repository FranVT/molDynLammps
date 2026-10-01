#=
    Script that compute the observables per simulation

    Energy and temperature are already in the fix file
=#

using DataFrames, CSV
using Statistics

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

    return simulations_dir, time_steps_range[test_filter], df_dat[test_filter,:]
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
    get_position_simulation(path::String)

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
    createqdom(qmax::Real, rq0::Real, dq0::Real, bin0::Real, numbin::Integer)
    
Create the reciprocal space of wave vectors
"""
function createqdom(l::Real, n_cp::Int64, qmax_0::Real)

    # Parameters
    l=2*l; #2*first(unique(df_set.L));    # Compute the length of the simulation box of the experiments
    #n_exp=nrow(df_set);    # Compute the amount of simulations
    #n_cp=first(unique(df_set.N_PP));    # Amount of central particles

    # Crea los dominios del vector de onda
    x_c = l;                 # longitud de la caja en x
    y_c = l;                 # longitud de la caja en y
    z_c = l;                 # longitud de la caja en z
    dq_0 = 2 * pi / x_c;      # Δq fundamental
    qmax = Int(floor(qmax_0 / dq_0));   # número entero de pasos hasta qmax0
    rq_0 = qmax * dq_0;       # valor máximo real de |q| usado
    bin_0 = dq_0;             # nuevo ancho de bin (bin0 original * dq0)
    n_bin = Int(floor(qmax * dq_0 / bin_0));  # número total de bines

    # Parametros para el factor de estructura
    n_tot_av = n_cp;   # Número total de partículas

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

    return qmean, n_bin
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

# Compute Observables from dump files 

# Select a time step to analyze
time_step_to_analyze = last.(time_steps_range);

# Create the file name
file_name_dump_to_analyse = string.("traj_assembly.",time_step_to_analyze,".dumpf");

# Create the path to the file
path_dump_to_analyse = joinpath.(simulations_dir,"traj",file_name_dump_to_analyse);

# Compute the structure factor for each simulation

# Get the length of each simulation
box_length = df_dat.L;

# Get the amount of central particles in the simulation
N_central = df_dat.N_PP; 

# Select one simulation
it_sim = 1;

    # Compute the domain
    q_mean, n_bin = createqdom(box_length[it_sim], N_central[it_sim], 2*pi);

    # Extract the positions
    r = get_position_simulation(path_dump_to_analyse[it_sim]);

    # Save memory space to save the structure factor and compute the mean 
    info = zeros(n_bin, 3)

    # Store the q domain
    info[:,1] = q_mean;

    # Store the structure factor
    # S_q and S_q/max(S_q)
    info[:, 2:3] = computeSq(n_bin, n_tot_av, qx_his, qy_his, qz_his, q_his, r);




#    j


# Create the wave vector domain
#


#=

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:time_heat,:time_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

=#



