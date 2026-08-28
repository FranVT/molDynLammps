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

    return (steps_heat_compute, steps_isothermal_compute)
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

    return (steps_heat_compute, steps_isothermal_compute)
end
 
"""
    create_paths_compute(path_traj_timestep::String, df_set::AbstractDataFrame, steps_heat_compute::Vector{Int64}, steps_isothermal_compute::Vector{Int64})

Create the paths to the dump files to analyze
"""
function create_paths_compute(path_traj_timestep::String, df_set::AbstractDataFrame, steps_heat_compute::Vector{Int64}, steps_isothermal_compute::Vector{Int64})
    # read the directory
    files_traj_simulation=readdir(path_traj_timestep);

    # Test if the simulation is complete
    aux = Int64(first(df_set.N_heat .+ df_set.N_isothermal)/first(df_set.N_dump)) + 1

    # for mid simulation analysis
    if length(files_traj_simulation) != aux
        (steps_heat_compute, steps_isothermal_compute) = correct_steps(files_traj_simulation,steps_heat_compute,steps_isothermal_compute)
    end

    # Complete steps domain
    steps_analyze = [steps_heat_compute; steps_isothermal_compute];

    # Get the names of the files
    files_compute = string.("traj_assembly.",steps_analyze,".dumpf");

    # Create the paths
    paths_compute = joinpath.(path_traj_timestep,files_compute);

    return paths_compute
end
#=
    Start script
=#

# Paths and directories
DIR_DATA = "/run/media/franvt/rogelio/DinMol/gel-batch-3-long/data/";
DIR_MAIN = pwd();
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

# Select one set
df_set = df_groups[1];

    # Define the amoun of steps to compute the structure factor
    n_compute = 10;

    # Get the ids
    simulation_id = df_set.id;

    # Compute the tiome steps to analyze
    (steps_heat_compute, steps_isothermal_compute) = get_steps_compute(df_set,n_compute)

    # Create the paths to the dump files
    paths_traj_timestep = joinpath.(df_set.dir,"traj")

    # Select one path
    path_traj_timestep = paths_traj_timestep[1];

    paths_compute = create_paths_compute(path_traj_timestep,df_set,steps_heat_compute,steps_isothermal_compute)


        # Get the number
        #id_time_step=first(unique(id_time_step));

