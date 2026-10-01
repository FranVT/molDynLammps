#=
    Script that compute the observables per simulation

    Energy and temperature are already in the fix file
=#

using DataFrames, CSV

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

    return simulations_dir, time_steps_range[test_filter]
    
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
simulations_dir, time_steps_range = directories_to_analyze(DIR_MAIN,FILE_DAT);

# Compute the structure factor of the final configuration
time_step_to_analyze = last.(time_steps_range);

# Create the file name
file_name_dump_to_analyse = string.("traj_assembly.",time_step_to_analyze,".dumpf");

# Create the path to the file
path_dump_to_analyse = joinpath.(simulations_dir,"traj",file_name_dump_to_analyse);



#=

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:time_heat,:time_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

=#



