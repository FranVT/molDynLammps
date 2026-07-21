#=
    Script to analyze the pores 
=#

using DataFrames, CSV
using Statistics

#=
    Functions
=#
"""
    get_paths_simulation(path_dumpf_simulation::String,N_steps_Sq::Integer)

Get N_steps_Sq time steps positions paths
"""
function get_paths_simulation(path_dumpf_simulation::String,n_steps::Integer)

    # Read the directory to get the time steps stored 
    files_traj_simulation=readdir(path_dumpf_simulation);
    
    # Stored timesteps
    time_step=[parse(Int, split(s, ".")[2]) for s in files_traj_simulation]
   
    # Sort the time steps from low to high
    time_step=sort(time_step);

    if n_steps == 1
        # Get the time steps to analyse last or first
        time_step_analyze=time_step[1];

        # Re-create the filenane
        files_time_step = [replace(FILE_DUMP, "*" => string(it)) for it in time_step_analyze];

        # Create the paths
        files_traj_simulation=joinpath.(path_dumpf_simulation,files_time_step);

        return files_traj_simulation
    else
        # Create an array with the ids of the elements to analyzed
        ids_files_simulation=range(1,length(files_traj_simulation),length=n_steps);

        # Ensure integer numbers
        ids_files_simulation=Int64.(floor.(ids_files_simulation));

        # Get the time steps to analyse
        time_step_analyze=time_step[ids_files_simulation];

        # Re-create the filenane
        files_time_step = [replace(FILE_DUMP, "*" => string(it)) for it in time_step_analyze];

        # Create the paths
        files_traj_simulation=joinpath.(path_dumpf_simulation,files_time_step);

        return files_traj_simulation
    end
end

"""
    get_dump(path::String)

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
    wrap(dr::Real, L::Real)

Compute distance tacking into account periodic boundary conditions
"""
function wrap(dr::Real, L::Real)
    return dr - round(dr / L) * L
end

"""
    probe_sphere(rad::Real, l_x::Real, l_y::Real, l_z::Real, R_CP::Real)

A function that select a random point within a simulation box.
That point is the center of a sphere which raidus increases until colides with a particle.
Then return the radius of the sphere at colission 
"""
function probe_sphere(position_timestep_simulation::Vector{Vector{Float64}},l_x::Real, l_y::Real, l_z::Real, delta_rad::Real, R_CP::Real)
    # Initialize values 
    no_overlap=true;    # Status of no overlap 
    rad=0.0;            # Sphere radius

    # Create a random sphere center within the simulation box
    center_rand=rand(3).*[l_x;l_y;l_z];

    # Start the cycle
    while no_overlap
        # Increase the size of the sphere
        rad+=delta_rad;

        # Compute the distance with all particles
        d_x=center_rand[1] .- position_timestep_simulation[1];
        d_y=center_rand[2] .- position_timestep_simulation[2];
        d_z=center_rand[3] .- position_timestep_simulation[3];

        # Consider PBC
        d_x=wrap.(d_x,l_x);
        d_y=wrap.(d_y,l_y);
        d_z=wrap.(d_z,l_z);

        # Compute the magnitude of the distance
        dist=sqrt.(sum([d_x d_y d_z].^2,dims=2));

        # Create a boolean test
        overlap_test=dist.<=(R_CP+rad);

        if sum(overlap_test) != 0
            no_overlap=false;
            return rad
            break
        end
    end
end

"""
    compute_pore_histogram(n_samples::Integer,l_x::Real, l_y::Real, l_z::Real, delta_rad::Real, R_CP::Real)

Function that probes spheres enough times to create a histogram in order to get a pore distribution
"""
function compute_pore_histogram(n_samples::Integer, position_timestep_simulation::Vector{Vector{Float64}}, l_x::Real, l_y::Real, l_z::Real, delta_rad::Real, R_CP::Real)
    
    # Initialize the variable 
    hist_pore=zeros(n_samples);

    # Initialize auxiliary variables
    aux_count=1;        # Counts
    aux_ave=0.0;        # Average

    for _ in 1:n_samples 
        rad=probe_sphere(position_timestep_simulation,l_x,l_y,l_z,delta_rad,R_CP)

        # Store the information at the counter
        hist_pore[aux_count]=rad;

        # Increase the counter
        aux_count+=1;

        # Update the average varaible
        aux_ave+=rad;
    end

    # Compute the mean radius
    mean_pore=aux_ave/aux_count;

    return (hist_pore, mean_pore)
end

"""
    save_pore_histogram(df_set::AbstractDataFrame, n_samples::Integer, n_steps::Integer, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol}, DIR_SAVE::String)

This function stores the histogram of pores of a set of simulations in a time domain
"""
function save_pore_histogram(df_set::AbstractDataFrame, n_samples::Integer, n_steps::Integer, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol}, DIR_SAVE::String)

    # Definition of parameters for the analysis
    R_CP=0.5*first(unique(df_set.r_CP));                           # Radius of the central particles
    delta_rad=0.2;                      # Increment of the radius
    l=2*first(unique(df_set.L));        # Compute the length of the simulation box of the experiments
    l_x = l;                            # Length of the box at x 
    l_y = l;                            # Length of the box at y 
    l_z = l;                            # Length of the box at z
    n_cp=first(unique(df_set.N_PP));    # Amount of central particles

# Get the directories of all experiments of the system 
    dir_set=String.(df_set.dir);

    # Create the paths to the files
    path_dumpf=joinpath.(dir_set,"traj");

    # Get all central particles position of all simulations at a given time domain 
    paths_dumpf_simulations=get_paths_simulation.(path_dumpf,n_steps);

    # Define a set of categories
    categories_total=[categories_system; categories_experiment];

    # Create a file name from the categories_total
    ids_set_info=[df_set[1, col] for col in categories_total];

    # Iterate per each time step in each simulation 
    for (it_sim,paths_dumpf_simulation) in enumerate(paths_dumpf_simulations)

        # Get the time step analyzed from the files
        ids_time_step=[parse(Int, match(r"traj_assembly\.(\d+)\.dumpf", s).captures[1]) for s in paths_dumpf_simulation];

        # Initialize a variable to store histograms
        histograms_pore=[zeros(n_samples) for _ in 1:n_steps];

        # Initialize a variable to mean pore 
        means_pore=zeros(n_steps);

        # Data frame of the dump
        df_dump_timestep=get_dump.(paths_dumpf_simulation);

        # Get the positions
        position_timestep_simulations=get_position_simulation.(df_dump_timestep);

        for (it,position_timestep_simulation) in enumerate(position_timestep_simulations)
            # Get the histogram and mean pore of a configuration
            (hist_pore, mean_pore)=compute_pore_histogram(n_samples,position_timestep_simulation,l_x,l_y,l_z,delta_rad,R_CP);

            # Store the histogram
            histograms_pore[it].=hist_pore;

            # Store the mean
            means_pore[it]=mean_pore;

            println(it," time step done of ",n_steps," time steps.")
        end # for each time step

        # Reshape the array
        histograms_pore_set=reduce(vcat,histograms_pore);

        # Compute the mean of the set
        mean_pore_set=mean(means_pore);

        # Store the information per timestep
        for it_time in eachindex(ids_time_step)
            df_to_store=DataFrame([repeat([ids_time_step[it_time]], n_samples*n_steps) repeat([mean_pore_set], n_samples*n_steps) histograms_pore_set],
                   [:timeStep, :mean_pore_set, :histogram_pore]);
            
            for (col, val) in zip(categories_total, ids_set_info)
                df_to_store[!, col] .= val 
            end

            # Create a file name from the ids 
            file_name=string("pore_analysis_",join(string.(ids_set_info)),"_step_",ids_time_step[it_time],"_simulation_",it_sim,".csv");

            # Save the information
            CSV.write(joinpath(DIR_SAVE, file_name), df_to_store)

            println(file_name," stored.")
        end # For each time step

    end # For each simulation

end

#=
    Script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
FILE_DUMP = "traj_assembly.*.dumpf";


categories_system=[:phi,:chi_4,:temp,:damp,:tstep];    # Select the categories that define a system
categories_experiment=[:N_heat,:N_isothermal];  # Create categories to select different experiments (Just in case)

# Read the dat file
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Group by system 
df_systems=groupby(df_dat,categories_system);

# Set the parameters
n_steps=2;                          # Amount of time steps to analyzed
n_samples=100000;                     # Amount of sphere samples to construct the histogram

# Save the mean of S(q) of a system given a set of experiments and a time domain
for df_system in df_systems
    # Group by experiments
    df_experiments=groupby(df_system,categories_experiment);

    # Save the mean of S(q) of an experiment given a set of simulation and a time domain 
    foreach(df_set->save_pore_histogram(df_set,n_samples,n_steps,categories_system,categories_experiment,DIR_SAVE) ,df_experiments)

    println("One system done")
end

