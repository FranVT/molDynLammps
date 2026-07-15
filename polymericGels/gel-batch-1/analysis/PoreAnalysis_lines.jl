# =========================================================================
# Program:      PoreAnalysis_lines 
# Description:  Calculate the distribution of pore lengths in a system
#               of spheres with a radius of 0.5, launching segments between pairs of
#               particles chosen at random and verifying that they do not intersect
#               any other particle.
#               Reads multiple dump files (LAMMPS dump) and writes a
#               histogram for each configuration.
# =========================================================================

using DataFrames, CSV
using Statistics, Random

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
        time_step_analyze=time_step[end];

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

    return reduce(hcat,[dump_filtered.x, dump_filtered.y, dump_filtered.z])
end

"""
    wrap(dr::Real, L::Real)

Compute distance tacking into account periodic boundary conditions
"""
function wrap(dr::Real, L::Real)
    return dr - round(dr / L) * L
end

"""
    detect_collision(position_simulation::Matrix{Float64}, R_CP::Float64, l_x::Float64, l_y::Float64, l_z::Float64, n_cp::Int64, bin_size::Float64, n_bins::Int64)

Create a histogram of the length of the pores 
"""
function detect_collision(position_simulation::Matrix{Float64}, R_CP::Float64, l_x::Float64, l_y::Float64, l_z::Float64, n_cp::Int64, bin_size::Float64, n_bins::Int64)
   
    # Allocate memory
    hist_pore_length = zeros(n_bins);

    # Start counter of crashes
    choles = 0;
    
    while choles <= n_samples
        # Select one random particle
        ind_i = rand(1:n_cp);

        # Select another random particle
        ind_j = ind_i;
        while ind_j == ind_i
            ind_j = rand(1:n_cp);
        end

        # Get positions
        pos_i = position_simulation[ind_i,:];
        pos_j = position_simulation[ind_j,:];

        # Compute the distance between particles with periodic boundary conditions
        v_ij = wrap.(pos_j .- pos_i , [l_x, l_y, l_z]);

    # Create two unit vectors
        v_1=2.0 .* randn(3) .- 1.0;
        u_1=v_1/sqrt(v_1'*v_1);

        v_2=2.0 .* randn(3) .- 1.0;
        u_2=v_2/sqrt(v_2'*v_2);

        # Create the direction of the line
        v_dir=u_2 .- u_1;

        # Move the vector distance into a random direction considering the radius of the particles 
        v_line = v_ij .+ (R_CP).*v_dir;

        # Compute unit vector
        d_line = sqrt(v_line'*v_line); # d_12
        u_line = v_line/d_line; # eu12

        id_bin = trunc(Int64,d_line/bin_size);

        if id_bin > n_bins 
            continue # Ignore if it is outside the domain
        end

# Check if there is a collision 

        # Auxiliary to identify collisions
        collision_1 = 0;

        # Check if it collides with intermediate spheres
        no_overlap = 0;     # Parameter

        # Run over all particles 
        for ind_part in 1:n_cp 
            # Pass over the reference particles
            if ind_part == ind_i || ind_part == ind_j
                continue
            end
    
            # Auxiliary to identify collisions
            collision_2=0;

            # Select one particle
            pos_k = position_simulation[ind_part];

        # Distance i-k 
            v_ik = wrap.(pos_k .- pos_i , [l_x, l_y, l_z]);

            # Displace to the surface of the central particle 
            vd_ik = v_ik .+ (R_CP)*u_1;

            # Compute distance
            dd_ik = sqrt(vd_ik'*vd_ik); # d1
            
            # Unit vector
            ud_ik = vd_ik/dd_ik; # r_ij

            # Wierd proyection
            test=ud_ik'*u_line; # eu1eu2

            if test<=0.0 
                # Collision
                collision_2 += 1
            else
                aux = test*dd_ik; # dp
                aux_ik = sqrt(dd_ik^2 - aux^2)
                if aux_ik < R_CP && aux < d_line
                    # Collision
                    no_overlap += 1
                    break # Break from the for
                end
            end

        # Distance j-k
            v_jk = wrap.(pos_k .- pos_j , [l_x, l_y, l_z]);

            # Displace to the surface of the central particle 
            vd_jk = v_jk .+ (R_CP)*u_2;

            # Compute distance
            dd_jk = sqrt(vd_jk'*vd_jk); # d1
            
            # Unit vector
            ud_jk = vd_jk/dd_jk; # r_ij

            # Wierd proyection
            test=ud_jk'*-u_line; # eu1eu2

            if test<=0.0 
                # Collision
                collision_2 += 1
            else
                aux = test*dd_jk; # dp
                aux_jk = sqrt(dd_jk^2 - aux^2)
                if aux_jk < R_CP && aux < d_line
                    # Collision
                    no_overlap += 1
                    break # Break from the for
                end
            end

            if collision_2 == 0
                collision_1 += 1
            end

        end # End of the for

        if no_overlap == 0 && collision_1 > 0
            hist_pore_length[id_bin] += 1
            choles += 1
        end

    end

    return hist_pore_length
end

"""
    store_pore_length_hist_experiment(df_set::AbstractDataFrame, n_steps::Int64, n_samples::Int64, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol}, DIR_SAVE::String)

Store pore length histogram per experiment
"""
function store_pore_length_hist_experiment(df_set::AbstractDataFrame, n_steps::Int64, n_samples::Int64, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol}, DIR_SAVE::String)
    
    # Definition of parameters for the analysis
    R_CP=0.5;                           # Radius of the central particles
    l_half=first(unique(df_set.L));        # Compute the length of the simulation box of the experiments
    l_half_x = l_half;
    l_half_y = l_half;
    l_half_z = l_half;
    l_x = 2*l_half_x;                            # Length of the box at x 
    l_y = 2*l_half_y;                            # Length of the box at y 
    l_z = 2*l_half_z;                            # Length of the box at z
    n_cp=Int64(first(unique(df_set.Npart)));    # Amount of central particles
    bin_size=0.2;       # Size of the bins. Is equal to the patches radius

    # Compute the amount of bins
    n_bins = trunc(Int64,sqrt(4*l_x^2 + l_y^2 + l_z^2)/bin_size) + 1;

# Get the directories of all experiments of the system 
    dir_set=String.(joinpath.(df_set.PARENT_DIR,df_set.dir));

    # Create the paths to the files
    path_dumpf=joinpath.(dir_set,"traj");

    # Get all central particles position of all simulations at a given time domain 
    paths_dumpf_simulations=get_paths_simulation.(path_dumpf,n_steps);

    # Create the domain for the bins
    pore_domain=range(start=0,length=n_bins,step=bin_size);

    # Define a set of categories
    categories_total=[categories_system; categories_experiment];

    # Create a file name from the categories_total
    ids_set_info=[df_set[1, col] for col in categories_total];

    # Iterate per each time step in each simulation 
    for (it_sim,paths_dumpf_simulation) in enumerate(paths_dumpf_simulations)

        # Get the time step analyzed from the files
        ids_time_step=[parse(Int, match(r"traj_assembly\.(\d+)\.dumpf", s).captures[1]) for s in paths_dumpf_simulation];

        # Allcoate memory
        hist_pore_length = zeros(n_steps,n_bins);    # Histogram of the lengths of the pores

        # Iterate thru all the desired time steps
        for (it,path_dumpf_simulation) in enumerate(paths_dumpf_simulation)
            # Get the dump
            df_dump_timestep=get_dump(path_dumpf_simulation);

            # Get the positions
            position_simulation=get_position_simulation(df_dump_timestep);

            hist_pore_length[it,:] = detect_collision(position_simulation,R_CP,l_x,l_y,l_z,n_cp,bin_size,n_bins);
            println(it," time step done of ",n_steps," time steps.")
        end # for of histogram per time step

        # Store the information per timestep
        for it_time in eachindex(ids_time_step)
            df_to_store=DataFrame([pore_domain, hist_pore_length[it_time,:]],[:poreDomain, :histLength]);
            for (col, val) in zip(categories_total, ids_set_info)
                df_to_store[!, col] .= val 
            end

            # Add the time step to the dataframe
            df_to_store[!,:timeStep] .= ids_time_step[it_time]

            # Create a file name from the ids 
            file_name=string("pore_length_histogram_",join(string.(ids_set_info)),"_step_",ids_time_step[it_time],"_simulation_",it_sim,".csv");

            # Save the information
            CSV.write(joinpath(DIR_SAVE, file_name), df_to_store)

            println(file_name," stored.")
        end # for of dave dataframes

    end # For per simulation

end # Of the function


#=
    Script
=#

# Set sed
Random.seed!(1234)

# Paths and directories
DIR_MAIN = pwd();
DAT_PATH   = joinpath(DIR_MAIN, "datFiles")
DIR_SAVE = joinpath(DIR_MAIN,"analyzedData");
FILE_DAT = "systemDatfiles.csv";
FILE_FIX = "system_assembly.fixf";
FILE_DUMP = "traj_assembly.*.dumpf";


categories_system=[:phi,Symbol("CL-Con"),:Temperature,:damp];    # Select the categories that define a system
categories_experiment=[:N_heat,:N_isot];  # Create categories to select different experiments (Just in case)

# Read the dat file
df_dat=CSV.read(joinpath(DAT_PATH,FILE_DAT), DataFrame);

# Group by system 
df_systems=groupby(df_dat,categories_system);

# Set the parameters
n_steps=2;                          # Amount of time steps to analyzed
n_samples=1000;                    # Amount of tries per simulation

# Save the mean of S(q) of a system given a set of experiments and a time domain
for df_system in df_systems
    # Group by experiments
    df_experiments=groupby(df_system,categories_experiment);

    # Store per each set of experiments
    foreach(df_set->store_pore_length_hist_experiment(df_set,n_steps,n_samples,categories_system,categories_experiment,DIR_SAVE), df_experiments)

    println("One system done")
end

