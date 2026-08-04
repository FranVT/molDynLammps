#   Start to analyse the topology of the stuff


using DataFrames, CSV
using NearestNeighbors, Graphs
using GLMakie, GraphMakie


#=
    FUNCTIONS
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
        # Get the time steps to analyse last
        time_step_analyze=time_step[end];

        # Re-create the filenane
        files_time_step = [replace(FILE_DUMP, "*" => string(it)) for it in time_step_analyze];

        # Create the paths
        files_traj_simulation=joinpath.(path_dumpf_simulation,files_time_step);

        return [files_traj_simulation]
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
    
    # Create the dataframe
    df = DataFrame(INFO, HEADERS)

    # Filter the data frame to only give the central particle information
    mask = (df.type .== 1) .| (df.type .== 2.0)
    dump_filtered = df[mask, :]
    
    return dump_filtered
end

"""
    get_position_simulation(dump::DataFrame)

Get the position of the central particles of a given dump
"""
function get_position_simulation(dump::DataFrame)
    return [dump.x'; dump.y'; dump.z']
end

"""
    dist_pbc(dr::Real, L::Real)

Compute distance tacking into account periodic boundary conditions
"""
function dist_pbc(dr::Real, L::Real)
    return dr - round(dr / L) * L
end

"""
    explore_chain(id_neigh::Int64,visited_id::Vector{Int64},id_to_pos::Dict{Int64, Vector{Float64}},ind_to_id::Dict{Int64, Int64},id_to_type::Dict{Int64, Int64},graph::SimpleGraph{Int64})

Explore nearest neighbors and add edges to a graph 
"""
function explore_chain(id_neigh::Int64,visited_id::Vector{Int64},id_to_pos::Dict{Int64, Vector{Float64}},id_to_ind::Dict{Int64, Int64},ind_to_id::Dict{Int64, Int64},id_to_type::Dict{Int64, Int64},graph::SimpleGraph{Int64},tree_pbc)
            # Add the id to the vissited
            push!(visited_id,id_neigh)

            # Get its position
            pos_id = id_to_pos[id_neigh];

            # Get the first nearest neighbors
            inds_neigh, dists_neigh = knn(tree_pbc,pos_id,3) # For monomers

            # Pass from index to id
            ids_neigh = map(s-> ind_to_id[s], inds_neigh); 

            # Filter by visited ids and stuff
            mask_ids = mapreduce(s->ids_neigh .!= s,.&,visited_id);
 
            # Filter by distance
            mask_distance = 1.4 .>= dists_neigh .> 1.0; 

            # Combine masks
            mask = mask_ids .& mask_distance;

            # Update the array with the filters
            inds_neigh_filter = inds_neigh[mask];
            ids_neigh_filter = ids_neigh[mask];
            dists_neigh_filter = dists_neigh[mask];

            # Know the types of the neighbors
            #type_neigh_filter = map(s-> id_to_type[s],ids_neigh_filter);

            # Check if there is a Crosslinker
            #if isempty(findall(==(1),type_neigh_filter)) == false
            #    println("Found a Crosslinker")
            #end

            # Add the connections
            #foreach(s-> add_edge!(graph,id_to_ind[id_neigh], id_to_ind[s]), ids_neigh_filter)

        for s in ids_neigh_filter
            # Skip the crosslinkers as neighbors
            if id_to_type[s] == 1  
                continue
            end
            
            # Add the bond with a monomer
            add_edge!(graph,id_to_ind[id_neigh], id_to_ind[s])
        end




    return (graph,visited_id,ids_neigh_filter)

end
 
"""
    explore_nodes(id_cl,graph,visited_id,tree_pbc)

Explore the chains around each node (Cross Linkers)
"""
function explore_nodes(id_cl,graph,visited_id,tree_pbc)

        # Add the id to the vissited
        push!(visited_id,id_cl)

        # Get its position
        pos_id = id_to_pos[id_cl];

        # Get the first nearest neighbors
        inds_neigh, dists_neigh = knn(tree_pbc,pos_id,6) # Get the 4 patches

        # Pass from index to id
        ids_neigh = map(s-> ind_to_id[s], inds_neigh); 

        # Filter by visited ids and stuff
        mask_ids = mapreduce(s->ids_neigh .!= s,.&,visited_id);
 
        # Filter by distance
        mask_distance = 1.4 .>= dists_neigh .> 1.0; 

        # Combine masks
        mask = mask_ids .& mask_distance;

        # Update the array with the filters
        inds_neigh_filter = inds_neigh[mask];
        ids_neigh_filter = ids_neigh[mask];
        dists_neigh_filter = dists_neigh[mask];

        # Know the types of the neighbors
        #type_neigh_filter = map(s-> id_to_type[s],ids_neigh_filter);

        for s in ids_neigh_filter
            # Skip the crosslinkers as neighbors
            if id_to_type[s] == 1  
                println("Warning: Two crosslinkers as neighbors") 
            end
            
            # Add the bond with a monomer
            add_edge!(graph,id_to_ind[id_cl], id_to_ind[s])
        end

        # Add the connections
        #foreach(s-> add_edge!(graph,id_to_ind[id_cl], id_to_ind[s]), ids_neigh_filter)

    # Search each neighbor
       # Inicializar la cola con los primeros vecinos
        to_explore_id = copy(ids_neigh_filter);
       
        while !isempty(to_explore_id);

            id_neigh = popfirst!(to_explore_id);

        #for id_neigh in ids_neigh_filter
            graph,visited_id,new_ids = explore_chain(id_neigh,visited_id,id_to_pos,id_to_ind,ind_to_id,id_to_type,graph,tree_pbc)
        #end

            new_ids = setdiff(new_ids,visited_id);

            append!(to_explore_id, new_ids)         

        end

        return (graph,visited_id)
end



#=
    SCRIPT 
=#

# Paths and directories
DIR_MAIN = dirname(pwd());
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
FILE_DUMP = "traj_assembly.*.dumpf";

# Critera to diffirientiate from system and experiment
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];    # Select the categories that define a system
categories_experiment=[:N_heat,:N_isothermal];  # Create categories to select different experiments (Just in case)

# Select the amount of time steps to analyze
n_steps = 1; # Implies the final configuration

# Read the dat file
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Group by system 
df_systems=groupby(df_dat,categories_system);

# Select one system
df_system = df_systems[1];

    # Group by experiments
    df_experiments=groupby(df_system,categories_experiment);

    # Select all simulations of one experiment
    df_set = df_experiments[1];

    N_cpart =  first(df_set.N_PP);
    l=first(unique(df_set.L));        # Compute the length of the simulation box of the experiments
    l_x = 2*l;                            # Length of the box at x 
    l_y = 2*l;                            # Length of the box at y 
    l_z = 2*l;                            # Length of the box at z

    # Get the directories of all simulations of the system 
    dir_set=String.(df_set.dir);

    # Create the paths to the files
    path_dumpf=joinpath.(dir_set,"traj");

    # Get all central particles position of all simulations at a given time domain 
    paths_dumpf_simulations=get_paths_simulation.(path_dumpf,n_steps);

    # Iterate per each time step in each simulation 
    #for (it_sim,paths_dumpf_simulation) in enumerate(paths_dumpf_simulations)
    it_sim = 1;
    paths_dumpf_simulation = paths_dumpf_simulations[1];


        # Get the time step analyzed from the files
        ids_time_step=[parse(Int, match(r"traj_assembly\.(\d+)\.dumpf", s).captures[1]) for s in paths_dumpf_simulation];

        # Data frame of the dump
        df_dump_timesteps=get_dump.(paths_dumpf_simulation);

        # Select one time step
        df_dump = df_dump_timesteps[1];

        # Parse the id floeat64 to Int64
        df_dump.id = Int64.(df_dump.id)

        # Same to particle types
        df_dump.type = Int64.(df_dump.type)

        # Parse the mol floeat64 to Int64
        df_dump.mol = Int64.(df_dump.mol)

        # Same to cluster id 
        df_dump."c_clusters" = Int64.(df_dump."c_clusters")

        # Sort the dataframe by id
        # IMPORTANT FOR POS_TO_ID DICTIONARY
        sort!(df_dump,:id);

        # Prepare Nerest Neighbor searh
        positions=get_position_simulation(df_dump)

        # Dictionary from ind to id
        ind_to_id =  Dict{Int64, Int64}();

        # fill the Dictionary
        foreach(s->ind_to_id[s] = df_dump.id[s], 1:N_cpart);

        # Dictionary from ind to id
        id_to_ind = Dict(zip(values(ind_to_id), keys(ind_to_id)));

        # Dictionary from id to type
        id_to_type = Dict(zip(df_dump.id,df_dump.type));

        # Dictionary from positions to id
        pos_to_id = Dict{Vector{Float64}, Int64}();

        # Create the dictirionary
        foreach(s->pos_to_id[positions[:,s]] = df_dump.id[s], 1:N_cpart);

        # Dictionary from id to positions
        id_to_pos = Dict(zip(values(pos_to_id), keys(pos_to_id)));

        # Create a tree
        tree=KDTree(positions)

        # Consider periodic boundary conditions
        tree_pbc=PeriodicTree(tree,[-l, -l, -l],[l, l, l])

        # Group by particle type
        df_by_types=groupby(df_dump,:type)

        # The initial nodes are going to be all the crosslinkers
        ids_cl = df_by_types[1].id[:];

        ids_mo = df_by_types[2].id[:];


        # Start the graph
        graph=SimpleGraph(N_cpart); # Based on index

        # Create an array with id visited
        visited_id=Array{Int64,1}();

        #explore_id=Array{Int64,1}();

# tests 

        # Start with one crosslinker
        #id_cl = ids_cl[1];

        # Create tha graphs with only the crosslinkers
        for id_cl in ids_cl
            global graph,visited_id = explore_nodes(id_cl,graph,visited_id,tree_pbc)
        end

        # Create a copy of the graph for further modification
        graph_cl = deepcopy(graph);

    # Clean the graph created by the crosslinkers
        # Find all particles with degree 0 (There are isolated)
        inds_isolated = findall(degree(graph_cl).==0)

        # Remove isolated particles
        rem_vertices!(graph_cl,inds_isolated)

        if (length(visited_id) - length(unique(visited_id))) != 0
            println("Three body particle bonds found")
        end


    # Continue creating monomers chains
        notvisited_ids=setdiff([ids_cl; ids_mo], unique(visited_id))
       
        # Finish the lonely polymer chains        
        to_explore_id = copy(notvisited_ids);
       
        while !isempty(to_explore_id);

            id_neigh = popfirst!(to_explore_id);

        #for id_neigh in ids_neigh_filter
            global graph,visited_id,new_ids = explore_chain(id_neigh,visited_id,id_to_pos,id_to_ind,ind_to_id,id_to_type,graph,tree_pbc)
        #end

            new_ids = setdiff(new_ids,visited_id);

            append!(to_explore_id, new_ids)         

        end

        # Create a list with the inds of particles in a cluster
        list_inds_clusters=connected_components(graph);

    # Remove clusters of one elements

        # Create a mask of one
        mask_cluster = length.(list_inds_clusters) .!= 1;

        list_inds_clusters=list_inds_clusters[mask_cluster];


        #for id_notvisited in notvisited_ids
        #    global graph,visited_id,new_ids = explore_chain(id_notvisited,visited_id,id_to_pos,id_to_ind,ind_to_id,id_to_type,graph,tree_pbc)
        #end

#explore_chain(id_neigh::Int64,visited_id::Vector{Int64},id_to_pos::Dict{Int64, Vector{Float64}},id_to_ind::Dict{Int64, Int64},ind_to_id::Dict{Int64, Int64},id_to_type::Dict{Int64, Int64},graph::SimpleGraph{Int64},tree_pbc)

     #end # For enumerate(paths_dumpf_simulations)




