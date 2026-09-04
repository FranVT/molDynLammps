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
    mask = (df.type .== 1) .| (df.type .== 2.0) .| (df.type .== 3.0) .| (df.type .== 4.0)
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
    get_patches(visited_id::Vector{Int64}, ind_to_id::Dict{Int64, Int64}, id_to_type::Dict{Int64, Int64}, pos_id::Vector{Float64}, tree_pbc)

Function that return the patches information of a crosslinker
"""
function get_patches(id_to_pos::Dict{Int64, Vector{Float64}}, ind_to_id::Dict{Int64, Int64}, id_to_type::Dict{Int64, Int64}, id::Int64, tree_pbc)
        # Get its position
        pos_id = id_to_pos[id];

        # Get the first nearest neighbors
        inds_neigh, dists_neigh = knn(tree_pbc,pos_id,5) # Get the 4 patches

        # Pass from index to id
        ids_neigh = map(s-> ind_to_id[s], inds_neigh); 

        # Get types
        types_neigh = map(s-> id_to_type[s], ids_neigh); 

        # Filter by type, only patches of crosslinkers
        mask_types = (types_neigh .== 4) .| (types_neigh .== 3)

        # Filter by distance of pathces
        mask_distance = 0.5 .>= dists_neigh .> 0.2; 

        # Combine masks
        mask = mask_distance .& mask_types;

        # Update the array with the filters
        inds_neigh_filter = inds_neigh[mask];
        ids_neigh_filter = ids_neigh[mask];
        dists_neigh_filter = dists_neigh[mask];

        return (inds_neigh_filter, ids_neigh_filter)
end


"""
    get_patch_patch

Functions that retunr patch information bonded with the patch given by pos_id
"""
function get_patch_patch(id_to_pos::Dict{Int64, Vector{Float64}}, ind_to_id::Dict{Int64, Int64}, id_to_type::Dict{Int64, Int64}, id::Int64, tree_pbc)
        # Get the position
        pos_id = id_to_pos[id];

        # Get the first nearest neighbors
        inds_neigh, dists_neigh = knn(tree_pbc,pos_id,4) # Get the 4 patches

        # Pass from index to id
        ids_neigh = map(s-> ind_to_id[s], inds_neigh); 

        # Get types
        types_neigh = map(s-> id_to_type[s], ids_neigh); 

        # Filter by type, only patches of monomers 
        mask_types = (types_neigh .== 4) .| (types_neigh .== 3);

        # Filter by distance of patches
        mask_distance = 0.6 .>= dists_neigh .> 0.3; 

        # Combine masks
        mask = mask_distance .& mask_types;

        # Update the array with the filters
        inds_neigh_filter = inds_neigh[mask];
        ids_neigh_filter = ids_neigh[mask];
        dists_neigh_filter = dists_neigh[mask];

        return (inds_neigh_filter, ids_neigh_filter)
end

"""
    get_cp(visited_id::Vector{Int64}, ind_to_id::Dict{Int64, Int64}, id_to_type::Dict{Int64, Int64}, pos_id::Vector{Float64}, tree_pbc)

Function that return the central particle neihgbors information 
"""
function get_cp(id_to_pos::Dict{Int64, Vector{Float64}}, ind_to_id::Dict{Int64, Int64}, id_to_type::Dict{Int64, Int64}, id::Int64, tree_pbc)
        # Get its position
        pos_id = id_to_pos[id];

        # Get the first nearest neighbors
        inds_neigh, dists_neigh = knn(tree_pbc,pos_id,3) # Get the 4 patches

        # Pass from index to id
        ids_neigh = map(s-> ind_to_id[s], inds_neigh); 

        # Get types
        types_neigh = map(s-> id_to_type[s], ids_neigh); 

        # Filter by type, only patches of crosslinkers
        mask_types = (types_neigh .== 1) .| (types_neigh .== 2)

        # Filter by distance of pathces
        mask_distance = 0.5 .>= dists_neigh .> 0.2; 

        # Combine masks
        mask = mask_distance .& mask_types;

        # Update the array with the filters
        inds_neigh_filter = inds_neigh[mask];
        ids_neigh_filter = ids_neigh[mask];
        dists_neigh_filter = dists_neigh[mask];

        return (inds_neigh_filter,ids_neigh_filter)
end

"""
    create_clusters(N_part::Int64, ids_central::Vector{Int64}, id_to_pos::Dict{Int64, Vector{Float64}}, ind_to_id::Dict{Int64, Int64}, id_to_type::Dict{Int64, Int64}, tree_pbc)

Function that return a graph with all connection between particles
"""
function create_clusters(N_part::Int64, ids_central::Vector{Int64}, tree_pbc)

    # Start the graph
    graph=SimpleGraph(N_part); # Based on index

    # Count the interaction that are between three patches.
    count_threebody=0;

    while !isempty(ids_central)
        # Start with one crosslinker
        id_central = popfirst!(ids_central);

    # cl -> patch
        # Get the index and ids of the patches of the crosslinker
        (inds_patch, ids_patch)=get_patches(id_to_pos,ind_to_id,id_to_type,id_central,tree_pbc)

        # Add bonds to the graph 
        foreach(s-> add_edge!(graph,id_to_ind[id_central], s), inds_patch);

        # Explore the patches of one central particle 
        for id_patch_explore in ids_patch
    # patch-patch
            # Get the index and ids  of the neighbors of the patch
            (inds_neigh, ids_neigh)=get_patch_patch(id_to_pos,ind_to_id,id_to_type,id_patch_explore,tree_pbc)

            if length(ids_neigh) > 1
                count_threebody += 1;
            end

            # Add bonds to the graph 
            foreach(s-> add_edge!(graph,id_to_ind[id_patch_explore], s), inds_neigh);

            # explore the patches to find central particles 
            for id_patch in ids_neigh
    # patch - central
                (ind_neigh, id_neigh)=get_cp(id_to_pos,ind_to_id,id_to_type,id_patch,tree_pbc)

                if length(id_neigh) == 0
                    ##println("Chain end")
                    continue
                elseif length(id_neigh) > 1
                    #println("Two central particles for one patch.")
                    ind_neigh = first(ind_neigh);
                    id_neigh = first(id_neigh);
                else
                    ind_neigh = first(ind_neigh);
                    id_neigh = first(id_neigh);
                end # if
           
                # Add the bond
                add_edge!(graph,id_to_ind[id_patch], ind_neigh)

            end # for patch - central
        end # for patch - patch 
    end # while central

    return (graph, div(count_threebody,3))

end

"""
    get_CL_cistances_euclidean(list_inds_clusters::Vector{Vector{Int64}}, l_x::Float64, l_y::Float64, l_z::Float64)

Compute the arithmetic distances between crosslinkers inside a cluster
"""
function get_CL_cistances_euclidean(list_inds_clusters::Vector{Vector{Int64}}, l_x::Float64, l_y::Float64, l_z::Float64)

        # Distance between crosslinkers
        distances = Array{Float64,1}();

        for cluster in list_inds_clusters
            # Create mask for central particles
            cluster_type = map(s->id_to_type[ind_to_id[s]],cluster)
            mask_type = cluster_type .== 1;

            # Filter thru the mask
            cluster_type = cluster[mask_type];

            while !isempty(cluster_type)

                # Select one particle and delete it from the list
                ind_i = popfirst!(cluster_type);

                # Get the position
                pos_i = id_to_pos[ind_to_id[ind_i]];
                pos = map(s->id_to_pos[ind_to_id[s]],cluster_type);

                # Compute the distances with the rest
                dist = map(s->dist_pbc.(pos_i.-s,[l_x; l_y; l_z]),pos);

                # Append the distances
                append!(distances,map(s-> sqrt(s'*s),dist));
            end # while 
        end # for

    return distances
end

"""
    create_hist_CL_distances(distances::Vector{Float64})

Create a histogram with the distances between cl
"""
function create_hist_CL_distances(distances::Vector{Float64})

    # Create a histogram of the distances between CL
        
        # Define a bin size
        bin_size = 1.4;

        # Maximum number of bins
        N_bins = ceil(Int64,ceil(div(maximum(distances),bin_size))) + 1;

        # Create the bins
        bins = bin_size.*(0:N_bins);

        # Compute bins id
        bins_id = floor.(Int64,distances./bin_size).+1;

        # Create the histogram
        hist_distances=map(s->count(==(s),bins_id),eachindex(bins));

        return (domain=bins,range=hist_distances)

end

"""
    explore_chain_cl(chain::Vector{Int64}, graph::SimpleGraph{Int64}, id_to_type::Dict{Int64, Int64}, ind_to_id::Dict{Int64, Int64})

Function that explores a chain until it finds another CL 
"""
function explore_chain_cl(chain::Vector{Int64}, graph::SimpleGraph{Int64}, id_to_type::Dict{Int64, Int64}, ind_to_id::Dict{Int64, Int64})

            while true 
    
                # Select the last item of the chain.
                # MUST BE A PATCH
                patch_ind = last(chain);

    # patch - patch
                # Find neighbors of the patch
                patch_patch_inds = all_neighbors(graph,patch_ind);

                # Filter from central particles
                # Map from inds to id to type
                patch_patch_type = map(s-> id_to_type[ind_to_id[s]],patch_patch_inds);

                # Create a mask
                mask_type1 = patch_patch_type .!= 1;
                mask_type2 = patch_patch_type .!= 2;
                mask_type = mask_type1 .& mask_type2;

                # Apply the mask
                patch_patch_inds = patch_patch_inds[mask_type];

                if length(patch_patch_inds) < 1
                    break # Chain end
                end

                # Select one patch
                #  HERE I DONT TAKE INTO ACCOUTN THREEBODY STUFF
                patch2_ind = patch_patch_inds[1];

                # Eliminate particles in the chain
                mask_repeat = mapreduce(s->patch2_ind .!= s,.&,chain);

                # Apply the mask
                patch2_ind = patch2_ind[mask_repeat];

                # Add the patch to the chain
                append!(chain,patch2_ind)
                    
    # patch - central
                # Find neighbors of the patch (central particle)
                patch_central_ind = all_neighbors(graph,last(chain));

                # Filter from patches 
                # Map from inds to id to type
                patch_central_type = map(s-> id_to_type[ind_to_id[s]],patch_central_ind);

                # Create a mask
                mask_type3 = patch_central_type .!= 3;
                mask_type4 = patch_central_type .!= 4;
                mask_type = mask_type3 .& mask_type4;

                # Apply the mask
                patch_central_ind = first(patch_central_ind[mask_type]);

                # Check if it is a Monomer or CL
                patch_central_type = first(patch_central_type[mask_type]);

                if patch_central_type == 1
                    #println("Is CL")
                    
                    # Eliminate particles in the chain
                    mask_repeat = mapreduce(s->patch_central_ind .!= s,.&,chain);

                    if mask_repeat == true 
                        
                        # Apply the mask
                        patch_central_ind = patch_central_ind[mask_repeat];

                        append!(chain,patch_central_ind)

                        break
                    
                    elseif mask_repeat == false
                        # No more CL
                        #println("No more CL")
                        break
                    end

                else
                    # Eliminate particles in the chain
                    mask_repeat = mapreduce(s->patch_central_ind .!= s,.&,chain);

                    # Apply the mask
                    patch_central_ind = patch_central_ind[mask_repeat];

                    append!(chain,patch_central_ind)
                end # if

    # central - patch
                # Select one central particle 
                central_ind = last(chain);

                # get the inds neighbors (patches)
                patches_inds = all_neighbors(graph,central_ind);

                # Filter from central particles
                # Map from inds to id to type
                patch_patch_type = map(s-> id_to_type[ind_to_id[s]],patches_inds);

                # Create a mask
                mask_type1 = patch_patch_type .!= 1;
                mask_type2 = patch_patch_type .!= 2;
                mask_type = mask_type1 .& mask_type2;

                # Apply the mask
                patches_inds = patches_inds[mask_type];

                # Eliminate particles in the chain
                mask_repeat = mapreduce(s->patches_inds .!= s,.&,chain);

                # Apply the mask
                patches_inds = patches_inds[mask_repeat];

                # Add the patch to the chain
                append!(chain,first(patches_inds))

            end # while

        return chain

end


"""
    compute_cl_cl_distances_same_chain(all_chains)
Compute the distances between crosslinkers in the same chain
"""
function compute_cl_cl_distances_same_chain(all_chains)
            # Reduce the array
            all_chains = reduce(vcat,all_chains);

            # Eliminate duplicates
            all_chains = unique_undirected(all_chains)

            # To store the distances between CL
            distances_cl_cl = Array{Float64,1}();

            # Iterate in all chains with the cl
            for chain_to_analyze in all_chains
                first_particle_id = ind_to_id[first(chain_to_analyze)];
                last_particle_id = ind_to_id[last(chain_to_analyze)];

                # check if it finishes with a CL
                if 1 == id_to_type[last_particle_id]
                    #println("Ended with a CL")

                    # Get the ids of the chain
                    chain_ids = map(s->ind_to_id[s],chain_to_analyze);

                    # Get the position of the chains
                    chain_positions = map(s->id_to_pos[s],chain_ids);

                    # Compute the distance between particles in the chain
                    chain_displacements = chain_positions[2:end] .- chain_positions[1:end-1];

                    # Get the distances
                    chain_distances = map(s->sqrt(s'*s),chain_displacements);

                    # Get the euclidean distance
                    euclidean_distance = id_to_pos[last_particle_id] .- id_to_pos[first_particle_id];

                    euclidean_distance = sqrt(euclidean_distance'*euclidean_distance);
                    
                    # Distance between CL
                    append!(distances_cl_cl,sum(chain_distances));
                end
            end

    return distances_cl_cl
end

"""
    unique_undirected(arrays)

Done by DeepSeek.
Para cada vector v:

Calcula rev = reverse(v).

Compara v y rev con < (orden lexicográfico, igual que comparar listas).

La canónica es el menor de los dos (canon = v < rev ? v : rev).

Convierte esa canónica a tupla (para que sea hashable y usable en un Set).

Si la tupla no está en seen, la agregamos y guardamos el vector original en el resultado.
"""
function unique_undirected(arrays)
    seen = Set{Tuple}()                     # guarda tuplas canónicas
    result = Vector{Vector{Int}}()          # resultado final

    for v in arrays
        # Asegurar que v es un Vector{Int} (por si es Any)
        v_int = Vector{Int}(v)
        rev = reverse(v_int)                # reverso del vector
        # Canon = el menor entre v y su reverso (orden lexicográfico)
        canon = v_int < rev ? v_int : rev
        canon_tuple = Tuple(canon)          # convertir a tupla para usar en Set

        if !(canon_tuple in seen)           # si no lo hemos visto antes
            push!(seen, canon_tuple)
            push!(result, v_int)            # guardamos el vector original
        end
    end
    return result
end


"""
    get_chains_cl_start(cluster_inds::Vector{Int64})

Get all the chains in a cluster starting from a crosslinker
"""
function get_chains_cl_start(cluster_inds::Vector{Int64})
        # Create mask for CrossLinkers 
        cluster_type = map(s->id_to_type[ind_to_id[s]],cluster_inds)
        mask_type = cluster_type .== 1;

        # Filter thru the mask
        cl_cluster = cluster_inds[mask_type];

        if length(cl_cluster) == 1
            #println("Only one Cl in the cluster")
            # continue or break
            return []
        else
            # To store all the chains
            all_chains = [[] for _ in eachindex(cl_cluster)]; 

            for (it,cl_explore) in enumerate(cl_cluster)

                # Explore a chain
                # Until it reaches a CL or
                # It ends

                # To store the 4 chains per CL
                chains_cl = [Vector{Int64}() for _ in 1:4];

                # Find the patches of the initial CL
                patches_inds = all_neighbors(graph,cl_explore);

                # Explore the 4 chains of one crosslinker
                for (it2,patch_ind) in enumerate(patches_inds)
                
                    # Start the chain
                    # IMPORTANT: CENTRAL THEN PATCH
                    chain = [cl_explore; patch_ind];

                    # Get the chain starting from one patch
                    chains_cl[it2] = explore_chain_cl(chain,graph,id_to_type,ind_to_id)

                end #for 

                # Store the chains
                all_chains[it] = chains_cl;

            end # for

            return all_chains
        end

end

"""
    get_cl_cl_distances(list_inds_clusters::Vector{Vector{Int64}})

Get the distances betwen cl going throuht the chain. For all chains in the cluster
"""
function get_cl_cl_distances(list_inds_clusters::Vector{Vector{Int64}})
        # Store the distances between cl following the chain
        distances_cl_cl = Array{Float64,1}();

        for cluster in list_inds_clusters

            # Going thru the monomer chain
            all_chains = get_chains_cl_start(cluster);

            # Check if there are more than one cl in the cluster
            if isempty(all_chains)
                # Implies only one cl in the cluster
                continue
            else 

                # Compute the distances between CL
                aux = compute_cl_cl_distances_same_chain(all_chains);

                # Store the results
                append!(distances_cl_cl,aux)
            end
        end

    return distances_cl_cl

end

"""
    quantify_loops(graph::SimpleGraph{Int64})
    
Compute the total loops and stuff
"""
function quantify_loops(graph::SimpleGraph{Int64})
        # Get the loops in the box
        loops_entire_box = cycle_basis(graph);

        # Get the types of the particles in the loops
        loops_box_type = [map(s->id_to_type[ind_to_id[s]],loops) for loops in loops_entire_box];

        # Count the size of the loops
        size_loops_type = length.(loops_box_type);

        # Mask of threebody
        mask_threebody_type = size_loops_type .== 3;

        # Apply the mask
        threebody_loops_type = loops_box_type[mask_threebody_type];

        # Mask for only patches
        patches_loops_type = [map(s-> s==4 || s==3,loops) for loops in threebody_loops_type];

        # Make sure that all of them are patches
        test_patch = first(unique(sum.(patches_loops_type)));

        if test_patch != 3
            @warn "A central particle has a threebody interaction."
        end
    
        count_threebody_loops = length(threebody_loops_type);

        # Loops between central particles
        loops_real_ind = loops_entire_box[.!mask_threebody_type];

        # Loops between central particles
        loops_real_type = loops_box_type[.!mask_threebody_type]; 

        # Amoun of loops
        N_loops = length(loops_real_type);

        # Amount of central particles in each loop
        N_size_real_loop = [sum((loop .== 1) .|| (loop .== 2)) for loop in loops_real_type];

        return (count_threebody_loops,N_loops,N_size_real_loop,loops_real_ind)

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

        # Get the number of particles to analyse
        N_part = nrow(df_dump);

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
        foreach(s->ind_to_id[s] = df_dump.id[s], 1:N_part);

        # Dictionary from ind to id
        id_to_ind = Dict(zip(values(ind_to_id), keys(ind_to_id)));

        # Dictionary from id to type
        id_to_type = Dict(zip(df_dump.id,df_dump.type));

        # Dictionary from positions to id
        pos_to_id = Dict{Vector{Float64}, Int64}();

        # Create the dictirionary
        foreach(s->pos_to_id[positions[:,s]] = df_dump.id[s], 1:N_part);

        # Dictionary from id to positions
        id_to_pos = Dict(zip(values(pos_to_id), keys(pos_to_id)));

        # Create a tree NearestNeighbors and distances
        tree=KDTree(positions)

        # Consider periodic boundary conditions
        tree_pbc=PeriodicTree(tree,[-l, -l, -l],[l, l, l])

        # Group by particle type
        df_by_types=groupby(df_dump,:type)

        # The initial nodes are going to be all the crosslinkers
        ids_cl = df_by_types[1].id[:];
        ids_mo = df_by_types[2].id[:];
        ids_pc = df_by_types[3].id[:];
        ids_pm = df_by_types[4].id[:];

        # To store ids to explore 
        ids_central=deepcopy([ids_cl; ids_mo]);    
#        ids_patches=Array{Int64,1}();

        # Create a graph with the position of the particles and cutoff distances of the potentials
        (graph, count_threebody) = create_clusters(N_part,ids_central,tree_pbc)


# Analysis of the graph

        # Create a list with the inds of particles in a cluster
        list_inds_clusters=connected_components(graph);

    # Quantify the amount of clusters 

        # Create a mask that select clusters bigger than one particle 
        mask_cluster = length.(list_inds_clusters) .> 3;

        list_inds_clusters=list_inds_clusters[mask_cluster];

        # Number of clusters in the system
        N_clusters = length(list_inds_clusters);

        # Biggest cluster
        Max_cluster = maximum(length.(list_inds_clusters));

        # Get the euclidean distance between CL
        euclidean_cl_cl = get_CL_cistances_euclidean(list_inds_clusters,l_x,l_y,l_z)

        # Create a histogram with the euclidean distances
        hist_dist_euclidean = create_hist_CL_distances(euclidean_cl_cl);
      
        distances_cl_cl = get_cl_cl_distances(list_inds_clusters);
        
        # Create a histogram with the euclidean distances
        hist_dist_chain = create_hist_CL_distances(distances_cl_cl);
 
    # Quantify loops and threebody indetractions
       (count_threebody_loops,N_loops,N_size_real_loop,loops_real_ind) = quantify_loops(graph) 




    # Quantify dangling ends
 

        # Count threebody interactions
#        unique_types_loops = unique.(loops_box_type);

        # Mask of types of particles
        #mask_n_types = length.(unique_types_loops).<=2


        #miau = 


        


#=
        it_cluster = 1;

        chains_in_cluster = get_chains_cl_start(list_inds_clusters[it_cluster]);

        all_chains = reduce(vcat,chains_in_cluster)

        # Eliminate duplicates
        all_chains = unique_undirected(all_chains)

        loops = cycle_basis(graph);
=#

     #end # For enumerate(paths_dumpf_simulations)




