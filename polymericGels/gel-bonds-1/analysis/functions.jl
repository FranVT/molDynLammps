"""
    Functions 
"""

function getDatFile(date)
"""
    Return a dataframe with the parameters of the simulation of a given directory.
"""

    DIR=getDir(date);
    filename="dataAssembly.dat";
    aux=split.(readlines(joinpath(DIR,filename)),",");

    header=aux[1][5:end];
    info=parse.(Float64,aux[2][5:end]);

    return DataFrame(info',header)
end

function structureFactor(vec_r,vec_K)
"""
    Get the factor of structure of a set of points
"""

    # Dot product between the wave vector and the position vector
    dot_kr=map(s->dot(vec_K,vec_r[s,:]),eachindex(vec_r[:,1]));
    I_k=sum(exp.(im.*dot_kr))

   return I_k
end



function createNeighborList(cluster,cluster_patch,L)
"""
    Create analysis of the cluster for each dump
    cluster:    A data frame with the information of the cluster 
    L:          The length of one side of the simulation box.
"""
    
    # Data frame with the ccluster information
    df=cluster;
    n = nrow(df);
    ids_mol=df.mol; # Id of the molecule

    # Agregar los vecinos
    df_patchs=cluster_patch;
    n_patchs = nrow(df_patchs);
    ids_mol_patchs=df_patchs.mol; # Id of the molecule

    # Crear una lista de vecinos

    # Alocar memoria
    df_patchs.neigh=[Float64[] for _ in 1:n_patchs];

    for it1 in 1:n_patchs
        # Get the reference particle
        pos1=[df_patchs.x[it1], df_patchs.y[it1], df_patchs.z[it1]];
    
        # Compute the distacnes with other particles and determine the neighbors
        for it2 in 1:n_patchs
            it1 == it2 && continue
            # Compute the distances
            pos2=[df_patchs.x[it2], df_patchs.y[it2], df_patchs.z[it2]];
            box_size=[L,L,L];
            dist=evaluate(PeriodicEuclidean(box_size),pos1,pos2);

            # Classify as a neighbor or not
            if (dist <= 0.6) && (ids_mol_patchs[it1] != ids_mol_patchs[it2])
                push!(df_patchs.neigh[it1],ids_mol_patchs[it2])
            end
        end
    end

    # Guardar los grados de cada partícula
    df_patchs.grado=length.(df_patchs.neigh)

    # Checar concistencia (Ver problemas con potencial de 3 cuerpos)
    df_patchs.inconsistente = df_patchs.grado.>1

    if nrow(df_patchs) == 0 # Implica que no tiene vecinos
        # Crear redPatch con las mismas moléculas que df y valores por defecto
        df.neigh=[[] for _ in 1:nrow(df)];
        df.grado=zeros(Int, nrow(df));
        df.inconsistente=zeros(Int, nrow(df));

        #redPatch = DataFrame(mol = df.mol,
        #                 neigh = [[] for _ in 1:nrow(df)],
        #                 grado = zeros(Int, nrow(df)),
        #                 inconsistente = zeros(Int, nrow(df)))
    else
        # Agrupar y combinar como antes
        redPatch = combine(groupby(df_patchs, :mol),
                       :neigh => (s -> [vcat(s...)]) => :neigh,
                       :grado => sum => :grado,
                       :inconsistente => maximum => :inconsistente)
        # Pasar la información de vecinos al Dataframe con partículas centrales
        df = leftjoin(df, redPatch, on = :mol)
    end

        return df
end

function getClusters(data,L)
"""
    Get a dataframe with a list of neighbors for each cluster in the system
"""
    
    clusters=map(collect(groupby(data,:c_clusters,sort=true))) do s
        s[(s.type.==1.0).|(s.type.==2.0),:]
    end |> x -> filter(!isempty, x)

    # Create a dataframe with only the patches. (From this df the neoghbor list is constructed)
    clusters_patchs=map(collect(groupby(data,:c_clusters,sort=true))) do s
        s[(s.type.==3.0).|(s.type.==4.0),:]
    end |> x -> filter(!isempty, x)

    # Agregar los vecinos
    clusters = map(s->createNeighborList(clusters[s],clusters_patchs[s],L),eachindex(clusters));

    return clusters

end

function createGraph(df)
"""
    Create a graph from a Dataframe with list of neighbors
"""

    # Crear las cadenas
    ady = Dict{Int, Vector{Int}}()
    tipo = Dict{Int, Int}()
    for row in eachrow(df)
        id = Int(row.mol)
        tipo[id] = Int(row.type)
        ady[id] = [Int(n) for n in row.neigh]
    end


    # Mapeo ID → índice
    ids = sort(unique(Int.(df.mol)))
    map_id = Dict(id => i for (i, id) in enumerate(ids))
    n = length(ids)
    g = SimpleGraph(n)

    for row in eachrow(df)
        u = map_id[Int(row.mol)]
        for v in row.neigh
            add_edge!(g, u, map_id[Int(v)])
        end
    end

    # Guardar tipo como propiedad de vértice (opcional)
    tipo_vertices = [tipo[ids[i]] for i in 1:n]
    # Asignar colores según tipo (CL=1, MO=2)
    nodecolor = [tipo_vertices[i] == 1.0 ? "red" : "blue" for i in 1:n]

    return (g,nodecolor)
end



function getDump(dir,file_name)
"""
    Get the data from a single dump file that stores one timeste information
"""
    data = split.(readlines(joinpath(dir,file_name))," ")[9:end];
    HEADERS=data[1][3:end];
    INFO=parse.(Float64,reduce(hcat,data[2:end]))';

    return DataFrame(INFO,HEADERS)
end

function extractFixScalar(path_system,file_name)
"""
    Function that extracts the information of fix files that stores global scalar values
"""
    aux=split.(readlines(joinpath(path_system,file_name))," ");
    header=aux[2][2:end];
    info=reduce(hcat,map(s->parse.(Float64,s),aux[3:end]));

    return DataFrame(info',header)
end

function getDir(date)
"""
    To get the directory of the simulation
"""

# Get the directories
MAIN_DIR=dirname(pwd());
DATA_DIR=joinpath(MAIN_DIR,"data");

# Read the directory where data is stored
sims=filter(isdir,readdir(DATA_DIR,join=true));

# For this script, evetually we are going to get assembly averages
id_c=date; # filter!(s->s==id_c,readdir(DATA_DIR))

# Get the idex for the directory
indx=findall(!isempty,findall.(id_c,sims));

# Get the directory
DIR=sims[indx][1];

    return DIR

end


