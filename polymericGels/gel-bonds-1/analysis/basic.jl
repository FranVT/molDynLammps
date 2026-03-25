"""
    Basic analysis script
"""

using DataFrames, CSV
using Statistics, StatsBase
using Distances, Graphs
using GraphMakie, LinearAlgebra
using BenchmarkTools

# https://juliagraphs.org/Graphs.jl/v1.5/
# https://graphsjl-docs.readthedocs.io/en/latest/

# Load the functions
include("functions.jl")
include("graphs-functions.jl") # Includes the graphical packages

# Selection of an specific simulation
date="2026-03-25-093253";
#"2026-03-24-163630";
#"2026-03-24-102830";
#"2026-03-20-112940";
#"2026-03-20-105611";


#"2026-03-19-130353"; # COMPUTE CLUSTER
#"2026-03-20-105611";
#"2026-03-19-130353";
path=getDir(date);
path_dump=joinpath(path,"traj");

# Read the data file
datConfig=getDatFile(path);

dt=datConfig."time-step"[1];

fix_analysis=false;

if fix_analysis == true
    # Extract the data from the file
    dataSystem=extractFixScalar(path,"system_assembly.fixf");

    # Create and save the graphics
    fig_Temp(dt,dataSystem,datConfig."Temperature"[1],path,date);
    fig_Eng(dt,dataSystem,path,date);
    fig_EngB(dt,dataSystem,path,date);
    fig_EngSys(dt,dataSystem,path,date);
    fig_EngPair(dt,dataSystem,path,date);

    println("Imagenes del sistema listas")
end


# Cluster analysis
id_files=(0:Int(datConfig."save-dump"[1]):Int(datConfig."N_heat"[1]+datConfig."N_isot"[1]));
file_names=string.("traj_assembly.",id_files,".dumpf");

L=2*datConfig.L[1]; # Length of the box

# Reducción de la cantidad de time steps a analizar
N = 2^8  # cantidad deseada de puntos
log_min = log(1.0);
log_max = log(length(file_names));
pts = unique(round.(Int, exp.(range(log_min, log_max, length=N))))
pts = pts[pts .<= 100001]  # asegura límite superior

println("Inicio de analisis de cluster")

#info_cluster=map(pts) do s
#    (N_clusters,cluster_Size)=clusterAnalysis(path_dump,file_names[s],L);   
#end;

function clusterAnalysis_opt(DIR,FILE_NAME,L)
"""
    Function that performs a quick cluster analysis.
    Return a dataframe with the position of the particles and neighbors and stuff
"""
    # Get one dataframe per cluster in the system with its neighbors
    clusters=getClusters_opt(DIR,FILE_NAME,L); # Just central particles. Patches have been descarted

    N_clusters=length(clusters);    # Amount of clusters in the system
    cluster_Size=nrow.(clusters);   # Amount of central particles in each cluster in the system

    # Get one graph for each cluster (strand length and loops)
    #graphs=map(s->createGraph(clusters[s]),eachindex(clusters))

    # Get the structure factor


    return (N_clusters,cluster_Size)

end 

function getClusters_opt(DIR,FILE_NAME,L)
"""
    Get a dataframe with a list of neighbors for each cluster in the system
    Se usa clusters_patchs para definir los vecinos: 
        Cuando dos patches con distinto id de mol están a una distancia menor de rc.
        Los vecinos son id mol.
    Los id mol de los vecinos de patches se transfiere a clusters, que son las partículas centrales.
"""

    df=getDump(DIR,FILE_NAME);

    clusters = view(df, findall((df.type .== 1.0) .| (df.type .== 2.0)), :);

    clusters_patchs = view(df, findall((df.type .== 3.0) .| (df.type .== 4.0)), :);

    # Agregar los vecinos
    #clusters = map(s->createNeighborList(clusters[s],clusters_patchs[s],L),eachindex(clusters));

    return (clusters,clusters_patchs)

end


#@btime getClusters_opt(path_dump,file_names[pts[89]],L); # 21.257 ms
#@btime getClusters(getDump(path_dump,file_names[pts[89]]),L); # 58.392 ms



tests=getClusters_opt(path_dump,file_names[pts[89]],L)
names(first(tests))


function createNeighborList_opt(df_central,df_patch,L)
"""
    Create analysis of the cluster for each dump
    cluster:    A data frame with the information of the cluster 
    L:          The length of one side of the simulation box.

    df: es el cluster de las partículas centrales.
"""
    
    n_central = nrow(df_central);
    mol_central=df_central.mol; # Id of the molecule

    # Extrear posiciones para optimizar cálculos
    x_patch=df_patch.x;
    y_patch=df_patch.y;
    z_patch=df_patch.z;
    mol_patch=df_patch.mol;
    n_patch = nrow(df_patch);
    
    # Preasignar memoria
    neigh_patch=[Float64[] for _ in 1:n_patch];

    # Falta agregar coso de calcular distancias por celdas

    println("Inicio ciclo 1")
    for it1 in 1:n_patch
        # Get the reference particle
        x_ref=x_patch[it1];
        y_ref=y_patch[it1];
        z_ref=z_patch[it1];
        mol_ref=mol_patch[it1];

        # Compute the distances with other particles and determine the neighbors
        for it2 in 1:n_patch
            it1 == it2 && continue # Avoid the same patch
            mol_comp=mol_patch[it2];
            mol_ref == mol_comp && continue # Avoid if the patches are in the same patchy particle

            
            # Compute distances considering periodic boundary contiions
            dx = abs(x_ref - x_patch[it2]);
            dx = min(dx, L - dx);
            
            dy = abs(y_ref - y_patch[it2]);
            dy = min(dy, L - dy);
            
            dz = abs(z_ref - z_patch[it2]);
            dz = min(dz, L - dz);
            
            dist2 = dx*dx + dy*dy + dz*dz;


            # Classify as a neighbor or not
            if dist2 <= 0.36 
                push!(neigh_patch[it1],mol_comp)
            end
        end
    end
    println("Fin ciclo 1")


    df_patch.neigh=neigh_patch;

    grado=length.(neigh_patch);

    # Guardar los grados de cada partícula
    df_patch.grado=grado;

    # Checar concistencia (Ver problemas con potencial de 3 cuerpos)
    df_patch.inconsistente = grado.>1

    if nrow(df_patch) == 0 # Implica que no tiene vecinos
        # Crear redPatch con las mismas moléculas que df y valores por defecto
        redPatch = DataFrame(mol = df_patch.mol,
                         neigh = [[] for _ in 1:nrow(df_patch)],
                         grado = zeros(Int, nrow(df_patch)),
                         inconsistente = zeros(Int, nrow(df_patch)))
    else
        # Agrupar y combinar como antes
        redPatch = combine(groupby(df_patch, :mol),
                       :neigh => (s -> [vcat(s...)]) => :neigh,
                       :grado => sum => :grado,
                       :inconsistente => maximum => :inconsistente)
    end

     # Pasar la información de vecinos al Dataframe con partículas centrales
     df = leftjoin(df_central, redPatch, on = :mol,makeunique=true)
    
    return df
end

#println("Obtención de vecinos")
@btime createNeighborList_opt(first(tests),last(tests),L) # 60 s

#@btime createNeighborList(first(tests),last(tests),L)

#println("Figuras de analisis de cluster")
#fig_NumClusters(dt,id_files[pts],first.(info_cluster),path,date)
#fig_MaxClusterPart(dt,id_files[pts],datConfig.Npart[1].-maximum.(last.(info_cluster)),path,date)


nothing



#=
# Get data
# Structure factor
l=1.5; # "wave length"
vec_K=(2*pi).*[1/l,1/l,1/l]; # vector
# Structure factor per cluster
I_k=map(s->structureFactor([clusters[s].x clusters[s].y clusters[s].z],vec_K),eachindex(clusters));

#Structure factor of all the system
I_K=sum(I_k);





=#


