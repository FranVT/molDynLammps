"""
    Basic analysis script
"""

using DataFrames, CSV
using Statistics, StatsBase
using Distances, Graphs
using GraphMakie

# https://juliagraphs.org/Graphs.jl/v1.5/
# https://graphsjl-docs.readthedocs.io/en/latest/

# Load the functions
include("functions.jl")
include("graphs-functions.jl") # Includes the graphical packages

# Parameter to select the system
T=0.05;
N_particles=5000;
phi=0.15;
CL_con=0.05;

# Selection of an specific simulation
# 0.050.150.055000-2026-03-06-175535
date="2026-03-06-175535";
#"2026-03-06-112337";
#"2026-03-05-124121";

#"2026-01-23-150058";
#"2026-01-22-165623";
#"2026-01-22-163251";
#"2026-01-22-154929";
#"2026-01-21-161335";
#"2026-01-21-133721";
#"2026-01-20-151755";
#"2026-01-20-143651";
#"2026-01-20-135923";
#"2026-01-20-121005"; 
#"2026-01-20-111704";

#0.050.30.05500-2026-01-20-121005


# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=DIR[1];

# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
data=extractFixScalar(DIR,FILE_NAME);

# Convert the array into a DataFrame
DATA=DataFrame(data[2]',data[1]);

"""
    P L O T S 
"""

damp=1;
dt=0.001;

# Create and save the graphics
fig_Temp(dt,DATA,T,DIR,id_c);
fig_Eng(dt,DATA,DIR,id_c);
fig_EngB(dt,DATA,DIR,id_c);
fig_EngSys(dt,DATA,DIR,id_c);
fig_EngPair(dt,DATA,DIR,id_c)


# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=joinpath(DIR[1],"traj");

# Filename with the simulation data
FILE_NAME="traj_assembly.9000000.dumpf";

# Get data
data=getDump(DIR,FILE_NAME);

clusters=map(collect(groupby(data,:c_clusters,sort=true))) do s
    s[(s.type.==1.0).|(s.type.==2.0),:]
end



# Start creating the list of neighbors to analuze distances and stuff
#pos=[clusters[1].x clusters[1].y clusters[1].z]

#sum(map(s->evaluate(Euclidean(),pos[1,:],pos[s,:]),2:length(pos[:,1])).<=1.2)


# Agregar los vecinos
n_clu=2;
df=clusters[n_clu];
n = nrow(df);
ids_mol=df.mol; # Id of the molecule

clusters_patchs=map(collect(groupby(data,:c_clusters,sort=true))) do s
    s[(s.type.==3.0).|(s.type.==4.0),:]
end

# Agregar los vecinos
df_patchs=clusters_patchs[n_clu];
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
        L=2*26.592409;
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

# Reducir el dataframe a id de molecula
redPatch=combine(groupby(df_patchs,:mol),
                 :neigh => (s->[vcat(s...)]) => :neigh,
                 :grado => sum => :grado,
                 :inconsistente => maximum => :inconsistente );

# Pasar la información de vecinos al Dataframe con partículas centrales
df = leftjoin(df, redPatch, on = :mol)

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

#draw(PNG("grafo.png", 10cm, 10cm), gplot(g, nodefillc=nodecolor))

# Visualizer el cluster en grafo
graphplot(g, node_color=nodecolor, node_size=15, edge_width=1)


if is_cyclic(g)
    # Regresa los ciclos con los ids
    cylces=cycle_basis(g)
end

function cadenas_desde(g, inicio, tipo)
    # inicio es índice del CL
    cadenas = []
    for nb in neighbors(g, inicio)
        path = [inicio, nb]
        while tipo_vertices[last(path)] == 0   # mientras sea MO
            siguientes = setdiff(neighbors(g, last(path)), [path[end-1]])
            isempty(siguientes) && break
            push!(path, only(siguientes))
        end
        push!(cadenas, path)
    end
    return cadenas
end

cad1=cadenas_desde(g, 4, 1)


# Crear una lista de vecinos

#=
# Alocar memoria
df.neigh=[Float64[] for _ in 1:n];

for it1 in 1:n
    # Get the reference particle
    pos1=[df.x[it1], df.y[it1], df.z[it1]];
    
    # Compute the distacnes with other particles and determine the neighbors
    for it2 in 1:n
        it1 == it2 && continue
        # Compute the distances
        pos2=[df.x[it2], df.y[it2], df.z[it2]];
        L=2*26.592409;
        box_size=[L,L,L];
        dist=evaluate(PeriodicEuclidean(box_size),pos1,pos2);

        # Classify as a neighbor or not
        if dist <= 1.4 && dist >= 1.0 
            push!(df.neigh[it1],ids_mol[it2])
        end
    end
end

# Modificar la lista de vecinos para asegurar simetría
# Evitar grafos dirigidos
neigh_dict = Dict(df.mol[i] => Set(df.neigh[i]) for i in 1:n)

# Agregar vecinos faltantes (i → j implica j → i)
for i in 1:n
    id_i = df.mol[i]
    for j in df.neigh[i]
        push!(neigh_dict[j], id_i)
    end
end

# Actualizar df.neigh con listas simétricas 
for i in 1:n
    df.neigh[i] = collect(neigh_dict[df.mol[i]])
end

# Guardar los grados de cada partícula
df.grado=length.(df.neigh)


# Checar concistencia (Ver problemas con potencial de 3 cuerpos)
df.inconsistente = [if type == 1.0
                        grado > 4 ? 1 : 0
                    elseif type == 2.0
                        grado > 2 ? 1 : 0
                    else
                        1  # tipo desconocido se marca como inconsistente
                    end for (type, grado) in zip(df.type, df.grado)]
=#




#=

# Clasificar los grados para poner los CL como "nodos especiales"
especial = Dict{Float64, Bool}()
for id in keys(adyacencia)
    if grados[id] == 1
        especial[id] = true   # extremo
    elseif df[df.id .== id,:].type[1] == 1.0 && grados[id] >= 3
        especial[id] = true   # CL 
    else
        especial[id] = false  # lineal (grado 2, o 4-patch con grado ≤2)
    end
end


# Conjunto de nodos lineales (grado 2 y no especiales)
lineales=Set(id for (id,e) in especial if !e && grados[id] == 2.0);



# Creación del grafo compelto
visitados=Set{Float64}();
cadenas=Vector{Vector{Float64}}();

# Creación del grafo simplificado
caminos = Vector{Vector{Int}}()
visitados_lineas = Set{Int}()  # nodos lineales ya procesados




# Explorar una línea desde un nodo inicial en una dirección
function explorar_linea(inicio, direccion_primer_paso)
    linea = Float64[inicio]
    actual = inicio
    siguiente = direccion_primer_paso
    while true
        push!(linea, siguiente)
        visitado_actual = siguiente
        # Obtener vecinos del siguiente
        vecinos_sig = adyacencia[siguiente]
        # Excluir el nodo anterior para no retroceder
        candidatos = filter(v -> v != actual, vecinos_sig)
        if length(candidatos) == 1
            # Continuamos en línea recta
            actual, siguiente = siguiente, candidatos[1]
        else
            # Terminamos porque encontramos una bifurcación o final
            break
        end
    end
    return linea
end

# -----------------------------------------------------------------------------------------
function recorrer_linea(inicio_especial, primer_lineal)
    camino = [inicio_especial, primer_lineal]
    actual = primer_lineal
    anterior = inicio_especial
    while true
        vecinos = adyacencia[actual]
        # El nodo lineal tiene grado 2, así que hay un solo vecino diferente al anterior
        siguiente = filter(v -> v != anterior, vecinos)[1]
        if especial[siguiente]
            push!(camino, siguiente)
            break
        elseif siguiente in lineales
            push!(camino, siguiente)
            push!(visitados_lineas, siguiente)
            anterior, actual = actual, siguiente
        else
            # No debería pasar, pero por si acaso
            break
        end
    end
    return camino
end


# Procesar cada nodo especial
for s in sort(collect(keys(adyacencia)))  # orden para evitar duplicados
    if especial[s]
        for v in adyacencia[s]
            if especial[v]
                # Conexión directa entre dos especiales
                if s < v
                    push!(caminos, [s, v])
                end
            elseif v in lineales && !(v in visitados_lineas)
                # Explorar línea desde s a través de v
                camino = recorrer_linea(s, v)
                t = last(camino)
                if s < t
                    push!(caminos, camino)
                end
                # Los nodos lineales ya se marcaron en recorrer_linea
            end
        end
    end
end


function visualizar_caminos(df, caminos, especial)
    coords = Dict(row.id => (row.x, row.y, row.z) for row in eachrow(df))
   
    fig=Figure()
    ax=Axis3(fig[1,1],xlabel="x", ylabel="y", zlabel="z")
    colores = cgrad(:tab10, length(caminos); categorical=true)
    
    for (i, camino) in enumerate(caminos)
        puntos = [coords[id] for id in camino]
        xs = [p[1] for p in puntos]
        ys = [p[2] for p in puntos]
        zs = [p[3] for p in puntos]
        
        # Dibujar línea
        lines!(ax,xs, ys, zs, color=colores[i], linewidth=3)
        
        # Nodos: tamaño según sean especiales o no
        sizes = [especial[id] ? 8 : 4 for id in camino]
        scatter!(ax,xs, ys, zs, color=colores[i], markersize=8)
    end
    
    display(fig)
end

visualizar_caminos(df, caminos, especial)

# ----------------------------------------------------------------------------------------------

# Primero, buscar nodos con grado != 2 como puntos de partida (extremos o bifurcaciones)
for (id, grado) in grados
    if grado != 2 && !(id in visitados)
        # Este nodo puede ser inicio de varias líneas (una por cada vecino)
        for vecino in adyacencia[id]
            if !(vecino in visitados)  # para no repetir
                linea = explorar_linea(id, vecino)
                # Marcar todos los nodos de la línea como visitados
                union!(visitados, linea)
                push!(cadenas, linea)
            end
        end
    end
end

# Ahora, buscar posibles ciclos (nodos con grado 2 no visitados)
for (id, grado) in grados
    if grado == 2 && !(id in visitados)
        # Iniciamos un ciclo desde este nodo
        # Elegimos un vecino cualquiera para comenzar
        vecino1 = adyacencia[id][1]
        linea = explorar_linea(id, vecino1)
        # Si la línea termina porque se encuentra un nodo con grado !=2, no es un ciclo puro.
        # Pero si todos los nodos en la línea tienen grado 2, entonces debe cerrarse.
        # En explorar_linea, la condición de parada es cuando el siguiente nodo tiene más de un vecino (sin contar el anterior)
        # En un ciclo puro, al final deberíamos volver al inicio. Comprobemos:
        ultimo = last(linea)
        vecinos_ultimo = adyacencia[ultimo]
        # Si el último tiene grado 2 y uno de sus vecinos (distinto del anterior) es el primer nodo, entonces es ciclo.
        if length(vecinos_ultimo) == 2
            anterior = linea[end-1]
            posible_retorno = filter(v -> v != anterior, vecinos_ultimo)[1]
            if posible_retorno == linea[1]
                # Es un ciclo: añadimos el primer nodo al final para indicar cierre
                push!(linea, linea[1])
            end
        end
        union!(visitados, linea)
        push!(cadenas, linea)
    end
end


# Visualización de las cadenas anteriores
    coords = Dict(row.id => (row.x, row.y, row.z) for row in eachrow(df))
    
    fig = Figure()
    ax = Axis3(fig[1,1])
    
    colores = cgrad(:tab20, length(cadenas); categorical=true)
    
    for (i, cadena) in enumerate(cadenas)
        puntos = [coords[id] for id in cadena]
        xs = [p[1] for p in puntos]
        ys = [p[2] for p in puntos]
        zs = [p[3] for p in puntos]
        
        # Línea
        lines!(ax, xs, ys, zs, color=colores[i], linewidth=3)
        # Nodos
        scatter!(ax, xs, ys, zs, color=colores[i], markersize=8)
    end
    
    fig
=#
# Simplificar las cadenas/red
















#=

# The amount of CL and MO
nrow.(typeCl1|>collect)





patches=data[(data.type.==4.0).|(data.type.==3.0),:];

# Count all same id of cl_clusters (Amount of clusters)
N_patchclusters=values(countmap(patches.c_clusters));

# Count the size of each cluster
cl_patches=countmap(N_patchclusters);
=#

