"""
    Basic analysis script
"""

using DataFrames, CSV
using Statistics, StatsBase

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


#=
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
df=clusters[1];
n = nrow(df);
ids=df.id;
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
        dist=evaluate(Euclidean(),pos1,pos2);

        # Classify as a neighbor or not
        if dist <= 1.2
            push!(df.neigh[it1],ids[it2])
        end
    end
end

=#

# El df ahora tiene una columna con vecinos
# Ahora a creare el grafo

# Crear un diccionario de adyacencia y calcular grados
adyacencia=Dict{Float64,Vector{Float64}}();
for row in eachrow(df)
    adyacencia[row.id]=row.neigh
end

# Verificar simetría, ya que no es un grafo dirigido
for (id, vecinos) in adyacencia
    for v in vecinos
        if !(id in get(adyacencia, v, []))
            # Agregar vecinos para no hacerlo dirigido
            push!(adyacencia[v], id)
        end
    end
end

# Se calcula el grado de los nodos para facilitar creación del grafo
# (Cuanto vecinos tiene)
grados=Dict(id=>length(vecinos) for (id, vecinos) in adyacencia);

# Creación del grafo

visitados=Set{Float64}();
cadenas=Vector{Vector{Float64}}();

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


    coords = Dict(row.id => (row.x, row.y, row.z) for row in eachrow(df))
    
    fig = Figure()
    ax = Axis3(fig[1,1])
    
    colores = cgrad(:tab10, length(cadenas); categorical=true)
    
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

#=

# The amount of CL and MO
nrow.(typeCl1|>collect)





patches=data[(data.type.==4.0).|(data.type.==3.0),:];

# Count all same id of cl_clusters (Amount of clusters)
N_patchclusters=values(countmap(patches.c_clusters));

# Count the size of each cluster
cl_patches=countmap(N_patchclusters);
=#

