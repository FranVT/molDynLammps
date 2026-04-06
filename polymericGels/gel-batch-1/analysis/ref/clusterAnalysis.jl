"""
    Vibe codeado con DEepseek
"""

using DataFrames
using Plots

"""
    simplificar_grafo_patchy(df::DataFrame)

Devuelve un vector de caminos (vectores de IDs) que representan la simplificación del grafo
de partículas patchy, donde los nodos de tipo "4patch" con grado ≥ 3 son considerados centrales
y los nodos con grado 1 son extremos. Los caminos conectan estos nodos especiales a través de
nodos lineales (grado 2 o 4patch con grado ≤ 2).

# Argumentos
- `df`: DataFrame con columnas:
    - `id`: identificador único del nodo
    - `x`, `y`, `z`: coordenadas
    - `tipo`: "2patch" o "4patch"
    - `neigh`: vector de IDs vecinos (debe ser simétrico, aunque se fuerza simetría)

# Retorna
- `caminos`: Vector{Vector{Int}} con los caminos simplificados.
- `especiales`: Diccionario que indica si un nodo es especial (true/false).
"""
function simplificar_grafo_patchy(df)
    # 1. Construir diccionario de adyacencia y grados
    ady = Dict{Int, Vector{Int}}()
    for row in eachrow(df)
        ady[row.id] = copy(row.neigh)  # copia para no modificar original
    end
    # Asegurar simetría (grafo no dirigido)
    for (u, vecinos) in ady
        for v in vecinos
            if !(u in get(ady, v, []))
                push!(ady[v], u)
            end
        end
    end

    grado = Dict(id => length(vecs) for (id, vecs) in ady)

    # 2. Clasificar nodos en especiales y lineales
    tipo = Dict(row.id => row.tipo for row in eachrow(df))
    especial = Dict{Int, Bool}()
    for id in keys(ady)
        if grado[id] == 1
            especial[id] = true   # extremo
        elseif tipo[id] == "4patch" && grado[id] >= 3
            especial[id] = true   # central
        else
            especial[id] = false  # lineal (incluye 4patch con grado ≤2 y 2patch con grado 2)
        end
    end

    # Conjunto de nodos lineales (para recorridos)
    lineales = Set(id for (id, e) in especial if !e)

    # 3. Almacenar caminos
    caminos = Vector{Vector{Int}}()
    visitados_lineas = Set{Int}()  # nodos lineales ya procesados

    # Función auxiliar: dado un nodo especial `inicio` y un vecino lineal `primer_lineal`,
    # recorre la línea hasta encontrar otro especial.
    function recorrer_hasta_especial(inicio, primer_lineal)
        camino = [inicio, primer_lineal]
        actual = primer_lineal
        anterior = inicio
        while true
            vecinos = ady[actual]
            # En un nodo lineal (grado 2) solo hay un vecino diferente al anterior
            siguientes = filter(v -> v != anterior, vecinos)
            if length(siguientes) != 1
                # Esto no debería ocurrir si el nodo es lineal (grado 2), pero por si acaso
                @warn "Nodo $actual con grado $(length(vecinos)) no es lineal puro" maxlog=1
                break
            end
            siguiente = siguientes[1]
            if especial[siguiente]
                push!(camino, siguiente)
                break
            elseif siguiente in lineales
                push!(camino, siguiente)
                push!(visitados_lineas, siguiente)
                anterior, actual = actual, siguiente
            else
                # No debería pasar
                break
            end
        end
        return camino
    end

    # Procesar cada nodo especial
    for s in sort(collect(keys(ady)))  # orden para consistencia
        if especial[s]
            for v in ady[s]
                if especial[v]
                    # Conexión directa entre dos especiales
                    if s < v  # evitar duplicados
                        push!(caminos, [s, v])
                    end
                elseif v in lineales && !(v in visitados_lineas)
                    # Explorar línea desde s a través de v
                    camino = recorrer_hasta_especial(s, v)
                    t = last(camino)
                    if s < t   # evitar duplicados (cada camino se guarda una vez)
                        push!(caminos, camino)
                    end
                    # Los nodos lineales ya se marcaron en la función
                end
            end
        end
    end

    # 4. Detectar ciclos de nodos lineales puros (todos grado 2, sin especiales)
    for id in keys(ady)
        if !especial[id] && grado[id] == 2 && !(id in visitados_lineas)
            # Iniciar un ciclo desde este nodo
            ciclo = Int[]
            actual = id
            anterior = 0
            # Elegir el primer vecino (cualquiera)
            vecinos = ady[actual]
            siguiente = vecinos[1]
            push!(ciclo, actual, siguiente)
            push!(visitados_lineas, actual, siguiente)
            anterior, actual = actual, siguiente
            while true
                vecinos_act = ady[actual]
                # El siguiente es el vecino diferente al anterior
                sig = filter(v -> v != anterior, vecinos_act)
                if sig == ciclo[1]  # cerramos el ciclo
                    break
                else
                    push!(ciclo, sig[1])
                    push!(visitados_lineas, sig[1])
                    anterior, actual = actual, sig[1]
                end
            end
            # Opcional: repetir el primer nodo al final para indicar cierre
            push!(ciclo, ciclo[1])
            push!(caminos, ciclo)
        end
    end

    return caminos, especial
end

"""
    visualizar_caminos(df, caminos, especiales; kwargs)

Dibuja los caminos en 3D. Los nodos especiales se muestran con marcador más grande.
"""
function visualizar_caminos(df, caminos, especiales; kwargs...)
    coords = Dict(row.id => (row.x, row.y, row.z) for row in eachrow(df))
    
    p = plot3d(legend=false, xlabel="x", ylabel="y", zlabel="z"; kwargs...)
    colores = palette(:tab10, length(caminos))
    
    for (i, camino) in enumerate(caminos)
        puntos = [coords[id] for id in camino]
        xs = [p[1] for p in puntos]
        ys = [p[2] for p in puntos]
        zs = [p[3] for p in puntos]
        
        # Línea
        plot3d!(xs, ys, zs, color=colores[i], linewidth=2, label="")
        
        # Nodos: tamaño según sean especiales o no
        sizes = [especiales[id] ? 8 : 4 for id in camino]
        scatter3d!(xs, ys, zs, color=colores[i], markersize=sizes, markerstrokewidth=0)
    end
    
    display(p)
    return p
end

# =====================
# EJEMPLO DE USO
# =====================
# Crear datos de prueba (simulando una red con partículas patchy)
# Aquí definimos un grafo pequeño: un central 4patch con grado 3, conectado a extremos y a una cadena lineal que forma un ciclo.
df_ejemplo1 = DataFrame(
    id = 1:9,
    x = [0.0, 1.0, 2.0, 1.0, 2.0, 3.0, 3.0, 4.0, 4.0],
    y = [0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0],
    z = zeros(9),
    tipo = ["2patch", "4patch", "2patch", "2patch", "2patch", "2patch", "2patch", "2patch", "2patch"],
    neigh = [
        [2],          # 1 (extremo)
        [1,3,4],      # 2 (central 4patch, grado 3)
        [2,5],        # 3 (lineal)
        [2,6],        # 4 (lineal)
        [3,7],        # 5 (lineal)
        [4,8],        # 6 (lineal)
        [5,9],        # 7 (lineal)
        [6,9],        # 8 (lineal)
        [7,8]         # 9 (lineal, junto con 7 y 8 forma un ciclo 5-7-9-8-6-4? mejor revisar: en realidad 7-9-8 es un triángulo? no, 7 conecta con 5 y 9; 8 conecta con 6 y 9; 9 conecta con 7 y 8. Esto forma un ciclo 7-9-8, pero 7 también conecta a 5, y 8 a 6, así que no es un ciclo puro de grado 2. Para un ciclo puro necesitamos que todos tengan grado 2. Mejor ajustamos para que el ciclo sea 7-9-8 y luego 7,8 no tengan otros vecinos. Corrijamos:
    ]
)

# Sobrescribimos para tener un ciclo lineal puro: nodos 7,8,9 forman un triángulo? pero grado 2 no puede ser triángulo, sería un ciclo de 3 nodos donde cada uno tiene grado 2, eso es posible si es un ciclo simple: 7-8-9-7. Pero entonces cada nodo tiene dos vecinos. Vamos a redefinir:
df_ejemplo2 = DataFrame(
    id = 1:8,
    x = [0,1,2,1,2,3,3,4],
    y = [0,0,0,1,1,0,1,0],
    z = zeros(8),
    tipo = ["2patch", "4patch", "2patch", "2patch", "2patch", "2patch", "2patch", "2patch"],
    neigh = [
        [2],          # 1
        [1,3,4],      # 2 (central)
        [2,5],        # 3
        [2,6],        # 4
        [3,7],        # 5
        [4,8],        # 6
        [5,8],        # 7 (conecta a 5 y 8)
        [6,7]         # 8 (conecta a 6 y 7) → esto crea un ciclo 5-7-8-6, pero 5 también conecta a 3, y 6 a 4, entonces no es puro. Para un ciclo puro necesitamos que 5 y 6 no tengan otras conexiones. Vamos a hacer un ciclo aparte con nodos 9,10,11.
    ]
)

# Mejor creamos dos DataFrames separados para pruebas, pero en el script final el usuario usará sus datos reales.
# Para no complicar, aquí va un ejemplo más simple con un ciclo lineal puro:
df_ejemplo3 = DataFrame(
    id = 1:7,
    x = [0,1,2,2,3,3,4],
    y = [0,0,0,1,0,1,0],
    z = zeros(7),
    tipo = ["2patch", "4patch", "2patch", "2patch", "2patch", "2patch", "2patch"],
    neigh = [
        [2],          # 1
        [1,3,4],      # 2 (central)
        [2,5],        # 3
        [2,6],        # 4
        [3,7],        # 5
        [4,7],        # 6
        [5,6]         # 7 (forma un ciclo 5-7-6 con 5 y 6 también conectados a 3 y 4 respectivamente, así que no es ciclo puro. Mejor rindamos)
    ]
)

