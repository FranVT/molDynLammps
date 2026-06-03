"""
    Script en el que se adapta el código de FORTRAN de Claudia a Julia
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll
using SplitApplyCombine


include("functions.jl")

function createqdom(qmax,rq0,dq0,bin0,numbin)
"""
    Create the reciprocal space of wave vectors
"""
    # Crear el espacio del vector recíproco para no definirlo en cada experimento
    qx=(-qmax:qmax).*dq0;
    qy=(-qmax:qmax).*dq0;
    qz=(-qmax:qmax).*dq0;

    qhis = zeros(numbin, 5);        # histograma: col1 = suma de |S(q)|^2, col2 = conteos de la misma magnitud
    qhis .= 0.0;

    qxhis = [[] for _ in 1:numbin]#zeros(numbin, 2)
    qyhis = [[] for _ in 1:numbin]#zeros(numbin, 2)
    qzhis = [[] for _ in 1:numbin]#zeros(numbin, 2)
    qhis = [[] for _ in 1:numbin]

    # Calcular la magnitud cuadrada
    for x in qx
        for y in qy
            for z in qz

                # Excluir el origen del espacio recíproco
                if x == 0 && y == 0 && z == 0
                    continue
                end

                q=sqrt(x^2 + y^2 + z^2) 

                # Si la magnitud es grande que el maximo se salta
                if q > rq0
                    continue
                end

                # Determinar el índice del bin (manejo especial del borde)
                sbin = floor(Int, q / bin0) + 1

                # Coso para trnasformar de un dominio continup a uno discreto 
                if q % bin0 == 0.0
                    sbin -= 1
                end

                # Si la magnitud del vector supera la magnitud de interés se lo salta
                if sbin > numbin
                    continue
                end

                append!(qhis[sbin],q)
                append!(qxhis[sbin],x)
                append!(qyhis[sbin],y)
                append!(qzhis[sbin],z)

                #qhis[sbin,1]=x
                #qhis[sbin,2]=y
                #qhis[sbin,3]=z
                #qhis[sbin,4]=q
                #qhis[sbin,5]+= 1
            end
        end
    end

    return (qxhis,qyhis,qzhis,qhis) 
end



# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

dat_DF=data_bySystem[1];

    # Parameters of the system relevant to the structure factor
    L=2*dat_DF.L[1];    # Longitud de la caja
    N_instants=2;       # Instantes temporales a analizar
    N_exp=nrow(dat_DF); # Cantidad de experimentos por sistema
    N_part=dat_DF.Npart[1];

    # Create time steps range
    aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));

    # Get the index for the time instants that we are interested to analyzed
    ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));

    # Get the time steps
    timeSteps=aux_timeStep[ind];

    # Create file names
    file_names=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in timeSteps];

    # Path to the dumps
    dump_paths=joinpath.(dat_DF.PARENT_DIR,dat_DF.dir,"traj");

# SELECT ONE EXPERIMENT
r_system = getPositions(file_names[2], dump_paths);


r = r_system[1]

# System information
lxn=L;
lyn=L;
lzn=L;
sigma=1;
pi2=2*pi;
ntot=Int(N_part); # Número total de partículas

# Definition of parameters

# --- Parámetros de la caja y de la rejilla en q ---
xc = lxn      # longitud de la caja en x
yc = lyn      # longitud de la caja en y
zc = lzn      # longitud de la caja en z
s = sigma     # factor de escala (sigma)
x2 = xc / 2.0 # mitad de la caja en x (no usada después, posible para centrar)
y2 = yc / 2.0 # mitad de la caja en y
z2 = zc / 2.0 # mitad de la caja en z
bin0 = 1;
qmax0 = 10;
ball=0;     # No se para que es esta variable


# Espaciado mínimo en el espacio recíproco y ancho de bin escalado
dq0 = pi2 / xc               # Δq fundamental
bin0 = dq0 * bin0            # nuevo ancho de bin (bin0 original * dq0)
qmax = Int(floor(qmax0 / dq0)) # número entero de pasos hasta qmax0
rq0 = qmax * dq0             # valor máximo real de |q| usado
numbin = Int(floor(qmax * dq0 / bin0)) + 1  # número total de bines

# Reserva de memoria
r = reduce(hcat,r)              # posiciones (columnas 1:3)
sfhis = zeros(numbin, 2)        # histograma: col1 = suma de |S(q)|^2, col2 = conteos de la misma magnitud
kr = zeros(ntot)                # producto escalar q·r para cada partícula
ntotav = ntot                   # Número total de particulas


# Aviso si alguna partícula está exactamente en el origen
for i in 1:ntot
    if r[i,1] == 0.0 && r[i,2] == 0.0 && r[i,3] == 0.0
        println("No data")
    end
end

## --- Cálculo del factor de estructura ---
sfhis .= 0.0 # Se inicializa la variable en la que se guardara el factor de estructura
for iq in -qmax:qmax
    qr1 = iq * dq0          # componente x del vector q
    for jq in -qmax:qmax
        qr2 = jq * dq0      # componente y
        for kq in -qmax:qmax
            
            # Excluir el origen del espacio recíproco
            if iq == 0 && jq == 0 && kq == 0
                continue
            end
            
            # Componente z del vector q
            qr3 = kq * dq0                     # componente z
           
            qr = [qr1, qr2, qr3]               # vector q
            rq = sqrt(qr1^2 + qr2^2 + qr3^2)   # magnitud |q|
            if rq > rq0
                continue
            end

            # Determinar el índice del bin (manejo especial del borde)
            sbin = floor(Int, rq / bin0) + 1

            # Coso para trnasformar de un dominio continup a uno discreto 
            if rq % bin0 == 0.0
                sbin -= 1
            end
           
            # Si la magnitud del vector supera la magnitud de interés se lo salta
            if sbin > numbin
                continue
            end

            # Producto escalar q·r para todas las partículas
            kr = r[:, 1:3] * qr     # vector de tamaño ntot

            # Sumas con pesos (partes real e imaginaria de la suma de exp(i q·r))
            sfsum1 = sum(cos.(kr / s))  # parte real se divide entre el diametro de la partícula
            sfsum2 = sum(sin.(kr / s))  # parte imaginaria

            # Acumular |S(q)|^2 / N y aumentar el conteo del bin
            sfhis[sbin, 1] += (sfsum1^2 + sfsum2^2) / ntotav
            sfhis[sbin, 2] += 1
        end
    end
end



(qxhis,qyhis,qzhis,qhis)=createqdom(qmax,rq0,dq0,bin0,numbin);
#qmag=sum(qx.^2 + qy.^2 + qz.^2)














# Promedio del factor de estructura para cada magnitud
sfhis[:, 1] = sfhis[:, 1] ./ sfhis[:, 2]

# Encontrar el índice del bin con el valor máximo de S(q)
kmax = argmax(sfhis[:, 1])

# Valor máximo del factor de estructura
smax = sfhis[kmax, 1]

# Almacenar la información del factor de estructura
info=zeros(numbin-1,3)

# Bucle sobre todos los bines excepto el último (Fortran: ij=1,numbin-1)
for ij in 1:numbin-1
    ki = ij * bin0                        # centro del bin ij
    sq = sfhis[ij, 1]                     # S(q) promedio
    sq_norm = sq / smax                   # S(q) / S_max
    info[ij,1]=ki
    info[ij,2]=sq
    info[ij,3]=sq_norm
end


# --- Cálculo del factor de estructura ---
sfhis .= 0.0 # Se inicializa la variable en la que se guardara el factor de estructura
for iq in -qmax:qmax
    qr1 = iq * dq0          # componente x del vector q
    for jq in -qmax:qmax
        qr2 = jq * dq0      # componente y
        for kq in -qmax:qmax
            
            # Excluir el origen del espacio recíproco
            if iq == 0 && jq == 0 && kq == 0
                continue
            end
            
            # Componente z del vector q
            qr3 = kq * dq0                     # componente z
           
            qr = [qr1, qr2, qr3]               # vector q
            rq = sqrt(qr1^2 + qr2^2 + qr3^2)   # magnitud |q|
            if rq > rq0
                continue
            end

            # Determinar el índice del bin (manejo especial del borde)
            sbin = floor(Int, rq / bin0) + 1

            # Coso para trnasformar de un dominio continup a uno discreto 
            if rq % bin0 == 0.0
                sbin -= 1
            end
           
            # Si la magnitud del vector supera la magnitud de interés se lo salta
            if sbin > numbin
                continue
            end

            # Producto escalar q·r para todas las partículas
            kr = r[:, 1:3] * qr     # vector de tamaño ntot

            # Sumas con pesos (partes real e imaginaria de la suma de exp(i q·r))
            sfsum1 = sum(cos.(kr / s))  # parte real se divide entre el diametro de la partícula
            sfsum2 = sum(sin.(kr / s))  # parte imaginaria

            # Acumular |S(q)|^2 / N y aumentar el conteo del bin
            sfhis[sbin, 1] += (sfsum1^2 + sfsum2^2) / ntotav
            sfhis[sbin, 2] += 1
        end
    end
end
 
# Promedio del factor de estructura para cada magnitud
sfhis[:, 1] = sfhis[:, 1] ./ sfhis[:, 2]

# Encontrar el índice del bin con el valor máximo de S(q)
kmax = argmax(sfhis[:, 1])

# Valor máximo del factor de estructura
smax = sfhis[kmax, 1]

# Almacenar la información del factor de estructura
info=zeros(numbin-1,3)

# Bucle sobre todos los bines excepto el último (Fortran: ij=1,numbin-1)
for ij in 1:numbin-1
    ki = ij * bin0                        # centro del bin ij
    sq = sfhis[ij, 1]                     # S(q) promedio
    sq_norm = sq / smax                   # S(q) / S_max
    info[ij,1]=ki
    info[ij,2]=sq
    info[ij,3]=sq_norm
end


 --- Cálculo del factor de estructura ---
sfhis .= 0.0 # Se inicializa la variable en la que se guardara el factor de estructura
for iq in -qmax:qmax
    qr1 = iq * dq0          # componente x del vector q
    for jq in -qmax:qmax
        qr2 = jq * dq0      # componente y
        for kq in -qmax:qmax
            
            # Excluir el origen del espacio recíproco
            if iq == 0 && jq == 0 && kq == 0
                continue
            end
            
            # Componente z del vector q
            qr3 = kq * dq0                     # componente z
           
            qr = [qr1, qr2, qr3]               # vector q
            rq = sqrt(qr1^2 + qr2^2 + qr3^2)   # magnitud |q|
            if rq > rq0
                continue
            end

            # Determinar el índice del bin (manejo especial del borde)
            sbin = floor(Int, rq / bin0) + 1

            # Coso para trnasformar de un dominio continup a uno discreto 
            if rq % bin0 == 0.0
                sbin -= 1
            end
           
            # Si la magnitud del vector supera la magnitud de interés se lo salta
            if sbin > numbin
                continue
            end

            # Producto escalar q·r para todas las partículas
            kr = r[:, 1:3] * qr     # vector de tamaño ntot

            # Sumas con pesos (partes real e imaginaria de la suma de exp(i q·r))
            sfsum1 = sum(cos.(kr / s))  # parte real se divide entre el diametro de la partícula
            sfsum2 = sum(sin.(kr / s))  # parte imaginaria

            # Acumular |S(q)|^2 / N y aumentar el conteo del bin
            sfhis[sbin, 1] += (sfsum1^2 + sfsum2^2) / ntotav
            sfhis[sbin, 2] += 1
        end
    end
end
 
# Promedio del factor de estructura para cada magnitud
sfhis[:, 1] = sfhis[:, 1] ./ sfhis[:, 2]

# Encontrar el índice del bin con el valor máximo de S(q)
kmax = argmax(sfhis[:, 1])

# Valor máximo del factor de estructura
smax = sfhis[kmax, 1]

# Almacenar la información del factor de estructura
info=zeros(numbin-1,3)

# Bucle sobre todos los bines excepto el último (Fortran: ij=1,numbin-1)
for ij in 1:numbin-1
    ki = ij * bin0                        # centro del bin ij
    sq = sfhis[ij, 1]                     # S(q) promedio
    sq_norm = sq / smax                   # S(q) / S_max
    info[ij,1]=ki
    info[ij,2]=sq
    info[ij,3]=sq_norm
end



