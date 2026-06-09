"""
    Script en el que se adapta el código de FORTRAN de Claudia a Julia
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll
using SplitApplyCombine


include("functions.jl")

function getpathfilesSq(dat_DF,N_instants)
"""
    Get the names files 
"""
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
     return (timeSteps,dump_paths,file_names)
end



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

            end
        end
    end

    return (qxhis,qyhis,qzhis,qhis) 
end


function computeSq(numbin,ntotav,qxhis,qyhis,qzhis,qhis,r)
"""
    Compute the structure factor of a set of positions
"""

    Sq=zeros(numbin,2);
    rho=[[] for _ in 1:numbin];
    for it_bin in 1:numbin
    
        # Seleccion de las componentes del vector de onda
        qx=qxhis[it_bin];
        qy=qyhis[it_bin];
        qz=qzhis[it_bin];
        
        # Calcular la densidad para cada bin
        for it_q in eachindex(qx)
            vq=[qx[it_q], qy[it_q], qz[it_q]];
            dp=r*vq; # Producto punto del vector para cada partícula
            rho_re=sum(cos.(dp));
            rho_im=sum(sin.(dp));
            sq=(rho_re^2 + rho_im^2)/ntotav;
            append!(rho[it_bin],sq);
        end
    end

    # Guardamos información
    # Valor esperado del factor de estructura
    Sq[:,1]=sum.(rho)./length.(qhis);
    smax=maximum(Sq[:,1])
    Sq[:,2]=Sq[:,1]/smax

    return Sq
end

function computeSqmean(file_names,timeStep,file_path,dump_paths,dat_DF,numbin,ntotav,qxhis,qyhis,qzhis,qhis,info)
"""
    Compute the avg structure factor of a system at a given time step
"""
     # Preparar informacion para el factor de estructura
    r_system = getPositions(file_names, dump_paths);

    # Calcular el factor de estructura para N experimentos del mismo sistema
    for it_exp in eachindex(r_system)
        r = reduce(hcat,r_system[it_exp]);
        info[it_exp][:,2:3]=computeSq(numbin,ntotav,qxhis,qyhis,qzhis,qhis,r);
        println("One experiment done")
    end

    # Se obtienen promedios
    info=mean(info);        # Compute the assembly mean of the Sq
    info[:,1]=mean.(qhis);  # Compute the mean of the wave vectors 

    # Paths and file names and stuff
    filename=string("structureFactorPBC",first(dat_DF.id),"time",timeStep,".csv");

    # Create the data frame
    df=DataFrame([repeat([timeStep],numbin) info],[:timeStep, :qmean, :Sqmean, :Sqmeannorm])

    # Save the data frame
    CSV.write(joinpath(file_path,filename),df)   
end

function computeAllTimeSqmean(dat_DF,N_instants,qmax0)
"""
    Compute all structure factor of all system of N time instants
"""

    # Parameters of the system relevant to the structure factor
    L=2*dat_DF.L[1];            # Longitud de la caja
    N_exp=nrow(dat_DF);         # Cantidad de experimentos por sistema
    N_part=dat_DF.Npart[1];     # Numero de particulas centrales

    # Crea los dominios del vector de onda
    xc = L;                 # longitud de la caja en x
    yc = L;                 # longitud de la caja en y
    zc = L;                 # longitud de la caja en z
    dq0 = 2*pi / xc;        # Δq fundamental
    qmax = Int(floor(qmax0 / dq0)) # número entero de pasos hasta qmax0
    rq0 = qmax * dq0;       # valor máximo real de |q| usado
    bin0 = dq0;             # nuevo ancho de bin (bin0 original * dq0)
    numbin = Int(floor(qmax * dq0 / bin0)) + 1;  # número total de bines

    # Parametros para el factor de estructura
    ntotav=Int(N_part); # Número total de partículas
    (qxhis,qyhis,qzhis,qhis)=createqdom(qmax,rq0,dq0,bin0,numbin);  # vector de onda
    file_path=joinpath(pwd(),"analyzedData");

    # Obtener los paths y los nombres de los archivos a analizar
    (timeSteps,dump_paths,file_names)=getpathfilesSq(dat_DF,N_instants);

    # Select one time step
    # Alocar memoria
    info=[zeros(numbin,3) for _ in 1:N_exp];

    for it_time in eachindex(file_names)
        computeSqmean(file_names[it_time],timeSteps[it_time],file_path,dump_paths,dat_DF,numbin,ntotav,qxhis,qyhis,qzhis,qhis,info)
        println("One time step done")
    end
    
    println("One system done")
end

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);


# Main Parameters of the analysis
N_systems=length(data_bySystem);
N_instants=25;               # Instantes temporales a analizar
qmax0 = 6;              # 3 es el min sin que cause problemas

for it_system in (1,2,3,5) 
    computeAllTimeSqmean(data_bySystem[it_system],N_instants,qmax0)
end

#computeAllTimeSqmean(data_bySystem[5],N_instants,qmax0)
