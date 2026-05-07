"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra
using Distances, Random
using GLMakie

Random.seed!(1234)

# Include auxiliary files
include("functions.jl")

    function computeDensity(pp)
    """
        Compute the density for the structure factor
        pp is the dot product.
        pp is an array. Each row represent the dot product of onde direction with a distance
    """
        return mapreduce(s->cos(s),+,pp)^2 + mapreduce(s->sin(s),+,pp)^2
    end

    function computeDistances(r,lambda_f)
        r_x=reshape(r[1],1,N_part);
        r_y=reshape(r[2],1,N_part);
        r_z=reshape(r[3],1,N_part);

        # Calcular la distancia entre partículas considerando condiciones periodicas de frontera
        dist_x=pairwise(PeriodicEuclidean(lambda_f),r_x,dims=2);
        dist_y=pairwise(PeriodicEuclidean(lambda_f),r_y,dims=2);
        dist_z=pairwise(PeriodicEuclidean(lambda_f),r_z,dims=2);

        # Reducir la cantidad de elementos para evaluar distancias
        dist_x=[dist_x[it...] for it in ids_utri];
        dist_y=[dist_y[it...] for it in ids_utri];
        dist_z=[dist_z[it...] for it in ids_utri];

        return [dist_x dist_y dist_z]
        #return (r[1],r[2],r[3])
    end


# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","experiments_dat.csv");
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

dat_files=CSV.read(DAT_PATH,DataFrame);

# Selection of the system by parameters
phi=0.05;
Temp=0.05;
N_part=5000.0;
CL_con=0.05;

# Se filtra el dataframe 
dat_DF = subset(dat_files,
    :phi => ByRow(==(phi)),
    :Temperature => ByRow(==(Temp)),
    :Npart => ByRow(==(N_part)),
    :"CL-Con" => ByRow(==(CL_con))
)

# Path to the dumps
dump_paths=joinpath.(dat_DF.PARENT_DIR,dat_DF.dir,"traj");

# Parametros para obtener el factor de estructura
N_qu=2^7; # EXPONENTE DEBE SER IMPAR Cantidad de direcciones
lambda_o=1; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=64; # Cantidad de magnitudes
N_instants=2;


# Seleccion de time instants
aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
ind=round.(Int, LinRange(1, length(aux_timeStep), N_instants));
aux_id=aux_timeStep[ind];

time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];

time_instant=time_instants[end];

    N_part=Int(dat_DF.Npart[1]);

    # Vector unitario del vector de onda
    N_phi=Int64(sqrt(div(N_qu,2)));
    N_theta=Int64(2*N_phi);

    theta=2*pi*range(0,1,length=N_theta);
    phi=pi*range(0,1,length=N_phi); 

    # Calcular la densidad promedio de cada magnitud
    q_min=2*pi/lambda_f;
    q_max=2*pi/lambda_o;
    q_dom=range(q_min,q_max,length=N_lambda);

    # Para todas las distancias
    ids_utri=reduce(vcat,[[[r,c] for c in r+1:N_part] for r in 1:N_part]);

    # Para todas las posiciones/factor de estructura
    # Crear los vectores unitarios del vector de onda
    q_x=mapreduce(ph->map(th->cos(th)*sin(ph),theta),vcat,phi); 
    q_y=mapreduce(ph->map(th->sin(th)*sin(ph),theta),vcat,phi);
    q_z=mapreduce(ph->map(th->cos(ph),theta),vcat,phi); 

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    # Vector de N elementos. Cada elemento es la posición de N partículas de los N_exp.
    r_exp=[getPosition(df) for df in dumps];

    # Allocate memory      
    S_q=[zeros(N_lambda) for _ in eachindex(r_exp)];

    for it_exp in eachindex(r_exp)
        # Seleccionar un experimento
        dist=computeDistances(r_exp[it_exp],lambda_f);
 
        for it_lambda in 1:N_lambda
            # Seleccionar una sola magnitud
            # Calcular el producto punto para todas las direcciones usando multiplicación matricial 
            pp=map(s->dist*((q_dom[it_lambda]).*[q_x[s],q_y[s],q_z[s]]),1:N_qu);
            println("Dot product done")

            # Calcular la densidad para una sola magnitud, dada N_qu direcciones
            rho_q=map(s->computeDensity(s),pp);
            println("Density done")

            # Obtener el factor de estructura para una sola magnitud y almacenar el resultado
            S_q[it_exp][it_lambda]=mean(rho_q)/N_part;

            println(it_lambda)
        end
        println("Experiment",it_exp," done")
    end



#=

fig=Figure()
ax=Axis(fig[1:1,1:1],
        limits=(first(q_dom),last(q_dom),0,5.0e7)
       )
vlines!(ax,2*pi/(1.2),linestyle=:dash,color=:blue)
vlines!(ax,2*pi/(2*1.2),linestyle=:dash,color=:blue)
vlines!(ax,2*pi/(3*1.2),linestyle=:dash,color=:blue)

#vlines!(ax,2*pi/(lambda_f),linestyle=:solid,color=:black)
vlines!(ax,2*pi/(0.5*lambda_f),linestyle=:solid,color=:black)
vlines!(ax,2*pi/(0.25*lambda_f),linestyle=:solid,color=:black)
vlines!(ax,2*pi/(0.125*lambda_f),linestyle=:solid,color=:black)
vlines!(ax,2*pi/(0.0625*lambda_f),linestyle=:solid,color=:black)

scatterlines!(ax,q_dom,S_q)

=#



































#=
    S_q=zeros(N_lambda);
    for it_q in eachindex(q_dom)
    
        # Calcular el producto punto. 
        # Misma magnitud. Distintas direcciones
        l=q_x[it_q].*dist_x' .+ q_y[it_q].*dist_y' .+ q_z[it_q].*dist_z';

        # Calcular la densidad
        A=map(id->mapreduce(s->cos(s),+,l[id,:]),1:N_qu);
        B=map(id->mapreduce(s->sin(s),+,l[id,:]),1:N_qu);
        rho_q=A.^2 .+ B.^2;

        # Obtener el factor de structura para una magnitud
        S_q[it_q]=mean(rho_q)/N_part;

        println(it_q)
    end
=#

#=
    S_q=zeros(N_lambda);
    for it_q in eachindex(q_dom)
        
        # Compute the density
        rho_q=computeDensity(q_x[it_q],q_y[it_q],q_z[it_q],dist_x,dist_y,dist_z);

        # Obtener el factor de structura para una magnitud
        S_q[it_q]=mean(rho_q)/N_part;

        println(it_q)
    end
=#














