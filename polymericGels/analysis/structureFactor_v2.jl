"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra
using Distances

# Include auxiliary files
include("functions.jl")

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
N_qu=2^5; # EXPONENTE DEBE SER IMPAR Cantidad de direcciones
lambda_o=1; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=32; # Cantidad de magnitudes
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

    theta=2*pi*rand(N_theta);
    phi=pi*rand(N_phi); 

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    # Vector de N elementos. Cada elemento es la posición de N partículas de los N_exp.
    r_exp=[getPosition(df) for df in dumps];

    # Extraer la posición y separar por componentes
    r=r_exp[1];
    r_x=reshape(r[1],1,N_part);
    r_y=reshape(r[2],1,N_part);
    r_z=reshape(r[3],1,N_part);
   
    # Calcular la distancia entre partículas considerando condiciones periodicas de frontera
    dist_x=pairwise(PeriodicEuclidean(lambda_f),r_x,dims=2);
    dist_y=pairwise(PeriodicEuclidean(lambda_f),r_y,dims=2);
    dist_z=pairwise(PeriodicEuclidean(lambda_f),r_z,dims=2);

    utri=triu!(trues(size(dist_x)));
    utri=utri;

    # Crear los vectores unitarios del vector de onda
    q_x=reshape([cos(th)*sin(ph) for th in theta, ph in phi],N_qu,1);
    q_y=reshape([sin(th)*sin(ph) for th in theta, ph in phi],N_qu,1);
    q_z=reshape([cos(ph) for th in theta, ph in phi],N_qu,1);

    # Calcular el producto punto para distintas direcciones
    pp=[q_x[it]*dist_x+q_y[it]*dist_y+q_z[it]*dist_z for it in eachindex(q_x)];
    
    # Get upper triangle
    pp=reduce(hcat,[filter(!iszero,it[utri]) for it in pp]);


    #pp_x=[vq*dist_x for vq in q_x];
    #pp_y=[vq*dist_y for vq in q_y];
    #pp_z=[vq*dist_z for vq in q_z];

    #pp=pp_x .+ pp_y .+ pp_z;

    q_min=2*pi/lambda_f;
    q_max=2*pi/lambda_o;
    q_dom=range(q_min,q_max,length=N_lambda);

    #A=[mean(sum.([2*cos.(q*it) for it in pp])) for q in q_dom];
    
    data=zeros(N_lambda);
    Threads.@threads for id in eachindex(q_dom)
        data[id]=mean(sum.([2*cos.(q_dom[id]*it) for it in pp]));
    end



    #B=[mean(sum.([sin.(it) for it in pp])) for q in q_dom];


    #pp=[(cos(th)*sin(ph).*dist_x) .+ (sin(th)*sin(ph).*dist_y) .+ (cos(ph).*dist_z) for th in theta, ph in phi];

#=
    # Calcular el producto punto
    pp=[(cos(th)*sin(ph).*dist_x) .+ (sin(th)*sin(ph).*dist_y) .+ (cos(ph).*dist_z) for th in theta, ph in phi];

    # Calcular la densidad promedio de cada magnitud
    q_min=2*pi/lambda_f;
    q_max=2*pi/lambda_o;
    q_dom=range(q_min,q_max,length=N_lambda);

    #density=[mean([sum(cos.(mag.*it) - I)^2 + sum(sin.(mag.*it))^2 for it in pp])/N_part for mag in q_dom];

    data=zeros(N_lambda);
    Threads.@threads for id in eachindex(q_dom)
        data[id]=mean([sum(cos.(q_dom[id].*it) - I)^2 + sum(sin.(q_dom[id].*it))^2 for it in pp])/N_part;
    end
=#
