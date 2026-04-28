"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra

# Include auxiliary files
include("functions.jl")

# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","experiments_dat.csv");
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

dat_files=CSV.read(DAT_PATH,DataFrame);

# Selection of the system by parameters
phi=0.01;
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
N_qu=2^7; # Cantidad de direcciones
lambda_o=0.2; # Limites del rango a explorar (Monomero)
lambda_f=2*dat_DF.L[1]; # Limites del rango a explorar (Tamaño de la caja)
N_lambda=2^8; # Cantidad de magnitudes

# Seleccion de time instants
aux_timeStep=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
ind=round.(Int, LinRange(1, length(aux_timeStep), 2));
aux_id=aux_timeStep[ind];

time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];



function computeDensity(theta,phi,lambda_o,lambda_f,N_lambda,r)
"""
    Function that computes the static structure factor de different wave vectors.
    Returns a vector with the following interpretation of the values:
    [row] -> [magnitude]}
"""

    # Calculo del producto punto
    dot_qr=[dotSpherical(theta[s],phi[s],r) for s in 1:length(theta)];

    # Evaluación de la densidad y promedio
    # [renglon x columna] -> [ mag x direccion ]
    q_min=2*pi/lambda_f;
    q_max=2*pi/lambda_o;
    q_dom=range(q_min,q_max,length=N_lambda);

    rho_q=[densityRhoQ(l,d) for l in q_dom, d in dot_qr];

    return [q_dom,reduce(vcat,mean(rho_q,dims=2))]

end

function structureFactor(theta,phi,lambda_o,lambda_f,N_lambda,r_exp)
"""
    Compute the static structure factor
"""
    S_q_exp=[computeDensity(theta,phi,lambda_o,lambda_f,N_lambda,r) for r in r_exp];

    return mean(reduce(hcat,S_q_exp),dims=2)/dat_DF.Npart[1];
end



function getTimeEvolSq(N_qu,dump_paths,time_instant)
"""
    Compute the time evolution of the structure factor
"""

    # Vector unitario del vector de onda
    theta=2*pi*rand(N_qu);
    phi=pi*rand(N_qu); 

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    r_exp=[getPosition(df) for df in dumps];

    return structureFactor(theta,phi,lambda_o,lambda_f,N_lambda,r_exp);

end



getTimeEvolSq(dump_paths,time_instants[end])





#Sq_t=[getTimeEvolSq(dump_paths,time_instant) for time_instant in time_instants];

















#Sq_df=DataFrame(;
#    timeStep = aux_id,
#    Sq       = reduce.(vcat, Sq_t),
#    lambda_o = fill(lambda_o, length(aux_id)),
#    lambda_f = fill(lambda_f, length(aux_id)),
#    id=unique(dat_DF.id)[1]
#)

#CSV.write(joinpath(SAVE_DIR,string("structureFactor",Sq_df.id[1],".csv")),Sq_df)
