#=
    Factor de estructura
=#

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

#Sq_t=[getTimeEvolSq(N_qu,dump_paths,time_instant) for time_instant in time_instants];
Sq_t=for time_instant in time_instants
       # Vector unitario del vector de onda
    N_phi=Int64(sqrt(div(N_qu,2)));
    N_theta=Int64(2*N_phi);

    theta=2*pi*rand(N_theta);
    phi=pi*rand(N_phi); 

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    r_exp=[getPosition(df) for df in dumps];

    #Sq=structureFactor(theta,phi,lambda_o,lambda_f,N_lambda,r_exp)
    #data=[computeDensity(theta,phi,lambda_o,lambda_f,N_lambda,r) for r in r_exp];
    data=for r in r_exp
        # Calculo del producto punto
        dot_qr=[dotSpherical(th,ph,r) for th in theta, ph in phi];

        # Evaluación de la densidad y promedio
        # [renglon x columna] -> [ mag x direccion ]
        q_min=2*pi/lambda_f;
        q_max=2*pi/lambda_o;
        q_dom=range(q_min,q_max,length=N_lambda);

        # Cálculo de la densidad
        rho_q=[sum(cos.(l*d))^2 + sum(sin.(l*d))^2 for l in q_dom, d in dot_qr];

        # Compute the avg of the different direction but same magnitude
        rho_q=reduce(vcat,mean(mean(rho_q,dims=3),dims=2));

    end

    data=reduce(hcat,data);
    q_domain=collect(first(unique(data[1,:])));

    # Compute assembly average and scale of 1/N_particles 
    Sq=reduce(vcat,mean(reduce(hcat,data[2,:]),dims=2))./length(r_exp[1]);

end

q_domain=[Sq_t[it][1] for it in eachindex(Sq_t)];
Sq=[Sq_t[it][2] for it in eachindex(Sq_t)];

