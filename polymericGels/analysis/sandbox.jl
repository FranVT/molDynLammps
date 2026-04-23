"""
    Script that transforms the fix file into a DataFrame.
    An average of observables of a set of simulations of the same system 
"""

using DataFrames, CSV
using Statistics, StatsBase
using UUIDs, LinearAlgebra

# Include auxiliary files
include("functions.jl")

function getDump(dir,file_name)
"""
    Get the data from a single dump file that stores one timeste information
"""
    data = split.(readlines(joinpath(dir,file_name))," ")[9:end];
    HEADERS=data[1][3:end];
    INFO=parse.(Float64,reduce(hcat,data[2:end]))';

    return DataFrame(INFO,HEADERS)
end

function dotSpherical(theta,phi,r)
"""
    Compute the dot product betwen a position and a unit vector r in psherical coordinates.
"""
    q_x=cos(theta)*sin(phi);
    q_y=sin(theta)*sin(phi);
    q_z=cos(phi);
    return q_x*r[1]+q_y*r[2]+q_z*r[3]
end

function densityRhoQ(q_mag,dot_qr)
"""
    Compute the squared of the absolute value of the density at the reciprocal space.
    |rho(r)|^2 = A(vec{q}cdotvec{r})^2 + B(vec{q}cdotvec{r}^2)
    A = sumcos(); B = sumsin()
"""
    return sum(cos.(q_mag*dot_qr))^2 + sum(sin.(q_mag*dot_qr))^2
end

function computeDensity(N_qu,lambda_o,lambda_f,N_lambda,r)
"""
    Function that computes the static structure factor de different wave vectors.
    Returns a vector with the following interpretation of the values:
    [row] -> [magnitude]}
"""

    # Vector unitario del vector de onda
    theta=2*pi*rand(N_qu);
    phi=pi*rand(N_qu); 

    # Calculo del producto punto
    dot_qr=[dotSpherical(theta[s],phi[s],r) for s in 1:N_qu];

    # Evaluación de la densidad y promedio
    # [renglon x columna] -> [ mag x direccion ]
    rho_q=[densityRhoQ(2*pi/l,d) for l in range(lambda_o,lambda_f,length=N_lambda), d in dot_qr];
    return mean(rho_q,dims=2)

end

function getPosition(dump)
"""
    Get the position of the central particles of a given dump
"""

    # Filtrar
    mask=(dump.type .==1) .| (dump.type .== 2.0);
    dump_filtered=dump[mask,:];

   return [dump_filtered.x,dump_filtered.y,dump_filtered.z]

end

function structureFactor(N_qu,lambda_o,lambda_f,N_lambda,r_exp)
"""
    Compute the static structure factor
"""
    S_q_exp=[computeDensity(N_qu,lambda_o,lambda_f,N_lambda,r) for r in r_exp];

    return mean(reduce(hcat,S_q_exp),dims=2)/dat_DF.Npart[1];
end

function getTimeEvolSq(dump_paths,time_instant)
"""
    Compute the time evolution of the structure factor
"""

    # Obtenemos los dumps de los N experimentos para un instante de tiempo
    dumps=[getDump(path,time_instant) for path in dump_paths];

    r_exp=[getPosition(df) for df in dumps];

    return structureFactor(N_qu,lambda_o,lambda_f,N_lambda,r_exp);

end

######
#   SCRIPT
######

# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","experiments_dat.csv");

dat_files=CSV.read(DAT_PATH,DataFrame);

# Selection of the system by parameters
phi=0.02;
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
N_qu=2^5; # Cantidad de direcciones
lambda_o=1; # Limites del rango a explorar
lambda_f=100; # Limites del rango a explorar
N_lambda=2^6; # Cantidad de magnitudes

# Seleccion de time instants
aux_id=Int.((0:dat_DF."save-dump"[1]:(dat_DF."N_heat"[1] + dat_DF."N_isot"[1])));
aux_id=aux_id[1:100:end];


time_instants=[replace("traj_assembly.*.dumpf", "*" => string(it)) for it in aux_id];

time_instant=time_instants[1];


S_q=getTimeEvolSq(dump_paths,time_instant);


Sq_t=[getTimeEvolSq(dump_paths,time_instant) for time_instant in time_instants];



#=

M_frames=100;   # Cada 100 unidades temporales

Sq_observable=[];

for (DIR_id,dir) in enumerate(dat_DF.dir[1]) # Cycle thru directories
    # Stored timesteps
    aux_id=Int.((0:dat_DF."save-dump"[DIR_id]:(dat_DF."N_heat"[DIR_id] + dat_DF."N_isot"[DIR_id])));

    # Reduce the amount of frames, because I no time for optimization
    ind_graph=floor.(Int,range(1,length(aux_id),length=M_frames));

    # Directory of the files
    TRAJ_DIR=joinpath(DATA_DIR,dat_DF.dir[DIR_id],"traj");

    for (it,ind) in enumerate(ind_graph) # Cycle thru dump files

        aux_dump=getDump(TRAJ_DIR,string("traj_assembly.",aux_id[ind],".dumpf"));

        # Extraer las posiciones del dump de las partículas centrales.
        mask=(aux_dump.type .==1) .| (aux_dump.type .== 2.0);
        df_filtered=aux_dump[mask,:];

        # Central particles positions
        rx=df_filtered.x;
        ry=df_filtered.y;
        rz=df_filtered.z;

        # Range of wave vector
        L=2*dat_DF.L[5];

        # Factor de estructura
        Sq_ensamble=structureFactor(L,rx,ry,rz,dat_DF.Npart[DIR_id]);
    
        push!(Sq_observable,Sq_ensamble)
    end
end

=#
