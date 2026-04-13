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

function getDatInfo(DF_DIR)
"""
    This function get the dat.csv of experiments
"""

    archivos = filter(str -> occursin("dat-", str), DF_DIR)
    df_final = DataFrame()  # vacío inicial
    for (i, archivo) in enumerate(archivos)
        df_temp = CSV.read(archivo, DataFrame)
        
        # IF there are is missing stuff
        for col in names(df_temp)
            if eltype(df_temp[!, col]) <: Union{Missing, Number}
                df_temp[!, col] = coalesce.(df_temp[!, col], 0.0)
            end
        end

        append!(df_final, df_temp)
        df_temp = nothing  # liberar referencia
    end
    return df_final
end


# Get directories 
MAIN_DIR=dirname(pwd());
DATA_DIR=joinpath(MAIN_DIR,"data");
INFO_DIR=joinpath(pwd(),"data_mod");
DF_DIR=filter(isfile,readdir(INFO_DIR,join=true));

dat_files=getDatInfo(DF_DIR);

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

DIR_id=5;

M_frames=1000;

# Stored timesteps
aux_id=Int.((0:dat_DF."save-dump"[DIR_id]:(dat_DF."N_heat"[DIR_id] + dat_DF."N_isot"[DIR_id])));

# Directory of the files
TRAJ_DIR=joinpath(DATA_DIR,dat_DF.dir[DIR_id],"traj");

aux_dump=getDump(TRAJ_DIR,string("traj_assembly.",M_frames,".dumpf"));

# Extraer las posiciones del dump de las partículas centrales.
mask=(aux_dump.type .==1) .| (aux_dump.type .== 2.0);
df_filtered=aux_dump[mask,:];

rx=df_filtered.x;
ry=df_filtered.y;
rz=df_filtered.z;
r=[rx ry rz];

# número de partículas en la simulación
N_part=dat_DF.Npart[DIR_id];

# Compute Structure factor

# Vector de onda
#=
    Se obtienen vectores de onda agrupados por su  magnitud
=#

tol=10;
grupos=Dict{Float64, Vector{Tuple{Float64,Float64,Float64}}}();

N=30;
Q_max=sqrt(3*(2*pi/2*dat_DF.L[5])^2);

for nx in -N:N, ny in -N:N, nz in -N:N
    qx=2*pi*nx/2*dat_DF.L[5];
    qy=2*pi*ny/2*dat_DF.L[5];
    qz=2*pi*nz/2*dat_DF.L[5];
    qq=qx^2 + qy^2 + qz^2;
    q=sqrt(qq);
    
    if q<Q_max || q == 0.0
        continue
    end

    clave=round(qq, digits=tol);
    push!(get!(grupos,clave,[]),(qx,qy,qz))
end

# Compute the structure factor per magnitude
Sq_observable=[];

for (clave, lista_q) in grupos
    q_mod=sqrt(clave);
    S_acumulado=0.0;
    for v_q in lista_q
        qr=v_q[1].*rx .+ v_q[2].*ry .+ v_q[3].*rz;
        S_q=(1/N_part)*(sum(cos.(qr))^2 + sum(sin.(qr))^2);
        S_acumulado += S_q;
    end
    S_promedio = S_acumulado / length(lista_q);
    push!(Sq_observable,(q_mod,S_promedio));
end

# Ordena de mayor a menor
Sq_observable=sort!(Sq_observable);





