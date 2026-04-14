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

function waveVectorSq(N,L)
"""
    Se calcula el vector de onda para el factor de estructura.
    Regresa un diccionario con los vectores agrupados por magnitud
"""
    tol=10;
    grupos=Dict{Float64, Vector{Tuple{Float64,Float64,Float64}}}();

    #N=ceil(Int,8/sqrt(27)*L); Ese esta cool, pero se tarda un chingo
    N=20; # Por mientrás, que no está optimizado
    Q_max=sqrt(3*(2*pi/L)^2);

    for nx in -N:N, ny in -N:N, nz in -N:N
        qx=2*pi*nx/L;
        qy=2*pi*ny/L;
        qz=2*pi*nz/L;
        qq=qx^2 + qy^2 + qz^2;
        q=sqrt(qq);
    
        if q<Q_max || q == 0.0
            continue
        end

        clave=round(qq, digits=tol);
        push!(get!(grupos,clave,[]),(qx,qy,qz))
    end

    return grupos
end

function structureFactor(grupos,rx,ry,rz,N_part)
"""
    Compute the structure factor given the wave vector in a dictionary and positions.
"""

    # Vector de onda
    grupos=waveVectorSq(N,L); 

    Sq_ensamble=[];

    for (clave, lista_q) in grupos
        q_mod=sqrt(clave);
        S_acumulado=0.0;
    
        # Compute structure factor per group
        for v_q in lista_q
            # Compute dot product
            qr=v_q[1].*rx .+ v_q[2].*ry .+ v_q[3].*rz;
        
            # Compute the structure factor
            S_q=(1/N_part)*(sum(cos.(qr))^2 + sum(sin.(qr))^2);
        
            # Guardar el factor de estructura de cada vector
            push!(Sq_ensamble,(q_mod,S_q));
        end
    end

    # Ordena de mayor a menor
    Sq_ensamble=sort!(Sq_ensamble);

    return Sq_ensamble

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

M_frames=100;

Sq_observable=[];

for (DIR_id,dir) in enumerate(dat_DF.dir) # Cycle thru directories
    # Stored timesteps
    aux_id=Int.((0:dat_DF."save-dump"[DIR_id]:(dat_DF."N_heat"[DIR_id] + dat_DF."N_isot"[DIR_id])));

    # REduce the amount of frames, because I no time for optimization
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

        # número de partículas en la simulación
        N_part=dat_DF.Npart[DIR_id];

        # Range of wave vector
        L=2*dat_DF.L[5];

        # Factor de estructura
        Sq_ensamble=structureFactor(N,L,rx,ry,rz,N_part);
    
        push!(Sq_observable,Sq_ensamble)
    end

end


