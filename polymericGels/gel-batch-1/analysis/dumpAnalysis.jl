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

M_frames=100;

save_observables=[zeros(Int,M_frames,2) for _ in eachindex(dat_DF.dir)];

for (DIR_id,dir) in enumerate(dat_DF.dir)

    # Stored timesteps
    aux_id=Int.((0:dat_DF."save-dump"[DIR_id]:(dat_DF."N_heat"[DIR_id] + dat_DF."N_isot"[DIR_id])));

    # REduce the amount of frames, because I no time for optimization
    ind_graph=floor.(Int,range(1,length(aux_id),length=M_frames));

    # Directory of the files
    TRAJ_DIR=joinpath(DATA_DIR,dat_DF.dir[DIR_id],"traj");

    # Alocate memory
    observables=zeros(Int,length(ind_graph),2);


    for (it,ind) in enumerate(ind_graph)

        aux_dump=getDump(TRAJ_DIR,string("traj_assembly.",aux_id[ind],".dumpf"));

        # Filtrar
        mask=(aux_dump.type .==1) .| (aux_dump.type .== 2.0);
        df_filtered=aux_dump[mask,:];

        # Contar frecuencias de c_clusters (más rápido que groupby+combine)
        counts = countmap(df_filtered.c_clusters);
        cluster_observable = [counts[cluster] for cluster in sort(collect(keys(counts)))];

        observables[it,1]=length(cluster_observable);
        observables[it,2]=maximum(cluster_observable);
    end

    save_observables[DIR_id]=observables; 

end

# Compute the average
avg_observables=ceil.(Int,dropdims(mean(cat(save_observables..., dims=3), dims=3), dims=3));

# Save the observables
file_name=string("clusterObservables_",dat_DF.id_ref[1],".csv");
CSV.write(joinpath(INFO_DIR,file_name),DataFrame(avg_observables,[:Nclustres, :Maxparticles]));

