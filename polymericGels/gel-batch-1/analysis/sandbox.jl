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

# Compute Structure factor


