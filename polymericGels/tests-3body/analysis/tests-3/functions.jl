"""
    Functions 
"""

# Include the functions for the force and potentials
include("../../src/tables/functions.jl")

function evaluateFswap(dist_12,dist_13,dist_23,DATA_dump_1,DATA_dump_2,DATA_dump_3)
"""
    Evaluate the Force to get analytical result
"""

f_swap1=mapreduce(s->forceSwap(1.0,1.0,1.0,1.0,0.4,dist_12[s],dist_13[s]),vcat,1:1:nrow(DATA_dump_1));
f_swap2=mapreduce(s->forceSwap(1.0,1.0,1.0,1.0,0.4,dist_12[s],dist_23[s]),vcat,1:1:nrow(DATA_dump_2));
f_swap3=mapreduce(s->forceSwap(1.0,1.0,1.0,1.0,0.4,dist_13[s],dist_23[s]),vcat,1:1:nrow(DATA_dump_3));

# Stuff of the table
F_1s12 = f_swap1[:,1] .- f_swap2[:,1];
F_1s13 = f_swap1[:,2] .- f_swap3[:,1];

F_2s21 = f_swap2[:,1] .- f_swap1[:,1];
F_2s23 = f_swap2[:,2] .- f_swap3[:,2];

F_3s31 = f_swap3[:,1] .- f_swap1[:,2];
F_3s32 = f_swap3[:,2] .- f_swap2[:,2];

# Project the forces into the x,y basis
f_swap1xy=(F_1s12).*vr_12 .+ (F_1s13).*vr_13;
f_swap2xy=(F_2s21).*(vr_12) .+ (F_2s23).*vr_23;
f_swap3xy=(F_3s31).*(vr_13) .+ (F_3s32).*(vr_23);


    return (f_swap1xy,f_swap2xy,f_swap3xy)

end

function getTable3b(dir,file_name)
"""
    Get the data from the table file of two particle interaction. 
"""
    data=split.(readlines(joinpath(dir,file_name))," ")[4:end];
    HEADERS=["n","r_ij","r_ik","theta","f_i1","f_i2","f_j1","f_j2","f_k1","f_k2","e"];
    INFO=reduce(hcat,map(s->parse.(Float64,s),data))';
    df=DataFrame(INFO,HEADERS);

    return df 
end

function getTable(dir,file_name)
"""
    Get the data from the table file of two particle interaction. 
"""
    data=split.(readlines(joinpath(dir,file_name))," ")[7:end];
    HEADERS=["n","r","u","f"];
    INFO=reduce(hcat,map(s->parse.(Float64,s),data))';
    df=DataFrame(INFO,HEADERS);

    return df 
end


function getDump(dir,file_name)
"""
    Get the data from a single dump file that stores one timeste information
"""
    it=parse(Int64,split(file_name,'.')[2]);   
    data = split.(readlines(joinpath(dir,file_name))," ")[9:end];
    HEADERS=data[1][3:end];
    INFO=parse.(Float64,reduce(hcat,data[2:end]))';
    df=DataFrame(INFO,HEADERS);
    df.TimeStep .= it;

    return df
end

function extractFixScalar(path_system,file_name)
"""
    Function that extracts the information of fix files that stores global scalar values
"""
    aux=split.(readlines(joinpath(path_system,file_name))," ");
    return (aux[2][2:end],reduce(hcat,map(s->parse.(Float64,s),aux[3:end])));
end


function getDir(date)
"""
    To get the directory of the simulation
"""

# Get the directories
MAIN_DIR=dirname(dirname(pwd()));
DATA_DIR=joinpath(MAIN_DIR,"data");

# Read the directory where data is stored
sims=filter(isdir,readdir(DATA_DIR,join=true));

# Get the idex for the directory
indx=findall(!isempty,findall.(date,sims));

# Get the directory
DIR=sims[indx];

    return DIR

end

