"""
    Functions 
"""

function Upatch(eps_pair,sig_p,r)
"""
    Auxiliary potential to create Swap Mechanism based in Patch-Patch interaction
"""
    if r < 1.5*sig_p 
        return 2*eps_pair*( ((sig_p^4)./((2).*r.^4)) .-1).*exp.((sig_p)./(r.-(1.5*sig_p)).+2)
    else
        return 0.0
    end
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


