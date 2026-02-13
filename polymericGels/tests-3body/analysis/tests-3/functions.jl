"""
    Functions 
"""

function ForcePatch_fd(eps_pair,sig_p,r)
"""
    Get the central finite difference given the value of the position and the function.
"""
    dh=1e-6;
    fo=Upatch(eps_pair,sig_p,r+dh)
    ff=Upatch(eps_pair,sig_p,r-dh)
    return -(1/(2*dh))*( fo - ff );
end

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

function U3(eps_pair,eps_3,sig_p,r)
"""
    Auxiliary potential to create Swap Mechanism based in Patch-Patch interaction
"""
    if r <= sig_p 
        return 1.0
    elseif r >= 1.5*sig_p
        return 0.0 
    else 
        return -(1/eps_3)*Upatch(eps_pair,sig_p,r)
    end
end

function SwapU(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_ik)
"""
    Potential for the swap mechanism
"""
    return w.*eps_jk.*U3(eps_ij,eps_jk,sig_p,r_ij).*U3(eps_ik,eps_jk,sig_p,r_ik)
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


