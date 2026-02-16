"""
    Functions 
"""

function Upatch(eps,Dp,r)
"""
    Interaction between patches

    Input parameters
    eps ---- Energy of the interaction
    Dp ----- Patch diameter
    r ------ Distance between patches

    Derived parameters:
    r_c ---- Cut-off distance
"""
    r_c=1.5*Dp;
    if r<r_c
        return 2*eps*( (1/2)*(Dp/r)^4 - 1 )*exp( (Dp/(r-r_c)) + 2 )
    else
        return 0.0
    end
end

function U3(eps_jk,eps,Dp,r)
"""
    Potencial U3, which is an auxiliary function for the swap potential

    Input parameters
    eps_jk -- Interaction energy between particle j and particle k
    eps ----- Interaction energy between particle i and particle j/k
    Dp ------ Patch diameter
    r ------- Distance between patches

    Derived parameters
    r_min --- Distance such that Upatch(eps,Dp,r_min) = -eps
    r_c ----- Cut-off distance of Upatch
"""
    r_min=Dp;
    r_c=1.5*Dp;
    if r>0 && r<=r_min
        return 1.0
    elseif r>r_min && r<r_c
        return -Upatch(eps,Dp,r)/eps_jk
    else
        return 0.0
    end
end

function Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij,r_ik)
"""
    Swap potential of 3 body interactions

    Input parameters
    w ------- Interpolate the limits of swapping (w=1 swap, w>>1 no swap)
    eps_jk -- Interaction energy between particle j and particle k
    eps_ik -- Interaction energy between particle i and particle k
    eps_ij -- Interaction energy between particle i and particle j
    Dp ------ Patch diameter
    r_ij ---- Distance between patch i and patch j
    r_ik ---- Distance between patch i and patch k
"""
    return w*eps_jk*U3(eps_jk,eps_ij,Dp,r_ij)*U3(eps_jk,eps_ik,Dp,r_ik)
end

function forcePatch(eps,Dp,r)
"""
    Magnitude and components of the force given by the patch-patch interaction potential

    Input parameters
    eps ---- Energy of the interaction
    Dp ----- Patch diameter
    r ------ Distance between patches

    Derived parameters:
    r_c ---- Cut-off distance

    Output
    comp --- Component of the force
    map ---- Magnitude of the force

"""
    r_c=1.5*Dp;
    if r>=r_c
        comp=0;
    else
        comp=(1/r - 2*(eps/r)*exp( Dp/(r-r_c) + 2 ) - ( Dp/(r-r_c)^2 ) )*Upatch(eps,Dp,r)
    end
    mag=sqrt(comp^2)
    return (comp,mag) 
end

function forceSwap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij,r_ik)
"""
    Magnitude and components of the force given by the swap 3 body interaction potential

    Input parameters
    w ------- Interpolate the limits of swapping (w=1 swap, w>>1 no swap)
    eps_jk -- Interaction energy between particle j and particle k
    eps_ik -- Interaction energy between particle i and particle k
    eps_ij -- Interaction energy between particle i and particle j
    Dp ------ Patch diameter
    r_ij ---- Distance between patch i and patch j
    r_ik ---- Distance between patch i and patch k

    Output
    comp_ij - ij component of the force
    comp_ik - ik component of the force
    mag ----- The norm of the force
"""
    comp_ij=w*eps_jk*U3(eps_jk,eps_ik,Dp,r_ik)*first(forcePatch(eps_ij,Dp,r_ij))
    comp_ik=w*eps_jk*U3(eps_jk,eps_ij,Dp,r_ij)*first(forcePatch(eps_ik,Dp,r_ik))
    mag=sqrt( comp_ij^2 + comp_ik^2 )
    return [comp_ij comp_ik mag]
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


