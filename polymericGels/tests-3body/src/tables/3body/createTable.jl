"""
    Script to tabulate the table of the 3body potential
"""

include("functions.jl")

function createTable(N,rmin,rmax,info,filename)
"""
    Create the table for threebody/table filename input as follows:
        ind rij rik th fi1 fi2 fj1 fj2 fk1 fk2 e
        ind rij rik th fi1 fi2 -fi1 0 -fi2 0 e
        ind rij rik th fij fik -fij 0 -fik 0 e
"""

    # Start to write the data file
    touch(filename); # Create the file

    # Edit the file
    open(filename,"w") do f
        write(f,"SEC1\n")
        write(f,string("N ",N," rmin ",rmin," rmax ",rmax,"\n\n"))
        map(t->write(f,rstrip(join(map(s->s*" ",string.(info[t]))))*"\n" ),eachindex(info))
    end

    nothing
end

function force(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_ik,th)
"""
    Compute the scalars for the proyection of the forces
"""
    th = deg2rad(th);
    force=forceSwap(w,eps_jk,eps_ij,eps_ik,sig_p,r_ij,r_ik);
    f_ij=force[1];
    f_ik=force[2];
    #r_jk = sqrt(r_ij^2+r_ik^2-2*r_ij*r_ik*cos(th));
 
    f_i1=f_ij; 
    f_i2=f_ik; 
   
    f_j1=-f_i1; 
    f_j2=0; 

    f_k1=-f_i2; 
    f_k2=-f_j2; 

    eng=Uswap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_ik) 
    #+ SwapU(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_jk) + SwapU(w,eps_ij,eps_ik,eps_jk,sig_p,r_ik,r_jk)
    #eng=round(eng/3,digits=2^7)

    return (f_i1,f_i2,f_j1,f_j2,f_k1,f_k2,eng)

end



## Parameters for the file

N = 100;

eps_ij = 1.0;
eps_ik = 1.0;
eps_jk = 1.0;
sig = 0.4;
rc = 1.5*sig;
rmin = sig/10;
rmax = 2*sig;
thi = 180/(4*N)
thf = 180 - thi;
w=1;

filename1 = string("swapMechTab1new_w",w,".table");
filename2 = string("swapMechTab2new_w",w,".table");

# Create the domains of evaluation according filename nessetities
th_dom = range(thi,thf,2*N);
r_dom = range(rmin,rmax,N);

doms1 = reduce(vcat,reverse.(Iterators.product(th_dom,r_dom,r_dom)|>collect));
doms2 = reduce(vcat,map(s-> reshape(reverse.(Iterators.product(th_dom,r_dom[s:end],r_dom[s])),2*N*(N-(s-1)),1) ,eachindex(r_dom)));

# Create tuples with the information
# For Patch_j = Patch_k

# For different spicies
docs1 =  map(eachindex(doms1)) do s
            (
                 s,
                 doms1[s]...,
                 force(w,eps_ij,eps_ik,eps_jk,sig,doms1[s]...)...
            )
        end;

# For identical species
docs2 =  map(eachindex(doms2)) do s
            (
                 s,
                 doms2[s]...,
                 force(w,eps_ij,eps_ik,eps_jk,sig,doms2[s]...)...
            )
        end;

createTable(N,rmin,rmax,docs1,filename1)
createTable(N,rmin,rmax,docs2,filename2)

nothing
