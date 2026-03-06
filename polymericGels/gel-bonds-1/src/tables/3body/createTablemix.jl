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
    U_patch + U_swap
    F_patch + F_swap
"""

    # General parameters
    th = deg2rad(th);
    r_jk = sqrt(r_ij^2+r_ik^2-2*r_ij*r_ik*cos(th));

    # Two particle interactions
    p_ij=first(forcePatch(eps_ij,0.4,r_ij));
    p_ik=first(forcePatch(eps_ik,0.4,r_ik));

    p_ji=-p_ij;
    p_jk=first(forcePatch(eps_jk,0.4,r_jk));

    p_ki=-p_ik;
    p_kj=-p_jk;

    eng_2body=Upatch(eps_ij,0.4,r_ij) + Upatch(eps_ik,0.4,r_ik) + Upatch(eps_jk,0.4,r_jk);

     
    # Three particle interaction
    f_1=forceSwap(w,eps_jk,eps_ij,eps_ik,sig_p,r_ij,r_ik);
    f_2=forceSwap(w,eps_ik,eps_ij,eps_jk,sig_p,r_ij,r_jk);
    f_3=forceSwap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ik,r_jk);

        a=f_1[1]; b=-f_2[1]; c=a+b;
        f_swap_i1=c;
        
        a=f_1[2]; b=-f_3[1]; c=a+b;
        f_swap_i2=c;
    
        f_swap_j1=-f_swap_i1; 

        a=f_2[2]; b=-f_3[2]; c=a+b;
        f_swap_j2=c; 

        f_swap_k1=-f_swap_i2; 
        f_swap_k2=-f_swap_j2; 

    eng_3body=Uswap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_ik) + Uswap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ij,r_jk) + Uswap(w,eps_ij,eps_ik,eps_jk,sig_p,r_ik,r_jk);
        
    # Addition of the two particle and three particle interactions
    f_i1=p_ij+f_swap_i1;
    f_i2=p_ik+f_swap_i2;
    
    f_j1=p_ji+f_swap_j1; 
    f_j2=p_jk+f_swap_j2; 

    f_k1=p_ki+f_swap_k1; 
    f_k2=p_kj+f_swap_k2; 

    eng=eng_2body+eng_3body
    
    return (f_i1,f_i2,f_j1,f_j2,f_k1,f_k2,eng)

end

## Parameters for the file

N = 2^7; # N=100=2^6.643856

eps_ij = 1.0;
eps_ik = 1.0;
eps_jk = 1.0;
sig = 0.4;
rc = 1.5*sig;
rmin = sig/2;
rmax = 1.5*sig;
thi = 180/(4*N)
thf = 180 - thi;
w=1;

filename1 = string("swapMechTab1mix_w",w,".table");
filename2 = string("swapMechTab2mix_w",w,".table");

# Create the domains of evaluation according filename nessetities
th_dom = range(thi,thf,2*N);
r_dom = range(rmin,rmax,N);

#doms1 = reduce(vcat,reverse.(Iterators.product(th_dom,r_dom,r_dom)|>collect));
#doms2 = reduce(vcat,map(s-> reshape(reverse.(Iterators.product(th_dom,r_dom[s:end],r_dom[s])),2*N*(N-(s-1)),1) ,eachindex(r_dom)));

doms_1=map(rij->map(rik->map(th->[rij rik th],th_dom),r_dom),r_dom);
doms_1=reduce(vcat,reduce(vcat,doms_1));

doms_2=map(rij->map(rik->map(th->[r_dom[rij] r_dom[rik] th],th_dom),rij:N),eachindex(r_dom));
doms_2=reduce(vcat,reduce(vcat,doms_2));



# Create tuples with the information
# For Patch_j = Patch_k

# For different spicies
docs1 =  map(eachindex(doms_1)) do s
            (
                 s,
                 doms_1[s]...,
                 force(w,eps_ij,eps_ik,eps_jk,sig,doms_1[s]...)...
            )
        end;

# For identical species
docs2 =  map(eachindex(doms_2)) do s
            (
                 s,
                 doms_2[s]...,
                 force(w,eps_ij,eps_ik,eps_jk,sig,doms_2[s]...)...
            )
        end;

createTable(N,rmin,rmax,docs1,filename1)
createTable(N,rmin,rmax,docs2,filename2)

nothing
