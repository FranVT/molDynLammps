"""
    Script to tabulate the table of the 3body potential
"""

include("parameters.jl")
include("functions.jl")
include("auxFunctions.jl")

filename1 = string("swapMechTab1new_w",w,".table");
filename2 = string("swapMechTab2new_w",w,".table");

# Create the domains of evaluation according filename nessetities
th_dom = range(thi,thf,2*N);
r_dom = range(rmin,rmax,N);

doms_1=map(rij->map(rik->map(th->[rij rik th],th_dom),r_dom),r_dom);
doms_1=reduce(vcat,reduce(vcat,doms_1));

doms_2=map(rij->map(rik->map(th->[r_dom[rij] r_dom[rik] th],th_dom),rij:N),eachindex(r_dom));
doms_2=reduce(vcat,reduce(vcat,doms_2));

# Create tuples with the information

# For different spicies
docs1 =  map(eachindex(doms_1)) do s
            (
                 s,
                 doms_1[s]...,
                 forceSwapTable(w,eps_ij,eps_ik,eps_jk,sig,doms_1[s]...)...
            )
        end;

# For identical species
docs2 =  map(eachindex(doms_2)) do s
            (
                 s,
                 doms_2[s]...,
                 forceSwapTable(w,eps_ij,eps_ik,eps_jk,sig,doms_2[s]...)...
            )
        end;

createSwapTable(N,rmin,rmax,docs1,filename1)
createSwapTable(N,rmin,rmax,docs2,filename2)

nothing
