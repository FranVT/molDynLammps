"""
    Script to create the table filename parameter for Patch-PAtch interaction
"""
filename = "pachTab.table";

include("parameters.jl")
include("function.jl")
include("auxFunctions.jl")

r_dom = range(rmin,rmax,length=N);

# Create the table
#info = map(s->(s,r_dom[s],Upatch(eps,sig,r_dom[s]),Fpatch(eps,sig,rc,r_dom[s])),eachindex(r_dom));

info = map(s->(s,r_dom[s],Upatch(eps,sig,r_dom[s]),first(forcePatch(eps,sig,r_dom[s]))),eachindex(r_dom));

# Start to write the data file
    touch(filename); # Create the file

    # Edit the file
    open(filename,"w") do f
        write(f,"DATE: 2025-02-09 UNITS: lj CONTRIBUTOR: Fco.\n\n\n")
        write(f,"POT\n")
#        write(f,string("N ",N," RSQ ",rmin," ",rmax,"\n\n"))
        write(f,string("N ",N,"\n\n"))
        map(t->write(f,rstrip(join(map(s->s*" ",string.(info[t]))))*"\n" ),eachindex(info))
    end;

