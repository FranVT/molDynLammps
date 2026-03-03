"""
    Script with auxiliary functions to create the tables
"""

function createSwapTable(N,rmin,rmax,info,filename)
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

function DiffEval(eps_pair,sig_p,r)
"""
    Get the central finite difference given the value of the position and the function.
"""
    dh=1e-6;
    fo=Upatch(eps_pair,sig_p,r+dh)
    ff=Upatch(eps_pair,sig_p,r-dh)
    return (1/(2*dh))*( fo - ff );
end


