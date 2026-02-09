"""
    Script to create the table filename parameter for Patch-PAtch interaction
"""

filename = "pachTab.table";

# Create the functions
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

function DiffEval(eps_pair,sig_p,r)
"""
    Get the central finite difference given the value of the position and the function.
"""
    dh=1e-6;
    fo=Upatch(eps_pair,sig_p,r+dh)
    ff=Upatch(eps_pair,sig_p,r-dh)
    return (1/(2*dh))*( fo - ff );
end



# Parameters
M=12;
N = 2^M;
sig = 0.4;
eps = 1;
rc=1.5*sig;
rmin = sig/10;
rmax = 2*sig;
r_dom = range(rmin,rmax,length=N);

# Create the table
#info = map(s->(s,r_dom[s],Upatch(eps,sig,r_dom[s]),Fpatch(eps,sig,rc,r_dom[s])),eachindex(r_dom));

info = map(s->(s,r_dom[s],Upatch(eps,sig,r_dom[s]),-DiffEval(eps,sig,r_dom[s])),eachindex(r_dom));


# Start to write the data file
    touch(filename); # Create the file

    # Edit the file
    open(filename,"w") do f
        write(f,"DATE: 2025-02-09 UNITS: lj CONTRIBUTOR: Fco.\n\n\n")
        write(f,"POT\n")
        write(f,string("N ",N," RSQ ",rmin," ",rmax,"\n\n"))
        map(t->write(f,rstrip(join(map(s->s*" ",string.(info[t]))))*"\n" ),eachindex(info))
    end;

fig = Figure();
ax = Axis(fig[1,1], limits=(nothing,nothing,-2,2));
scatter!(ax, map(s->info[s][2],eachindex(info)),map(s->info[s][3],eachindex(info)) ,markersize=3);

