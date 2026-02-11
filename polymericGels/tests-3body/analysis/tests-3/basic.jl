"""
    Basic analysis
"""

using DataFrames, CSV
#using Plots, LaTeXStrings, Plots.PlotMeasures
#gr()
using Statistics
using GLMakie, LaTeXStrings, Typst_jll

# Load the functions
include("functions.jl")

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
    Get the data from a single dump file that stores one timeste information
"""
    data=split.(readlines(joinpath(dir,file_name))," ")[7:end];
    HEADERS=["n","r","u","f"];
    INFO=reduce(hcat,map(s->parse.(Float64,s),data))';
    df=DataFrame(INFO,HEADERS);

    return df 
end





# Selection of an specific simulation
date="2026-02-10-114759";

# 2026-02-09-104056
#"2026-02-05-140904";



# Get the directory of the desire system
DIR=getDir(date);
DIR=DIR[1];

table1=getTable(DIR,"pachTab.table");


# Activate extract info
act=1;

# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
if act == 1
    data=extractFixScalar(DIR,FILE_NAME);
end

# Convert the array into a DataFrame
DATA_fix=DataFrame(data[2]',data[1]);

# Get the directory of the desire system
DIR=joinpath(DIR,"traj");

# Get data
if act == 1
    DATA_dump=map(s->getDump(DIR,s),readdir(DIR));
end

# Select parameters to filter dataframe (patches)
# type 3: PA
# type 4: PB


# Parameters
dt=0.001;


"""
#    Reducing the DATA_dump
"""

# Concatenate the data frames into one
l=sort(reduce(vcat,DATA_dump),:TimeStep);

# Create Individual data frames
DATA_dump_1 = l[l.id.==1.0,:];
DATA_dump_2 = l[l.id.==2.0,:];

"""
#    SANDBOX
"""

# these are relative to 1 CSS px
inch = 96;
pt = 4/3;
cm = inch / 2.54;

L=1;

"""
#    Energy and that stuff
"""

time=dt.*DATA_fix.TimeStep;

e_lim=2;

fig_U=Figure(size = (15cm, 18.75cm));
tl_sz=0.55cm;
ot_sz=0.35cm;

clbr=:managua10;

ax1=Axis(fig_U[1,1],
    title=latexstring("\\mathrm{Potential~energy}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"t~[\tau]",
    ylabel=L"U~[J/\epsilon]",
    titlesize=tl_sz,
    xticklabelsize=ot_sz,
    yticklabelsize=ot_sz,
    xlabelsize=tl_sz,
    ylabelsize=tl_sz,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,-e_lim,e_lim), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )
ax2=Axis(fig_U[2,1],
    title=latexstring("\\mathrm{Distance~between~particles}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"t~[\tau]",
    ylabel=L"d~[x/D]",
    titlesize=tl_sz,
    xticklabelsize=ot_sz,
    yticklabelsize=ot_sz,
    xlabelsize=tl_sz,
    ylabelsize=tl_sz,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )
ax3=Axis(fig_U[3,1],
    title=latexstring("\\mathrm{Potential~energy}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"t~[\tau]",
    ylabel=L"U~[J/\epsilon]",
    titlesize=tl_sz,
    xticklabelsize=ot_sz,
    yticklabelsize=ot_sz,
    xlabelsize=tl_sz,
    ylabelsize=tl_sz,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,-e_lim,e_lim), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

dist = DATA_dump_1.x .- DATA_dump_2.x;
scatterlines!(ax2,time,dist,markersize=3)

scatterlines!(ax1,time,DATA_dump_1.c_potAtom,markersize=3,label=L"2")
scatterlines!(ax3,time,DATA_dump_2.c_potAtom,markersize=3,label=L"3")

Legend(fig_U[1,2],ax1,
      L"\mathrm{id}",
     labelsize=0.5cm)
Legend(fig_U[3,2],ax3,
      L"\mathrm{id}",
     labelsize=0.5cm)



fig_Udist=Figure(size = (12cm, 10cm));

clbr=:managua10;

ax=Axis(fig_Udist[1,1],
    title=latexstring("\\mathrm{Potential~energy}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"d~[x/D]",
    ylabel=L"U~[J/\epsilon]",
    titlesize=tl_sz,
    xticklabelsize=ot_sz,
    yticklabelsize=ot_sz,
    xlabelsize=tl_sz,
    ylabelsize=tl_sz,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(0.25,nothing,-e_lim,e_lim), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

N=100;
rmin = 0.4/1000;
rmax = 2*0.4;
r_dom = range(rmin,rmax,length=N);
lines!(ax,r_dom,map(r->Upatch(1,0.4,r),r_dom),label=L"U_{\mathrm{patch}}",linestyle=:solid,color=:black)
lines!(ax,table1.r,table1.u,label=L"U_{\mathrm{table}}",linewidth=5)

lines!(ax,dist,DATA_dump_2.c_potAtom,label=L"U_{\mathrm{atom}}")
lines!(ax,dist,DATA_fix.c_patchPair,label=L"U_{\mathrm{pe}}")

Legend(fig_Udist[1,2],ax,
      L"\mathrm{Labels}",
     labelsize=0.5cm)





# Total Energy
fig_E=Figure(size = (18.75cm, 15cm));

clbr=:managua10;

ax=Axis(fig_E[1:1,1:1],
    title=latexstring("\\mathrm{Total~energy}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"t~[\tau]",
    ylabel=L"U~[J/\epsilon]",
    titlesize=1cm,
    xticklabelsize=0.5cm,
    yticklabelsize=0.5cm,
    xlabelsize=1cm,
    ylabelsize=1cm,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,-e_lim,e_lim), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

scatterlines!(ax,time,DATA_fix.c_ep .+ DATA_fix.c_ek,label=L"\mathrm{U} + \mathrm{K}",markersize=3)
scatterlines!(ax,time,DATA_fix.c_ep,label=L"\mathrm{U}",markersize=3)
scatterlines!(ax,time,DATA_fix.c_ek,label=L"\mathrm{K}",markersize=3)

#hlines!(ax,mean(DATA_fix.c_ep .+ DATA_fix.c_ek))

Legend(fig_E[1,2],ax,
      L"\mathrm{Energy}",
     labelsize=0.5cm)










"""
    Position
"""

fig=Figure(size = (18.75cm, 15cm));

ms=3;

clbr=:managua10;

ax_1=Axis(fig[1,1],
    title=latexstring("\\mathrm{Time~evolution}"),
    #subtitle=latexstring(subtitle),
    ylabel=L"x~[r/D]",
#    xlabel=L"\mathrm{Time~step}",
    titlesize=1cm,
    xticklabelsize=0.5cm,
    yticklabelsize=0.5cm,
    xlabelsize=1cm,
    ylabelsize=1cm,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

ax_2=Axis(fig[2,1],
    #subtitle=latexstring(subtitle),
    ylabel=L"y~[r/D]",
    xlabel=L"t~[\tau]",
    titlesize=1cm,
    xticklabelsize=0.5cm,
    yticklabelsize=0.5cm,
    xlabelsize=1cm,
    ylabelsize=1cm,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

linkxaxes!(ax_1,ax_2)


scatterlines!(ax_1,time,l[l.id.==1.0,:].x,label=L"2",markersize=ms)
scatterlines!(ax_1,time,l[l.id.==2.0,:].x,label=L"3",markersize=ms)

scatterlines!(ax_2,time,l[l.id.==1.0,:].y,label=L"2",markersize=ms)
scatterlines!(ax_2,time,l[l.id.==2.0,:].y,label=L"3",markersize=ms)



Legend(fig[1:2,2],ax_2,
      L"\mathrm{id}",
     labelsize=0.5cm)

save("fig_comp.png", fig, px_per_unit = 300/inch)
save("fig_pot.png", fig_U, px_per_unit = 300/inch)
save("fig_potUComp.png", fig_Udist, px_per_unit = 300/inch)
#save("fig_dist_pot.png", fig_Udist, px_per_unit = 300/inch)

