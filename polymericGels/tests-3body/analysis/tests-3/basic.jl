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

# Selection of an specific simulation
date="2026-02-11-102800";
#"2026-02-11-130434";
#"2026-02-11-125220";
#

# Get the directory of the desire system
DIR=getDir(date);
DIR=DIR[1];

#table1=getTable(DIR,"pachTab.table");


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
DATA_dump_3 = l[l.id.==3.0,:];

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

tl_sz=0.55cm;
ot_sz=0.35cm;



fig_U=Figure(size = (15cm, 18.75cm));
# Create twin axis: Energy and distance

# First axis
ax1_eng=Axis(fig_U[1,1],
             title=latexstring("\\mathrm{Particle~1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax1_dist=Axis(fig_U[1,1],
             #title=latexstring("\\mathrm{Potential~energy}"),
             #xlabel=L"t~[\tau]",
             yticks = (0.0:0.2:2.0),
             ylabel=L"d~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             yaxisposition=:right
   )
hidespines!(ax1_dist)
hidexdecorations!(ax1_dist)

# Second axis
ax2_eng=Axis(fig_U[2,1],
             title=latexstring("\\mathrm{Particle~2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_dist=Axis(fig_U[2,1],
             #title=latexstring("\\mathrm{Potential~energy}"),
             #xlabel=L"t~[\tau]",
             yticks = (0.0:0.2:2.0),
             ylabel=L"d~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             yaxisposition=:right
   )
hidespines!(ax2_dist)
hidexdecorations!(ax2_dist)

# Third axis
ax3_eng=Axis(fig_U[3,1],
             title=latexstring("\\mathrm{Particle~3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_dist=Axis(fig_U[3,1],
             #title=latexstring("\\mathrm{Potential~energy}"),
             #xlabel=L"t~[\tau]",
             yticks = (0.0:0.2:2.0),
             ylabel=L"d~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             yaxisposition=:right
   )
hidespines!(ax3_dist)
hidexdecorations!(ax3_dist)


dist_12=sqrt.((DATA_dump_1.x .- DATA_dump_2.x).^2 .+ (DATA_dump_1.y .- DATA_dump_2.y).^2);
dist_13=sqrt.((DATA_dump_1.x .- DATA_dump_3.x).^2 .+ (DATA_dump_1.y .- DATA_dump_3.y).^2);
dist_23=sqrt.((DATA_dump_2.x .- DATA_dump_3.x).^2 .+ (DATA_dump_2.y .- DATA_dump_3.y).^2);



sl1_d12=scatterlines!(ax1_dist,time,dist_12,markersize=3,color=1,colormap=:tab10,colorrange=(1,10))
sl1_d13=scatterlines!(ax1_dist,time,dist_13,markersize=3,color=2,colormap=:tab10,colorrange=(1,10))
sl1_eg1=scatterlines!(ax1_eng,time,DATA_dump_1.c_potAtom,markersize=3,color=3, colormap=:tab10,colorrange=(1,10))

sl2_d12=scatterlines!(ax2_dist,time,dist_12,markersize=3,color=1,colormap=:tab10,colorrange=(1,10))
sl2_d23=scatterlines!(ax2_dist,time,dist_23,markersize=3,color=2,colormap=:tab10,colorrange=(1,10))
sl2_eg2=scatterlines!(ax2_eng,time,DATA_dump_2.c_potAtom,markersize=3,color=3, colormap=:tab10,colorrange=(1,10))

sl3_d12=scatterlines!(ax3_dist,time,dist_13,markersize=3,color=1,colormap=:tab10,colorrange=(1,10))
sl3_d23=scatterlines!(ax3_dist,time,dist_23,markersize=3,color=2,colormap=:tab10,colorrange=(1,10))
sl3_eg3=scatterlines!(ax3_eng,time,DATA_dump_3.c_potAtom,markersize=3,color=3, colormap=:tab10,colorrange=(1,10))


Legend(fig_U[1,2],
       [sl1_d12,sl1_d13,sl1_eg1],
       [L"d_{12}", L"d_{13}", L"U_{\mathrm{atom}}"],
     labelsize=0.5cm)

Legend(fig_U[2,2],
       [sl2_d12,sl2_d23,sl2_eg2],
       [L"d_{21}", L"d_{23}", L"U_{\mathrm{atom}}"],
     labelsize=0.5cm)

Legend(fig_U[3,2],
       [sl3_d12,sl3_d23,sl3_eg3],
       [L"d_{31}", L"d_{32}", L"U_{\mathrm{atom}}"],
     labelsize=0.5cm)


save("fig_Uatom.png", fig_U, px_per_unit = 300/inch)


"""
    Energy of the fix
"""

fig_Uf=Figure(size = (18.75cm, 15cm));
# Create twin axis: Energy and distance

# First axis
ax1_eng=Axis(fig_Uf[1,1],
             title=latexstring("\\mathrm{Potential~energy}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )

time_fix=DATA_fix.TimeStep.*0.001;
time_dump=DATA_dump_1.TimeStep.*0.001;

aux=combine(groupby(l, :TimeStep), :c_potAtom => sum => :suma_c_pot);

scatter!(ax1_eng,time_fix,DATA_fix.c_patchPair,label=L"U_{\mathrm{fix}}",color=1,colormap=:tab10,colorrange=(1,10))
lines!(ax1_eng,time_dump,DATA_dump_1.c_potAtom,label=L"U_{\mathrm{atom}1}",color=2,colormap=:tab10,colorrange=(1,10))
lines!(ax1_eng,time_dump,DATA_dump_2.c_potAtom,label=L"U_{\mathrm{atom}2}",color=3,colormap=:tab10,colorrange=(1,10))
lines!(ax1_eng,time_dump,DATA_dump_3.c_potAtom,label=L"U_{\mathrm{atom}3}",color=4,colormap=:tab10,colorrange=(1,10))
lines!(ax1_eng,time_dump,aux.suma_c_pot,label=L"\sum U_{\mathrm{atom}}",color=:black)



# jj = 

Legend(fig_Uf[1,2],ax1_eng,
      L"\mathrm{Labels}",
     labelsize=0.5cm)

save("fig_Uatom-fix.png", fig_Uf, px_per_unit = 300/inch)


"""
    Energy from function evaluated at the distances given by the simulation 
"""
fig_UF=Figure(size = (18.75cm, 15cm));
# Create twin axis: Energy and distance
ax2_eng=Axis(fig_UF[1,2],
             title=latexstring("\\mathrm{Comparisson~of~fix~and~function}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_dist=Axis(fig_UF[1,2],
             #title=latexstring("\\mathrm{Potential~energy}"),
             #xlabel=L"t~[\tau]",
             yticks = (0.0:0.2:2.0),
             ylabel=L"d~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             yaxisposition=:right
   )
hidespines!(ax2_dist)
hidexdecorations!(ax2_dist)

scatter!(ax2_eng,time_fix,DATA_fix.c_patchPair,label=L"U_{\mathrm{fix}}",color=7,colormap=:tab10,colorrange=(1,10))
lines!(ax2_eng,time_dump,map(r->Upatch(1,0.4,r),dist_12),label=L"U_{\mathrm{patch}}(d_{12})",color=4,colormap=:tab10,colorrange=(1,10))
lines!(ax2_eng,time_dump,map(r->Upatch(1,0.4,r),dist_13),label=L"U_{\mathrm{patch}}(d_{13})",color=5,colormap=:tab10,colorrange=(1,10))
lines!(ax2_eng,time_dump,map(r->Upatch(1,0.4,r),dist_23),label=L"U_{\mathrm{patch}}(d_{23})",color=6,colormap=:tab10,colorrange=(1,10))
lines!(ax2_eng,time_dump,map(r->Upatch(1,0.4,r),dist_12).+map(r->Upatch(1,0.4,r),dist_13).+map(r->Upatch(1,0.4,r),dist_23),label=L"\sum U_{\mathrm{patch}}(r_{ij})",color=:black)


lines!(ax2_dist,time_dump,dist_12,color=1,colormap=:tab10,colorrange=(1,10),label=L"d_{12}")
lines!(ax2_dist,time_dump,dist_13,color=2,colormap=:tab10,colorrange=(1,10),label=L"d_{13}")
lines!(ax2_dist,time_dump,dist_23,color=3,colormap=:tab10,colorrange=(1,10),label=L"d_{23}")

Legend(fig_UF[1,1],ax2_eng,
      L"\mathrm{Labels}",
     labelsize=0.5cm)
Legend(fig_UF[1,3],ax2_dist,
      L"\mathrm{Labels}",
     labelsize=0.5cm)

save("fig_Ufix-func.png", fig_UF, px_per_unit = 300/inch)


#=

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
=#
