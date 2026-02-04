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
date="2026-02-04-155132";
#"2026-02-04-102545";
#2026-02-04-123638
#"2026-02-04-102545";
#"2026-02-02-155528";


# Get the directory of the desire system
DIR=getDir(date);
DIR=DIR[1];

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
    SANDBOX
"""

# these are relative to 1 CSS px
inch = 96;
pt = 4/3;
cm = inch / 2.54;

# Concatenate the data frames into one
l=reduce(vcat,DATA_dump);

L=1;

"""
    Energy and that stuff
"""

time=dt.*DATA_fix.TimeStep;

U_2=l[l.id.==2.0,:].c_potAtom;
U_3=l[l.id.==3.0,:].c_potAtom;

fig_U=Figure(size = (18.75cm, 15cm));

clbr=:managua10;

ax=Axis(fig_U[1:1,1:1],
    title=latexstring("\\mathrm{Potential~energy}"),
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
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

scatterlines!(ax,time,U_2,markersize=3,label=L"2")
scatterlines!(ax,time,U_3,markersize=3,label=L"3")

scatterlines!(ax,time,DATA_fix.c_patchPair,markersize=3,label=L"\mathrm{Patch}")
scatterlines!(ax,time,DATA_fix.c_swapPair,markersize=3,label=L"\mathrm{Swap}")
scatterlines!(ax,time,DATA_fix.c_ep,markersize=3,label=L"\mathrm{System}")

#hlines!(ax,[T])

Legend(fig_U[1,2],ax,
      L"\mathrm{id}",
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
    limits=(nothing,nothing,nothing,nothing), 
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

ms=5;

# Get the positions of one particle
# Plot
fig_pos=Figure(size = (18.75cm, 15cm));

clbr=:managua10;

ax=Axis(fig_pos[1:1,1:1],
    title=latexstring("\\mathrm{Position~of~particles}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"x~[x/D]",
    ylabel=L"y~[y/D]",
    titlesize=1cm,
    xticklabelsize=0.5cm,
    yticklabelsize=0.5cm,
    xlabelsize=1cm,
    ylabelsize=1cm,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(-L,L,-L,L), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

scatterlines!(ax,l[l.id.==2.0,:].x,l[l.id.==2.0,:].y,markersize=ms,label=L"2")
scatterlines!(ax,l[l.id.==3.0,:].x,l[l.id.==3.0,:].y,markersize=ms,label=L"3")
scatterlines!(ax,l[l.id.==4.0,:].x,l[l.id.==4.0,:].y,markersize=ms,label=L"4")

#hlines!(ax,[T])

Legend(fig_pos[1,2],ax,
      L"\mathrm{id}",
     labelsize=0.5cm)

# Time evolution of the components of the position of each particle
# Plot

fig=Figure(size = (18.75cm, 15cm));


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


scatterlines!(ax_1,time,l[l.id.==2.0,:].x,label=L"2",markersize=ms)
scatterlines!(ax_1,time,l[l.id.==3.0,:].x,label=L"3",markersize=ms)

scatterlines!(ax_2,time,l[l.id.==2.0,:].y,label=L"2",markersize=ms)
scatterlines!(ax_2,time,l[l.id.==3.0,:].y,label=L"3",markersize=ms)



Legend(fig[1:2,2],ax_2,
      L"\mathrm{id}",
     labelsize=0.5cm)

