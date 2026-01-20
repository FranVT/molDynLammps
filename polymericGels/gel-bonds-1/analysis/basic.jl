"""
    Basic analysis script
"""

using DataFrames, CSV
#using Plots, LaTeXStrings, Plots.PlotMeasures
#gr()
using Statistics
using GLMakie, LaTeXStrings, Typst_jll

# Load the functions
include("functions.jl")


# Parameter to select the system
T=0.05;
N_particles=500;
phi=0.3;
CL_con=0.05;

# Selection of an specific simulation
date="2026-01-20-111704";

# Get the directory of the desire system
(first(DIR),id_c)=getDir(T,N_particles,phi,CL_con,date);

# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
data=extractFixScalar(DIR,FILE_NAME);

# Convert the array into a DataFrame
DATA=DataFrame(data[2]',data[1]);

"""
    P L O T S 
"""

# these are relative to 1 CSS px
inch = 96
pt = 4/3
cm = inch / 2.54


damp=1;
dt=0.001;

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Temperature}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"\mathrm{Time~units}~[LJ]",
    ylabel=L"\langle T \rangle~[kBT/\epsilon]",
    titlesize=24.0f0,
    subtitlesize=20.0f0,
    xticklabelsize=18.0f0,
    yticklabelsize=18.0f0,
    xlabelsize=22.0f0,
    ylabelsize=22.0f0,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(0,nothing,0,nothing), 
    yticks = 0:0.01:1.5*T
    #xticks=domain
   )

plot!(ax,dt.*DATA.TimeStep,DATA."c_t")
hlines!(ax,[T])

save(joinpath(DIR,string("temp-",id_c,".png")), fig, px_per_unit = 300/inch)
