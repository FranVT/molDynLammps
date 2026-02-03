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
date=;
#"2026-02-03-170509";
#"2026-02-03-164806";
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


"""
    SANDBOX
"""

# these are relative to 1 CSS px
inch = 96;
pt = 4/3;
cm = inch / 2.54;

# Concatenate the data frames into one
l=reduce(vcat,DATA_dump);

# Get the positions of one particle
# Plot
fig_pos=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig_pos[1:1,1:1],
    title=latexstring("\\mathrm{Position~of~particles}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"x~[x/D]",
    ylabel=L"y~[y/D]",
    titlesize=24.0f0,
    subtitlesize=20.0f0,
    xticklabelsize=18.0f0,
    yticklabelsize=18.0f0,
    xlabelsize=22.0f0,
    ylabelsize=22.0f0,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

scatterlines!(ax,l[l.id.==2.0,:].x,l[l.id.==2.0,:].y,label=L"2")
scatterlines!(ax,l[l.id.==3.0,:].x,l[l.id.==3.0,:].y,label=L"3")
scatterlines!(ax,l[l.id.==4.0,:].x,l[l.id.==4.0,:].y,label=L"4")

#hlines!(ax,[T])

Legend(fig_pos[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

# Time evolution of the components of the position of each particle
# Plot
fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax_1=Axis(fig[1,1],
    title=latexstring("\\mathrm{Time~evolution~Position~of~particles}"),
    #subtitle=latexstring(subtitle),
    ylabel=L"x|y~[r/D]",
    xlabel=L"\mathrm{Time~step}",
    titlesize=24.0f0,
    subtitlesize=20.0f0,
    xticklabelsize=18.0f0,
    yticklabelsize=18.0f0,
    xlabelsize=22.0f0,
    ylabelsize=22.0f0,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

ax_2=Axis(fig[2,1],
    #subtitle=latexstring(subtitle),
    ylabel=L"x|y~[r/D]",
    xlabel=L"\mathrm{Time~step}",
    titlesize=24.0f0,
    subtitlesize=20.0f0,
    xticklabelsize=18.0f0,
    yticklabelsize=18.0f0,
    xlabelsize=22.0f0,
    ylabelsize=22.0f0,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

ax_3=Axis(fig[3,1],
    #subtitle=latexstring(subtitle),
    ylabel=L"x|y~[r/D]",
    xlabel=L"\mathrm{Time~step}",
    titlesize=24.0f0,
    subtitlesize=20.0f0,
    xticklabelsize=18.0f0,
    yticklabelsize=18.0f0,
    xlabelsize=22.0f0,
    ylabelsize=22.0f0,
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(nothing,nothing,nothing,nothing), 
    #yticks = 0:0.01:1.5*T
    #xticks=domain
   )

scatterlines!(ax_1,l[l.id.==2.0,:].x,label=L"2.x",markersize=10)
scatterlines!(ax_1,l[l.id.==2.0,:].y,label=L"2.y",markersize=10)

scatterlines!(ax_2,l[l.id.==3.0,:].x,label=L"3.x",markersize=10)
scatterlines!(ax_2,l[l.id.==3.0,:].y,label=L"3.y",markersize=10)

scatterlines!(ax_3,l[l.id.==4.0,:].x,label=L"4.x",markersize=10)
scatterlines!(ax_3,l[l.id.==4.0,:].y,label=L"4.y",markersize=10)

Legend(fig[1:3,2],ax_3,
      L"\mathrm{Legend}",
     labelsize=12pt)

