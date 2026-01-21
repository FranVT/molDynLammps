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
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);

# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
data=extractFixScalar(DIR...,FILE_NAME);

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

function fig_Temp(cm, pt, dt,DATA,T)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Temperature}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
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

    return fig

end

function fig_Eng(cm, pt, dt,DATA)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle U \rangle~[J/\epsilon]",
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

plot!(ax,dt.*DATA.TimeStep,DATA."c_ep",label=L"\mathrm{ep}")
plot!(ax,dt.*DATA.TimeStep,DATA."c_ek",label=L"\mathrm{ek}")

#hlines!(ax,[T])

Legend(fig[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

    return fig

end



fig_eng=fig_Eng(cm, pt, dt,DATA)

#fig_temp=fig_Temp(cm, pt, dt,DATA,T)

#save(joinpath(DIR,string("temp-",id_c,".png")), fig, px_per_unit = 300/inch)

function fig_EngB(cm, pt, dt,DATA)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy~Bonds}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle U \rangle~[J/\epsilon]",
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

plot!(ax,dt.*DATA.TimeStep,DATA."v_eB",label=L"\mathrm{eB}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eA",label=L"\mathrm{eA}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eM",label=L"\mathrm{eM}")


#hlines!(ax,[T])

Legend(fig[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

    return fig

end


fig_engB=fig_EngB(cm, pt, dt,DATA)


function fig_EngSys(cm, pt, dt,DATA)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy~System}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle U \rangle~[J/\epsilon]",
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

plot!(ax,dt.*DATA.TimeStep,DATA."v_eT",label=L"\mathrm{eT}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_ec",label=L"\mathrm{ecouple}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eC",label=L"\mathrm{eConserve}")

#hlines!(ax,[T])

Legend(fig[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

    return fig

end


fig_engSys=fig_EngSys(cm, pt, dt,DATA)


