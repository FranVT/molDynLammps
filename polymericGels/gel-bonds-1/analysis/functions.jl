"""
    Functions 
"""

function getDump(dir,file_name)
"""
    Get the data from a single dump file that stores one timeste information
"""
    data = split.(readlines(joinpath(dir,file_name))," ")[9:end];
    HEADERS=data[1][3:end];
    INFO=parse.(Float64,reduce(hcat,data[2:end]))';

    return DataFrame(INFO,HEADERS)
end

function fig_EngPair(cm, pt, dt,DATA)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy~Interactions}"),
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

plot!(ax,dt.*DATA.TimeStep,DATA."c_wcaPair",label=L"\mathrm{WCA}")
plot!(ax,dt.*DATA.TimeStep,DATA."c_patchPair",label=L"\mathrm{patch}")
plot!(ax,dt.*DATA.TimeStep,DATA."c_swapPair",label=L"\mathrm{swap}")

#hlines!(ax,[T])

Legend(fig[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

    return fig

end



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

function fig_EngB(cm, pt, dt,DATA)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy~Bonds}"),
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
   )

plot!(ax,dt.*DATA.TimeStep,DATA."v_eB",label=L"\mathrm{eB}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eA",label=L"\mathrm{eA}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eM",label=L"\mathrm{eM}")


Legend(fig[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

    return fig

end



function fig_Eng(cm, pt, dt,DATA)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy}"),
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
   )

plot!(ax,dt.*DATA.TimeStep,DATA."c_ep",label=L"\mathrm{ep}")
plot!(ax,dt.*DATA.TimeStep,DATA."c_ek",label=L"\mathrm{ek}")


Legend(fig[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

    return fig

end

function fig_Temp(cm, pt, dt,DATA,T)

fig=Figure(size = (15cm, 12cm), fontsize = 12pt);

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Temperature}"),
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
   )

plot!(ax,dt.*DATA.TimeStep,DATA."c_t")
hlines!(ax,[T])

    return fig

end

function extractFixScalar(path_system,file_name)
"""
    Function that extracts the information of fix files that stores global scalar values
"""
    aux=split.(readlines(joinpath(path_system,file_name))," ");
    return (aux[2][2:end],reduce(hcat,map(s->parse.(Float64,s),aux[3:end])));
end



function getDir(T,N_particles,phi,CL_con,date)
"""
    To get the directory of the simulation
"""

# Get the directories
MAIN_DIR=dirname(pwd());
DATA_DIR=joinpath(MAIN_DIR,"data");

# Read the directory where data is stored
sims=filter(isdir,readdir(DATA_DIR,join=true));

# Creation of the id
id=string(T,phi,CL_con,N_particles);

# For this script, evetually we are going to get assembly averages
id_c=string(T,phi,CL_con,N_particles,"-",date); # filter!(s->s==id_c,readdir(DATA_DIR))

# Get the idex for the directory
indx=findall(!isempty,findall.(id_c,sims));

# Get the directory
DIR=sims[indx];

    return DIR, id_c

end


