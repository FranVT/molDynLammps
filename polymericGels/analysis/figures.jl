"""
    Scripts to create figures
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll


include("functions.jl")

# Include the theme of the plots and usefull functions
include("functions_graphs.jl")

# Main directory
MAIN_DIR=pwd();

# Directory of the data already processed
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:id];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Get the information of all systems of the fix files 
dfs_fix = map(s->get_fixInfo(data_bySystem[s],SAVE_DIR),eachindex(data_bySystem));

# Store the figure of the energy 
#figure_fixEnergy(dfs_fix,dat_files)

# Get the information of all systems of the structure factor files
dfs_sf = map(s->get_sfInfo(data_bySystem[s],SAVE_DIR),eachindex(data_bySystem));

# ids
ids = map(s->unique(data_bySystem[s].id)[1],eachindex(data_bySystem));


# Select one system
it_sys=5;
system=dfs_sf[it_sys];

# Get the time steps of the S(q)
time_aux=map(s->unique(s.timeStep)[1],system);

# Sort the time steps 
time_domain=sort(time_aux);

# Get the index from the sorted order to the original file
ind_sort=map(s->findall(x->x==s,time_aux)[1],time_domain);

# Create the graph of the time evolution of the structure factor

# Time domain label
dt=0.001;
time_domain=(dt).*time_domain;

# Labels
to=first(time_domain);
tf=last(time_domain);


# Prepare the color code
color_ref=maximum(time_domain);
color_norm=time_domain./color_ref;
color_min=minimum(color_norm);
color_max=maximum(color_norm);

    fig = Figure()
    ax_plot = Axis(fig[1:1, 1:1],
                   xlabel = L"|\vec{q}|~[\frac{1}{x^*}]",
                   ylabel = L"S(|\vec{q}|)",
                   xminorticksvisible = true,
                   xminorgridvisible = true,
                   limits = (nothing, nothing, nothing, nothing),
                   xscale = log10,
                   yscale = log10
                  )

    # Plot the potential energy for each system
    for (it_sort,it_og) in enumerate(ind_sort)
        scatterlines!(ax_plot, system[it_og].qmean, system[it_og].Sqmean,
                      color = color_norm[it_sort],
                      colorrange = (color_min, color_max)
                     )
    end

    # Legends in terms of the packing fraction
    Colorbar(fig[1, 2], label = L"\tau", colormap = :viridis, limits = (to, tf))

    save(string("fig_Sq_",ids[it_sys],"_timeDomain.png"), fig, px_per_unit = 300 / inch)


#figureCompareSq(dat_files)
#figureCompareSl(dat_files)

#Legend(fig[5,4],ax_plot,L"\phi")
#=

# Grafica del factor de estructura
MAIN_DIR=pwd();
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# File name of interested data
pattern=r"^structureFactor.*\.csv$";

# Create the paths to the files
file_names=filter(f->occursin(pattern,f),readdir(SAVE_DIR))

file_paths=joinpath.(SAVE_DIR,file_names);

ids=unique(map(s->match(r"C(.*?)t", s).captures[1],file_names))

file_paths=[file_paths[occursin.(it,file_paths)] for it in ids]


# Graficas por sistema
dict_phi=Dict(dat_files.id.=>dat_files.phi);
dict_CL=Dict(dat_files.id.=>dat_files."CL-Con");
dict_T=Dict(dat_files.id.=>dat_files.Temperature);

# Select one percentage
id_selec=4;

Sq=[CSV.read(file,DataFrame) for file in file_paths[id_selec]];

timeDomain=convert.(BigInt,first.([s.timeStep for s in Sq]));
t_max=maximum(timeDomain);
t_min=minimum(timeDomain);

timeDomain_norm=timeDomain./t_max;



fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"|\vec{q}|=\frac{2\pi}{L}|\vec{n}|",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(nothing,nothing,nothing,nothing),
        xscale=log10,
        yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)

    #vlines!(ax,1.2,linestyle=:dash,color=:grey)
    #vlines!(ax,5*1.2,linestyle=:dash,color=:grey)

    for it_sys in eachindex(Sq) 
        scatterlines!(ax,Sq[it_sys].qmean,Sq[it_sys].Sqmeannorm,
                      color=timeDomain_norm[it_sys],
                      colorrange=(t_min/t_max,1)
                      #label=latexstring(100*phis[it_sys])
              )
    end



    #Legend(fig[2,4],ax,L"\phi~\%",merge=true)

    p1=plot!(ax_f,[0],[-1],
             label=latexstring(100*dict_phi[ids[id_selec]])
            )
    p1.visible = false
    Legend(fig[5,4],ax_f,L"\phi")
    Colorbar(fig[1:3,4],label=L"\mathrm{Time}",colormap=:viridis,limits=(t_min/1000,t_max/1000))

save(string("fig_Sq_phi_",100*dict_phi[ids[id_selec]],"_timeseries.png"), fig, px_per_unit = 300/inch)

=#


#=

# Dictionaries to map id to info.
dict_phi=Dict(dat_files.id.=>dat_files.phi);
dict_CL=Dict(dat_files.id.=>dat_files."CL-Con");
dict_T=Dict(dat_files.id.=>dat_files.Temperature);


# Get the information
Sq=[CSV.read(file,DataFrame) for file in file_paths];

Sq_1=Sq[1:2:end];
Sq_2=Sq[2:2:end];


phis=map(s->first(data_bySystem[s].phi),(1,2,3,5))





println("Se inicia a graficar")
fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"|\vec{q}|=\frac{2\pi}{L}|\vec{n}|",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(nothing,nothing,nothing,nothing),
        xscale=log10,
        yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)

    #vlines!(ax,1.2,linestyle=:dash,color=:grey)
    #vlines!(ax,5*1.2,linestyle=:dash,color=:grey)

    for it_sys in eachindex(ids) 
        scatterlines!(ax,Sq_2[it_sys].qmean,Sq_2[it_sys].Sqmeannorm,
                      label=latexstring(100*phis[it_sys])
              )
    end

    Legend(fig[2,4],ax,L"\phi~\%",merge=true)

    p1=plot!(ax_f,[0],[-1],
             label=L"t_f"
            )
    p1.visible = false
    Legend(fig[4,4],ax_f,L"\mathrm{Time}")

    
save(string("fig_SqcompPBC_phi_q_tf10e6_loglog.png"), fig, px_per_unit = 300/inch)
=#

#=
    selec=1;

    for it in 1:N_instants
        lines!(ax,qdomain_plot[selec][it],Sq_plot[selec][it],
            label=latexstring(0.001*time_domains[selec][it]))
    
        p1=plot!(ax_f,[0],[-1],
            label=latexstring(100*dict_CL[unique(Sq[it].id)[1]])
            )
        p2=plot!(ax_f2,[0],[-1],
            label=latexstring(dict_T[unique(Sq[it].id)[1]])
            )
        p3=plot!(ax_f3,[0],[-1],
            label=latexstring(100*dict_phi[unique(Sq[selec].id)[1]])
            )

        p1.visible = false
        p2.visible = false
        p3.visible = false


    end

    Legend(fig[1:2,4],ax,L"\mathrm{Time}")
    Legend(fig[3,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[4,4],ax_f2,L"\mathrm{T}",merge=true)
    Legend(fig[5,4],ax_f3,L"\phi",merge=true)
=#


