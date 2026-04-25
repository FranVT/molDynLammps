"""
    Script to make graphs
"""

using CSV, DataFrames
using GLMakie, LaTeXStrings, Typst_jll

# these are relative to 1 CSS px
inch = 96;
pt = 4/3;
cm = inch / 2.54;

# Condigure the theme for ALL figures
set_theme!(
    backgroundcolor = :white,
    fontsize = 16pt,                     # tamaño de letra base
    Axis = (
        titlesize = 20pt,                 # tamaño del título del eje
        xlabelsize = 14pt,                # tamaño etiqueta eje x
        ylabelsize = 14pt,                 # tamaño etiqueta eje y
        xticklabelsize = 12pt,             # tamaño números eje x
        yticklabelsize = 12pt,              # tamaño números eje y
        xgridstyle = :dash,               # estilo de la cuadrícula
        ygridstyle = :dash,
        spinewidth = 1.5pt,
    ),
    Legend = (
        labelsize = 10pt,                   # tamaño texto leyenda
        framewidth = 1.5pt,
    ),
    Colorbar = (
        labelsize = 14pt,
        ticklabelsize = 12pt,
    ),
    Figure = (
        size = (15cm, 12cm)                 # tamaño de la figura (ancho, alto)
    )
)

#=
    Script
=#

include("functions_graphs.jl")

# Get directories 
MAIN_DIR=pwd();
DAT_PATH=joinpath(MAIN_DIR,"datFiles","experiments_dat.csv");
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# Get the dat dataframes
dat_files=CSV.read(DAT_PATH,DataFrame);

# Dictionaries to map id to info.
dict_phi=Dict(dat_files.id.=>dat_files.phi);
dict_CL=Dict(dat_files.id.=>dat_files."CL-Con");
dict_T=Dict(dat_files.id.=>dat_files.Temperature);

# File name of interested data
filename_Sq="structureFactor*.csv";

# Create the paths to the files
file_paths=[joinpath(SAVE_DIR,replace(filename_Sq, "*" => it)) for it in unique(dat_files.id)];

# Get the information
Sq=[CSV.read(file,DataFrame) for file in file_paths];

# Parameters to graph
N_sims=length(Sq);
N_instants=nrow(Sq[1]);

# Prepare the data for the graphs

time_domains=[Sq[i].timeStep for i in eachindex(Sq)];
Sq_plot=[collect(eval(Meta.parse(s)) for s in Sq[it].Sq) for it in eachindex(Sq)];

l_domain=[range(2*pi/first(Sq[N].lambda_o),2*pi/first(Sq[N].lambda_f),length(Sq_plot[N][1]))*first(Sq[N].lambda_f) for N in eachindex(Sq)];


id_exp=5;

taus = 0.001 * time_domains[id_exp]          # vector de tiempos
cmap = :viridis

# Normaliza entre 0 y 1
norm_taus = (taus .- minimum(taus)) ./ (maximum(taus) - minimum(taus))
colors = cgrad(cmap)[norm_taus]

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])



    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"|\vec{q}|L=\frac{2\pi}{\lambda}L",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(nothing,nothing,0,nothing),
        #xscale=log10,
        #yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)
 

    for it in eachindex(time_domains[id_exp])

        q_domain=collect(l_domain[id_exp]);
        range=Sq_plot[id_exp][it];

        lines!(ax,q_domain,range,
               color=colors[it],
               label="" #latexstring(0.001*time_domains[id_exp][it])
           )
    
        p1=plot!(ax_f,[0],[-1],
            label=latexstring(100*dict_CL[unique(Sq[id_exp].id)[1]])
            )
        p2=plot!(ax_f2,[0],[-1],
            label=latexstring(dict_T[unique(Sq[id_exp].id)[1]])
            )
        p3=plot!(ax_f3,[0],[-1],
              label=latexstring(100*dict_phi[unique(Sq[id_exp].id)[1]])
            )

    p1.visible = false
    p2.visible = false
    p3.visible = false


    end

    #Legend(fig[1:3,4],ax,L"\tau")
    Colorbar(fig[1:3, 4], colormap=cmap,
         limits = (minimum(taus), maximum(taus)),
         label = L"\tau")

    Legend(fig[4,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[5,4],ax_f2,L"\mathrm{T}",merge=true)
    Legend(fig[6,4],ax_f3,L"\phi~\%",merge=true)






function structureFactorGraph()
id=5;

# Scketch of the graph

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])



    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"|\vec{q}|L=\frac{2\pi}{\lambda}L",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(nothing,nothing,0,nothing),
        xscale=log10,
        #yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)
 

    for it in eachindex(Sq)
        Sq_plot=[Float64.(M) for M in eval(Meta.parse(Sq[it].Sq[id]))];
        l_o=first(unique(Sq[it].lambda_o));
        l_f=first(unique(Sq[it].lambda_f));
        dom=range(l_o,l_f,length=length(Sq_plot));
        

        data=[[(2*pi/dom[it])*l_f,Sq_plot[it]]  for it in eachindex(dom)];
        data=reduce(hcat,data);



        lines!(ax,data[1,:],data[2,:],
            label=latexstring(100*dict_phi[unique(Sq[it].id)[1]]))
    
        p1=plot!(ax_f,[0],[-1],
            label=latexstring(100*dict_CL[unique(Sq[it].id)[1]])
            )
        p2=plot!(ax_f2,[0],[-1],
            label=latexstring(dict_T[unique(Sq[it].id)[1]])
            )
        p3=plot!(ax_f3,[0],[-1],
              label=latexstring(0.001*Sq[it].timeStep[id])
            )

    p1.visible = false
    p2.visible = false
    p3.visible = false


    end

    Legend(fig[1:3,4],ax,L"\phi~\%")
    
    Legend(fig[4,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[5,4],ax_f2,L"\mathrm{T}",merge=true)
    Legend(fig[6,4],ax_f3,L"\mathrm{Time}",merge=true)
end



