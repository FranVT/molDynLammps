"""
    Grphs jld plots
"""

using GLMakie, LaTeXStrings, Typst_jll
using Peaks, JLD
#, Forecast

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

# Get data

#info = load(joinpath(pwd(),"Sq_exps.jld"));
#data=info["Sq"];
phi=[1,2,3,4,5];

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

L=map(s->unique(data_bySystem[s].L)[1],1:5);

lambda=map(s->(L[s])./data[s][:,1],1:5);

id_time=3;

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])

    #y_max=10;

    ax=Axis(fig[1:6,1:3],
        #title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"|\vec{q}|=\frac{2\pi}{\lambda}=\frac{2\pi}{L}|\vec{n}|",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(-2.5,25,nothing,nothing),
        #xscale=log10,
        #yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)

    map(s->lines!(ax,(2*pi)./lambda[s],data[s][:,id_time],label=latexstring(phi[s])),eachindex(data))

    p1=plot!(ax_f,[0],[-1],
            label=L"t_f"
    
    )
    p1.visible = false

    Legend(fig[4,4],ax,L"\phi~\%",merge=true)
    Legend(fig[1:2,4],ax_f,L"\mathrm{Time}")

