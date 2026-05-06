"""
    Script para graficar arvhivos jld
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

d = load(joinpath(pwd(),"Sq_phi5.jld2"));
Sq_phi5_tf=d["Sq"];

d = load(joinpath(pwd(),"Sq_phi5_to.jld2"));
Sq_phi5_to=d["Sq"];

L=37.75791;
r_c=1.2;

q_min=2*pi/L;
q_max=2*pi/1;
qdomain=range(q_min,q_max,length=length(Sq_phi5_tf));

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])



    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"|\vec{q}|=\frac{2\pi}{\lambda}",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(nothing,nothing,nothing,nothing),
        #xscale=log10,
        #yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)

    #vlines!(ax,2*pi,linestyle=:solid)
    #for p in peaks_q
    #    vlines!(ax,p,linestyle=:solid)
    #end

    #vlines!(ax,(2*pi),linestyle=:solid)
    #vlines!(ax,(2*pi)/(1.2),linestyle=:solid)
    #vlines!(ax,(2*pi)/(3*1.2),linestyle=:solid)

    vlines!(ax,2*pi/(L),linestyle=:solid,color=:black)
    vlines!(ax,2*pi/(0.5*L),linestyle=:solid,color=:black)
    vlines!(ax,2*pi/(0.25*L),linestyle=:solid,color=:black)
    vlines!(ax,2*pi/(0.125*L),linestyle=:solid,color=:black)
    vlines!(ax,2*pi/(0.0625*L),linestyle=:solid,color=:black)

    vlines!(ax,2*pi/(1.2),linestyle=:dash,color=:blue)
    vlines!(ax,2*pi/(2*1.2),linestyle=:dash,color=:blue)
    vlines!(ax,2*pi/(3*1.2),linestyle=:dash,color=:blue)




        lines!(ax,qdomain,Sq_phi5_tf,label=L"t_f")
        lines!(ax,qdomain,Sq_phi5_to,label=L"t_o")
    

   p1=plot!(ax_f,[0],[-1],
            label=latexstring(5)
            )
    p1.visible = false
    Legend(fig[1:3,4],ax,L"\mathrm{Time}")
   Legend(fig[4,4],ax_f,L"\phi~\%",merge=true)


    #Legend(fig[4,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    #Legend(fig[5,4],ax_f2,L"\mathrm{T}",merge=true)
    #Legend(fig[6,4],ax_f3,L"\mathrm{Time}",merge=true)




