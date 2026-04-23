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

# PAra los labels
phi=unique(dat_files.phi);
temp=unique(dat_files.Temperature);
clCon=unique(dat_files."CL-Con");


# File name of interested data
filename_Sq="structureFactor*.csv";

# Create the paths to the files
file_paths=[joinpath(SAVE_DIR,replace(filename_Sq, "*" => it)) for it in unique(dat_files.id)];

# Get the information
Sq=[CSV.read(file,DataFrame) for file in file_paths];

# Scketch of the graph

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:3,1:1],
        title=latexstring("\\mathrm{Structure factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"\mathrm{1/d}~[D]",
        ylabel=L"U~[S(q)]",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(0,nothing,0,nothing),
        #xscale=log10,
        #yscale=log10
    )

    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)

    linkyaxes!(ax,ax_f,ax_f2)
    linkxaxes!(ax,ax_f,ax_f2)


    for it in eachindex(Sq)
        Sq[it].Sq
    end



    for it in 1:eachindex(Sq)
        plot!(ax_f,[-1],[0],
            label=latexstring(100*dat_DF."CL-Con"[it])
            )
        plot!(ax_f2,[-1],[0],
            label=latexstring(dat_DF."Temperature"[it])
            )
    end

    Legend(fig[1,2],ax,L"\phi~\%")
    Legend(fig[2,2],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[3,2],ax_f2,L"\mathrm{T}",merge=true)




