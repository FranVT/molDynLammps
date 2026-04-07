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
        labelsize = 14pt,                   # tamaño texto leyenda
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
    Functions
=#

function getDatInfo(DF_DIR)
"""
    This function get the dat.csv of experiments
"""

    archivos = filter(str -> occursin("dat-", str), DF_DIR)
    df_final = DataFrame()  # vacío inicial
    for (i, archivo) in enumerate(archivos)
        df_temp = CSV.read(archivo, DataFrame)
        unique!(df_temp,2)
        append!(df_final, df_temp)
        df_temp = nothing  # liberar referencia
    end
    return df_final
end

function getInfoSystem(DIR,files)
"""
    This functions returns the system information given a dat dataframe
"""
    archivos = joinpath.(DIR,files)
    df_final = [DataFrame() for i in eachindex(archivos)];  # vacío inicial

    for (i, archivo) in enumerate(archivos)
        df_final[i] = CSV.read(archivo, DataFrame)
    end
    return df_final
end

function potentialEnergyFig(dat_DF,system_DF)
"""
    Function that stores a comparisson of potential energy
"""
    y_min=round(minimum(map(s->minimum(s.c_ep),system_DF)),sigdigits=1);

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:3,1:1],
        title=latexstring("\\mathrm{Potential~Energy}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"\mathrm{Time~units}~[\tau^*]",
        ylabel=L"U~[\mathrm{J}/\epsilon]",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(0,nothing,y_min,nothing),
        #xscale=log10,
        #yscale=log10
    )

    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)

    linkyaxes!(ax,ax_f,ax_f2)
    linkxaxes!(ax,ax_f,ax_f2)

    for it in 1:nrow(dat_DF)
        plot!(ax,dat_DF."time-step"[it].*system_DF[it].TimeStep,system_DF[it].c_ep,
            label=latexstring(100*dat_DF.phi[it]),
            )
        plot!(ax_f,[-1],[0],
            label=latexstring(100*dat_DF."CL-Con"[it])
            )
        plot!(ax_f2,[-1],[0],
            label=latexstring(dat_DF."Temperature"[it])
            )
    end

    Legend(fig[1,2],ax,L"\phi~\%")
    Legend(fig[2,2],ax_f,L"\mathrm{CL}~\%")
    Legend(fig[3,2],ax_f2,L"\mathrm{T}")

    save(joinpath(pwd(),"ep.png"), fig, px_per_unit = 300/inch)
end

function patchEnergyFig(dat_DF,system_DF)
"""
    Function that stores a comparisson of potential energy
"""
    y_min=round(minimum(map(s->minimum(s.c_patchPair),system_DF)),sigdigits=1,RoundDown);

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:3,1:1],
        title=latexstring("\\mathrm{Patch~pair~Energy}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"\mathrm{Time~units}~[\tau^*]",
        ylabel=L"U~[\mathrm{J}/\epsilon]",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(0,nothing,y_min,nothing),
        #xscale=log10,
        #yscale=log10
    )

    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)

    linkyaxes!(ax,ax_f,ax_f2)
    linkxaxes!(ax,ax_f,ax_f2)

    for it in 1:nrow(dat_DF)
        plot!(ax,dat_DF."time-step"[it].*system_DF[it].TimeStep,system_DF[it].c_patchPair,
            label=latexstring(100*dat_DF.phi[it]),
            )
        plot!(ax_f,[-1],[0],
            label=latexstring(100*dat_DF."CL-Con"[it])
            )
        plot!(ax_f2,[-1],[0],
            label=latexstring(dat_DF."Temperature"[it])
            )
    end

    Legend(fig[1,2],ax,L"\phi~\%")
    Legend(fig[2,2],ax_f,L"\mathrm{CL}~\%")
    Legend(fig[3,2],ax_f2,L"\mathrm{T}")

    save(joinpath(pwd(),"patchyEnergy.png"), fig, px_per_unit = 300/inch)
end
#=
    Script
=#

# Get directories 
MAIN_DIR=dirname(pwd());
INFO_DIR=joinpath(pwd(),"data_mod");
DF_DIR=filter(isfile,readdir(INFO_DIR,join=true));

# Get the dat dataframes
dat_DF=getDatInfo(DF_DIR);

# Get the info of each system
system_DF=getInfoSystem(INFO_DIR,dat_DF.file_system);

#=
    Figures
=#

potentialEnergyFig(dat_DF,system_DF)
patchEnergyFig(dat_DF,system_DF)

