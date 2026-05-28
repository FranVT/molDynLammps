"""
    Script that execute all the analysis from a given system
"""

using DataFrames, CSV
using Statistics, StatsBase
using LinearAlgebra , Random
using GLMakie, LaTeXStrings, Typst_jll

Random.seed!(1234)

# Include auxiliary files
include("functions.jl")
include("functions_graphs.jl")

# Activate functions
up=0;
Sq=0;
Sq_PBC=0;


# Update the data file with new systems
if up ==1 
    createDatFiles();
end

# Extraer la información del dat file
dat_files=extractDatFiles();

# Selección de categorias
categories=[:phi];

# Creación de los subdataframes por sistema
data_bySystem=groupby(dat_files,categories);

# Analyze the structure factor per system
if Sq_PBC == 1
    storeAllSq(data_bySystem,[dumpAnalysis(dat_DF) for dat_DF in data_bySystem])
end


#figureCompareSq(dat_files)
#figureCompareSl(dat_files)

# Grafica del factor de estructura
MAIN_DIR=pwd();
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# File name of interested data
filename_Sq="structureFactorPBC*.csv";

# Create the paths to the files
file_paths=[joinpath(SAVE_DIR,replace(filename_Sq, "*" => it)) for it in unique(dat_files.id)];

# Dictionaries to map id to info.
dict_phi=Dict(dat_files.id.=>dat_files.phi);
dict_CL=Dict(dat_files.id.=>dat_files."CL-Con");
dict_T=Dict(dat_files.id.=>dat_files.Temperature);

# Get the information
Sq=[CSV.read(file,DataFrame) for file in file_paths];

#=
println("Se inicia a graficar")
fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"|\vec{q}|",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(0,2,0,25),
        #xscale=log10,
        #yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)

    #vlines!(ax,1.2,linestyle=:dash,color=:grey)
    #vlines!(ax,5*1.2,linestyle=:dash,color=:grey)

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

save(string("fig_SqcompPBC_phi_",100*dict_phi[unique(Sq[selec].id)[1]],"_zoom.png"), fig, px_per_unit = 300/inch)

=#
