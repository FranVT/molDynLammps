"""
    Script to make graphs
"""

using CSV, DataFrames
using GLMakie, LaTeXStrings, Typst_jll
using Peaks
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


#[collect(eval(Meta.parse(s)) for s in Sq[it].Sq) for it in eachindex(Sq)];

time_domains=[Sq[i].timeStep for i in eachindex(Sq)];
qdomain_plot=[collect(eval(Meta.parse(s)) for s in Sq[it].q_domain) for it in eachindex(Sq)];
Sq_plot=[collect(eval(Meta.parse(s)) for s in Sq[it].Sq) for it in eachindex(Sq)];

id_time=1;

qdomain_tf=[qdomain_plot[id_exp][id_time] for id_exp in eachindex(Sq)];
Sq_tf=[Sq_plot[id_exp][id_time] for id_exp in eachindex(Sq)];

#peaks_Sq_tf=[findmaxima(s) for s in Sq_tf];
#peaks_q=qdomain_tf[1][peaks_Sq_tf[1].indices];


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
        limits=(nothing,nothing,0,5),
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

    #vlines!(ax,(2*pi)/(1.2*0.5),linestyle=:solid)


    #    annotation!(-200, 0, 0, 0, path = Ann.Paths.Line(), text = "Line()")

    annotation!(ax, pi,4, 2*pi/1.2,4,
        text = L"\frac{2\pi}{1.2}",
        path = Ann.Paths.Arc(0.3),
        style = Ann.Styles.LineArrow(),
        labelspace = :data
        )
    
    annotation!(ax, 2*pi,4, 2*pi/(1.2*0.5),4,
        text = L"\frac{2\pi}{0.24}",
        path = Ann.Paths.Arc(0.3),
        style = Ann.Styles.LineArrow(),
        labelspace = :data
        )



    for it in eachindex(Sq_tf)
        lines!(ax,qdomain_tf[it],Sq_tf[it],
            label=latexstring(100*dict_phi[unique(Sq[it].id)[1]]))
    
        p1=plot!(ax_f,[0],[-1],
            label=latexstring(100*dict_CL[unique(Sq[it].id)[1]])
            )
        p2=plot!(ax_f2,[0],[-1],
            label=latexstring(dict_T[unique(Sq[it].id)[1]])
            )
        p3=plot!(ax_f3,[0],[-1],
            label=latexstring(0.001*Sq[it].timeStep[id_time])
            )

    p1.visible = false
    p2.visible = false
    p3.visible = false


    end

    Legend(fig[1:3,4],ax,L"\phi~\%")
    
    Legend(fig[4,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[5,4],ax_f2,L"\mathrm{T}",merge=true)
    Legend(fig[6,4],ax_f3,L"\mathrm{Time}",merge=true)
















#=
Sq2analyze_tf=Sq_plot[id_exp][end];

peaks_Sq_tf=[findmaxima(s) for s in Sq_tf];
peaks_q_tf=[l_domain[s][peaks_Sq_tf[s][1]] for s in eachindex(Sq_tf)];

#peaks_q_tf=[peaks_q_tf[s][1:id_min] for s in eachindex(Sq_tf)];
#peaks_q_tf=reduce(vcat,mean(reduce(hcat,peaks_q_tf),dims=2));





#(ind_peaks_tf,peaks_tf,~)=findmaxima(Sq2analyze_tf);


cmap = :nipy_spectral #:viridis

# Normaliza entre 0 y 1
norm_peaks = (eachindex(Sq_tf) .- minimum(eachindex(Sq_tf))) ./ (maximum(eachindex(Sq_tf)) - minimum(eachindex(Sq_tf)))
colors = cgrad(cmap)[norm_peaks]

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

    #for it in eachindex(peaks_q_tf) 
    #    vlines!(ax,peaks_q_tf[it],linestyle=:dash,color=(:black,0.25))
    #end

    vlines!(ax,2*pi,linestyle=:solid)
    vlines!(ax,(2*pi)*(1/1.1),linestyle=:solid)
    vlines!(ax,(2*pi)*(1/1.2),linestyle=:solid)
    vlines!(ax,(2*pi)*(1/1.3),linestyle=:solid)
    vlines!(ax,(2*pi)*(1/1.4),linestyle=:solid)
    vlines!(ax,(2*pi)*(1/1.5),linestyle=:solid)
    vlines!(ax,(2*pi)*(1/2),linestyle=:solid)
    vlines!(ax,(2*pi)*(1/3),linestyle=:solid)





    for it in eachindex(Sq_tf)

        q_domain=collect(l_domain[id_exp]);
        range=Sq_tf[it];

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
              label=L"t_f"
            )

    p1.visible = false
    p2.visible = false
    p3.visible = false


    end

    #Legend(fig[1:3,4],ax,L"\tau")
    Colorbar(fig[1:3, 4], colormap=cmap,
         limits = (1,5),
         label = L"\phi\%")

    Legend(fig[4,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[5,4],ax_f2,L"\mathrm{T}",merge=true)
    Legend(fig[6,4],ax_f3,L"\mathrm{Time}",merge=true)


=#


#=
id_exp=5;
Sq2analyze_to=Sq_plot[id_exp][1];
Sq2analyze_tf=Sq_plot[id_exp][end];

(ind_peaks_to,peaks_to,~)=findmaxima(Sq2analyze_to);
(ind_peaks_tf,peaks_tf,~)=findmaxima(Sq2analyze_tf);

peaks_timeDomain=[[] for _ in eachindex(time_domains[id_exp])];

for it in eachindex(time_domains[id_exp])
    peaks_timeDomain[it]=Sq_plot[id_exp][it][ind_peaks_to];
end    

#q_peaks=reduce(hcat,q_peaks);
peaks_timeDomain=reduce(hcat,peaks_timeDomain);



#peaks_sq=[findmaxima(Sq_plot[id_exp][it_time]) for it_time in eachindex(Sq_plot[id_exp])]

# Plot the time evolution of the peaks

cmap = :nipy_spectral #:viridis

# Normaliza entre 0 y 1
norm_peaks = (peaks_to .- minimum(peaks_to)) ./ (maximum(peaks_to) - minimum(peaks_to))
colors = cgrad(cmap)[norm_peaks]

    fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Peaks~evolution}"),
        #subtitle=latexstring(subtitle),
        xlabel=latexstring("\\mathrm{Peaks~evolution}"),
        ylabel=L"S(q)_{\mathrm{peak}}",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(10^(2.5),10^(4.1),0,25),
        xscale=log10,
        #yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)

    for it in eachindex(ind_peaks_to)

        scatterlines!(ax,(0.001).*time_domains[1],peaks_timeDomain[it,:],color=colors[it])
   
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
         limits = (minimum(peaks_to), maximum(peaks_to)),
         label = L"|\vec{q}|L")

    Legend(fig[4,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[5,4],ax_f2,L"\mathrm{T}",merge=true)
    Legend(fig[6,4],ax_f3,L"\phi~\%",merge=true)
=#



function timeSeries()
#id_exp=5;

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

    for it in ind_peaks_to
        vlines!(ax,l_domain[id_exp][it],linestyle=:dash,color=:black)
    end

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

    for it in eachindex(ind_peaks_to)
        scatter!(ax,l_domain[id_exp][ind_peaks_to[it]],peaks_to[it],color=:black)
    end
    for it in eachindex(ind_peaks_tf)
        scatter!(ax,l_domain[id_exp][ind_peaks_tf[it]],peaks_tf[it],color=:red)
    end

    #Legend(fig[1:3,4],ax,L"\tau")
    Colorbar(fig[1:3, 4], colormap=cmap,
         limits = (minimum(taus), maximum(taus)),
         label = L"\tau")

    Legend(fig[4,4],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[5,4],ax_f2,L"\mathrm{T}",merge=true)
    Legend(fig[6,4],ax_f3,L"\phi~\%",merge=true)
end





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



