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


####
#   INFORMATION ORGANIZATION
####

function classify_files(SAVE_DIR, id_str)
    # Filtrar archivos CSV que contengan el id (escapando puntos para seguridad)
    all_files = readdir(SAVE_DIR)
    
    # Get all the files with the id of the system
    matched = filter(f -> occursin(id_str, f) && endswith(f, ".csv"), all_files)
    
    # Separar por prefijo
    fix_files = filter(f -> startswith(f, "fix_avg_"), matched)
    sf_files  = filter(f -> startswith(f, "structureFactorPBC"), matched)
    
    return (fix_files, sf_files)
end

"""
    Get the information of the fix files
"""
function get_fixInfo(df,SAVE_DIR)
    # Get the id of the system
    id_str = string(first(unique(df[:, :id])));

    # Get the list of all files
    fix_list, sf_list = classify_files(SAVE_DIR, id_str);

    # Create the paths to the files
    fix_paths = joinpath.(SAVE_DIR, fix_list);
    sf_paths  = joinpath.(SAVE_DIR, sf_list);

    # Read the information of the fix files
    df_fix=CSV.read(fix_paths[1],DataFrame);

    # Select the categories needed to do the plots
    categories_plot=[:id,:TimeStep,:c_t,:c_ep,:c_ek];

    # Select the information to plot
    return df_fix[:,categories_plot]
end












function figureCompareSl(dat_files)


# Grafica del factor de estructura
MAIN_DIR=pwd();
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# File name of interested data
filename_Sq="structureFactor*.csv";

# Create the paths to the files
file_paths=[joinpath(SAVE_DIR,replace(filename_Sq, "*" => it)) for it in unique(dat_files.id)];

# Dictionaries to map id to info.
dict_phi=Dict(dat_files.id.=>dat_files.phi);
dict_CL=Dict(dat_files.id.=>dat_files."CL-Con");
dict_T=Dict(dat_files.id.=>dat_files.Temperature);

# Get the information
Sq=[CSV.read(file,DataFrame) for file in file_paths];

time_domains=[Sq[i].timeStep for i in eachindex(Sq)];

# Parameters to graph
N_sims=length(Sq);
N_instants=nrow(Sq[1]);

qdomain_plot=[collect(eval(Meta.parse(s)) for s in Sq[it].q_domain) for it in eachindex(Sq)];
Sq_plot=[collect(eval(Meta.parse(s)) for s in Sq[it].Sq) for it in eachindex(Sq)];

id_time=1; #length(Sq_plot[1]);

qdomain_tf=[qdomain_plot[id_exp][id_time] for id_exp in eachindex(Sq)];
Sq_tf=[Sq_plot[id_exp][id_time] for id_exp in eachindex(Sq)];

fig=Figure()
    ax_f=Axis(fig[1:1,1:1])
    ax_f2=Axis(fig[1:1,1:1])
    ax_f3=Axis(fig[1:1,1:1])

    ax=Axis(fig[1:6,1:3],
        title=latexstring("\\mathrm{Structure~factor}"),
        #subtitle=latexstring(subtitle),
        xlabel=L"\lambda",
        ylabel=L"S(q)",
        xminorticksvisible=true,
        xminorgridvisible=true,
        limits=(0,10,0,10),
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


    for it in eachindex(Sq_tf)
        lines!(ax,(2*pi)./qdomain_tf[it],Sq_tf[it],
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

save(string("fig_Slcomp_t_",time_domains[1][id_time],".png"), fig, px_per_unit = 300/inch)

end

function figureCompareSq(dat_files)

# Grafica del factor de estructura
MAIN_DIR=pwd();
SAVE_DIR=joinpath(MAIN_DIR,"analyzedData");

# File name of interested data
filename_Sq="structureFactor*.csv";

# Create the paths to the files
file_paths=[joinpath(SAVE_DIR,replace(filename_Sq, "*" => it)) for it in unique(dat_files.id)];

# Dictionaries to map id to info.
dict_phi=Dict(dat_files.id.=>dat_files.phi);
dict_CL=Dict(dat_files.id.=>dat_files."CL-Con");
dict_T=Dict(dat_files.id.=>dat_files.Temperature);

# Get the information
Sq=[CSV.read(file,DataFrame) for file in file_paths];

time_domains=[Sq[i].timeStep for i in eachindex(Sq)];

# Parameters to graph
N_sims=length(Sq);
N_instants=nrow(Sq[1]);

qdomain_plot=[collect(eval(Meta.parse(s)) for s in Sq[it].q_domain) for it in eachindex(Sq)];
Sq_plot=[collect(eval(Meta.parse(s)) for s in Sq[it].Sq) for it in eachindex(Sq)];

id_time=length(Sq_plot[1]);

qdomain_tf=[qdomain_plot[id_exp][id_time] for id_exp in eachindex(Sq)];
Sq_tf=[Sq_plot[id_exp][id_time] for id_exp in eachindex(Sq)];

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
        limits=(nothing,nothing,0,10),
        #xscale=log10,
        #yscale=log10
    )
    hidespines!(ax_f)
    hidedecorations!(ax_f)
    hidespines!(ax_f2)
    hidedecorations!(ax_f2)
    hidespines!(ax_f3)
    hidedecorations!(ax_f3)

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

save(string("fig_Sqcomp_t_",time_domains[1][id_time],".png"), fig, px_per_unit = 300/inch)

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
    Legend(fig[2,2],ax_f,L"\mathrm{CL}~\%",merge=true)
    Legend(fig[3,2],ax_f2,L"\mathrm{T}",merge=true)

    save(joinpath(pwd(),"ep.png"), fig, px_per_unit = 300/inch)
end


