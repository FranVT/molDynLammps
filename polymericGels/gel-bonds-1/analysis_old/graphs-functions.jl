"""
    Script with the functions to create graphics
"""

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

function fig_EngPair(dt,DATA,DIR,id_c)

fig=Figure();

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy~Interactions}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle U \rangle~[J/\epsilon]",
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(0,nothing,nothing,nothing) 
   )

plot!(ax,dt.*DATA.TimeStep,DATA."c_wcaPair",label=L"\mathrm{WCA}")
plot!(ax,dt.*DATA.TimeStep,DATA."c_patchPair",label=L"\mathrm{patch}")
plot!(ax,dt.*DATA.TimeStep,DATA."c_swapPair",label=L"\mathrm{swap}")

#hlines!(ax,[T])

Legend(fig[1,2],ax,
      L"\mathrm{Legend}"
     )

save(joinpath(DIR,string("ePair-",id_c,".png")), fig, px_per_unit = 300/inch)

    return fig

end



function fig_EngSys(dt,DATA,DIR,id_c)

fig=Figure();

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy~System}"),
    #subtitle=latexstring(subtitle),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle U \rangle~[J/\epsilon]",
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(0,nothing,nothing,nothing)
   )

plot!(ax,dt.*DATA.TimeStep,DATA."v_eT",label=L"\mathrm{eT}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_ec",label=L"\mathrm{ecouple}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eC",label=L"\mathrm{eConserve}")

#hlines!(ax,[T])

Legend(fig[1,2],ax,
      L"\mathrm{Legend}")

save(joinpath(DIR,string("eSys-",id_c,".png")), fig, px_per_unit = 300/inch)

    return fig

end

function fig_EngB(dt,DATA,DIR,id_c)

fig=Figure();

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy~Bonds}"),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle U \rangle~[J/\epsilon]",
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(0,nothing,nothing,nothing) 
   )

plot!(ax,dt.*DATA.TimeStep,DATA."v_eB",label=L"\mathrm{eB}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eA",label=L"\mathrm{eA}")
plot!(ax,dt.*DATA.TimeStep,DATA."v_eM",label=L"\mathrm{eM}")


Legend(fig[1,2],ax,
      L"\mathrm{Legend}")

save(joinpath(DIR,string("eB-",id_c,".png")), fig, px_per_unit = 300/inch)

    return fig

end

function fig_Eng(dt,DATA,DIR,id_c)

fig=Figure();

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Energy}"),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle U \rangle~[J/\epsilon]",
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(0,nothing,nothing,nothing) 
   )

plot!(ax,dt.*DATA.TimeStep,DATA."c_ep",label=L"\mathrm{ep}")
plot!(ax,dt.*DATA.TimeStep,DATA."c_ek",label=L"\mathrm{ek}")


Legend(fig[1,2],ax,
      L"\mathrm{Legend}",
     labelsize=12pt)

save(joinpath(DIR,string("epk-",id_c,".png")), fig, px_per_unit = 300/inch)

    return fig

end

function fig_Temp(dt,DATA,T,DIR,id_c)

fig=Figure();

clbr=:managua10;

ax=Axis(fig[1:1,1:1],
    title=latexstring("\\mathrm{Temperature}"),
    xlabel=L"\mathrm{Time~units}~[\tau^*]",
    ylabel=L"\langle T \rangle~[kBT/\epsilon]",
    xminorticksvisible=true,
    xminorgridvisible=true,
    limits=(0,nothing,0,1.5*T), 
    yticks = 0:0.01:1.5*T
   )

plot!(ax,dt.*DATA.TimeStep,DATA."c_t")
hlines!(ax,[T])

save(joinpath(DIR,string("temp-",id_c,".png")), fig, px_per_unit = 300/inch)

    return fig

end


