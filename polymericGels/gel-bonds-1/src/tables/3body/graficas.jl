"""
    Script para crear graficas de las funciones de los potenciales
"""

using GLMakie, LaTeXStrings, Typst_jll
#using DataFrames

# Parameters for the figure size
# these are relative to 1 CSS px
inch = 96;
pt = 4/3;
cm = inch / 2.54;

# Include the functions
include("functions.jl")

# Parameters for the interactions
Dp=0.4;
eps_ij=1.0;
eps_ik=1.0;
eps_jk=1.0;
w=1.0;

# Domain
N=2^9;
r_min=Dp/4;
r_max=2*Dp;
r_ij=(r_min:(r_max-r_min)/(N-1):r_max);
r_ik=(r_min:(r_max-r_min)/(N-1):r_max);

# Range: Patch-patch interaction
Upatch_ij=map(r->Upatch(eps_ij,Dp,r),r_ij);
feval=map(r->forcePatch(eps_ij,Dp,r),r_ij);
Fpatch_ijComp=first.(feval);
Fpatch_ijMag=last.(feval);

# Ranges: 3 body interaction
r_ijVal=last.(Iterators.partition(eachindex(r_ij),2^4));
Uswap_eval=map(s->map(r->Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij[s],r),r_ik),r_ijVal);

Fswap_eval=map(r_ijVal) do it 
    (Fswap_ij,Fswap_ik,Fswap_mag)=map(r->forceSwap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij[it],r),r_ik)
end

Fswap_eval=reduce.(vcat,Fswap_eval);

#;

"""
    Script to create the graphics/figures
"""
tl_sz=0.6cm;
ot_sz=0.4cm;

fig_Comp=Figure(size = (18.75cm, 15cm));
ax1=Axis(fig_Comp[1,1],
             title=latexstring("\\mathrm{Potentials}"),
             xlabel=L"r_{ik}~[r/Dp]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             limits=(nothing,nothing,-1.5,1.5),
             xticks = r_min:0.2:r_max
   )
ax2=Axis(fig_Comp[2,1],
             title=latexstring("\\mathrm{Magnitude~of~the~forces}"),
             xlabel=L"r_{ik}~[r/Dp]",
             ylabel=L"|\vec{F}|~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             limits=(nothing,nothing,-1,20),
             xticks = r_min:0.2:r_max
   )

# Potentials

pot_patch=scatterlines!(ax1,r_ij,Upatch_ij,label=L"\mathrm{Patch}")

# Definir un mapa de color (puedes cambiarlo, ej: :thermal, :plasma, etc.)
cmap = :thermal

for (idx, g_eval) in enumerate(Uswap_eval)
    r_val = r_ij[r_ijVal[idx]]
    # Normalizar el valor de r_ij al intervalo [0,1] para el mapeo de color
    color_norm = (r_val - r_min) / (r_max - r_min)
    lines!(ax1, r_ik, g_eval,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1),
           linewidth = 1, label=L"\mathrm{Swap}")
end

# Forces

f_mag=scatterlines!(ax2,r_ij,Fpatch_ijMag,label=L"\mathrm{Patch}")

for (idx, g_eval) in enumerate(Fswap_eval)
    r_val = r_ij[r_ijVal[idx]]
    # Normalizar el valor de r_ij al intervalo [0,1] para el mapeo de color
    color_norm = (r_val - r_min) / (r_max - r_min)
    g_ij=g_eval[:,1]
    g_ik=g_eval[:,2]
    g_mag=g_eval[:,3]

    lines!(ax2, r_ik, g_mag,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1), label=L"\mathrm{Swap}"
          )
end

Colorbar(fig_Comp[1:2,2],
         limits = (r_min, r_max),
         colormap = cmap,
         label = L"r_{ij}")   # etiqueta de la barra

axislegend(ax1,
       merge = true, 
       unique = true,
       labelsize=0.5cm,position=:rb)

axislegend(ax2,
       merge = true, 
       unique = true,
       labelsize=0.5cm,position=:rt)


#axislegend(ax, [sc1, sc2], ["One", "Two"], "Selected Dots", position = :rb,
#    orientation = :horizontal)
save("fig_Comp.png", fig_Comp, px_per_unit = 300/inch)

"""
    Components of the swap potencial
"""

fig_comp=Figure(size = (18.75cm, 15cm));
ax1=Axis(fig_comp[1,1],
             title=latexstring("\\hat{r}_{ij}"),
             xlabel=L"r_{ik}~[r/Dp]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             limits=(nothing,nothing,-20,20),
             xticks = r_min:0.2:r_max
   )
ax2=Axis(fig_comp[2,1],
             title=latexstring("\\hat{r}_{ik}"),
             xlabel=L"r_{ik}~[r/Dp]",
             ylabel=L"|\vec{F}|~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             limits=(nothing,nothing,-20,20),
             xticks = r_min:0.2:r_max
   )

for (idx, g_eval) in enumerate(Fswap_eval)
    r_val = r_ij[r_ijVal[idx]]
    # Normalizar el valor de r_ij al intervalo [0,1] para el mapeo de color
    color_norm = (r_val - r_min) / (r_max - r_min)
    g_ij=g_eval[:,1]
    g_ik=g_eval[:,2]
    g_mag=g_eval[:,3]

    lines!(ax1, r_ik, g_ij,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1)#, label=L"\mathrm{Swap}"
          )
end

for (idx, g_eval) in enumerate(Fswap_eval)
    r_val = r_ij[r_ijVal[idx]]
    # Normalizar el valor de r_ij al intervalo [0,1] para el mapeo de color
    color_norm = (r_val - r_min) / (r_max - r_min)
    g_ij=g_eval[:,1]
    g_ik=g_eval[:,2]
    g_mag=g_eval[:,3]

    lines!(ax2, r_ik, g_ik,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1)#, label=L"\mathrm{Swap}"
          )
end

Colorbar(fig_comp[1:2,2],
         limits = (r_min, r_max),
         colormap = cmap,
         label = L"r_{ij}")   # etiqueta de la barra

save("fig_comp.png", fig_comp, px_per_unit = 300/inch)


