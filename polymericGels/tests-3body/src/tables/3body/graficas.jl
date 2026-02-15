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
tl_sz=0.55cm;
ot_sz=0.35cm;

# Patch
fig_Patch=Figure(size = (17cm, 18.75cm));
ax1_U=Axis(fig_Patch[1,1],
             title=latexstring("\\mathrm{Patch-Patch~Interaction}"),
             xlabel=L"r~[r/Dp]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax1_F=Axis(fig_Patch[1,1],
             yaxisposition = :right, 
             ylabel=L"F~[N^*]",
             limits=(nothing,nothing,-2,20)
   )
hidespines!(ax1_F)
hidexdecorations!(ax1_F)

linkyaxes!(ax1_U, ax1_F)

pot=lines!(ax1_U,r_ij,Upatch_ij,color=:black)

f_comp=lines!(ax1_F,r_ij,Fpatch_ijComp)
f_mag=lines!(ax1_F,r_ij,Fpatch_ijMag)

Legend(fig_Patch[1,2],
       [pot,f_comp,f_mag],
       [L"U_{\mathrm{patch}}(r_{ij})", L"f_{\mathrm{patch}}(r_{ij})\hat{r}_{ij}", L"|\vec{F}_{\mathrm{patch}}(r_{ij})|"],
     labelsize=0.5cm)


# Swap interaction
fig_Swap=Figure(size = (17cm, 18.75cm));
ax1_U=Axis(fig_Swap[1,1],
             title=latexstring("\\mathrm{Swap~Interaction}"),
             xlabel=L"r~[r_{ik}/Dp]",
             ylabel=L"U_{\mathrm{swap}}(r_{ij},r_{ik})~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             limits=(nothing,nothing,-0.5,1.5)
   )
ax1_F=Axis(fig_Swap[1,1],
             yaxisposition = :right, 
             ylabel=L"F~[N^*]",
             limits=(nothing,nothing,-0.5,20)
   )

ax1_Fc=Axis(fig_Swap[2,1],
             yaxisposition = :right, 
             ylabel=L"F~[N^*]",
             limits=(nothing,nothing,-50,50)
   )
linkyaxes!(ax1_U, ax1_F)



# Definir un mapa de color (puedes cambiarlo, ej: :thermal, :plasma, etc.)
cmap = :thermal

for (idx, g_eval) in enumerate(Uswap_eval)
    r_val = r_ij[r_ijVal[idx]]
    # Normalizar el valor de r_ij al intervalo [0,1] para el mapeo de color
    color_norm = (r_val - r_min) / (r_max - r_min)
    lines!(ax1_U, r_ik, g_eval,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1),
           linewidth = 1)
end

Colorbar(fig_Swap[1:2,2],
         limits = (r_min, r_max),
         colormap = cmap,
         label = L"r_{ij}")   # etiqueta de la barra

#cmap=:viridis

for (idx, g_eval) in enumerate(Fswap_eval)
    r_val = r_ij[r_ijVal[idx]]
    # Normalizar el valor de r_ij al intervalo [0,1] para el mapeo de color
    color_norm = (r_val - r_min) / (r_max - r_min)
    g_ij=g_eval[:,1]
    g_ik=g_eval[:,2]
    g_mag=g_eval[:,3]

    lines!(ax1_F, r_ik, g_mag,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1),
           linestyle=:dash
          )

    lines!(ax1_Fc, r_ik, g_ij,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1),
           linestyle=:dash
          )

    lines!(ax1_Fc, r_ik, g_ik,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1),
           linestyle=:dot
          )
end





