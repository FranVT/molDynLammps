"""
    Script to test the summatory stuff for the swap potential
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

# Assuming gamma = pi 
# r_jk^2 = r_ij^2+r_ik^2-2*r_ij*r_ik*cos(gamma)
r_jk=r_ij.^2 .+ r_ik.^2 .- 2*r_ij.*r_ik.*cos(pi);



# Ranges: 3 body interaction
r_ijVal=last.(Iterators.partition(eachindex(r_ij),2^4));

Uswap_eval1=map(s->map(r->Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij[s],r),r_ik),r_ijVal);
Uswap_eval2=map(s->map(r->Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij[s],r),r_jk),r_ijVal);
Uswap_eval3=map(s->map(r->Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ik[s],r),r_jk),r_ijVal);

# Add the potentials

Uswap_eval=map(s->Uswap_eval1[s]+Uswap_eval2[s]+Uswap_eval3[s],eachindex(Uswap_eval1))

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
             limits=(nothing,nothing,-0.5,3.5),
             xticks = r_min:0.2:r_max
   )
#=
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
=#

# Potentials

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
#=
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
=#

Colorbar(fig_Comp[1:2,2],
         limits = (r_min, r_max),
         colormap = cmap,
         label = L"r_{ij}")   # etiqueta de la barra

#=
axislegend(ax1,
       merge = true, 
       unique = true,
       labelsize=0.5cm,position=:rb)

axislegend(ax2,
       merge = true, 
       unique = true,
       labelsize=0.5cm,position=:rt)
=#

