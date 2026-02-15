"""
    Script para crear graficas de las funciones de los potenciales
"""

using GLMakie, LaTeXStrings, Typst_jll

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

# Ranges
Upatch_ij=map(r->Upatch(eps_ij,Dp,r),r_ij);
feval=map(r->forcePatch(eps_ij,Dp,r),r_ij);
Fpatch_ijComp=first.(feval);
Fpatch_ijMag=last.(feval);



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

lines!(ax1_U,r_ij,Upatch_ij,color=:black)

lines!(ax1_F,r_ij,Fpatch_ijComp)
lines!(ax1_F,r_ij,Fpatch_ijMag)

# Swap interaction
fig_Swap=Figure(size = (17cm, 18.75cm));
ax1_U=Axis(fig_Swap[1,1],
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
