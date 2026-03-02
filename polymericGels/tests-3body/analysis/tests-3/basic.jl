"""
    Basic analysis
"""

using DataFrames, CSV
#using Plots, LaTeXStrings, Plots.PlotMeasures
#gr()
using Statistics
using GLMakie, LaTeXStrings, Typst_jll

# these are relative to 1 CSS px
inch = 96;
pt = 4/3;
cm = inch / 2.54;

# Load the functions
include("functions.jl")


# 2026-02-23-112958

# Selection of an specific simulation
date="2026-03-02-093124";

#"2026-02-27-163258";

#"2026-02-27-155337";

#"2026-02-27-122735";

#"2026-02-27-151246";
#"2026-02-27-144232";
#
#"2026-02-27-121522";
#"2026-02-27-120321";

#"2026-02-27-101516";
#"2026-02-27-100222";
#"2026-02-27-091917";

#"2026-02-26-142919";
#"2026-02-26-142254";
#"2026-02-26-123205";
#"2026-02-26-122533";


#"2026-02-23-112958";
#"2026-02-26-110302";

#"2026-02-25-144129";
#"2026-02-25-143337";
#"2026-02-25-142453";
#"2026-02-25-140610";
#"2026-02-25-135928";
#"2026-02-25-135336";
#"2026-02-25-132053";
#"2026-02-25-115459";
#"2026-02-25-114530";
#"2026-02-25-113617";
#"2026-02-25-113000";
#"2026-02-25-110028";
#"2026-02-25-105552";
#"2026-02-25-101906";

#"2026-02-23-111821";
#"2026-02-23-120646";
# Config2: 

# Config1: "2026-02-23-111821";


#"2026-02-17-121845";
#"2026-02-17-121845";

#"2026-02-17-120748";
#"2026-02-17-111013";

#"2026-02-13-124253";
#"2026-02-13-152655";

#"2026-02-12-155912";
#"2026-02-12-154430";
#"2026-02-12-152030";
#"2026-02-12-145030";

# Get the directory of the desire system
DIR=getDir(date);
DIR=DIR[1];


# Activate extract info
act=1;

# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
if act == 1
    data=extractFixScalar(DIR,"system_assembly.fixf");
    # Convert the array into a DataFrame
    DATA_fix=DataFrame(data[2]',data[1]);

    # Get table
    tableSwap=getTable3b(DIR,"swapMechTab2new_w1.table");
    swapTabPlot=Iterators.partition(tableSwap.e,200)|>collect;

    # Obtener la evolución de r_ik
    id=1:200:nrow(tableSwap) #first.(Iterators.partition(1:1:length(tableSwap.theta),200)|>collect);

    # Obtener cada rango de r_ik
    table_pot1=tableSwap[id, [:r_ij, :r_ik, :e]];

    # Separa los rangos de r_ik por cada valor de r_ij
    table_pot1=groupby(table_pot1, :r_ij);

    # Get the directory of the desire system
    DIR=joinpath(DIR,"traj");
    DATA_dump=map(s->getDump(DIR,s),readdir(DIR));

    """
        Reducing the DATA_dump
    """

    # Concatenate the data frames into one
    l=sort(reduce(vcat,DATA_dump),:TimeStep);

    # Create Individual data frames
    DATA_dump_1 = l[l.id.==1.0,:];
    DATA_dump_2 = l[l.id.==2.0,:];
    DATA_dump_3 = l[l.id.==3.0,:];

end

"""
    Usefull stuff for the graphs
"""

time_fix=DATA_fix.TimeStep.*0.001;
time_dump=DATA_dump_1.TimeStep.*0.001;


# Compute the unitary vectors r_ik and r_jk and express them in the cartisian basis

# The components of the vectors in terms of x,y,z basis
deltaX_12=DATA_dump_2.x .- DATA_dump_1.x;
deltaY_12=DATA_dump_2.y .- DATA_dump_1.y;

deltaX_13=DATA_dump_3.x .- DATA_dump_1.x;
deltaY_13=DATA_dump_3.y .- DATA_dump_1.y;

deltaX_23=DATA_dump_3.x .- DATA_dump_2.x;
deltaY_23=DATA_dump_3.y .- DATA_dump_2.y;

# The norms
dist_12=sqrt.(deltaX_12.^2 .+ deltaY_12.^2);
dist_13=sqrt.(deltaX_13.^2 .+ deltaY_13.^2);
dist_23=sqrt.(deltaX_23.^2 .+ deltaY_23.^2);

# Unitary vectors
vr_12=[deltaX_12 deltaY_12]./dist_12;
vr_13=[deltaX_13 deltaY_13]./dist_13;
vr_23=[deltaX_23 deltaY_23]./dist_23;

# Compute the magnitude of the force
DATA_dump_1[!, :norma_f] = sqrt.(DATA_dump_1.fx.^2 .+ DATA_dump_1.fy.^2 .+ DATA_dump_1.fz.^2)
DATA_dump_2[!, :norma_f] = sqrt.(DATA_dump_2.fx.^2 .+ DATA_dump_2.fy.^2 .+ DATA_dump_2.fz.^2)
DATA_dump_3[!, :norma_f] = sqrt.(DATA_dump_3.fx.^2 .+ DATA_dump_3.fy.^2 .+ DATA_dump_3.fz.^2)

F_systemL=DATA_dump_1.norma_f .+ DATA_dump_2.norma_f .+DATA_dump_3.norma_f;

# Evaluate the forces

# Evaluation and projection of the forces from the theorical definition
f_patch12=map(r->last(forcePatch(1.0,0.4,r)),dist_12);
f_patch13=map(r->last(forcePatch(1.0,0.4,r)),dist_13);
f_patch23=map(r->last(forcePatch(1.0,0.4,r)),dist_23);

f_patch12xy=f_patch12.*vr_12;
f_patch13xy=f_patch13.*vr_13; 
f_patch23xy=f_patch23.*vr_23; 

# Evaluate the Swap function


function evaluateFswap(r12,r13,r23)
    r_c=1.5*0.4;

    if r12<r_c && r13<r_c && r23<r_c
        f_1=forceSwap(1.0,1.0,1.0,1.0,0.4,r12,r13)
        f_2=forceSwap(1.0,1.0,1.0,1.0,0.4,r12,r23)
        f_3=forceSwap(1.0,1.0,1.0,1.0,0.4,r12,r13)
end
end


f_swap1=mapreduce(s->forceSwap(1.0,1.0,1.0,1.0,0.4,dist_12[s],dist_13[s]),vcat,1:1:nrow(DATA_dump_1));
f_swap2=mapreduce(s->forceSwap(1.0,1.0,1.0,1.0,0.4,dist_12[s],dist_23[s]),vcat,1:1:nrow(DATA_dump_2));
f_swap3=mapreduce(s->forceSwap(1.0,1.0,1.0,1.0,0.4,dist_13[s],dist_23[s]),vcat,1:1:nrow(DATA_dump_3));

# Stuff of the table
F_1s12 = f_swap1[:,1] .- f_swap2[:,1];
F_1s13 = f_swap1[:,2] .- f_swap3[:,1];

F_2s21 = f_swap2[:,1] .- f_swap1[:,1];
F_2s23 = f_swap2[:,2] .- f_swap3[:,2];

F_3s31 = f_swap3[:,1] .- f_swap1[:,2];
F_3s32 = f_swap3[:,2] .- f_swap2[:,2];

# Project the forces into the x,y basis
f_swap1xy=(F_1s12).*vr_12 .+ (F_1s13).*vr_13;
f_swap2xy=(F_2s21).*(vr_12) .+ (F_2s23).*vr_23;
f_swap3xy=(F_3s31).*(vr_13) .+ (F_3s32).*(vr_23);

#=
f_swap1xy=f_swap1[:,1].*vr_12 .+ f_swap1[:,2].*vr_13;
f_swap2xy=f_swap2[:,1].*(-vr_12) .+ f_swap2[:,2].*vr_23;
f_swap3xy=f_swap3[:,1].*(-vr_13) .+ f_swap3[:,2].*(-vr_23);
=#

f_part1=f_patch12xy .+ f_patch13xy .+ f_swap1xy;
f_part2=-f_patch12xy .+ f_patch23xy .+ f_swap2xy;
f_part3=-f_patch13xy .- f_patch23xy .+ f_swap3xy;

f_part1N=sqrt.(reduce(vcat,sum(f_part1.^2,dims=2)));
f_part2N=sqrt.(reduce(vcat,sum(f_part2.^2,dims=2)));
f_part3N=sqrt.(reduce(vcat,sum(f_part3.^2,dims=2)));

F_systemA=f_part1N .+ f_part2N .+ f_part3N;

# Graphics

tl_sz=0.55cm;
ot_sz=0.35cm;

"""
    Full components of the swap force analitical
"""

fig_Fswap=Figure(size=(18.75cm,18.75cm));
ax1_Fp=Axis(fig_Fswap[1,1],
             title=latexstring("\\vec{F}_{\\mathrm{patch}~1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_Fp=Axis(fig_Fswap[2,1],
             title=latexstring("\\vec{F}_{\\mathrm{patch}~2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_Fp=Axis(fig_Fswap[3,1],
             title=latexstring("\\vec{F}_{\\mathrm{patch}~3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )

ax1_Fs=Axis(fig_Fswap[1,2],
             title=latexstring("\\vec{F}_{\\mathrm{swap}~1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_Fs=Axis(fig_Fswap[2,2],
             title=latexstring("\\vec{F}_{\\mathrm{swap}~2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_Fs=Axis(fig_Fswap[3,2],
             title=latexstring("\\vec{F}_{\\mathrm{swap}~3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )


ax1_Ft=Axis(fig_Fswap[1,3],
             title=latexstring("\\vec{F}_{1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_Ft=Axis(fig_Fswap[2,3],
             title=latexstring("\\vec{F}_{2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_Ft=Axis(fig_Fswap[3,3],
             title=latexstring("\\vec{F}_{3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )


linkyaxes!(ax1_Fp, ax2_Fp, ax3_Fp, ax1_Ft, ax2_Ft, ax3_Ft, ax1_Fs, ax2_Fs, ax3_Fs) #, ax3_F, ax1_Fa, ax2_Fa, ax3_Fa)

lines!(ax1_Fp,time_dump,f_patch12xy[:,1].+f_patch13xy[:,1],
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax1_Fp,time_dump,f_patch12xy[:,2].+f_patch13xy[:,2],
        color=1,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )

lines!(ax1_Fs,time_dump,f_swap1xy[:,1],
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax1_Fs,time_dump,f_swap1xy[:,2],
        color=1,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )

lines!(ax1_Ft,time_dump,f_part1[:,1],
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax1_Ft,time_dump,f_part1[:,2],
        color=1,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )


lines!(ax2_Fp,time_dump,-f_patch12xy[:,1] .+ f_patch23xy[:,1],
        color=2,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax2_Fp,time_dump,-f_patch12xy[:,2] .+ f_patch23xy[:,2],
        color=2,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )

lines!(ax2_Fs,time_dump,f_swap2xy[:,1],
        color=2,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax2_Fs,time_dump,f_swap2xy[:,2],
        color=2,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )

lines!(ax2_Ft,time_dump,f_part2[:,1],
        color=2,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax2_Ft,time_dump,f_part2[:,2],
        color=2,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )


lines!(ax3_Fp,time_dump,-f_patch13xy[:,1] .- f_patch23xy[:,1],
        color=3,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax3_Fp,time_dump,-f_patch13xy[:,2] .- f_patch23xy[:,2],
        color=3,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )

lines!(ax3_Fs,time_dump,f_swap3xy[:,1],
        color=3,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax3_Fs,time_dump,f_swap3xy[:,2],
        color=3,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )

lines!(ax3_Ft,time_dump,f_part3[:,1],
        color=3,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax3_Ft,time_dump,f_part3[:,2],
        color=3,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )


#=
lines!(ax2_F,time_dump,f_swap1[:,2],
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax2_F,time_dump,f_swap2[:,2],
        color=2,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax2_F,time_dump,f_swap3[:,2],
        color=3,
        colormap=:tab10,
        colorrange=(1,10)
       )


lines!(ax1_Fxy,time_dump,f_swap1xy[:,1],
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax1_Fxy,time_dump,f_swap2xy[:,1],
        color=2,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax1_Fxy,time_dump,f_swap3xy[:,1],
        color=3,
        colormap=:tab10,
        colorrange=(1,10)
       )

lines!(ax2_Fxy,time_dump,f_swap1xy[:,2],
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax2_Fxy,time_dump,f_swap2xy[:,2],
        color=2,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax2_Fxy,time_dump,f_swap3xy[:,2],
        color=3,
        colormap=:tab10,
        colorrange=(1,10)
       )
=#

save("fig_ForceSwapcomp.png", fig_Fswap, px_per_unit = 300/inch)

"""
    Plot the total force of  
"""


"""
    Plot the components
"""
fig_F=Figure(size = (17cm, 18.75cm));

# LAMMPS
ax1_F=Axis(fig_F[1,1],
             title=latexstring("\\mathrm{LAMMPS~Particle~1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_F=Axis(fig_F[2,1],
             title=latexstring("\\mathrm{LAMMPS~Particle~2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_F=Axis(fig_F[3,1],
             title=latexstring("\\mathrm{LAMMPS~Particle~3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax4_F=Axis(fig_F[4,1],
           title=latexstring("\\sum|f_{i}|~\\mathrm{LAMMPS}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )


ax1_Fa=Axis(fig_F[1,2],
             title=latexstring("\\mathrm{Analytical~Particle~1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_Fa=Axis(fig_F[2,2],
             title=latexstring("\\mathrm{Analytical~Particle~2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_Fa=Axis(fig_F[3,2],
             title=latexstring("\\mathrm{Analytical~Particle~3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )


ax4_Fa=Axis(fig_F[4,2],
            title=latexstring("\\sum|f_{i}|~\\mathrm{Analytical}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )


linkyaxes!(ax1_F, ax2_F, ax3_F, ax1_Fa, ax2_Fa, ax3_Fa)
linkyaxes!(ax4_F, ax4_Fa)

lines!(ax1_F,time_dump,DATA_dump_1.fx,
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax1_F,time_dump,DATA_dump_1.fy,
        color=1,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )


lines!(ax2_F,time_dump,DATA_dump_2.fx,
       color=2,
       colormap=:tab10,
       colorrange=(1,10)
      )
lines!(ax2_F,time_dump,DATA_dump_2.fy,
       color=2,
       colormap=:tab10,
       colorrange=(1,10),
       linestyle=:dash
      )

lines!(ax3_F,time_dump,DATA_dump_3.fx,
       color=3,
       colormap=:tab10,
       colorrange=(1,10)
      )
lines!(ax3_F,time_dump,DATA_dump_3.fy,
       color=3,
       colormap=:tab10,
       colorrange=(1,10),
       linestyle=:dash
      )


lines!(ax1_Fa,time_dump,f_part1[:,1],
        color=1,
        colormap=:tab10,
        colorrange=(1,10)
       )
lines!(ax1_Fa,time_dump,f_part1[:,2],
        color=1,
        colormap=:tab10,
        colorrange=(1,10),
        linestyle=:dash
       )


lines!(ax2_Fa,time_dump,f_part2[:,1],
       color=2,
       colormap=:tab10,
       colorrange=(1,10)
      )
lines!(ax2_Fa,time_dump,f_part2[:,2],
       color=2,
       colormap=:tab10,
       colorrange=(1,10),
       linestyle=:dash
      )

lines!(ax3_Fa,time_dump,f_part3[:,1],
       color=3,
       colormap=:tab10,
       colorrange=(1,10)
      )
lines!(ax3_Fa,time_dump,f_part3[:,2],
       color=3,
       colormap=:tab10,
       colorrange=(1,10),
       linestyle=:dash
      )



# Norm 
lines!(ax4_F,time_dump,F_systemL,
        color=:black
       )

# Force system
lines!(ax4_Fa,time_dump,F_systemA,
        color=:black
       )

save("fig_FComp.png", fig_F, px_per_unit = 300/inch)

# Components of the force
fig_Fcomp=Figure(size = (17cm, 18.75cm));
ax1_F=Axis(fig_Fcomp[1,1],
             title=latexstring("\\mathrm{Force~components~Particle~1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_F=Axis(fig_Fcomp[2,1],
             title=latexstring("\\mathrm{Force~components~Particle~2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_F=Axis(fig_Fcomp[3,1],
             title=latexstring("\\mathrm{Force~components~Particle~3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
linkyaxes!(ax1_F, ax2_F, ax3_F)


#lines!(ax1_F,time_dump,f_func1,linestyle=:dash,label=L"\sum f_{1i}(r)",linewidth=4,alpha=0.5,color=:black)
lines!(ax1_F,time_dump,f_part1[:,1],label=L"f_{xt}",
        linewidth=5,
        color=1,
        colormap=:viridis,
        colorrange=(1,10)
       )
lines!(ax1_F,time_dump,f_part1[:,2],label=L"f_{xt}",
       linewidth=5,
       color=2,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax1_F,time_dump,DATA_dump_1.fx,label=L"f_x",
       linewidth=2.5,
       linestyle=:dash,
       color=10,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax1_F,time_dump,DATA_dump_1.fy,label=L"f_y",
       linewidth=2.5,
       linestyle=:dash,
       color=9,
       colormap=:viridis,
       colorrange=(1,10)
      )

#lines!(ax2_F,time_dump,f_func2,linestyle=:dash,label=L"\sum f_{2i}(r)",linewidth=4,alpha=0.5,color=:black)
lines!(ax2_F,time_dump,f_part2[:,1],label=L"f_{xt}",
       linewidth=5,
       color=1,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax2_F,time_dump,f_part2[:,2],label=L"f_{xt}",
       linewidth=5,
       color=2,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax2_F,time_dump,DATA_dump_2.fx,label=L"f_x",
       linewidth=2.5,
       linestyle=:dash,
       color=10,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax2_F,time_dump,DATA_dump_2.fy,label=L"f_y",
       linewidth=2.5,
       linestyle=:dash,
       color=9,
       colormap=:viridis,
       colorrange=(1,10)
      )



#lines!(ax3_F,time_dump,f_func3,linestyle=:dash,label=L"\sum f_{3i}(r)",linewidth=4,alpha=0.5,color=:black)
lines!(ax3_F,time_dump,f_part3[:,1],label=L"f_{xt}",
       linewidth=5,
       color=1,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax3_F,time_dump,f_part3[:,2],label=L"f_{yt}",
       linewidth=5,
       color=2,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax3_F,time_dump,DATA_dump_3.fx,label=L"f_x",
       linewidth=2.5,
       linestyle=:dash,
       color=10,
       colormap=:viridis,
       colorrange=(1,10)
      )
lines!(ax3_F,time_dump,DATA_dump_3.fy,label=L"f_y",
       linewidth=2.5,
       linestyle=:dash,
       color=9,
       colormap=:viridis,
       colorrange=(1,10)
      )



Legend(fig_Fcomp[1,2],ax1_F,
      L"\mathrm{Labels}",
     labelsize=0.5cm)
Legend(fig_Fcomp[2,2],ax2_F,
      L"\mathrm{Labels}",
     labelsize=0.5cm)
Legend(fig_Fcomp[3,2],ax3_F,
      L"\mathrm{Labels}",
     labelsize=0.5cm)



save("fig_ForceComp.png", fig_Fcomp, px_per_unit = 300/inch)


#=
"""
    Plot the Forces
"""
fig_Force=Figure(size = (17cm, 18.75cm));
ax1_F=Axis(fig_Force[1,1],
             title=latexstring("\\mathrm{Vectors}"),
             xlabel=L"x~[r/D_p]",
             ylabel=L"y~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )


lines!(ax1_F,DATA_dump_1.x,DATA_dump_1.y)
lines!(ax1_F,DATA_dump_2.x,DATA_dump_2.y)
lines!(ax1_F,DATA_dump_3.x,DATA_dump_3.y)

arrows2d!(ax1_F,DATA_dump_1.x,DATA_dump_1.y,DATA_dump_1.fx./DATA_dump_1.norma_f,DATA_dump_1.fy./DATA_dump_1.norma_f,lengthscale=0.2)
arrows2d!(ax1_F,DATA_dump_2.x,DATA_dump_2.y,DATA_dump_2.fx./DATA_dump_2.norma_f,DATA_dump_2.fy./DATA_dump_2.norma_f,lengthscale=0.2)
arrows2d!(ax1_F,DATA_dump_3.x,DATA_dump_3.y,DATA_dump_3.fx./DATA_dump_3.norma_f,DATA_dump_3.fy./DATA_dump_3.norma_f,lengthscale=0.2)
=#


"""
    Plot the threebody potential
"""

# Definir un mapa de color (puedes cambiarlo, ej: :thermal, :plasma, etc.)
cmap = :tokyo

fig_swapPot=Figure(size = (17cm, 18.75cm));
ax1_eng=Axis(fig_swapPot[1,1],
             title=latexstring("\\mathrm{Swap~potential~from~table}"),
             xlabel=L"r_{ik}~[r/D_p]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )

ax2_eng=Axis(fig_swapPot[2,1],
             title=latexstring("\\mathrm{Swap~potential~from~function}"),
             xlabel=L"r_{ik}~[r/D_p]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )

r_ij_values=unique(tableSwap.r_ij);
r_min, r_max = extrema(r_ij_values)

for (idx, df) in enumerate(table_pot1)
    r_val = r_ij_values[idx]
    # Normalizar el valor de r_ij al intervalo [0,1] para el mapeo de color
    color_norm = (r_val - r_min) / (r_max - r_min)
    scatterlines!(ax1_eng, df.r_ik, df.e,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1),
           markersize = 1.5)
    potTeo=map(s->Uswap(1,1,1,1,0.4,r_val,s),df.r_ik)

    scatterlines!(ax2_eng,df.r_ik,potTeo,
           color = color_norm,
           colormap = cmap,
           colorrange = (0,1),
           markersize=1.5
            )
end

# Añadir la barra de color (leyenda) a la derecha del eje
Colorbar(fig_swapPot[1,2],
         limits = (r_min, r_max),
         colormap = cmap,
         label = L"r_{ij}")   # etiqueta de la barra

# Añadir la barra de color (leyenda) a la derecha del eje
Colorbar(fig_swapPot[2,2],
         limits = (r_min, r_max),
         colormap = cmap,
         label = L"r_{ij}")   # etiqueta de la barra

save("fig_swapPotTable.png", fig_swapPot, px_per_unit = 300/inch)


# Parameters
dt=0.001;



"""
#    SANDBOX
"""

L=1;



"""
#    Energy and that stuff
"""

time=dt.*DATA_fix.TimeStep;

e_lim=2;



fig_U=Figure(size = (15cm, 18.75cm));
# Create twin axis: Energy and distance

# First axis
ax1_eng=Axis(fig_U[1,1],
             title=latexstring("\\mathrm{Particle~1}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax1_dist=Axis(fig_U[1,1],
             #title=latexstring("\\mathrm{Potential~energy}"),
             #xlabel=L"t~[\tau]",
             yticks = (0.0:0.2:2.0),
             ylabel=L"d~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             yaxisposition=:right
   )
hidespines!(ax1_dist)
hidexdecorations!(ax1_dist)

# Second axis
ax2_eng=Axis(fig_U[2,1],
             title=latexstring("\\mathrm{Particle~2}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_dist=Axis(fig_U[2,1],
             #title=latexstring("\\mathrm{Potential~energy}"),
             #xlabel=L"t~[\tau]",
             yticks = (0.0:0.2:2.0),
             ylabel=L"d~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             yaxisposition=:right
   )
hidespines!(ax2_dist)
hidexdecorations!(ax2_dist)

# Third axis
ax3_eng=Axis(fig_U[3,1],
             title=latexstring("\\mathrm{Particle~3}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax3_dist=Axis(fig_U[3,1],
             #title=latexstring("\\mathrm{Potential~energy}"),
             #xlabel=L"t~[\tau]",
             yticks = (0.0:0.2:2.0),
             ylabel=L"d~[r/D_p]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true,
             yaxisposition=:right
   )
hidespines!(ax3_dist)
hidexdecorations!(ax3_dist)




sl1_d12=lines!(ax1_dist,time,dist_12,linestyle=:dash,color=1,colormap=:tab10,colorrange=(1,10))
sl1_d13=lines!(ax1_dist,time,dist_13,linestyle=:dash,color=2,colormap=:tab10,colorrange=(1,10))
sl1_eg1=scatterlines!(ax1_eng,time,DATA_dump_1.c_potAtom,markersize=3,color=3, colormap=:tab10,colorrange=(1,10))

sl2_d12=lines!(ax2_dist,time,dist_12,linestyle=:dash,color=1,colormap=:tab10,colorrange=(1,10))
sl2_d23=lines!(ax2_dist,time,dist_23,linestyle=:dash,color=2,colormap=:tab10,colorrange=(1,10))
sl2_eg2=scatterlines!(ax2_eng,time,DATA_dump_2.c_potAtom,markersize=3,color=3, colormap=:tab10,colorrange=(1,10))

sl3_d12=lines!(ax3_dist,time,dist_13,linestyle=:dash,color=1,colormap=:tab10,colorrange=(1,10))
sl3_d23=lines!(ax3_dist,time,dist_23,linestyle=:dash,color=2,colormap=:tab10,colorrange=(1,10))
sl3_eg3=scatterlines!(ax3_eng,time,DATA_dump_3.c_potAtom,markersize=3,color=3, colormap=:tab10,colorrange=(1,10))


Legend(fig_U[1,2],
       [sl1_d12,sl1_d13,sl1_eg1],
       [L"d_{12}", L"d_{13}", L"U_{\mathrm{atom}}"],
     labelsize=0.5cm)

Legend(fig_U[2,2],
       [sl2_d12,sl2_d23,sl2_eg2],
       [L"d_{21}", L"d_{23}", L"U_{\mathrm{atom}}"],
     labelsize=0.5cm)

Legend(fig_U[3,2],
       [sl3_d12,sl3_d23,sl3_eg3],
       [L"d_{31}", L"d_{32}", L"U_{\mathrm{atom}}"],
     labelsize=0.5cm)


save("fig_Uatom.png", fig_U, px_per_unit = 300/inch)


"""
    Energy of the fix
"""

fig_Uf=Figure(size = (18.75cm, 15cm));
# Create twin axis: Energy and distance

# First axis
ax1_eng=Axis(fig_Uf[1,1],
             title=latexstring("\\mathrm{Potential~energy}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )

aux=combine(groupby(l, :TimeStep), :c_potAtom => sum => :suma_c_pot);

scatterlines!(ax1_eng,time_fix,DATA_fix.c_patchPair,label=L"U_{\mathrm{fixPatch}}",color=1,colormap=:tab10,colorrange=(1,10),markersize=1.5)
scatterlines!(ax1_eng,time_fix,DATA_fix.c_swapPair,label=L"U_{\mathrm{swap}}",color=2,colormap=:tab10,colorrange=(1,10),markersize=1.5)
scatterlines!(ax1_eng,time_fix,DATA_fix.c_ep,label=L"U_{\mathrm{ep}}",color=3,colormap=:tab10,colorrange=(1,10),markersize=1.5)


scatterlines!(ax1_eng,time_dump,DATA_dump_1.c_potAtom,label=L"U_{\mathrm{atom}1}",color=2,colormap=:viridis,colorrange=(1,10),markersize=1.5)
scatterlines!(ax1_eng,time_dump,DATA_dump_2.c_potAtom,label=L"U_{\mathrm{atom}2}",color=4,colormap=:viridis,colorrange=(1,10),markersize=1.5)
scatterlines!(ax1_eng,time_dump,DATA_dump_3.c_potAtom,label=L"U_{\mathrm{atom}3}",color=6,colormap=:viridis,colorrange=(1,10),markersize=1.5)
lines!(ax1_eng,time_dump,aux.suma_c_pot,label=L"\sum U_{\mathrm{atom}}",color=:black,linestyle=:dash)



# jj = 

Legend(fig_Uf[1,2],ax1_eng,
      L"\mathrm{Labels}",
     labelsize=0.5cm)

save("fig_Uatom-fix.png", fig_Uf, px_per_unit = 300/inch)


"""
    Energy from function evaluated at the distances given by the simulation 
"""
fig_UF=Figure(size = (18.75cm, 15cm));
# Create twin axis: Energy and distance
ax2_eng=Axis(fig_UF[1,1],
             title=latexstring("\\mathrm{Comparisson~of~fix~and~function}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )

"""
    Evaluate the swap potential
"""
swapEval1=map(s->Uswap(1.0,1.0,1.0,1.0,0.4,dist_12[s],dist_13[s]),eachindex(time_fix));
swapEval2=map(s->Uswap(1.0,1.0,1.0,1.0,0.4,dist_12[s],dist_23[s]),eachindex(time_fix));
swapEval3=map(s->Uswap(1.0,1.0,1.0,1.0,0.4,dist_13[s],dist_23[s]),eachindex(time_fix));

swapEvalsum=swapEval1 .+ swapEval2 .+ swapEval3;

patchEval1=map(r->Upatch(1,0.4,r),dist_12);
patchEval2=map(r->Upatch(1,0.4,r),dist_13);
patchEval3=map(r->Upatch(1,0.4,r),dist_23);

patchEvalsum= patchEval1 .+ patchEval2 .+ patchEval3

scatter!(ax2_eng,time_fix,DATA_fix.c_patchPair,label=L"U_{\mathrm{fixPatch}}",color=1,colormap=:tab10,colorrange=(1,10))
scatter!(ax2_eng,time_fix,DATA_fix.c_swapPair,label=L"U_{\mathrm{swap}}",color=2,colormap=:tab10,colorrange=(1,10))
scatter!(ax2_eng,time_fix,DATA_fix.c_ep,label=L"U_{\mathrm{ep}}",color=3,colormap=:tab10,colorrange=(1,10))



#lines!(ax2_eng,time_fix,swapEval1,label=L"U_{\mathrm{swap1}}",color=1,colormap=:darkrainbow,colorrange=(1,5))
#lines!(ax2_eng,time_fix,swapEval2,label=L"U_{\mathrm{swap2}}",color=2,colormap=:darkrainbow,colorrange=(1,5))
#lines!(ax2_eng,time_fix,swapEval3,label=L"U_{\mathrm{swap3}}",color=3,colormap=:darkrainbow,colorrange=(1,5))
#lines!(ax2_eng,time_fix,swapEvalsum,label=L"\sum U_{\mathrm{swap}2,3}",color=5,colormap=:darkrainbow,colorrange=(1,5))

#lines!(ax2_eng,time_dump,,label=L"U_{\mathrm{patch}}(d_{12})",color=4,colormap=:tab10,colorrange=(1,10))
#lines!(ax2_eng,time_dump,,label=L"U_{\mathrm{patch}}(d_{13})",color=5,colormap=:tab10,colorrange=(1,10))
#lines!(ax2_eng,time_dump,map(r->Upatch(1,0.4,r),dist_23),label=L"U_{\mathrm{patch}}(d_{23})",color=6,colormap=:tab10,colorrange=(1,10))

lines!(ax2_eng,time_fix,swapEvalsum,label=L"\sum U_{\mathrm{swap}}",color=:black,linestyle=:dash)
lines!(ax2_eng,time_dump,patchEvalsum,label=L"\sum U_{\mathrm{patch}}",color=:black,linestyle=:dot)
lines!(ax2_eng,time_dump,patchEvalsum .+ swapEvalsum,label=L"\sum U_{\mathrm{patch}} + \sum U_{\mathrm{swap}}",color=:black)




Legend(fig_UF[1,2],ax2_eng,
      L"\mathrm{Labels}",
     labelsize=0.5cm)

save("fig_Ufix-func.png", fig_UF, px_per_unit = 300/inch)


"""
    Numeric derivative
"""

fig_UForce=Figure(size = (18.75cm, 15cm));
# Create twin axis: Energy and distance
ax1_eng=Axis(fig_UForce[1,1],
             title=latexstring("\\mathrm{Total~potencial~energy~of~the~system}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_eng=Axis(fig_UForce[2,1],
             title=latexstring("\\mathrm{Derivative~of~the~total~potencial~energy~of~the~system}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )



forcel=(-1).*(DATA_fix.c_ep[2:end] .- DATA_fix.c_ep[1:end-1])./(time_fix[2:end] .- time_fix[1:end-1]);


scatterlines!(ax1_eng,time_fix[1:end-1],DATA_fix.c_ep[1:end-1],label=L"U_{\mathrm{ep}}",color=3,colormap=:tab10,colorrange=(1,10),markersize=2.5)

scatterlines!(ax2_eng,time_fix[1:end-1],forcel,markersize=2.5,color=3,colormap=:tab10,colorrange=(1,10))

save("fig_UforceDErivative.png", fig_UForce, px_per_unit = 300/inch)



ind_aux1=2001:length(time_fix)-2;
ind_aux=2001:length(forcel)-1;
fig_UForce=Figure(size = (18.75cm, 15cm));
# Create twin axis: Energy and distance
ax1_eng=Axis(fig_UForce[1,1],
             title=latexstring("\\mathrm{Total~potencial~energy~of~the~system}"),
             xlabel=L"t~[\tau]",
             ylabel=L"U~[J/\epsilon]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )
ax2_eng=Axis(fig_UForce[2,1],
             title=latexstring("\\mathrm{Derivative~of~the~total~potencial~energy~of~the~system}"),
             xlabel=L"t~[\tau]",
             ylabel=L"F~[N^*]",
             titlesize=tl_sz,
             xticklabelsize=ot_sz,
             yticklabelsize=ot_sz,
             xlabelsize=tl_sz,
             ylabelsize=tl_sz,
             xminorticksvisible=true,
             xminorgridvisible=true
   )

scatterlines!(ax1_eng,time_fix[ind_aux1],DATA_fix.c_ep[ind_aux],label=L"U_{\mathrm{ep}}",color=3,colormap=:tab10,colorrange=(1,10),markersize=2.5)

scatterlines!(ax2_eng,time_fix[ind_aux1],forcel[ind_aux],markersize=2.5,color=3,colormap=:tab10,colorrange=(1,10))






save("fig_UForceZoom.png", fig_UForce, px_per_unit = 300/inch)


