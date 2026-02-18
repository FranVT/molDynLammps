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
r_jk=(r_min:(r_max-r_min)/(N-1):r_max);



# Ranges: 3 body interaction
r_ijVal=last.(Iterators.partition(eachindex(r_ij),2^4));

Uswap_eval1=map(s->map(r->Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij[s],r),r_ik),r_ijVal);
Uswap_eval2=map(s->map(r->Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ij[s],r),r_jk),r_ijVal);
Uswap_eval3=map(s->map(r->Uswap(w,eps_jk,eps_ij,eps_ik,Dp,r_ik[s],r),r_jk),r_ijVal);



