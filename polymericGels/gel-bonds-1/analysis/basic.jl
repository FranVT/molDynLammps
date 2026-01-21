"""
    Basic analysis script
"""

using DataFrames, CSV
#using Plots, LaTeXStrings, Plots.PlotMeasures
#gr()
using Statistics
using GLMakie, LaTeXStrings, Typst_jll

# Load the functions
include("functions.jl")


# Parameter to select the system
T=0.05;
N_particles=500;
phi=0.3;
CL_con=0.05;

# Selection of an specific simulation
date="2026-01-21-161335";
#"2026-01-21-133721";
#"2026-01-20-151755";
#"2026-01-20-143651";
#"2026-01-20-135923";
#"2026-01-20-121005"; 
#"2026-01-20-111704";

#0.050.30.05500-2026-01-20-121005

# Get the directory of the desire system
(DIR,id_c)=getDir(T,N_particles,phi,CL_con,date);
DIR=DIR[1];


# Filename with the simulation data
FILE_NAME="system_assembly.fixf";

# Extract the data from the file
data=extractFixScalar(DIR,FILE_NAME);

# Convert the array into a DataFrame
DATA=DataFrame(data[2]',data[1]);

"""
    P L O T S 
"""

# these are relative to 1 CSS px
inch = 96
pt = 4/3
cm = inch / 2.54


damp=1;
dt=0.001;

fig_temp=fig_Temp(cm, pt, dt,DATA,T);
fig_eng=fig_Eng(cm, pt, dt,DATA);
fig_engB=fig_EngB(cm, pt, dt,DATA);
fig_engSys=fig_EngSys(cm, pt, dt,DATA);
fig_engPair=fig_EngPair(cm, pt, dt,DATA)

save(joinpath(DIR,string("temp-",id_c,".png")), fig_temp, px_per_unit = 300/inch)
save(joinpath(DIR,string("epk-",id_c,".png")), fig_eng, px_per_unit = 300/inch)
save(joinpath(DIR,string("eB-",id_c,".png")), fig_engB, px_per_unit = 300/inch)
save(joinpath(DIR,string("eSys-",id_c,".png")), fig_engSys, px_per_unit = 300/inch)
save(joinpath(DIR,string("ePair-",id_c,".png")), fig_engPair, px_per_unit = 300/inch)




