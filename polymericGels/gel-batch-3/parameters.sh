#!/bin/bash

# Control parameters
export phi_o=0.01           # Inital packing fraction 
export phi_f=0.05           # Final packing fraction
export phi_delta=0.01       # Delta of packing fraction simulate
export N_PP=5               # Total amount of patchy particles
export N_exp=1              # Amount of experiments per each configuration
export chi_4o=0.1           # Initial Crosslinker concentration
export chi_4f=0.1           # Final Crosslinker concentration
export chi_4delta=0.1       # Delta of Crosslinker concentration

# System parameters
export temp=0.05            # Temperature of thermal bath
export damp=1.0             # Damp, related with viscosity
export N_pCL=1              # Amount of crosslinkers particles Initite value
export N_pMO=1              # Amount of monomers Initite value
export L=1                  # Length of the simulation box Initite value

# Physics parameters
export m_CP=1.0             # Mass of central particles
export m_PT=1.0             # Mass of the patches
export r_cWCA=1.12          # Cutoff distance of WCA potential
export r_cPP=0.6            # Cutoff distance of patch-patch interaction
export r_CP=1.0             # Distance of interaction between central particles
export r_bondP=0.4          # Distance of bond
export theta_PA=109.4712    # Angle between patches of monomers
export theta_PB=180.0       # Angle between patches of monomers
export E_CP=1.0             # Energy of interaction between central particles
export E_bondP=100.0        # Bond energy between central particles and patches

# Analysis parameters
export N_heat=5000          # Time steps for going from 0 to temp
export N_isothermal=10000   # Time steps for isothermal process
export N_save=1000          # Create temporal averages
export N_dump=2000          # Frequency for saving dump files

# Lammps parameters
export N_3body=128          # Points of evaluations for 3body table
export r_noInter=0.5        # Distance of no interaction zero potential lammps
export seed1=1              # Seed for random initial position
export seed2=2              # Seed for random initial position
export seed3=3              # Seed for the thermal bath
export tstep=0.001          # Time step

# File names and stuff
export file1_name="system_assembly.fixf"          # Name for fix file
export file2_name="stress_assembly.fixf"       # Name for stress tensor file
export file3_name="traj_assembly.*.dumpf"       # Name for dump files
export file4_name="data.hydrogel"       # Name for the last configuration and restart

# Extra
export cs=6                 # Significant digits
export Vol_MO=0.5376050428  # Volume of the monomer
export Vol_CL=0.5516113101  # Volume of the CrossLinker
