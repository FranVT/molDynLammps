#!/bin/bash
# Archivo de configuración de parámetros para LAMMPS

export r_cWCA=1.12          # Cutoff distance of WCA potential
export N_3body=128          # Points of evaluations for 3body table
export L=5.0                # Length of the simulation box
export m_CP=1.0             # Mass of central particles
export m_PT=1.0             # Mass of the patches
export E_bondP=100.0        # Bond energy between central particles and patches
export r_bondP=0.4          # Distance of bond
export theta_PA=109.4712    # Angle between patches of monomers
export theta_PB=180.0       # Angle between patches of monomers
export r_noInter=0.5        # Distance of no interaction zero potential lammps
export r_CP=1.0             # Distance of interaction between central particles
export E_CP=1.0             # Energy of interaction between central particles
export r_cPP=0.6            # Cutoff distance of patch-patch interaction
export N_pCL=2              # Amount of crosslinkers particles
export N_pMO=2              # Amount of monomers
export seed1=1              # Seed for random initial position
export seed2=2              # Seed for random initial position
export N_save=1000          # Create temporal averages
export N_dump=2000          # Frequency for saving dump files
export tstep=0.001          # Time step
export temp=0.05            # Temperature of thermal bath
export damp=1.0             # Damp, related with viscosity
export seed3=3              # Seed for the thermal bath
export N_heat=5000          # Time steps for going from 0 to temp
export N_isothermal=10000   # Time steps for isothermal process

