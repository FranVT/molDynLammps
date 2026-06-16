#!/bin/bash


# Main Directories
dir_home="/home/franvt/GitRepos/molDynLammps/polymericGels/gel-batch-3";
dir_src="$dir_home/src";
dir_sim="$dir_src/sim";
dir_data="$dir_home/data";

# Create the id for the simulation
id="$(date +%F-%H%M%S)";
#"${T}${phi}${CL_con}${N_particles}-$(date +%F-%H%M%S)";

# Directory to save the informatio# Directory to save the informationn
dir_save="$dir_data/$id";

# Create directory to save data
mkdir "$dir_save";
mkdir "$dir_save/traj";

# Parametros
r_cWCA=1.12     # Cutoff distance of WCA potential 
N_3body=128     # Points of evaluations for 3body table
L=5.0             # Length of the simulation box
m_CP=1.0          # Mass of central particles
m_PT=1.0          # Mass of the patches

E_bondP=100.0      # Bond energy between central particles and patches
r_bondP=0.4         # Distance of bond
theta_PA=109.4712          # Angle between patches of monomers
theta_PB=180.0          # Angle between patches of monomers
r_noInter=0.5       # Distance of no interacion zero potential lammps
r_CP=1.0            # Distance of interaction between central particles
E_CP=1.0            # Energy of interaction between central particles
r_cPP=0.6           # Cutt of distance of patch-patch interaction

N_pCL=2                 # Amount of crosslinkers particles
N_pMO=2                 # Amount of monomers 
seed1=1                 # Seed for random initial position
seed2=2                 # Seed for random initial position

N_save=1000               # Create temporal averages
N_dump=2000             # Frequency for saving dump files
#Dir="/home/franvt/GitRepos/molDynLammps/polymericGels/gel-bath-3/data"              # Directory to store the information
file1_name="system_assembly.fixf"          # Name for fix file
file2_name="stress_assembly.fixf"       # Name for stress tensor file
file3_name="traj_assembly.*.dumpf"       # Name for dump files

tstep=0.001         # Time step
temp=0.05           # Temperature of thermal bath
damp=1.0            # Damp, related with viscosity
seed3=3             # Seed for the thermal bath 
N_heat=5000      # Time steps for going from 0 to temp
N_isothermal=10000    # Time steps for isothermal process 

file4_name="data.hydrogel"       # Name for the last configuration and restart


# Ejecutar la simulacion
mpirun -np 8 lmp -sf omp -in src_lamp/in.assembly.lmp -var $r_cWCA r_cWCA -var $N_3body N_3body -var $L L -var $m_CP m_CP -var $m_PT m_PT -var $E_bondP E_bondP -var $r_bondP r_bondP -var $theta_PA theta_PA -var $theta_PB theta_PB -var $r_noInter r_noInter -var $r_CP r_CP -var $E_CP E_CP -var $r_cPP r_cPP -var $N_pCL N_pCL -var $N_pMO N_pMO -var $seed1 seed1 -var $seed2 seed2 -var $N_save N_save -var $N_dump N_dump -var $Dir dir_save -var $file1_name file1_name -var $file2_name file2_name -var $file3_name file3_name -var $tstep tstep -var $temp temp -var $damp damp -var $seed3 seed3 -var $N_heat N_heat -var $N_isothermal N_isothermal -var $file4_name file4_name

