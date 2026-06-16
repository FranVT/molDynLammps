#!/bin/bash

# Main Directories
dir_home="/home/franvt/GitRepos/molDynLammps/polymericGels/gel-batch-3";
dir_src="$dir_home/src";
dir_sim="$dir_src/sim";
dir_data="$dir_home/data";

# Create the id for the simulation
id="$(date +%F-%H%M%S)";

# Directory to save the informatio# Directory to save the informationn
dir_save="$dir_data/$id";

# Create directory to save data
mkdir "$dir_save";
mkdir "$dir_save/traj";

# Parameters
source parameters.sh

# Ejecutar la simulacion
mpirun -np 8 lmp -sf omp -in src_lmp/in.assembly.lmp -var r_cWCA $r_cWCA -var N_3body $N_3body -var L $L -var m_CP $m_CP -var m_PT $m_PT -var E_bondP $E_bondP -var r_bondP $r_bondP -var theta_PA $theta_PA -var theta_PB $theta_PB -var r_noInter $r_noInter -var r_CP $r_CP -var E_CP $E_CP -var r_cPP $r_cPP -var N_pCL $N_pCL -var N_pMO $N_pMO -var seed1 $seed1 -var seed2 $seed2 -var N_save $N_save -var N_dump $N_dump -var Dir $dir_save -var file1_name $file1_name -var file2_name $file2_name -var file3_name $file3_name -var tstep $tstep -var temp $temp -var damp $damp -var seed3 $seed3 -var N_heat $N_heat -var N_isothermal $N_isothermal -var file4_name $file4_name
