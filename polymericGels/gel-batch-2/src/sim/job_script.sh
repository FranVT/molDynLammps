#!/bin/bash

phi=$1
Nexp=$2

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Configuración para el cluster
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#SBATCH --job-name=sims
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

. /etc/profile.d/modules.sh	# Inicializa el sistema de módulos de enetorno en la terminal

# Cargar modulos
module load gcc/11.5.0 			# Carga el compilador
module load openmpi/4.1.6-gcc-11.5.0 	# Poder paralelizar simulaciones
module load lammps/gcc/22jul25		# Cargar Lammps version 22jul25

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Script
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load the parameters file
source load_parameters.sh system.parameters

# Main Directories
dir_home="/mnt/data/cferreiro/fvazquez/experiments/gel-batch-2";
dir_src="$dir_home/src";
dir_sim="$dir_src/sim";
dir_data="$dir_home/data";

# Create the id for the simulation
id="$(date +%F-%H%M%S)_${Nexp}";

# Directory to save the information Directory to save the informationn
dir_save="$dir_data/$id";

echo -e "\n$dir_save\n"

        # Create the config fle
bash createConfigFile_assembly.sh $dir_sim $id $phi $Nexp

        # Create directory to save data
mkdir "$dir_save";
mkdir "$dir_save/traj";

        # Create the data file for the assembly simulation
bash  createDatFile_assembly.sh $dir_save $dir_src $id

        # Load the config file
source load_parameters.sh assembly$id.config 

env OMP_NUM_THREADS=1 lmp -sf omp -in in.assembly.lmp -var Dir "$dir_save" -var L "$L" -var seed1 "$seed1" -var seed2 "$seed2" -var NCL "$N_CL" -var NMO "$N_MO" -var Nsave "$Nsave" -var file1_name "${files_name[0]}" -var NsaveStress "$NsaveStress" -var file2_name "${files_name[1]}" -var Ndump "$Ndump" -var file3_name "${files_name[2]}" -var tstep "$dt" -var temp "$T" -var damp "$damp" -var seed3 "$seed3" -var stepsheat "$steps_heat" -var steps "$steps_isot" -var file4_name "${files_name[3]}"
