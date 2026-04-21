#!/bin/bash

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Configuración para el cluster
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#SBATCH --job-name=sims
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1

. /etc/profile.d/modules.sh	# Inicializa el sistema de módulos de enetorno en la terminal

# Cargar modulos
module load gcc/11.5.0 			# Carga el compilador
module load openmpi/4.1.6-gcc-11.5.0 	# Poder paralelizar simulaciones
module load lammps/gcc/22jul25		# Cargar Lammps version 22jul25

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Script
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Main Directories
dir_home="/mnt/data/cferreiro/fvazquez/experiments/gel-batch-2";
dir_src="$dir_home/src";
dir_sim="$dir_src/sim";
dir_data="$dir_home/data";

# Load the parameters file
source load_parameters.sh system.parameters

for var_phi in $(seq $phi_o $phi_d $phi_f);
do
    for var_Nexp in $(seq 1 1 $NexpT);
    do
        # Create the id for the simulation
        id="$(date +%F-%H%M%S)_${var_Nexp}";

        # Directory to save the information Directory to save the informationn
        dir_save="$dir_data/$id";

        # Create the config fle
        bash createConfigFile_assembly.sh $dir_sim $id $var_phi $var_Nexp

        # Create directory to save data
        mkdir "$dir_save";
        mkdir "$dir_save/traj";

        # Create the data file for the assembly simulation
        bash  createDatFile_assembly.sh $dir_save $dir_src $id

        # Load the config file
        source load_parameters.sh assembly$id.config 

        mpirun -np 4 env OMP_NUM_THREADS=1 lmp -sf omp -in in.assembly.lmp -var temp $T -var damp $damp -var L $L -var NCL $N_CL -var NMO $N_MO -var seed1 $seed1 -var seed2 $seed2 -var seed3 $seed3 -var tstep $dt -var Nsave $Nsave -var NsaveStress $NsaveStress -var Ndump $Ndump -var steps $steps_isot -var stepsheat $steps_heat -var Dir $dir_save -var file1_name ${files_name[0]} -var file2_name ${files_name[1]} -var file3_name ${files_name[2]} -var file4_name ${files_name[3]};

    done
done
