: '
    Script that runs the assembly procotol with given parameters.
'

#!/bin/bash

# Main Directories
dir_home="/home/franvt/GitRepos/molDynLammps/polymericGels/gel-bonds-1";
dir_src="$dir_home/src";
dir_sim="$dir_src/sim";
dir_data="$dir_home/data";
dir_save="$dir_data/$id";


# Load the parameters file
source load_parameters.sh system.parameters

# Create the id for the simulation
id="${T}${phi}${CL_con}${N_particles}-$(date +%F-%H%M%S)";


# Create the config fle
bash createConfigFile_assembly.sh $dir_src $id

# Compute values for dependent variables
#N_CL=$(echo "scale=0; $CL_con * $N_particles" | bc);
#N_CL=${N_CL%.*};
#N_MO=$(( $N_particles - $N_CL ));
#Vol_MO=$(echo "scale=$cs; $Vol_MO1 * $N_MO" | bc);         # Vol of N f=2 patchy particles
#Vol_CL=$(echo "scale=$cs; $Vol_CL1 * $N_CL" | bc);          # Vol of N f=4 patchy particles 
#Vol_Totg=$(echo "scale=$cs; $Vol_MO + $Vol_CL" | bc);       # Total volume of a mixture of N f=2 and M f=4 patchy particles
#Vol_Tot=$(echo "scale=$cs; $Vol_Totg / $phi" | bc);
#L_real=$(echo "scale=$cs; e( (1/3) * l($Vol_Tot) )" | bc -l );
#L=$(echo "scale=$cs; $L_real / 2" | bc);
#aux2=$(echo "scale=$cs; 1 / $damp" | bc);
#aux=$(echo "scale=$cs; 1 / $dt" | bc);
#Nsave=$(echo "scale=0; $aux2 * $aux" | bc);
#Nsave=${Nsave%.*};
#NsaveStress=$(echo "scale=0; $aux2 * $aux" | bc);
#NsaveStress=${NsaveStress%.*};
#Ndump=$(echo "scale=0; $aux" | bc);
#Ndump=${Ndump%.*};

# Create the data file for the assembly simulation
bash  createDatFile_assembly.sh $dir_system $dir_src $id


# Assign values to file management variables

#mkdir "$dir_save";
#mkdir "$dir_save/traj";
files_name=("system_assembly.fixf" "stress_assembly.fixf" "clust_assembly.fixf" "profiles_assembly.fixf" "traj_assembly.*.dumpf" "data.hydrogel" "system_shear.fixf" "stress_shear.fixf" "clust_shear.fixf" "profiles_shear.fixf" "traj_shear.*.dumpf" "strain-*.restart")

# Create the data file for the assembly simulation
bash create-file_dataAssembly.sh $dir_save $dir_src $var_ccL $id

#mpirun -np 8 lmp -sf omp -in in.assembly.lmp -log $log_name -var temp $T -var damp $damp -var L $L -var NCL $N_CL -var NMO $N_MO -var seed1 $seed1 -var seed2 $seed2 -var seed3 $seed3 -var tstep $dt -var Nsave $Nsave -var NsaveStress $NsaveStress -var Ndump $Ndump -var steps $steps_isot -var stepsheat $steps_heat -var Dir $dir_save -var file1_name ${files_name[0]} -var file2_name ${files_name[1]} -var file3_name ${files_name[2]} -var file4_name ${files_name[3]} -var file5_name ${files_name[4]} -var file6_name ${files_name[5]};
