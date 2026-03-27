: '
    Script that runs the assembly procotol with given parameters.
'

#!/bin/bash

# Main Directories
dir_home="/home/franvt/GitRepos/molDynLammps/polymericGels/gel-bonds-1";
dir_src="$dir_home/src";
dir_sim="$dir_src/sim";
dir_data="$dir_home/data";


# Load the parameters file
source load_parameters.sh system.parameters

# Create the id for the simulation
id="$(date +%F-%H%M%S)";
#"${T}${phi}${CL_con}${N_particles}-$(date +%F-%H%M%S)";

# Directory to save the informatio# Directory to save the informationn
dir_save="$dir_data/$id";

# Create the config fle
bash createConfigFile_assembly.sh $dir_sim $id

# Create directory to save data
mkdir "$dir_save";
mkdir "$dir_save/traj";

# Create the data file for the assembly simulation
bash  createDatFile_assembly.sh $dir_save $dir_src $id

# Load the config file
source load_parameters.sh assembly$id.config 

mpirun -np 8 lmp -sf omp -in in.assembly.lmp -var temp $T -var damp $damp -var L $L -var NCL $N_CL -var NMO $N_MO -var seed1 $seed1 -var seed2 $seed2 -var seed3 $seed3 -var tstep $dt -var Nsave $Nsave -var NsaveStress $NsaveStress -var Ndump $Ndump -var steps $steps_isot -var stepsheat $steps_heat -var Dir $dir_save -var file1_name ${files_name[0]} -var file2_name ${files_name[1]} -var file3_name ${files_name[2]} -var file4_name ${files_name[3]};
