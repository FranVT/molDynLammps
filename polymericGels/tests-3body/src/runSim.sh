: '
    Script that runs the assembly procotol with given parameters.
'

#!/bin/bash

# Assign values to the independent variables
T=0.05;
NPA=3;
NPB=3;
L=1;

# Assign computational parameters
damp=1;
dt=0.001;
seed1=3234;
seed2=6321;
seed3=2010;

steps_isot=50000;
steps_heat=30000;

log_name="log-3btest.log";

# Assign significant decimals
cs=6;

# Compute values for dependent variables
aux2=$(echo "scale=$cs; 1 / $damp" | bc);
aux=$(echo "scale=$cs; 1 / $dt" | bc);
Nsave=$(echo "scale=0; $aux2 * $aux" | bc);
Nsave=${Nsave%.*};
NsaveStress=$(echo "scale=0; $aux2 * $aux" | bc);
NsaveStress=${NsaveStress%.*};
Ndump=$(echo "scale=0; $aux" | bc);
Ndump=${Ndump%.*};

# Assign values to file management variables
dir_data="/home/franvt/GitRepos/molDynLammps/polymericGels/tests-3body/data";
id="${T}${NPA}${NPB}${damp}-$(date +%F-%H%M%S)";
dir_save="$dir_data/$id";
mkdir "$dir_save";
mkdir "$dir_save/traj";
files_name=("system_assembly.fixf" "stress_assembly.fixf" "clust_assembly.fixf" "profiles_assembly.fixf" "traj_assembly.*.dumpf" "data.hydrogel" "system_shear.fixf" "stress_shear.fixf" "clust_shear.fixf" "profiles_shear.fixf" "traj_shear.*.dumpf" "strain-*.restart")

mpirun -np 8 lmp -sf omp -in in.assembly.lmp -log $log_name -var temp $T -var damp $damp -var L $L -var NPA $NPA -var NPB $NPB -var seed1 $seed1 -var seed2 $seed2 -var seed3 $seed3 -var tstep $dt -var Nsave $Nsave -var NsaveStress $NsaveStress -var Ndump $Ndump -var steps $steps_isot -var stepsheat $steps_heat -var Dir $dir_save -var file1_name ${files_name[0]} -var file2_name ${files_name[1]} -var file3_name ${files_name[2]} -var file4_name ${files_name[3]} -var file5_name ${files_name[4]} -var file6_name ${files_name[5]};
