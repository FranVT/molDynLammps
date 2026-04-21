#!/bin/bash

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
# Script
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Main Directories
dir_home="/mnt/data/cferreiro/fvazquez/experiments/gel-batch-2";
dir_src="$dir_home/src";
dir_sim="$dir_src/sim";
dir_data="$dir_home/data";

# Remove old config files
rm assembly*.config

# Remove old slurm files
rm slurm*.out


# Load the parameters file
source load_parameters.sh system.parameters

for var_phi in $(seq $phi_o $phi_d $phi_f);
do
    for var_Nexp in $(seq 1 1 $NexpT);
    do
        sbatch "$dir_sim/job_script.sh" "$var_phi" "$var_Nexp"
        sleep 2
    done
done
