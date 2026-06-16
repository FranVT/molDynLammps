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

# Write variables to the file
for phi in "${VAR_NAMES[@]}"; do
  # Check if the variable is set (exists and not empty)
  if [ -n "${!var}" ]; then
    echo "$var=${!var}" >> "$OUTPUT_FILE"
  else
    echo "Warning: $var is not set. Skipping." >&2
  fi
done


for var_phi in $(seq $phi_o $phi_d $phi_f);
do
    for var_Nexp in $(seq 1 1 $NexpT);
    do
        # Create the id for the simulation
        id="$(date +%F-%H%M%S)_${var_Nexp}";

        # Directory to save the information Directory to save the informationn
        dir_save="$dir_data/$id";
    

    done
done


# Compute amount of patchy particles of each type
N_pCL=$(echo "scale=0; $chi_4 * $N_PP" | bc);
N_pCL=${N_pCL%.*};
N_pMO=$(( $N_PP - $N_pCL ));

# Box size given a concentration and number of particles
VolT_MO=$(echo "scale=$cs; $Vol_MO * $N_MO" | bc);         # Vol of N f=2 patchy particles
VolT_CL=$(echo "scale=$cs; $Vol_CL * $N_pCL" | bc);          # Vol of N f=4 patchy particles 
Vol_Totg=$(echo "scale=$cs; $VolT_MO + $VolT_CL" | bc);       # Total volume of a mixture of N f=2 and M f=4 patchy particles
Vol_Tot=$(echo "scale=$cs; $Vol_Totg / $phi" | bc);
L_real=$(echo "scale=$cs; e( (1/3) * l($Vol_Tot) )" | bc -l );
L=$(echo "scale=$cs; $L_real / 2" | bc);


# ------------------------------------------------------------------
# 2. Overwrite the parameters given the cyle 
# ------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --L) export L="$2"; shift 2 ;;
        --temp) export temp="$2"; shift 2 ;;
        --N_pCL) export N_pCL="$2"; shift 2 ;;
        --N_pMO) export N_pMO="$2"; shift 2 ;;
        --E_CP) export E_CP="$2"; shift 2 ;;
        --r_cPP) export r_cPP="$2"; shift 2 ;;
        --seed1) export seed1="$2"; shift 2 ;;
        --seed2) export seed2="$2"; shift 2 ;;
        --seed3) export seed3="$2"; shift 2 ;;
        *) echo "Unkwon parameter: $1"; exit 1 ;;
    esac
done

# ------------------------------------------------------------------
# 3. Confirm the main parameters 
# ------------------------------------------------------------------
echo "Parámetros activos para esta simulación:"
echo "   L       = $L"
echo "   temp    = $temp"
echo "   N_pCL   = $N_pCL"
echo "   N_pMO   = $N_pMO"
echo "   E_CP    = $E_CP"
echo "   r_cPP   = $r_cPP"
echo "   seed1   = $seed1"
echo "   seed2   = $seed2"
echo "   seed3   = $seed3"



# Ejecutar la simulacion
#mpirun -np 8 lmp -sf omp -in src_lmp/in.assembly.lmp -var r_cWCA $r_cWCA -var N_3body $N_3body -var L $L -var m_CP $m_CP -var m_PT $m_PT -var E_bondP $E_bondP -var r_bondP $r_bondP -var theta_PA $theta_PA -var theta_PB $theta_PB -var r_noInter $r_noInter -var r_CP $r_CP -var E_CP $E_CP -var r_cPP $r_cPP -var N_pCL $N_pCL -var N_pMO $N_pMO -var seed1 $seed1 -var seed2 $seed2 -var N_save $N_save -var N_dump $N_dump -var Dir $dir_save -var file1_name $file1_name -var file2_name $file2_name -var file3_name $file3_name -var tstep $tstep -var temp $temp -var damp $damp -var seed3 $seed3 -var N_heat $N_heat -var N_isothermal $N_isothermal -var file4_name $file4_name
