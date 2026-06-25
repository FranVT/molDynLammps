#!/bin/bash

# ------------------------------------------------------------------
# 1. Directorios principales
# ------------------------------------------------------------------
dir_home="/home/franvt/GitRepos/molDynLammps/polymericGels/gel-batch-3"
dir_src="$dir_home/src_lmp"
dir_data="$dir_home/data"

# ------------------------------------------------------------------
# 2. Cargar parámetros por defecto
# ------------------------------------------------------------------
source parameters.sh
source utils.sh   # <--- asegúrate de que la ruta sea correcta

# ------------------------------------------------------------------
# 4. Barrido de parámetros
# ------------------------------------------------------------------
for phi in $(seq $phi_o $phi_delta $phi_f); do
    for chi_4 in $(seq $chi_4o $chi_4delta $chi_4f); do
        # --- 4d. Ejecutar N_exp réplicas ---
        for (( exp=1; exp<=N_exp; exp++ )); do
            srun -p interactive --pty bash "$dir_home/job_script.sh" "$phi" "$chi_4" "$exp"
            # Pequeña pausa para evitar colisiones en fechas
            sleep 1
        done
    done
done

echo "Todas las simulaciones han enviado"



