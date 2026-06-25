#!/bin/bash

# ------------------------------------------------------------------
# 1. Directorios principales
# ------------------------------------------------------------------
dir_home="/home/franvt/GitRepos/molDynLammps/polymericGels/gel-batch-3"
dir_src="$dir_home/src"
dir_sim="$dir_src/sim"
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

        # --- 4a. Calcular número de partículas de cada tipo ---
        N_pCL=$(echo "scale=0; $chi_4 * $N_PP" | bc)
        N_pCL=${N_pCL%.*}          # eliminar decimales (entero)
        N_pMO=$(( N_PP - N_pCL ))

        # --- 4b. Calcular tamaño de caja L ---
        VolT_MO=$(echo "scale=$cs; $Vol_MO * $N_pMO" | bc)
        VolT_CL=$(echo "scale=$cs; $Vol_CL * $N_pCL" | bc)
        Vol_Totg=$(echo "scale=$cs; $VolT_MO + $VolT_CL" | bc)
        Vol_Tot=$(echo "scale=$cs; $Vol_Totg / $phi" | bc)
        L_real=$(echo "scale=$cs; e( (1/3) * l($Vol_Tot) )" | bc -l)
        L=$(echo "scale=$cs; $L_real / 2" | bc)

        # Exportar variables para que estén disponibles en LAMMPS
        export L N_pCL N_pMO

        # --- 4d. Ejecutar N_exp réplicas ---
        for (( exp=1; exp<=N_exp; exp++ )); do

        # --- 4c. Mostrar combinación actual ---
        echo "================================================"
        echo "Nexp  = $exp"
        echo "phi   = $phi"
        echo "chi_4 = $chi_4"
        echo "N_pCL = $N_pCL"
        echo "N_pMO = $N_pMO"
        echo "L     = $L"
        echo "================================================"

            # Generar semillas únicas
            generate_seeds $exp

            # Crear identificador único (fecha + número de experimento)
            id="$(date +%F-%H%M%S)_${exp}_phi${phi}_chi${chi_4}"
            dir_save="$dir_data/$id"

            # Crear directorios
            mkdir -p "$dir_save/traj"

            # Crea el config
            write_params "$dir_save/params.log" \
                "phi=$phi" \
                "chi_4=$chi_4" \
                "N_pCL=$N_pCL" \
                "N_pMO=$N_pMO" \
                "L=$L" \
                "seed1=$seed1" \
                "seed2=$seed2" \
                "seed3=$seed3" 

            # ------------------------------------------------------------------
            # 4e. Ejecutar LAMMPS con los parámetros actualizados
            # ------------------------------------------------------------------
            # Ejecutar la simulacion
            mpirun -np 8 lmp -sf omp -in src_lmp/in.assembly.lmp -var r_cWCA $r_cWCA -var N_3body $N_3body -var L $L -var m_CP $m_CP -var m_PT $m_PT -var E_bondP $E_bondP -var r_bondP $r_bondP -var theta_PA $theta_PA -var theta_PB $theta_PB -var r_noInter $r_noInter -var r_CP $r_CP -var E_CP $E_CP -var r_cPP $r_cPP -var N_pCL $N_pCL -var N_pMO $N_pMO -var seed1 $seed1 -var seed2 $seed2 -var N_save $N_save -var N_dump $N_dump -var Dir $dir_save -var file1_name $file1_name -var file2_name $file2_name -var file3_name $file3_name -var tstep $tstep -var temp $temp -var damp $damp -var seed3 $seed3 -var N_heat $N_heat -var N_isothermal $N_isothermal -var file4_name $file4_name 
            
            # Pequeña pausa para evitar colisiones en fechas
            sleep 1

        done
    done
done

echo "Todas las simulaciones han finalizado."


