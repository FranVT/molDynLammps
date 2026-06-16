#!/bin/bash

# Función que guarda todas las variables exportadas que comienzan con un prefijo,
# o bien todas las variables de entorno que definiste en parameters.sh.
# Aquí guardamos todas las variables que empiezan con mayúscula (opcional),
# o podemos pasar una lista explícita.

save_params() {
    local dir_save="$1"
    local log_file="$dir_save/params.log"

    # Opción A: Guardar todas las variables que contienen caracteres
    # que normalmente usas (mayúsculas, minúsculas, guiones bajos, números).
    # Pero cuidado: esto puede guardar muchas variables del sistema.
    # Mejor limitarnos a un conjunto de variables conocidas.
    
    # ------------------------------------------------------------------
    # Opción B (Recomendada): Definir explícitamente qué variables guardar.
    # Mantienes el control y evitas guardar variables del sistema.
    # ------------------------------------------------------------------
    local vars_to_save=(
        "phi_o" "phi_f" "phi_delta"
        "chi_4o" "chi_4f" "chi_4delta"
        "N_PP" "N_exp"
        "temp" "damp"
        "m_CP" "m_PT"
        "r_cWCA" "r_cPP" "r_CP" "r_bondP"
        "theta_PA" "theta_PB"
        "E_CP" "E_bondP"
        "N_heat" "N_isothermal" "N_save" "N_dump"
        "N_3body" "r_noInter" "seed1" "seed2" "seed3" "tstep"
        "file1_name" "file2_name" "file3_name" "file4_name"
        "cs" "Vol_MO" "Vol_CL"
        # Ahora las variables que se modifican dinámicamente
        "phi" "chi_4" "N_pCL" "N_pMO" "L"
    )

    # Escribir el encabezado
    echo "# Parámetros usados en esta simulación" > "$log_file"
    echo "# Fecha: $(date)" >> "$log_file"
    echo "# Directorio: $dir_save" >> "$log_file"
    echo "# ------------------------------------------------------" >> "$log_file"

    # Recorrer la lista y escribir cada variable con su valor actual
    for var in "${vars_to_save[@]}"; do
        # Verificar que la variable está definida (no vacía)
        if [ -n "${!var}" ]; then
            printf "%-15s = %s\n" "$var" "${!var}" >> "$log_file"
        else
            printf "%-15s = %s\n" "$var" "(no definida)" >> "$log_file"
        fi
    done

    # Añadir un comentario adicional si quieres
    echo "# ------------------------------------------------------" >> "$log_file"
    echo "# Todos los parámetros tomados de parameters.sh y actualizados" >> "$log_file"
}


