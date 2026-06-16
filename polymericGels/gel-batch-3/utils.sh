#!/bin/bash
# utils.sh – Funciones auxiliares para las simulaciones

# ------------------------------------------------------------------
# Función para escribir todos los parámetros en un archivo de log
# Uso: write_params <archivo_salida> [par_extra1=valor1 ...]
# ------------------------------------------------------------------
write_params() {
    local output_file="$1"
    shift

    # Lista de todas las variables a registrar
    local all_vars=(
        phi_o phi_f phi_delta chi_4o chi_4f chi_4delta N_PP N_exp
        temp damp N_pCL N_pMO L
        m_CP m_PT r_cWCA r_cPP r_CP r_bondP theta_PA theta_PB E_CP E_bondP
        N_heat N_isothermal N_save N_dump
        N_3body r_noInter seed1 seed2 seed3 tstep
        file1_name file2_name file3_name file4_name
        cs Vol_MO Vol_CL
    )

    declare -A params

    # Valores por defecto
    for var in "${all_vars[@]}"; do
        if [ -n "${!var}" ]; then
            params["$var"]="${!var}"
        fi
    done

    # Sobrescritura con valores específicos de la simulación
    for pair in "$@"; do
        IFS='=' read -r key value <<< "$pair"
        params["$key"]="$value"
    done

    # Escribir archivo sin comentarios, solo líneas "clave = valor"
    {
        for key in $(printf '%s\n' "${!params[@]}" | sort); do
            echo "$key = ${params[$key]}"
        done
    } > "$output_file"
}


# ------------------------------------------------------------------
# 3. Función para generar semillas únicas basadas en fecha y contador
# ------------------------------------------------------------------
generate_seeds() {
    local exp_num=$1
    local seed_base=$(( exp_num + 10#$(date +%S) + 10#$(date +%M) ))
    seed1=$(( seed_base + 1 ))
    seed2=$(( seed_base + 2 ))
    seed3=$(( seed_base + 3 ))
    export seed1 seed2 seed3
}

