"""
    File with all the functions needed for multiple.jl script
"""

function getDat(path)
"""
    Creates a dataframe from the dat file of the experiment.
"""
    file_path=joinpath(path,"dataAssembly.dat");
    aux=split.(readlines(file_path),",");

    df_aux=DataFrame();
    for (col, val) in zip(aux[1], aux[2])
        # Convertir a Float64 si es posible, sino mantener como String
        parsed_val = tryparse(Float64, val)
        if parsed_val !== nothing
            df_aux[!, col] = [parsed_val]
        else
            df_aux[!, col] = [val]
        end
    end

    df_aux[!, :dir] = [path];

    return df_aux
end


