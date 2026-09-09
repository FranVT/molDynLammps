# To read the dataframe to store


using DataFrames, CSV

#=
    Functions
=#

"""
    convert_from_an(s::Any)

Convierte un `String` a un tipo de datos Julia adecuado, detectando automáticamente si
representa un número, booleano, `nothing`, `missing`, arreglo, tupla o diccionario.

La conversión sigue este orden:
1. `nothing`, `missing`, `true`, `false` literales.
2. Números enteros y de punto flotante.
3. Expresiones estructuradas (arrays `[...]`, diccionarios `Dict(...)`, tuplas `(...)`) usando `eval(Meta.parse())`.
4. Si nada funciona, devuelve el `String` original.

!!! warning "Seguridad"
    El uso de `eval` en datos no confiables puede ejecutar código arbitrario.
    Esta función solo evalúa expresiones que comienzan con `[`, `Dict(` o `(`, pero
    sigue siendo recomendable usarla solo con datos generados internamente.

# Ejemplos
```julia
convert_from_string("123")         # 123 (Int64)
convert_from_string("3.14")        # 3.14 (Float64)
convert_from_string("[1,2,3]")     # [1, 2, 3] (Vector{Int64})
convert_from_string("Dict(1=>2)")  # Dict(1=>2) (Dict{Int64,Int64})
convert_from_string("nothing")     # nothing
convert_from_string("texto")       # "texto" (String)

Done by deepseek.
"""
function convert_from_any(s::Any)
  
    #println(typeof(s))

    if typeof(s) == Int64
        return s
    end

    if typeof(s) == Float64
        return s
    end

    s = strip(s)
    # 1. Casos especiales
    if s == "nothing"
        return nothing
    elseif s == "missing"
        return missing
    elseif s == "true"
        return true
    elseif s == "false"
        return false
    end

    # 2. Intentar parsear como número (Int, luego Float)
    try
        return parse(Int, s)
    catch
        try
            return parse(Float64, s)
        catch
        end
    end

    # 3. Intentar evaluar como expresión Julia (arrays, dicts, tuplas, etc.)
    #    Solo si el string parece una expresión estructurada.
    #    Para evitar eval peligroso, limitamos a casos que empiecen con '[', 'Dict(', etc.
    if startswith(s, '[') && endswith(s, ']')
        try
            expr = Meta.parse(s)
            return eval(expr)
        catch
        end
    elseif startswith(s, "Dict(") && endswith(s, ')')
        try
            expr = Meta.parse(s)
            return eval(expr)
        catch
        end
    elseif startswith(s, '(') && endswith(s, ')')
        try
            expr = Meta.parse(s)
            return eval(expr)
        catch
        end
    end


    # 3. Intentar evaluar como expresión Julia
    #    Esto cubre arrays, dicts, tuplas, y también Float64[] (arreglos vacíos tipados)
    try
        expr = Meta.parse(s)
        return eval(expr)
    catch
        @warn "No conversion done"
        return s
    end

    

    # 4. Si nada funcionó, devolver el string original
    #return s
end

"""
convert_dataframe(df::DataFrame)

Transform each column into computational information
"""
function convert_dataframe(df::DataFrame)
# Create the new dataframe to store the data
df_new = DataFrame();

# Suponiendo que df es tu DataFrame
for col in names(df)
    # Aplica la conversión a cada elemento de la columna
    df_new[!, col] = convert_from_any.(df[!, col])
end

return df_new
end

"""
    extract_connectivity_analysis(DIR_DATA::String)
Get the averages
"""
function extract_connectivity_analysis(DIR_DATA::String)

    # Read the directory 
    files=readdir(DIR_DATA);

    # Get only those of the structure factor
    files=filter(s -> occursin("connectivity_analysis_", s), files);

    # Read the files
    df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

    # Create one dataframe
    df_files = reduce(vcat,df_files)
  
    # Convert the information into computational data
    df_data = convert_dataframe(df_files)

    return  df_data

end



#=
    Script
=#

# Paths and directories
DIR_MAIN = dirname(pwd());
DIR_DATA = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
DIR_SAVE = joinpath(DIR_MAIN,"figures");

# Paths and directories
df_group=extract_connectivity_analysis(DIR_DATA);

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:time_heat,:time_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# Group by experiment
df_experiments = groupby(df_group,categories_experiment);

# Select one experiment
df_experiment = df_experiments[1];

    # Group by systems
    df_systems = groupby(df_experiment,categories_system);

    # Select one system
    df_system = df_systems[1];

        # Group by time instant
        df_time = groupby(df_system,:time_instant)

