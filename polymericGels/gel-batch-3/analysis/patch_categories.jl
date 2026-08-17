#=
    Script to fix categories missing and stuff
=#

#=
    Functions
=#

"""
    get_group_files(pattern::string,files::Vector{String})

Function that return and array of arrays of names of files of the same experiment at different time steps given an array of files names
"""
function get_group_files(patron::Regex,files::Vector{String})

    # Read by experiments
    grupos = Dict{Tuple{String, Int}, Vector{String}}()
    for f in files 
        m = match(patron, f)
        if m !== nothing
            sistema = m.captures[1]   # cadena con los parámetros
            step = parse(Int, m.captures[2])
            clave = (sistema, step)
            # Agregar al grupo correspondiente
            if haskey(grupos, clave)
                push!(grupos[clave], f)
            else
                grupos[clave] = [f]
            end
        else
            @warn "Nombre no coincide con el patrón: $f"
        end
    end

    return grupos
end



using DataFrames, CSV

#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_DATA = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
DIR_SAVE = joinpath(DIR_MAIN,"figures");

# Read the dat file
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Read the directory 
files=readdir(DIR_DATA);

# Get only those of the structure factor
files=filter(s -> occursin("structure_factor_", s), files);

# Read the files
df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

# The idea is that all df have the exact sames keys

symbols_per_simulation = names.(df_files);

# Get the unique keys
symbols = unique(symbols_per_simulation);

if length(symbols) == 1
    println("Same amount of keys")
end

# Get all the keys
all_symbols = union(symbols...);

# Amount of keys. For first test of similarity
n_symbols = length(all_symbols);

# Find the simulation that does not have the same keys
diff_symbols_per_simulation = map(k_p_s->mapreduce(s->k_p_s.== s,+,all_symbols),symbols_per_simulation);

# Known wich are equal and what are different
mask_diff = sum.(diff_symbols_per_simulation) .!= n_symbols;

# Get the files that have different keys
files_diff = files[mask_diff];

# Prepare for reading the information
pattern = r"structure_factor_(.+)_step_(\d+).csv"

# Get the group the simulations by expriment and time domain
groups_diff = get_group_files(pattern,files_diff)

# Get id and time step
keys_diff = collect(keys(groups_diff));

# Select one simulation
key_simulation = keys_diff[1];

# get the id and the time step
(id_simulation,time_step_simulation) = key_simulation;


