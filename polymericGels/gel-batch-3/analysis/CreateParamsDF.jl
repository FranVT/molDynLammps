module CreateParamsDF

#=
    Script to create dat files from the information stored 
=#

using DataFrames, CSV

"""
    parse_value(v::String)

Function that converts a string into float, int or string automatically

# Argument
- `v::String`: String to be transformed.

# Examples
julia> parse_value("3.14")
3.14

julia> parse_value.(["3.14","2"])
2-element Vector{Real}:
 3.14
 2

"""
function parse_value(v::String)
    if (i = tryparse(Int, v)) !== nothing
        return i
    elseif (f = tryparse(Float64, v)) !== nothing
        return f
    else
        return v
    end
end

"""
    create_df_params(path_file::String)

Read the parameter file `params.log` from the path `path_file` and return a DataFrame with a row containing all the parameters as columns. Additionally, add the columns `:id` (name of the containing folder) and `:dir` (parent directory).

# Arguments
- `path_file::String`: absolute path to the file `params.log`.

# Example
julia> df = create_df_params("/ruta/experimento/params.log")
1×44 DataFrame

"""
function create_df_params(path_file::String)
    # Read the file
    info_params=readlines(path_file);

    # Split the information into keys and values
    info_split=split.(info_params," = ");

    # Get the keys and values as a Vector{String}
    info_keys=String.(first.(info_split));
    info_values=String.(last.(info_split));

    # Transform the values into usefull data
    info_values_parse=parse_value.(info_values);

    # Create a dataframe with the information
    df_params=DataFrame(Dict(info_keys .=> info_values_parse));

    # Split the path top get the id
    id_system=String(split(path_file,"/")[end-1]);

    # Add the id to the dataframe
    df_params[!,:id].=[id_system];

    # Split the path to create the dir
    dir_experiment=joinpath(split(path_file,"/")[1:end-1]);

    # Add the dir
    df_params[!,:dir].=[dir_experiment];

    return df_params

end

"""
    collect_params_experiments(root_dir::String,log_filename::String="parapms.log") -> DataFrame

Recursively traverse the root_dir directory looking for all the log_filename files in the subfolders, and combine their DataFrames into one. log_filename in the subfolders, and combines their DataFrames into one.

# Arguments
- root_dir::String: root directory that contains the folders of each experiment.
- log_filename::String: name of the parameters file (default "params.log").
"""
function collect_experiments_metadata(root_dir::String, log_filename::String="params.log")

    # Get subdirectories of root_dir
    subdirs=filter(isdir,readdir(root_dir, join = true));

    # Create paths to the params.log files
    file_paths=joinpath.(subdirs, log_filename);

    # Filter to get only the existing files
    existing_paths=filter(isfile, file_paths);

    # Create the data frame for each file 
    dfs = [create_df_params(p) for p in existing_paths]
    return reduce(vcat, dfs)
end

# Export the functions that can be used in the public API
export create_df_params, collect_experiments_metadata

end # Module
