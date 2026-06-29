#=
    Script to create dat files from the information stored 
=#

using DataFrames, CSV


#=
    Functions
=#
"""
    parse_value(v::String)

Function that converts a string into float, int or string automatically
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
    create_dfParams(path_file::String)
Function that returns a dataframe of the params.log file given the path.
"""
function create_dfParams(path_file::String)
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

#=
    Start of the Script
=#

# Directory where the data is
dir_data="/run/media/franvt/rogelio/DinMol/gel-batch-3/data/";

# This directory
dir_main=pwd();

# Read the directory and find subdirectories
dirs_info=filter(isdir, readdir(dir_data, join = true));

# Name of the file
FileLog="params.log";

# Create the path to the data
path_params=joinpath.(dir_data,dirs_info,FileLog);

# Read the params.log file
df_params=[create_dfParams(s) for s in path_params];

# Create a combine dataframe
df_dat=reduce(vcat,df_params);

# Define a file name for the data file
FileDat="dat.csv";

# File path of dataframe data 
file_path=joinpath(dir_main,FileDat);

# Save the dataframe of the metadata of the experiments 
CSV.write(file_path,df_dat)



# Search for files in each path
#files_info=map(s->filter(isfile, readdir(s, join = true)),path_info);

