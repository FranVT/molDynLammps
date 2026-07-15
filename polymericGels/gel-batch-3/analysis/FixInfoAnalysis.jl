#=
    Script to analyze the infomration from the fix files
=#

using DataFrames, CSV
using Statistics

#=
    Functions
=#
"""
    extract_fix_scalar(path_file::String)
Function that extracts the information of fix files that stores global scalar values
"""
function extract_fix_scalar(path_file::String)
    aux = split.(readlines(path_file), " ")
    header = aux[2][2:end]
   
    aux_info = map(s -> parse.(Float64, s), aux[3:end])
    info = reduce(hcat, aux_info)

    return (header,info)
end

"""
    compute_mean_set_fixf(df_set::DataFrame, FILE_FIX::String="system_assembly.fixf")

Compute the mean of the values stored in the fix files of a set of the same same experiments
 
# Arguments
- `df_set::DataFrame`: DataFrame with all metada data of the set
- `FILE_FIX::String`: File name
"""
function compute_mean_set_fixf(df_set::AbstractDataFrame, FILE_FIX::String="system_assembly.fixf")

    # Get the directories of all experiments of the system 
    dir_set=String.(df_set.dir);

    # Create the paths to the files
    path_fixf=joinpath.(dir_set,FILE_FIX);

    # Extract all the information from the fix files
    info_fixf=extract_fix_scalar.(path_fixf)

    # Get the data to compute the mean 
    data_info_fix = last.(info_fixf);

    # Get the minimum time step 
    min_time_step = minimum(last.(size.(data_info_fix)))

    # Format to create the array to compute the mean 
    data_info_fix = [data[:,1:min_time_step] for data in data_info_fix];

    # Compute the average of the set of experiments
    mean_set=first(mean(data_info_fix,dims=1));

    # Get the headers
    headers_set=String.(first.(info_fixf)[1]);

    # Create a dataframe of the average
    df_set_info=DataFrame(mean_set',headers_set);

    # Add the amount of experiments
    df_set_info[!,:N_exp].=[nrow(df_set)];

    return df_set_info
end

"""
    save_mean_fix_analysis(df_set::AbstractDataFrame, FILE_FIX::String, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol})

Store in analyzed_data the mean of a set of experiments
"""
function save_mean_fix_analysis(df_set::AbstractDataFrame, DIR_SAVE::String, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol})
    # Create a dataframe with the information of the fix files
    df_set_info=compute_mean_set_fixf(df_set);

    # Define a set of categories
    categories_total=[categories_system; categories_experiment];

    # Create a file name from the categories_total
    ids_set_info=[df_set[1, col] for col in categories_total];

    # Add the values of the categories to the dataframe 
    for (col, val) in zip(categories_total, ids_set_info)
        df_set_info[!, col] .= val 
    end

    # Create a file name from the ids 
    file_name=string("fix_mean_",join(string.(ids_set_info)),".csv");

    # Save the data frame
    CSV.write(joinpath(DIR_SAVE,file_name),df_set_info);
  
    return nothing
end
#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
FILE_FIX = "system_assembly.fixf";

# Read the dat file
df_dat=CSV.read(joinpath(DIR_MAIN,FILE_DAT), DataFrame);

# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp];

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# Group by system 
df_systems=groupby(df_dat,categories_system);

# Save the mean of all experiments
for df_system in df_systems
    # Group by experiments
    df_experiments=groupby(df_system,categories_experiment);

    # Save the dataframes per set of experiments
    for df_set in df_experiments
        save_mean_fix_analysis(df_set,DIR_SAVE,categories_system,categories_experiment)
    end
end
