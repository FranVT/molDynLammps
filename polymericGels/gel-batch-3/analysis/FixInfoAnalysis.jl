#=
    Script to analyze the infomration from the fix files
=#

using DataFrames, CSV
using Statistics

#=
    Functions
=#
"""
    getDump(path::String)

Get the data from a single dump file that stores one timeste information
"""
function get_dump(path::String)
    data = split.(readlines(path), " ")[9:end]
    HEADERS = data[1][3:end]
    INFO = parse.(Float64, reduce(hcat, data[2:end]))'

    return DataFrame(INFO, HEADERS)
end



"""
    extract_fix_scalar(path_file::String)
Function that extracts the information of fix files that stores global scalar values
"""
function extract_fix_scalar(path_file::String)
    aux = split.(readlines(path_file), " ")
    header = aux[2][2:end]
   
    aux_info = map(s -> parse.(Float64, s), aux[3:end])
    info = Matrix(reduce(hcat, aux_info)')

    return (header,info)
end

"""
    promedios_acumulados(arrays)

Compute the mean of an arrays with different amount of elements
"""
function promedios_acumulados(arrays)
    
    # Get the maximum aomun of rows
    max_len = maximum(length, arrays)
    
    # To store the avergae
    promedios = Float64[]

    # Go into each row
    for i in 1:max_len
        # Start the addition
        suma = 0.0

        # Count the amount of elements
        count = 0

        # Explore each array
        for arr in arrays

            # Check if the length is correct
            if length(arr) >= i
                # Add to the mean
                suma += arr[i]
                count += 1
            end
        end
        
        # Store the mean
        push!(promedios, count > 0 ? suma / count : NaN)
    end
    return promedios
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

    # headers
    headers=first(unique(first.(info_fixf)));

    # info
    info = last.(info_fixf);

    # Amount of simulations
    n_sims = nrow(df_set); 

    # Number of columns
    n_cols = length(headers);

    # Reshape the stuff
    info_avg = Vector{Vector{Float64}}();

    # Compute the avergae per each column
    for it_col in 1:n_cols
        aux = [];

        for it_sim in 1:length(info) 
            append!(aux,[info[it_sim][:,it_col]])
        end
            
        # Compute the mean
        mean_mod = promedios_acumulados(aux);

        # Append the mean to the final array
        append!(info_avg,[mean_mod])
    end

    # Reduce the dimensions of the array
    info_avg = reduce(hcat,info_avg);

    # Create a dataframe of the average
    df_set_info=DataFrame(info_avg,headers);

    # Add the amount of experiments
    #df_set_info[!,:N_exp].=[nrow(df_set)];

    return df_set_info #data_info_fix #df_set_info
end

"""
    save_mean_fix_analysis(df_set::AbstractDataFrame, FILE_FIX::String, categories_system::Vector{Symbol}, categories_experiment::Vector{Symbol})

Store in analyzed_data the mean of a set of experiments
"""
function save_mean_fix_analysis(meta_data_simulations::AbstractDataFrame, DIR_SAVE::String, categories_id::Vector{Symbol})
        # Compute the assemble average of the fix info
        df_set_info = compute_mean_set_fixf(meta_data_simulations);

        # Create a file name from the categories_total
        ids_set_info=[meta_data_simulations[1, col] for col in categories_id];

        # Add the values of the categories to the dataframe 
        for (col, val) in zip(categories_id, ids_set_info)
            df_set_info[!, col] .= val; 
        end

        # Create a file name from the ids 
        file_name=string("fix_mean_",join(string.(ids_set_info)),".csv");

        # Save the data frame
        CSV.write(joinpath(DIR_SAVE,file_name),df_set_info);
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
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:N_heat,:N_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# Group by experiments
meta_data_experiments=groupby(df_dat,categories_experiment);

# Save the infomration of average assembles
for meta_data_experiment in meta_data_experiments

    # Group by system
    meta_data_systems = groupby(meta_data_experiment,categories_system);

    # Go through all systems and experiments 
    for meta_data_simulations in meta_data_systems
        # Compute the assemble average and store the data
        save_mean_fix_analysis(meta_data_simulations,DIR_SAVE,categories_id)
    end # end systems
end # for experiments

