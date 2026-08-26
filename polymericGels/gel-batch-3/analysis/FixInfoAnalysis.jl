#=
    Script to analyze the infomration from the fix files
=#

using DataFrames, CSV
using Statistics, LsqFit

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
    df_set_info[!,:N_exp].=[nrow(df_set)];

    return df_set_info #data_info_fix #df_set_info
end

"""
    create_energy_fit(data_system::AbstractDataFrame)

Function that returns the parameters of the fit eng(t) = p[1]/.t^(0.5)+p[2]
"""
function create_energy_fit(data_system::AbstractDataFrame)
    # Define the model to fit the energy
    model(s,p) = (p[1])./s.^(0.5) .+ p[2];

        # Make sure the time evolution
        sort!(data_system,[:TimeStep]);

        # Extract the time domain and potential energy
        time_domain = (data_system.tstep.*data_system.TimeStep)|>collect;
        pot_eng = data_system.c_ep|>collect;

        # Set intial values for the fit
        p_initial = [1.0, sum(pot_eng)/length(pot_eng)];

        # Fit the data
        fit = curve_fit(model, time_domain, pot_eng, p_initial);

        # Get the parameters
        p_final = fit.param|>collect;

        # Get the standard error
        se = standard_errors(fit)|>collect;

        # Get the margin of error
        margin_of_error = margin_error(fit)|>collect;

        # Store the parameters
        data_system[!, :fit_ep_parameters] = [p_final for _ in 1:nrow(data_system)]
        data_system[!, :fit_ep_margin_error] = [margin_of_error for _ in 1:nrow(data_system)]
        data_system[!, :fit_ep_standard_error] = [se for _ in 1:nrow(data_system)]


    return data_system
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

        # Create the fit of the energy
        df_set_info = create_energy_fit(df_set_info)

        # Create a file name from the ids 
        file_name=string("fix_mean_",join(string.(ids_set_info)),".csv");

        # Save the data frame
        CSV.write(joinpath(DIR_SAVE,file_name),df_set_info);
end

"""
    store_avg()
Function that store the avg of the simulations
"""
function store_avg(meta_data::AbstractDataFrame, categories_experiment::Vector{Symbol}, categories_system::Vector{Symbol}, categories_id::Vector{Symbol}, DIR_SAVE::String)

# Group by experiments
meta_data_experiments = groupby(meta_data,categories_experiment);

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


end

# Función para parsear un string como "[1.0, 2.0]" a Vector{Float64}
function parse_vector(str::String)
    # Quitar corchetes y split por comas
    cleaned = replace(str, r"\[|\]" => "")  # elimina [ y ]
    if isempty(cleaned)
        return Float64[]
    else
        return parse.(Float64, split(cleaned, ","))
    end
end

"""
    extract_fix_avg(DIR_DATA::String)
Get the averages
"""
function extract_fix_avg(DIR_DATA::String)

    # Read the directory 
    files=readdir(DIR_DATA);

    # Get only those of the structure factor
    files=filter(s -> occursin("fix_mean_", s), files);

    # Read the files
    df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

    # Create one dataframe
    df_files = reduce(vcat,df_files)

    # PArse the vectors of strings into vectors of floats
    df_files[!, :fit_ep_parameters] = parse_vector.(df_files[!, :fit_ep_parameters])
    df_files[!, :fit_ep_margin_error] = parse_vector.(df_files[!, :fit_ep_margin_error])
    df_files[!, :fit_ep_standard_error] = parse_vector.(df_files[!, :fit_ep_standard_error])
    
    return  df_files   

end
#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_SAVE = joinpath(DIR_MAIN,"analyzed_data");
DIR_DATA = joinpath(DIR_MAIN,"analyzed_data");
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

# Store the avg.
store_avg(df_dat,categories_experiment,categories_system,categories_id,DIR_SAVE)

# Open the data
df_group=extract_fix_avg(DIR_DATA)

# Group by experiments
data_per_experiment = groupby(df_group,categories_experiment);

# Select one experiment
data_experiment = data_per_experiment[1];

#=
    # Group by systems
    data_per_system = groupby(data_experiment,categories_system);

    # To store the parameter of the function
    systems_parameters = Vector{Vector{Float64}}();

    # To store the margin error
    systems_margin_error = Vector{Vector{Float64}}();

    # To store standard error
    systems_standard_error = Vector{Vector{Float64}}();

    # Define the model to fit the energy
    model(s,p) = (p[1])./s.^(0.5) .+ p[2];

    # Select one system
#    data_system = data_per_system[1];

    for (it,data_system) in enumerate(data_per_system)
        # Make sure the time evolution
        sort!(data_system,[:TimeStep]);

        # Extract the time domain and potential energy
        time_domain = (data_system.tstep.*data_system.TimeStep)|>collect;
        pot_eng = data_system.c_ep|>collect;

        # Set intial values for the fit
        p_initial = [1.0, sum(pot_eng)/length(pot_eng)];

        # Fit the data
        fit = curve_fit(model, time_domain, pot_eng, p_initial);

        # Get the parameters
        p_final = fit.param;

        # Get the standard error
        se = standard_errors(fit);

        # Get the margin of error
        margin_of_error = margin_error(fit);

        # Store the parameters
        data_per_system[it][!,:fit_parameters].=[p_final]
        data_per_system[it][!,:fit_margin_error].=[margin_of_error]
        data_per_system[it][!,:fit_standard_error].=[se]

    end
=#

