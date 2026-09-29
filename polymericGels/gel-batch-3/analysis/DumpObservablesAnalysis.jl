#=
    Script to analyse the observables computed by the dump info.
=#

using DataFrames, CSV
using GLMakie
using Statistics, LsqFit

#=
    Functions 
=#

"""
    extract_Sq_observable(DIR_DATA::String)
Get the averages
"""
function extract_Sq_observable(DIR_DATA::String)

    # Read the directory 
    files=readdir(DIR_DATA);

    # Get only those of the structure factor
    files=filter(s -> occursin("structure_factor_", s), files);

    # Read the files
    df_files=[CSV.read(joinpath(DIR_DATA,file), DataFrame) for file in files];

    # Create one dataframe
    df_files = reduce(vcat,df_files)
   
    return  df_files   

end

"""
    derivada_no_uniforme(x::AbstractVector, y::AbstractVector)

Por DeepSeek
"""
function derivada_no_uniforme(x::AbstractVector, y::AbstractVector)
    n = length(x)
    @assert n == length(y) "x e y deben tener la misma longitud"
    @assert n ≥ 3 "Se necesitan al menos 3 puntos"

    dy = zeros(eltype(y), n)

    # --- Primer punto: hacia adelante (3 puntos) ---
    h1 = x[2] - x[1]
    h2 = x[3] - x[2]
    dy[1] = y[1] * (-(2h1 + h2) / (h1 * (h1 + h2))) +
            y[2] * ((h1 + h2) / (h1 * h2)) +
            y[3] * (-h1 / (h2 * (h1 + h2)))

    # --- Puntos interiores: centrada (3 puntos) ---
    @inbounds for i in 2:n-1
        ha = x[i]   - x[i-1]   # espaciado hacia atrás
        hb = x[i+1] - x[i]     # espaciado hacia adelante
        dy[i] = y[i-1] * (-hb / (ha * (ha + hb))) +
                y[i]   * ((hb - ha) / (ha * hb)) +
                y[i+1] * ( ha / (hb * (ha + hb)))
    end

    # --- Último punto: hacia atrás (3 puntos) ---
    h1 = x[n-1] - x[n-2]
    h2 = x[n]   - x[n-1]
    dy[n] = y[n-2] * ( h2 / (h1 * (h1 + h2))) +
            y[n-1] * (-(h1 + h2) / (h1 * h2)) +
            y[n]   * ((h1 + 2h2) / (h2 * (h1 + h2)))

    return dy
end

"""
    find_spatial_regions_Sq(Sq_range::Vector{Float64})

Function that returns the index of the max and minums of the second derivative of Sq at a log scale
"""
function find_spatial_regions_Sq(Sq_range::Vector{Float64}, q_domain::Vector{Float64}; w=2)
                
    # Smooth the range
    Sq_range_smooth = [mean(Sq_range[max(1,i-w÷2):min(end,i+w÷2)]) for i in eachindex(Sq_range)];

    # Transform the domain and range into the log scale
    q_log_domain = log.(10,q_domain);
    Sq_log_range = log.(10,Sq_range_smooth);
    
    # the idea is to compute the derivative of the range.
    # When the derivative surpaes a trashhold, the interval is defined.
    derivative_Sq = derivada_no_uniforme(q_domain, Sq_range_smooth);
    
    # Smooth derivative
    derivative_Sq_smooth = [mean(derivative_Sq[max(1,i-w÷2):min(end,i+w÷2)]) for i in eachindex(derivative_Sq)];

    # Compute the second derivative
    second_derivative_Sq_smooth = derivada_no_uniforme(q_domain, derivative_Sq_smooth);
                
    # Change of variable (Lazziness)
    d2 = deepcopy(second_derivative_Sq_smooth);

    # Threshold
    umbral = 1; 

    # picos locales en |d2|
    picos_ind = Int[]
    for i in 2:length(d2)-1
        if abs(d2[i]) > abs(d2[i-1]) && abs(d2[i]) > abs(d2[i+1]) && abs(d2[i]) > umbral
            push!(picos_ind, i)  # índice en q
        end
    end

    return picos_ind

end

"""
    compute_the_fit(ind_peaks::Vector{Int64}, q_domain::Vector{Float64}, Sq_range::Vector{Float64})

Function that gives the values of m and b for the Sq at log scale in the fractal region
"""
function compute_the_fit(ind_peaks::Vector{Int64}, q_domain::Vector{Float64}, Sq_range::Vector{Float64})
   
    # Transform the domain and range into the log scale
    q_log_domain = log.(10,q_domain);
    Sq_log_range = log.(10,Sq_range_smooth);

    # Get the cut near the particle size
    q_fractal = q_domain[first(ind_peaks)];
    q_particle = 2*pi*0.7; # 70% of the particle size

    # Get the index at the middle
    ind_network = q_fractal .< q_domain .< q_particle

    # Select the region for the linear fit
    q_network = q_log_domain[ind_network]
    Sq_network = Sq_log_range[ind_network]

    # Create the fit
    model(t,p) = p[1].*t.+p[2]

    # Set intial values for the fit
    p_initial = [-1.0, 0.0];

    p_lower = [-Inf, -Inf];
    p_upper = [Inf, Inf];

    # Fit the data
    fit = curve_fit(model, q_network, Sq_network, p_initial; lower=p_lower, upper=p_upper);

    # Get the parameters
    params_final = fit.param|>collect;

    return params_final

end


function eval_model(t,p)
    return p[1].*t.+p[2]
end





#=
    Start the script
=#

# Paths and directories
DIR_MAIN = pwd();
DIR_DATA = joinpath(DIR_MAIN,"analyzed_data");
FILE_DAT = "dat.csv";
DIR_SAVE = joinpath(DIR_MAIN,"figures");

# Combine the dataframes
df_group=extract_Sq_observable(DIR_DATA);

# Add the time instant
df_group[!,:time] = df_group.timeStep.*df_group.tstep;


# NEED TO BE THE SAME AS THE FixInfoAnalysis.jl
# Select the categories that define a system
categories_system=[:phi,:chi_4,:temp,:damp,:tstep];

# Create categories to select different experiments (Just in case)
categories_experiment=[:time_heat,:time_isothermal];

# For id
categories_id = [categories_system; categories_experiment];

# Group by experiments
data_per_experiment = groupby(df_group,categories_experiment);

# Select one experiment
data_experiment = data_per_experiment[1];

    # Group by system
    data_per_system = groupby(data_experiment,categories_system);

    # Select one system
    data_system = data_per_system[1];

        # Group by time instant
        data_per_time = groupby(data_system,:time);

        # Select one time instant
        data_time = data_per_time[end];

            # Group by simulation
            data_per_simulation = groupby(data_time,:Nsim);

            # Select one simulation
            data_simulation = data_per_simulation[1];

                # Get the q domain
                q_domain = data_simulation.q_mean[1:end-1];

                # Get the Sq_mean
                Sq_range = data_simulation.Sq_mean[1:end-1];

                # Get the index for the spatial domains
                ind_peaks = find_spatial_regions_Sq(Sq_range::Vector{Float64}, q_domain::Vector{Float64};

                # Perform the fit 
                params_fit = compute_the_fit(ind_peaks,q_domain,Sq_range)


    cortes = [1; ind_peaks; length(q_domain)]



    
f = Figure()

ax1 = Axis(f[1, 1], yticklabelcolor = :blue)
ax2 = Axis(f[1, 1], yticklabelcolor = :red, yaxisposition = :right)
hidespines!(ax2)
hidexdecorations!(ax2)

lines!(ax1,q_log_domain,eval_model(q_log_domain,p_final),linestyle=:dash, color =:black,
       linewidth=2.5)
scatterlines!(ax1, q_log_domain, Sq_log_range, color = :blue)
scatterlines!(ax1, q_log_domain, log.(10,Sq_range), color = :grey)
scatterlines!(ax2, q_log_domain, derivative_Sq_smooth, color = (:red,0.5))
scatterlines!(ax2, q_log_domain, second_derivative_Sq_smooth, color = :red)

vlines!(ax1,q_log_domain[cortes])
vlines!(ax1,log(10,2*pi*0.7))

display(f)
