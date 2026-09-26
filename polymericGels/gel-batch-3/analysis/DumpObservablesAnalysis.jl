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

                # Smoothe the range
                w=3
                Sq_range_smooth = [mean(Sq_range[max(1,i-w÷2):min(end,i+w÷2)]) for i in eachindex(Sq_range)];

                # Transform the domain and range into the log scale
                q_log_domain = log.(10,q_domain);
                Sq_log_range = log.(10,Sq_range_smooth);
    
                # the idea is to compute the derivative of the range.
                # When the derivative surpaes a trashhold, the interval is defined.
                derivative_Sq = derivada_no_uniforme(q_domain, Sq_range_smooth);
                #diff(Sq_range)./diff(q_domain);

                # Smoothiung the derivative
                #w = 8;   # window of the mean
                derivative_Sq_smooth = [mean(derivative_Sq[max(1,i-w÷2):min(end,i+w÷2)]) for i in eachindex(derivative_Sq)];

                # Compute the second derivative
                second_derivative_Sq_smooth = derivada_no_uniforme(q_domain, derivative_Sq_smooth);
                #diff(derivative_Sq_smooth)./diff(q_domain[1:end-1]);
                d2 = deepcopy(derivative_Sq_smooth);

# --- 3. PELT con coste Normal y penalización por defecto (log(n)) ---
#    Nota: PELT espera el vector de datos, no (x, y).
#cps, costo = PELT(y_suave, Normal(:?, 1.0))

    umbral = 1; 

                 # picos locales en |d2|
    picos = Int[]
    for i in 2:length(d2)-1
        if abs(d2[i]) > abs(d2[i-1]) && abs(d2[i]) > abs(d2[i+1]) && abs(d2[i]) > umbral
            push!(picos, i)  # índice en q
        end
    end
    cortes = [1; picos; length(q_domain)]

    # Get the cut near the particle size
    q_fractal = q_domain[first(picos)];
    q_particle = 2*pi*0.7; # 70% of the particle size

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
    p_final = fit.param|>collect;

function eval_model(t,p)

    return p[1].*t.+p[2]
end



    
                # Smoothiung the derivative
                #w = div(length(Sq_log_range),10);   # window of the mean
                #second_derivative_Sq_log_smooth = [mean(second_derivative_Sq_log_smooth[max(1,i-w÷2):min(end,i+w÷2)]) for i in eachindex(second_derivative_Sq_log_smooth)];

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

f
