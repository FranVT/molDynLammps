#=
    Script to analyse the observables computed by the dump info.
=#

using DataFrames, CSV
using GLMakie, LaTeXStrings
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
function find_spatial_regions_Sq(Sq_range::Vector{Float64}, q_domain::Vector{Float64})
               
w=2
    # Smooth the range
    Sq_range_smooth = [mean(Sq_range[max(1,i-w÷2):min(end,i+w÷2)]) for i in eachindex(Sq_range)];

    # Transform the domain and range into the log scale
    #q_log_domain = log.(10,q_domain);
    #Sq_log_range = log.(10,Sq_range_smooth);
    
    # the idea is to compute the derivative of the range.
    # When the derivative surpaes a trashhold, the interval is defined.
    Sq_prime = derivada_no_uniforme(q_domain, Sq_range_smooth);
    
    # Smooth derivative
    Sq_prime_smooth = [mean(Sq_prime[max(1,i-w÷2):min(end,i+w÷2)]) for i in eachindex(Sq_prime)];

    # Compute the second derivative
    Sq_dprime_smooth = derivada_no_uniforme(q_domain, Sq_prime_smooth);
                
    # Change of variable (Lazziness)
    d2 = deepcopy(Sq_dprime_smooth);

    # Threshold
    umbral = 1; 

    # picos locales en |d2|
    picos_ind = Int[]
    for i in 2:length(d2)-1
        if abs(d2[i]) > abs(d2[i-1]) && abs(d2[i]) > abs(d2[i+1]) && abs(d2[i]) > umbral
            push!(picos_ind, i)  # índice en q
        end
    end

    return picos_ind, Sq_prime_smooth, Sq_dprime_smooth

end

"""
    compute_the_fit(ind_peaks::Vector{Int64}, q_domain::Vector{Float64}, Sq_range::Vector{Float64})

Function that gives the values of m and b for the Sq at log scale in the fractal region
"""
function compute_the_fit(ind_peaks::Vector{Int64}, q_domain::Vector{Float64}, Sq_range::Vector{Float64})
   
    # Transform the domain and range into the log scale
    q_log_domain = log.(10,q_domain);
    Sq_log_range = log.(10,Sq_range);

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

function eval_model_log(t,p)
    return p[1].*t.+p[2]
end

function eval_model_linear(t,p)
    return (10).^(p[2]).*t.^(p[1])
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

            # Allocate for the mean
            Sq_mean = zeros(length(q_domain));

            # Compute the mean 
            for aux in data_per_simulation
                Sq_mean[:] += aux.Sq_mean[1:end-1]
            end
            Sq_mean = Sq_mean./length(data_per_simulation);

                # Get the index for the spatial domains
                ind_peaks, Sq_prime_smooth, Sq_dprime_smooth = find_spatial_regions_Sq(Sq_mean,q_domain);

                # Perform the fit 
                params_fit = compute_the_fit(ind_peaks,q_domain,Sq_mean)





    
fig = Figure()

    # Prepare the ticks
    n_ticks = 10;
    q_aux_ticks = q_domain;
    l_domain = 2*pi./q_aux_ticks;
    ind_range = floor.(Int64,(10).^(range(log(10,1),log(10,length(q_aux_ticks)),length=n_ticks)));
    q_positions = round.(q_aux_ticks[ind_range],digits=2);
    q_ticks = latexstring.(q_positions);
    l_ticks = latexstring.(round.(l_domain[ind_range],digits=2));

    # --- Define tick positions (in q-space) and their top labels (λ = 2π/q) ---
    ax_bottom = Axis(fig[1:4, 1:5],
                         xlabel = L"|\vec{q}|",
                         ylabel = L"\mathrm{Intensity}",
                         xticks = (q_positions, q_ticks),
                             xscale = log10,
                             yscale = log10,
                             xticklabelrotation = pi/4
                            )

    # --- Top axis: wavelength λ ---
    ax_top = Axis(fig[1:4, 1:5],
                          xaxisposition = :top,
                          yaxisposition = :right,

        # Place ticks at the same data coordinates (q values),
        # but display the corresponding λ labels.
                          xticks = (q_positions, l_ticks),
                          xlabel = L"\mathrm{Wavelength}",

        # Spines: show only the top spine
                          topspinevisible = true,
                          bottomspinevisible = false,
                          leftspinevisible = false,
                          rightspinevisible = false,

                          xgridvisible = false,
                          #ygridvisible = false,

        # Hide all y‑axis decorations on the top axis
                          yticks = ([], []),
                          ylabelvisible = false,
                          ygridvisible = false,
                          yticklabelsvisible = false,
                             xscale = log10,
                             yscale = log10,
                             xticklabelrotation = pi/4
                         )

    # Synchronise limits and zoom/pan behaviour
    linkaxes!(ax_bottom, ax_top)


    # Transform the domain and range into the log scale
    q_log_domain = log.(10,q_domain);
    Sq_log_range = log.(10,Sq_range);

    cortes = [ind_peaks[1:end-1]]

    q_cut = [q_domain[ind_peaks[1:end-1]]; 2*pi*0.7];

 

#ax1 = Axis(fig[1, 1], yticklabelcolor = :blue)
ax2 = Axis(fig[1:4, 1:5], yticklabelcolor = :red, yaxisposition = :right)
hidespines!(ax2)
hidexdecorations!(ax2)

q_reg_low = q_cut[1:end-1];
q_reg_high = q_cut[2:end];

vspan!(ax_bottom,[q_reg_low[1]],[q_reg_high[1]], color = (:dodgerblue,0.5))
vspan!(ax_bottom,[q_reg_low[2]],[q_reg_high[2]], color = (:orange,0.5))


scatterlines!(ax2, q_log_domain, Sq_prime_smooth, color = (:red,0.5))
scatterlines!(ax2, q_log_domain, Sq_dprime_smooth, color = :red)


    lines!(ax_bottom,q_domain,eval_model_linear(q_domain,params_fit),linestyle=:dash, color =:black,
       linewidth=2.5)

    # Add the line to the plot
    scatterlines!(ax_bottom, q_domain, Sq_mean,
        #linestyle = label_systems[it_system],
        linewidth = 4
       )

#scatterlines!(ax1, q_log_domain, Sq_log_range, color = :blue)
#scatterlines!(ax1, q_log_domain, log.(10,Sq_range), color = :grey)
#vlines!(ax_bottom,q_domain[ind_peaks])
#vlines!(ax_bottom,2*pi*0.7)




display(fig)
