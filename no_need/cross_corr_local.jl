using FileIO, JLD2, Plots, StatsBase, Statistics, Random, Measures

ENV["GKSwstype"] = "100"

# Sliding window average function
function sliding_window_average_matrix(data, window_size, step_size)
    n_neurons, total_time = size(data)
    output_time = ceil(Int, (total_time - window_size) / step_size) + 1
    averaged_data = Array{Float64}(undef, n_neurons, output_time)
    
    for i in 1:n_neurons
        output_col = 1
        for t in 1:step_size:(total_time - window_size + 1)
            start_idx = t
            end_idx = t + window_size - 1
            averaged_data[i, output_col] = mean(data[i, start_idx:end_idx])
            output_col += 1
        end
    end
    return averaged_data
end

# Load the data
function load_data(directory)
    @load joinpath(directory, "E_input.jld2") E_input
    @load joinpath(directory, "I_input.jld2") I_input

    E_input = sliding_window_average_matrix(E_input, 350, 15)
    I_input = sliding_window_average_matrix(I_input, 350, 15)

    return E_input, I_input
end

# Compute cross-correlation using the StatsBase crosscor function
function calculate_cross_correlation(input1, input2, max_lag)
    crosscor(input1, input2, -max_lag:max_lag)
end

# Custom function to compute cross-correlation across random pairs
function compute_cross_correlation(E_input::Matrix, I_input::Matrix, tau_range::UnitRange{Int}=-30:30, num_pairs::Int=50)
    Ncells, Nsteps = size(E_input)
    num_pairs = min(round(Int, Ncells / 1.5), num_pairs)  # Control the number of pairs

    avg_correlations = Dict{Int, Float64}()

    for tau in tau_range
        if abs(tau) >= Nsteps
            avg_correlations[tau] = 0
        else
            total_correlation = 0.0
            for i = 1:num_pairs
                e_index, i_index = randperm(Ncells)[1:2]  # Distinct indices

                if tau > 0
                    e_data = E_input[e_index, 1:end-tau]
                    i_data = I_input[i_index, tau+1:end]
                else
                    e_data = E_input[e_index, -tau+1:end]
                    i_data = I_input[i_index, 1:end+tau]
                end

                total_correlation += cor(e_data, i_data)
            end
            avg_correlations[tau] = total_correlation / num_pairs
        end
    end

    return avg_correlations
end

function plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I, cross_corr_I_E, output_path)
    # Check if the dictionaries are empty
    if isempty(cross_corr_E_E)
        error("The correlation dictionaries are empty.")
    end
    
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]
    
    
    p = plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright, size=(1200, 800), margin=25mm)
    plot!(p, sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(p, sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(p, sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    

    xlabel!(p, "Tau")
    ylabel!(p, "Correlation")
    title!(p, "Cross-correlation")
    savefig(p, output_path)
end

# Main execution block
function main()
    println("Start loading")
    dir_name_event = "/Users/dracoxu/Research/data/5.4-real"
    dir_name_non_event = "/Users/dracoxu/Research/data/0.0-real"
    
    # Load event data
    E_input_event, I_input_event = load_data(dir_name_event)
    println("Loaded event data")
    time_indices = Int(round(1000/15)):Int(round(2200/15))
    E_input_non_event = E_input_non_event[:, time_indices]
    I_input_non_event = I_input_non_event[:, time_indices]

    # Compute cross-correlations for event
    cross_corr_E_E = compute_cross_correlation(E_input_event, E_input_event)

    cross_corr_I_I = compute_cross_correlation(I_input_event, I_input_event)

    cross_corr_E_I = compute_cross_correlation(E_input_event, I_input_event)

    cross_corr_I_E = compute_cross_correlation(I_input_event, E_input_event)

    plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I, cross_corr_I_E, "cross_event.png")

    E_input_non_event, I_input_non_event = load_data(dir_name_non_event)
    println("Loaded event data")

    # Compute cross-correlations for event
    cross_corr_E_E = compute_cross_correlation(E_input_non_event, E_input_non_event)

    cross_corr_I_I = compute_cross_correlation(I_input_non_event, I_input_non_event)

    cross_corr_E_I = compute_cross_correlation(E_input_non_event, I_input_non_event)

    cross_corr_I_E = compute_cross_correlation(I_input_non_event, E_input_non_event)

    plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I, cross_corr_I_E, "cross_non_event.png")
    
    
end

main()
