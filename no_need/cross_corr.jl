using FileIO, JLD2, Plots, StatsBase, Statistics

ENV["GKSwstype"] = "100"

# Define the function to calculate cross-correlation
function calculate_cross_correlation(input1, input2, max_lag)
    # Assuming input1 and input2 are time series data for two different neurons or neuron groups
    return crosscor(input1, input2, lag = -max_lag:max_lag)
end

function sliding_window_average_matrix(data, window_size, step_size)
    n_neurons, total_time = size(data)
    # Calculate the number of output points based on window size and step size
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

step_size = 15 
window_size = 350

# Load the data
function load_data(directory)

    @load joinpath(directory, "E_input.jld2") E_input
    @load joinpath(directory, "I_input.jld2") I_input

    E_input = sliding_window_average_matrix(E_input, window_size,step_size)
    I_input = sliding_window_average_matrix(I_input, window_size,step_size)

    return E_input, I_input
end

# Here's where we'd call the functions defined above and use the data
dir_name_event = "/gpfs/data/doiron-lab/draco/results_corr/5.4-real"
dir_name_non_event = "/gpfs/data/doiron-lab/draco/results_corr/0.0-real"
println("start loading")
# Load event data
E_input_event, I_input_event = load_data(dir_name_event)
println("loaded event data")
time_indices = Int(round(5800/step_size)):Int(round(7000/step_size))
E_input_event = E_input_event[:, time_indices]
I_input_event = I_input_event[:, time_indices]

# Load non-event data
println("start loading nonevent data")
E_input_non_event, I_input_non_event = load_data(dir_name_non_event)
time_indices = Int(round(1000/step_size)):Int(round(2200/step_size))
E_input_non_event = E_input_non_event[:, time_indices]
I_input_non_event = I_input_non_event[:, time_indices]
println("loaded event data")

println("start calulating")
# Set the maximum lag for cross-correlation
max_lag = 4 # Adjust this according to your needs

# Compute cross-correlations for event
cross_correlation_EE_event = calculate_cross_correlation(E_input_event, E_input_event, max_lag)
cross_correlation_EI_event = calculate_cross_correlation(E_input_event, I_input_event, max_lag)
cross_correlation_IE_event = calculate_cross_correlation(I_input_event, E_input_event, max_lag)
cross_correlation_II_event = calculate_cross_correlation(I_input_event, I_input_event, max_lag)

# Compute cross-correlations for non-event
cross_correlation_EE_non_event = calculate_cross_correlation(E_input_non_event, E_input_non_event, max_lag)
cross_correlation_EI_non_event = calculate_cross_correlation(E_input_non_event, I_input_non_event, max_lag)
cross_correlation_IE_non_event = calculate_cross_correlation(I_input_non_event, E_input_non_event, max_lag)
cross_correlation_II_non_event = calculate_cross_correlation(I_input_non_event, I_input_non_event, max_lag)

# Function to create the plot based on the calculated cross-correlations
function create_plot(title, cross_correlation_EE, cross_correlation_EI, cross_correlation_IE, cross_correlation_II)
    lags = -max_lag:max_lag
    plot(lags, [cross_correlation_EE, cross_correlation_EI, cross_correlation_IE, cross_correlation_II], label=["E-E" "E-I" "I-E" "I-I"], title=title, xlabel="Tau", ylabel="Correlation")
end

# Create and save the plots
plot_event = create_plot("Cross-correlation (Event)", cross_correlation_EE_event, cross_correlation_EI_event, cross_correlation_IE_event, cross_correlation_II_event)
plot_non_event = create_plot("Cross-correlation (Non-event)", cross_correlation_EE_non_event, cross_correlation_EI_non_event, cross_correlation_IE_non_event, cross_correlation_II_non_event)

savefig(plot_event, "cross_correlation_event.png")
savefig(plot_non_event, "cross_correlation_non_event.png")
