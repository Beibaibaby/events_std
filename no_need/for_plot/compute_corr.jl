using FileIO, JLD2, Plots, StatsBase, Statistics

ENV["GKSwstype"] = "100"
println("starting")
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

# Load data from the non-event directory
dir_name = "/gpfs/data/doiron-lab/draco/results_corr/0.0-real"
@load joinpath(dir_name, "E_input.jld2") E_input
@load joinpath(dir_name, "I_input.jld2") I_input
#@load joinpath(dir_name, "v_history.jld2") v_history


# Assuming E_input and I_input are 2D arrays where each row corresponds to a neuron's input data
window_size = 350  # Window size of 100 data points
E_input = sliding_window_average_matrix(E_input, window_size,step_size)
I_input = sliding_window_average_matrix(I_input, window_size,step_size)

E_input_non_event = E_input
I_input_non_event = I_input
#v_history_non_event = v_history


# Load data from the event directory
dir_name_event = "/gpfs/data/doiron-lab/draco/results_corr/5.4-real"
@load joinpath(dir_name_event, "E_input.jld2") E_input
@load joinpath(dir_name_event, "I_input.jld2") I_input
#@load joinpath(dir_name, "v_history.jld2") v_history
@load joinpath(dir_name_event, "times.jld2") times
@load joinpath(dir_name_event, "ns.jld2") ns


E_input = sliding_window_average_matrix(E_input, window_size,step_size)
I_input = sliding_window_average_matrix(I_input, window_size,step_size)

E_input_event = E_input
I_input_event = I_input
#v_history_event = v_history






plot(size=(1000, 600), xlim=(0, 2500), ylim=(0, 4000), xlabel="Time", ylabel="Neuron")


for ci in 1:4000
    
    vals = times[ci, 1:ns[ci]]
    y = fill(ci, length(vals))
    scatter!(vals, y, markersize=1, color=:black, marker=:circle, legend=false)
    

end

savefig("output.png")

# Parameters
neurons_size = 4000
max_time = 2500  # This should be set to the maximum time you want to consider.

# Initialize a sparse matrix for efficiency, given most values will be 0
spike_matrix = zeros(neurons_size, max_time)

# Populate the spike matrix
for neuron in 1:neurons_size
    for j in 1:ns[neuron]
        spike_time = Int(round(times[neuron, j]+1))
        
        if spike_time <= max_time
            spike_matrix[neuron, spike_time] = 1
        end
    end
end



# Time window of interest
#time_indices = 5800:6800  # corresponds to 100 ms to 200 ms
time_indices = Int(round(1000/step_size)):Int(round(2200/step_size))

function average_correlation(input_data, neuron_indices, num_pairs=5000)
    correlations = Float64[]
    for _ = 1:num_pairs
        pair = sample(neuron_indices, 2, replace=false)
        cor = Statistics.cor(input_data[pair[1], time_indices], input_data[pair[2], time_indices])
        push!(correlations, cor)
    end
    mean(correlations)
end
function load_data_from_jld2(file_path, data_key)
    # Open the .jld2 file and read the data associated with the given key
    data = jldopen(file_path, "r") do file
        return file[data_key]
    end

    return data
end

# Assuming both sets have the same neurons to skip
A = load_data_from_jld2("/gpfs/data/doiron-lab/draco/results_150_old/d_ee=0.24+f_ie=0.0+d_ie=0.24+2024-03-20_12-52-39/top_n_e_neurons_noise.jld2", "top_n_e_neurons")#[1155, 1236, 995, 1417, 140, 208, 512, 775, 883, 969, 1122, 1398, 1950, 2086, 2545, 3552, 3998, 21, 98, 117, 302, 389, 392, 430, 494, 547, 643, 681, 808, 893, 965, 1006, 1184, 1187, 1229, 1401, 1647, 1715, 1735, 1933, 2176, 2183, 2292, 2298, 2379, 2392, 2438, 2465, 2509, 2609, 2672, 2698, 2724, 2891, 2892, 3074, 3092, 3112, 3222, 3310, 3376, 3407, 3429, 3435, 3497, 3910, 12, 69, 100, 114, 136, 144, 162, 223, 238, 303, 349, 358, 369, 385, 411, 426, 445, 451, 504, 563, 580, 595, 646, 711, 776, 786, 800, 849, 876, 879, 897, 940, 1056, 1065, 1081, 1088, 1099, 1135, 1180, 1200, 1201, 1205, 1262, 1300, 1322, 1349, 1443, 1488, 1489, 1491, 1587, 1624, 1691, 1713, 1732, 1749, 1753, 1755, 1769, 1852, 1868, 1935, 1973, 1997, 2008, 2018, 2047, 2077, 2147, 2167, 2193, 2196, 2220, 2222, 2252, 2263, 2301, 2420, 2440, 2488, 2503, 2513, 2588, 2708, 2711, 2740, 2747, 2811, 2843, 2853, 2859, 2886, 3049, 3086, 3167, 3169, 3186, 3208, 3260, 3265, 3271, 3345, 3349, 3350, 3483, 3484, 3502, 3546, 3559, 3572, 3580, 3625, 3631, 3644, 3685, 3689, 3719, 3732, 3754, 3756, 3757, 3764, 3769, 3770, 3841, 3866, 3916, 3922, 3948, 3965, 1, 10, 13, 42]

skip_indices = Set(A)
E_neurons = setdiff(1:4000, skip_indices)
I_neurons = setdiff(4001:5000, skip_indices)
# Compute correlations for non-event data
E_E_correlation_non_event = average_correlation(E_input_non_event, E_neurons)
I_I_correlation_non_event = average_correlation(I_input_non_event, E_neurons)
Total_Total_correlation_non_event = average_correlation(E_input_non_event+I_input_non_event, E_neurons)


time_indices = Int(round(5800/step_size)):Int(round(7000/step_size))

function average_correlation(input_data, neuron_indices, num_pairs=15000)
    correlations = Float64[]
    for _ = 1:num_pairs
        pair = sample(neuron_indices, 2, replace=false)
        cor = Statistics.cor(input_data[pair[1], time_indices], input_data[pair[2], time_indices])
        push!(correlations, cor)
    end
    mean(correlations)
end
# Compute correlations for event data
E_E_correlation_event = average_correlation(E_input_event, E_neurons)
I_I_correlation_event = average_correlation(I_input_event, E_neurons)
Total_Total_correlation_event = average_correlation(E_input_event+I_input_event, E_neurons)
# Define the correlation values and labels
correlation_values_non_event = [E_E_correlation_non_event, Total_Total_correlation_non_event, I_I_correlation_non_event]
correlation_values_event = [E_E_correlation_event, Total_Total_correlation_event, I_I_correlation_event]
labels = ["EE", "Composite", "II"]

using Plots

# Sample data
labels = ["EE", "Composite", "II"]
correlation_values_non_event = [0.25, 0.1, 0.3]
correlation_values_event = [1.0, 0.3, 1.0]

# Initial plot for "Non-event"
plot(labels, correlation_values_non_event, 
     label="Non-event", 
     seriestype=:scatterpath, 
     color=:black, 
     ylabel="Correlation", 
     legend=:top, 
     markerstrokecolor=:black,
     markersize=12,
     ylim=(0.0, 1.05), 
     size=(600,500),
     grid=false,
     xtickfont=font(18),
     ytickfont=font(18),
     legendfontsize=18,
     titlefontsize=18,
     xlabelfontsize=18,
     ylabelfontsize=18,
     lw=2,
     dti=500, # Remove legend title border
     legendcolor=:white,  # Legend background color
     legendstroke=:transparent, # Legend border
     foreground_color_legend=nothing,
     )  # Remove legend border

# Add plot for "Event"
plot!(labels, correlation_values_event, 
      label="Event", 
      seriestype=:scatterpath, 
      color=:grey,
      markerstrokecolor=:grey, 
      markersize=12,
      legend=:top,
      lw=2,
      dti=500)  # Set legend text colors

# Display the plot
savefig("/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper/Neuron_Correlations_comparison.png")
savefig("/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper/Neuron_Correlations_comparison.svg")


# Cell selection
cell_indices = [5, 10]  # Ensure these are not in skip_indices and adjust as necessary

# Extract voltage data for selected time_indices and cells
#v_non_event_cell1 = v_history_non_event[cell_indices[1], time_indices]
#v_non_event_cell2 = v_history_non_event[cell_indices[2], time_indices]
#v_event_cell1 = v_history_event[cell_indices[1], time_indices]
#v_event_cell2 = v_history_event[cell_indices[2], time_indices]

# Create subplots
#p1 = plot(time_indices, v_non_event_cell1, title="Non-event Cell 1", label="Voltage", color=:blue)
#p2 = plot(time_indices, v_non_event_cell2, title="Non-event Cell 2", label="Voltage", color=:blue)
#p3 = plot(time_indices, v_event_cell1, title="Event Cell 1", label="Voltage", color=:red)
#p4 = plot(time_indices, v_event_cell2, title="Event Cell 2", label="Voltage", color=:red)

# Combine into a single figure
#plot(p1, p2, p3, p4, layout=(2, 2), legend=false, xlabel="Time Index", ylabel="Voltage")

# Save the figure
#savefig("Voltage_Traces_comparison.png")

function average_input_correlation(E_input, I_input, neuron_indices, time_indices)
    correlations = Float64[]
    for neuron in neuron_indices
        cor = Statistics.cor(E_input[neuron, time_indices], I_input[neuron, time_indices])
        push!(correlations, cor)
    end
    mean(correlations)
end


# Define neuron indices, ensuring they are not skipped
neuron_indices = setdiff(1:size(E_input_non_event, 1), skip_indices)  # Adjust as per neuron index range and skip list

# Define time indices for correlation analysis
time_window_start = 1000  # Adjust start time as needed
time_window_end = 2200    # Adjust end time as needed
time_indices = Int(round(time_window_start/step_size)):Int(round(time_window_end/step_size))

E_I_correlation = average_input_correlation(E_input_non_event, I_input_non_event, neuron_indices, time_indices)

println("Correlation E-I nonevent: $E_I_correlation")
