using Statistics
using Plots
using Dates
using JSON3
using Random
using JLD2
using Distributed
using SharedArrays
using Measures
using Distributions
using FileIO
using StatsBase

# Set the GKS terminal for off-screen rendering to avoid running messages
ENV["GKSwstype"] = "100"
println("starting")


# Define the directory and file paths
directory = "/gpfs/data/doiron-lab/draco/results_nn/exp_same_10549930/30"
e_rate_file = joinpath(directory, "e_rate.jld2")
i_rate_file = joinpath(directory, "i_rate.jld2")
times_file = joinpath(directory, "times.jld2")
ns_file = joinpath(directory, "ns.jld2")
record_kick_file = joinpath(directory, "record_kick.jld2")

# Load the data from the JLD2 files
e_rate = load(e_rate_file, "e_rate")
i_rate = load(i_rate_file, "i_rate")
times = load(times_file, "times")
ns = load(ns_file, "ns")
record_kick = load(record_kick_file, "record_kick")
println(record_kick)
# Ensure e_rate and i_rate have the same length
n_steps = length(e_rate)

# Define step size and window size
window_size = 25
step_size = 0.1

# Generate time values
time_values = [i * step_size + window_size for i in 1:n_steps]

# Define custom colors
blue_color = RGB(55/255, 120/255, 201/255)
red_color = RGB(228/255, 26/255, 28/255)

# Define smoothing function
function moving_average(data, window_size)
    return [mean(data[i:i+window_size-1]) for i in 1:(length(data)-window_size+1)]
end

# Smooth the data
window_size_smoothing = 50  # Adjust the window size for smoothing as needed
e_rate_smooth = moving_average(e_rate, window_size_smoothing)
i_rate_smooth = moving_average(i_rate, window_size_smoothing)

# Adjust time values to match the smoothed data
time_values_smooth = time_values[1:length(e_rate_smooth)]

# Filter the data to include only the desired time range (0 ms to 2500 ms)
start_time = 0
end_time = 2500

filtered_indices = (start_time .<= time_values_smooth .<= end_time)
time_values_filtered = time_values_smooth[filtered_indices]
e_rate_filtered = e_rate_smooth[filtered_indices]
i_rate_filtered = i_rate_smooth[filtered_indices]

# Plot size and margin settings
plot_size = (1000, 700)
plot_margin = 5mm

# Create the smoothed firing rates plot
p1 = plot(time_values_filtered, i_rate_filtered,
    xlabel = "Time (ms)",
    ylabel = "Firing rate (Hz)",
    label = "Inhibitory",
    lw = 2,
    linecolor = blue_color,
    xlim=(0, 2500),
    xtickfont = font(12),
    ytickfont = font(12),
    legendfontsize = 12,
    titlefontsize = 12,
    xlabelfontsize = 12,
    ylabelfontsize = 12,
    grid = false,
    ylim = (0, 15),  # Set y-axis limits
    dpi=500,
    foreground_color_legend=nothing,
    left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin)

plot!(p1, time_values_filtered, e_rate_filtered,
    label = "Excitatory",
    lw = 2,
    linecolor = red_color,
    foreground_color_legend=nothing,
    grid = false)

# Create the raster plot
p2 = plot(size=(1000, 400), xlim=(0, 2500), ylim=(0, 4000), xlabel="Time (ms)", ylabel="Neuron", dpi=500,
         left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin,
         xtickfont=font(12), ytickfont=font(12), legendfontsize=12, titlefontsize=12, xlabelfontsize=12, ylabelfontsize=12)

for ci in 1:4000
    vals = times[ci, 1:ns[ci]]
    y = fill(ci, length(vals))
    scatter!(p2, vals, y, markersize=1, color=:black, marker=:circle, legend=false)
end

# Combine the plots into a single figure with subplots
combined_plot = plot(p2, p1, layout = @layout([a; b]), size=plot_size, link=:x, dpi=500)

# Define the directory to save the plot
directory_for_plot = "/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper"

# Save the combined plot as SVG, PDF, and PNG
savefig(combined_plot, joinpath(directory_for_plot, "combined_same.svg"))
savefig(combined_plot, joinpath(directory_for_plot, "combined_same.pdf"))
savefig(combined_plot, joinpath(directory_for_plot, "combined_same.png"))

println("Combined plot saved.")


# Define the time range
start_time = 1182.9
end_time = 1262.9

# Initialize an array to store spike counts for each neuron
spike_counts = zeros(Int, 4000)

# Count the spikes within the specified time range for each neuron
for ci in 1:4000
    vals = times[ci, 1:ns[ci]]
    spike_counts[ci] = sum((start_time .<= vals) .& (vals .<= end_time))
end

# Calculate the number of neurons with zero spike counts
zero_count_neurons = sum(spike_counts .== 0)
println("Number of neurons with zero spike counts: ", zero_count_neurons)

# Sort the spike counts in descending order and keep the indices
sorted_indices = sortperm(spike_counts, rev=true)
sorted_spike_counts = spike_counts[sorted_indices]

# Plot the histogram of sorted spikes
bar(1:4000, sorted_spike_counts,
    xlabel="Neuron (sorted)",
    ylabel="Spike Count",
    title="Spike Counts per Neuron (1182.9 ms to 1262.9 ms)",
    legend=false,
    dpi=500,
    xtickfont=font(12),
    ytickfont=font(12),
    xlabelfontsize=12,
    ylabelfontsize=12,
    grid=false)

# Save the plot as SVG and PNG
savefig("spike_counts_histogram_sorted.png")