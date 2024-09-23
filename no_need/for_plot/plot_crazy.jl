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

# Set the GKS terminal for off-screen rendering to avoid running messages
ENV["GKSwstype"] = "100"

# Define the directory and file paths
directory = "/gpfs/data/doiron-lab/draco/results_crazy/exp_crazy_11269599/5"
e_rate_file = joinpath(directory, "e_rate.jld2")
i_rate_file = joinpath(directory, "i_rate.jld2")

# Load the data from the JLD2 files
e_rate = load(e_rate_file, "e_rate")
i_rate = load(i_rate_file, "i_rate")

# Ensure e_rate and i_rate have the same length
n_steps = length(e_rate)

# Define step size and window size
window_size = 100
step_size = 0.1

# Generate time values
time_values = [i * step_size + window_size for i in 1:n_steps]

# Define custom colors
blue_color = RGB(55/255, 120/255, 201/255)
red_color = RGB(228/255, 26/255, 28/255)
green_color = RGB(0/255, 201/255, 20/255)
grey_color = RGB(145/255, 55/255, 19/255)

# Define smoothing function
function moving_average(data, window_size)
    return [mean(data[i:i+window_size-1]) for i in 1:(length(data)-window_size+1)]
end

# Smooth the data
window_size_smoothing = 80  # Adjust the window size for smoothing as needed
e_rate_smooth = moving_average(e_rate, window_size_smoothing)
i_rate_smooth = moving_average(i_rate, window_size_smoothing)

# Adjust time values to match the smoothed data
time_values_smooth = time_values[1:length(e_rate_smooth)]

# Filter the data to include only the desired time range (0 ms to 2500 ms)
start_time = 0
end_time = 2000

filtered_indices = (start_time .<= time_values_smooth .<= end_time)
time_values_filtered = time_values_smooth[filtered_indices]
e_rate_filtered = e_rate_smooth[filtered_indices]
i_rate_filtered = i_rate_smooth[filtered_indices]

# Shift the time values to start from zero
time_values_filtered .-= time_values_filtered[1]

# Plot size and margin settings
plot_size = (1200, 600)
plot_margin = 10mm

# Function to create a rectangle shape
rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])

# Generate the plot
p2 = plot(time_values_filtered, e_rate_filtered,
    xlabel = "Time (ms)",
    ylabel = "Firing rate (Hz)",
    label = "Excitatory",
    lw = 3,
    linecolor = red_color,
    size = plot_size,
    xtickfont = font(18),
    ytickfont = font(18),
    legendfontsize = 16,
    titlefontsize = 16,
    xlabelfontsize = 18,
    ylabelfontsize = 18,
    left_margin = plot_margin,
    bottom_margin = plot_margin,
    grid = false,
    ylim = (-5, 200),  # Set y-axis limits
    dpi=500,
    legend = false,
    foreground_color_legend=nothing)

plot!(p2, time_values_filtered, i_rate_filtered,
    label = "Inhibitory",
    lw = 3,
    linecolor = blue_color,
    grid = false,
    foreground_color_legend=nothing)

# Adding shaded regions
plot!(p2, rectangle(20, 205, 100, -5), fillcolor=green_color, fillalpha=0.1, seriestype=:shape, legend=false, linecolor=nothing)
plot!(p2, rectangle(100, 205, 1000, -5), fillcolor=grey_color, fillalpha=0.1, seriestype=:shape, legend=false, linecolor=nothing)

# Save the plot
directory_for_plot = "/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper"
savefig(p2, joinpath(directory_for_plot, "plot_crazy.svg"))
savefig(p2, joinpath(directory_for_plot, "plot_crazy.pdf"))
savefig(p2, joinpath(directory_for_plot, "plot_crazy.png"))