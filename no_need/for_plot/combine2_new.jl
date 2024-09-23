using Statistics
using Plots
using JLD2
using Measures
ENV["GKSwstype"] = "100"


# Function to load and process data
function load_and_process_data(directory, e_rate_file_name, i_rate_file_name, start_time, end_time, window_size_smoothing=50, step_size=0.1, window_size=25)
    e_rate_file = joinpath(directory, e_rate_file_name)
    i_rate_file = joinpath(directory, i_rate_file_name)

    # Load the data from the JLD2 files
    e_rate = load(e_rate_file, "e_rate")
    i_rate = load(i_rate_file, "i_rate")

    # Ensure e_rate and i_rate have the same length
    n_steps = length(e_rate)

    # Generate time values
    time_values = [i * step_size + window_size for i in 1:n_steps]

    # Define smoothing function
    function moving_average(data, window_size)
        return [mean(data[i:i+window_size-1]) for i in 1:(length(data)-window_size+1)]
    end

    # Smooth the data
    e_rate_smooth = moving_average(e_rate, window_size_smoothing)
    i_rate_smooth = moving_average(i_rate, window_size_smoothing)

    # Adjust time values to match the smoothed data
    time_values_smooth = time_values[1:length(e_rate_smooth)]

    # Filter the data to include only the desired time range
    filtered_indices = (start_time .<= time_values_smooth .<= end_time)
    time_values_filtered = time_values_smooth[filtered_indices]
    e_rate_filtered = e_rate_smooth[filtered_indices]
    i_rate_filtered = i_rate_smooth[filtered_indices]

    return time_values_filtered, e_rate_filtered, i_rate_filtered
end

# Plot settings
plot_size = (400, 350)
plot_margin = 5mm
blue_color = RGB(55/255, 120/255, 201/255)
red_color = RGB(228/255, 26/255, 28/255)
green_color = RGB(0/255, 201/255,20/255)

# Load and process data for each condition
#dir1 = "/gpfs/data/doiron-lab/draco/results_nn/exp_flex_small_10336866_true/16_true"
dir1 ="/gpfs/data/doiron-lab/draco/results_nn/exp_flex_small_10567909_0.37/36"
start_time1 = 750
end_time1 = 1750
time_values1, e_rate1, i_rate1 = load_and_process_data(dir1, "e_rate.jld2", "i_rate.jld2", start_time1, end_time1)
kick_file1 = joinpath(dir1, "record_kick.jld2")
kick_times1 = load(kick_file1, "record_kick") / 10  # Convert to milliseconds

dir2 = "/gpfs/data/doiron-lab/draco/results_nn/exp_fixed_10329654/19_good/"
start_time2 = 700
end_time2 = 1700
time_values2, e_rate2, i_rate2 = load_and_process_data(dir2, "e_rate.jld2", "i_rate.jld2", start_time2, end_time2)

#dir3 = "/gpfs/data/doiron-lab/draco/results_nn/exp_same_10336916/37_good"
dir3="/gpfs/data/doiron-lab/draco/results_nn/exp_same_10549930/31"
start_time3 = 1300
end_time3 = 2300
time_values3, e_rate3, i_rate3 = load_and_process_data(dir3, "e_rate.jld2", "i_rate.jld2", start_time3, end_time3)
kick_file3 = joinpath(dir3, "record_kick.jld2")
kick_times3 = load(kick_file3, "record_kick") / 10  # Convert to milliseconds

# Create plots
p1 = plot(time_values1, i_rate1,
    size = plot_size,
    ylabel = "Firing rate (Hz)",
    label = "Inhibitory",
    lw = 2,
    linecolor = blue_color,
    xlim=(750, 1750),
    ylim=(0, 10),
    xtickfont = font(12),
    ytickfont = font(12),
    legend = false,
    grid = false,
    xaxis=false,
    foreground_color_legend=nothing,
    left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin)
    
    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

    for kick_time1 in kick_times1

        #plot!(p1, rectangle(80,4000,kick_time1,0), fillcolor=green_color, fillalpha=0.1, seriestype=:shape, legend=false, linecolor=nothing)
    end

plot!(p1, time_values1, e_rate1,
    label = "Excitatory",
    lw = 2,
    linecolor = red_color,
    grid = false,
    xaxis=false,
    foreground_color_legend=nothing)

p2 = plot(time_values2, i_rate2,
    size = plot_size,
    ylabel = "",
    label = "Inhibitory",
    lw = 2,
    linecolor = blue_color,
    xlim=(700, 1700),
    ylim=(0, 10),
    xtickfont = font(12),
    ytickfont = font(12),
    legend = false,
    grid = false,
    xaxis=false,
    yaxis=false,
    foreground_color_legend=nothing,
    left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin)

    plot!(p2, time_values2, e_rate2,
    label = "Excitatory",
    lw = 2,
    linecolor = red_color,
    grid = false,
    xaxis=false,
    yaxis=false,
    foreground_color_legend=nothing)

p3 = plot(time_values3, i_rate3,
    size = plot_size,
    ylabel = "",
    label = "Inhibitory",
    lw = 2,
    linecolor = blue_color,
    xlim=(1300, 2300),
    ylim=(0, 10),
    xtickfont = font(12),
    ytickfont = font(12),
    legendfontsize = 12,
    titlefontsize = 12,
    xlabelfontsize = 12,
    ylabelfontsize = 12,
    legend = :topright,
    grid = false,
    xaxis=false,
    yaxis=false,
    foreground_color_legend=nothing,
    left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin)
    

    plot!(p3, time_values3, e_rate3,
    label = "Excitatory",
    lw = 2,
    linecolor = red_color,
    grid = false,
    xaxis=false,
    yaxis=false,
    foreground_color_legend=nothing)

    for kick_time3 in kick_times3
        #plot!(p3, rectangle(80,4000,kick_time3,0), fillcolor=green_color, fillalpha=0.1, seriestype=:shape, label=nothing, linecolor=nothing)
    end


# Combine plots into a single figure with subplots
combined_plot = plot(p1, p3, layout = @layout([a c]), size=(2*plot_size[1], plot_size[2]), link=:y,dpi=500)

# Define the directory to save the plot
directory_for_plot = "/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper"

# Save the combined plot as SVG, PDF, and PNG
savefig(combined_plot, joinpath(directory_for_plot, "combined_2_nos.svg"))
savefig(combined_plot, joinpath(directory_for_plot, "combined_2_nos.pdf"))
savefig(combined_plot, joinpath(directory_for_plot, "combined_2_nos.png"))

println("Combined plot saved.")