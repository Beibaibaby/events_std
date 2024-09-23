using JLD2
using FilePathsBase
using Plots
include("sim_6.jl") 
using JLD2
using Measures
using Statistics

main_dir = "/gpfs/data/doiron-lab/draco/results_corr/exp_3421191"
# Modified function to average raw activity and compute standard deviation



# ... [previous parts of the script] ...

function average_raw_activity_and_std(main_dir::String, file_name::String, max_length::Int64)
    total_activity = Vector{Float64}()
    squared_activity = Vector{Float64}()
    count = 0

    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            file_path = joinpath(sub_dir, file_name)
            if isfile(file_path)
                raw_activity_data = load(file_path)
                key = file_name == "e_rate_raw_after_peak.jld2" ? "e_rate_raw_after_peak" : "i_rate_raw_after_peak"
                raw_activity = raw_activity_data[key]

                # Filter out excessively long trials
                filtered_activity = filter(a -> length(a) <= max_length, raw_activity)

                for activity in filtered_activity
                    # Resize arrays if necessary and initialize new elements
                    if length(total_activity) < length(activity)
                        resize!(total_activity, length(activity))
                        total_activity[(end - length(activity) + 1):end] .= 0.0
                        resize!(squared_activity, length(activity))
                        squared_activity[(end - length(activity) + 1):end] .= 0.0
                    end

                    # Sum up activities and squared activities
                    for i in 1:length(activity)
                        total_activity[i] += activity[i]
                        squared_activity[i] += activity[i]^2
                    end
                    count += 1
                end
            end
        end
    end

    mean_activity = count > 0 ? total_activity ./ count : total_activity
    variance_activity = count > 0 ? (squared_activity ./ count) .- mean_activity.^2 : squared_activity
    std_activity = sqrt.(variance_activity)

    return mean_activity, std_activity
end

# ... [rest of the script] ...



# Function to plot average raw activities with error ribbons
function plot_avg_raw_activity(time_step::Int, avg_e_activity, std_e_activity, avg_i_activity, std_i_activity, output_file::String)
    time_axis = 0:time_step:(length(avg_e_activity) - 1) * time_step
    
    left_margin = 20mm
    p = plot(size=(800, 400),left_margin=left_margin)

    ribbon_alpha = 0.1
    # Plot E rates with error ribbon
    plot!(p, time_axis, avg_e_activity, ribbon = std_e_activity, label="E Rates", color=:red, fillalpha=ribbon_alpha)

    # Plot I rates with error ribbon
    plot!(p, time_axis, avg_i_activity, ribbon = std_i_activity, label="I Rates", color=:blue, fillalpha=ribbon_alpha)

    xlabel!(p, "Time (ms)")
    ylabel!(p, "Average Rate")
    title!(p, "Average Neural Activity Over Time")

    savefig(p, output_file)
end

# Main script execution

file_e = "e_rate_raw_after_peak.jld2"
file_i = "i_rate_raw_after_peak.jld2"
max_trial_length = 52 # Define buffer_before and win_buff as needed
time_step = 5 # Define your time step size in ms

# Average raw activities and compute standard deviations
avg_e_rate_raw_after_peak, std_e_rate_raw_after_peak = average_raw_activity_and_std(main_dir, file_e, max_trial_length)
avg_i_rate_raw_after_peak, std_i_rate_raw_after_peak = average_raw_activity_and_std(main_dir, file_i, max_trial_length)

# Save the averaged raw activities
@save joinpath(main_dir, "avg_e_rate_raw_after_peak.jld2") avg_e_rate_raw_after_peak
@save joinpath(main_dir, "std_e_rate_raw_after_peak.jld2") std_e_rate_raw_after_peak
@save joinpath(main_dir, "avg_i_rate_raw_after_peak.jld2") avg_i_rate_raw_after_peak
@save joinpath(main_dir, "std_i_rate_raw_after_peak.jld2") std_i_rate_raw_after_peak

# Plotting and saving the plot
plot_file = joinpath(main_dir, "average_raw_activity_plot.png")
plot_avg_raw_activity(time_step, avg_e_rate_raw_after_peak, std_e_rate_raw_after_peak, avg_i_rate_raw_after_peak, std_i_rate_raw_after_peak, plot_file)


