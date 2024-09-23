using Plots
using JLD2, FileIO
using Statistics  # for mean
ENV["GKSwstype"] = "100"  # Use the off-screen GKS terminal to aviod running messages
using Measures

function load_and_plot_activity(main_dir::String, output_file::String)
    println("Loading data...")

    # Define file paths
    avg_e_path = joinpath(main_dir, "avg_e_rate_raw_after_peak.jld2")
    std_e_path = joinpath(main_dir, "std_e_rate_raw_after_peak.jld2")
    avg_i_path = joinpath(main_dir, "avg_i_rate_raw_after_peak.jld2")
    std_i_path = joinpath(main_dir, "std_i_rate_raw_after_peak.jld2")

    # Load data
    avg_e_activity = load(avg_e_path, "avg_e_rate_raw_after_peak")
    std_e_activity = load(std_e_path, "std_e_rate_raw_after_peak")
    avg_i_activity = load(avg_i_path, "avg_i_rate_raw_after_peak")
    std_i_activity = load(std_i_path, "std_i_rate_raw_after_peak")

    println("Data loaded successfully, now processing...")

    # Bin the data into chunks of 10 and compute the mean for each bin
    function bin_data(data, bin_size)
        binned_data = [mean(data[i:min(i+bin_size-1, end)]) for i in 1:bin_size:length(data)]
        return binned_data
    end

    bin_size = 50  # size of each bin (10 data points, representing 1ms each)
    avg_e_activity_binned = bin_data(avg_e_activity, bin_size)
    std_e_activity_binned = bin_data(std_e_activity, bin_size)
    avg_i_activity_binned = bin_data(avg_i_activity, bin_size)
    std_i_activity_binned = bin_data(std_i_activity, bin_size)

    # Define the time step in seconds (each bin represents 1ms)
    time_step = 0.005  # 1ms

    # Create time axis in seconds
    total_time_ms = length(avg_e_activity_binned) * time_step  # Total time in ms
    time_axis = 0:time_step:(total_time_ms - time_step)

    # Select data within 200ms to 700ms
    select_time = (time_axis .>= 0.2) .& (time_axis .<= 0.7)
    time_axis_selected = time_axis[select_time]
    avg_e_selected = avg_e_activity_binned[select_time]
    std_e_selected = std_e_activity_binned[select_time]
    avg_i_selected = avg_i_activity_binned[select_time]
    std_i_selected = std_i_activity_binned[select_time]

    left_margin = 20mm
    p = plot(size=(800, 400), left_margin=left_margin, ylim=(0, 15), bottom_margin=5mm)

    ribbon_alpha = 0.1
    # Plot I rates with error ribbon
    plot!(p, time_axis_selected, avg_i_selected, ribbon = std_i_selected, label="I Rates", color=:blue, fillalpha=ribbon_alpha)
    plot!(p, time_axis_selected, avg_e_selected, ribbon = std_e_selected, label="E Rates", color=:red, fillalpha=ribbon_alpha)
    
    # Increase the size of x and y tick labels
    plot!(p, xtickfontsize=20, ytickfontsize=20)
    
    # Set the labels with increased size and add the title
    xlabel!(p, "Time (s)", xguidefontsize=20)
    ylabel!(p, "Average Rate", yguidefontsize=20)
    title!(p, "Average Neural Activity (200ms - 700ms)")
    
    # Save the figure to the specified output file
    savefig(p, output_file)
    println("Plot saved to $output_file")
end
# Example usage
main_dir = "/gpfs/data/doiron-lab/draco/results_nn/exp_9099088"
output_file = joinpath(main_dir, "avg_activity_plot.png")
load_and_plot_activity(main_dir, output_file)
