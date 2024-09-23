using JLD2, FileIO, Plots, Statistics, Measures


function plot_rates_with_stats(e_rates, i_rates, labels, output_file)
    # Filtering out values larger than 30
    e_rates_filtered = [filter(rate -> rate <= 25, rates) for rates in e_rates]
    i_rates_filtered = [filter(rate -> rate <= 25, rates) for rates in i_rates]

    e_means = [mean(rates) for rates in e_rates_filtered]
    e_stds = [std(rates) for rates in e_rates_filtered]
    i_means = [mean(rates) for rates in i_rates_filtered]
    i_stds = [std(rates) for rates in i_rates_filtered]

    plot_size = (1000, 450)
    left_margin = 10mm
    mean_marker_size = 10

    p = plot(size=plot_size, left_margin=left_margin, legend=false)

    # Define the x-coordinates for E and I data points
    x_coords_e = 1:2:(2 * length(labels))
    x_coords_i = 2:2:(2 * length(labels))

    for (i, (e_rate, i_rate)) in enumerate(zip(e_rates_filtered, i_rates_filtered))
        scatter!(p, fill(x_coords_e[i], length(e_rate)), e_rate, label=false, color=:pink, alpha=0.5)
        scatter!(p, fill(x_coords_i[i], length(i_rate)), i_rate, label=false, color=:lightblue, alpha=0.5)
        
        scatter!(p, [x_coords_e[i]], [e_means[i]], yerr=e_stds[i], label=false, color=:red, markershape=:cross, markersize=mean_marker_size)
        scatter!(p, [x_coords_i[i]], [i_means[i]], yerr=i_stds[i], label=false, color=:blue, markershape=:cross, markersize=mean_marker_size)
    end

    # Set x-ticks
    x_ticks_labels = [label * "-" * neuron_type for label in labels for neuron_type in ["E", "I"]]
    xticks!(p, 1:(2 * length(labels)), x_ticks_labels)

    ylabel!(p, "Rate")
    xlabel!(p, "Neuron Type and Condition")
    savefig(p, output_file)
end

exp_paths = ["/gpfs/data/doiron-lab/draco/results_corr/exp_3417121",
             "/gpfs/data/doiron-lab/draco/results_corr/exp_3417309",
             "/gpfs/data/doiron-lab/draco/results_corr/exp_3417363"]

labels = ["Hypothesized", "No Facilitation", "No Plasticity"]

e_rates_all = []
i_rates_all = []
for exp_path in exp_paths
    e_data = load(joinpath(exp_path, "merged_e_rates.jld2"))["e_rates_merged"]
    i_data = load(joinpath(exp_path, "merged_i_rates.jld2"))["i_rates_merged"]
    push!(e_rates_all, e_data)
    push!(i_rates_all, i_data)
end

output_file = "/gpfs/data/doiron-lab/draco/results_corr/merged_rates_plot.png"
plot_rates_with_stats(e_rates_all, i_rates_all, labels, output_file)
