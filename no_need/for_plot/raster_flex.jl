using FileIO
using JLD2
using Plots
using StatsBase
using Statistics
using Measures

ENV["GKSwstype"] = "100"
println("starting")

# Define the directory and file paths
directory = "/gpfs/data/doiron-lab/draco/results_nn/exp_flex_10329684/17_good"
times_file = joinpath(directory, "times.jld2")
ns_file = joinpath(directory, "ns.jld2")

# Load the data from the JLD2 files
times = load(times_file, "times")
ns = load(ns_file, "ns")

# Plot size and margin settings
plot_size = (800, 400)
plot_margin = 5mm

# Generate the raster plot
p = plot(size=plot_size, xlim=(0, 2500), ylim=(0, 4000), xlabel="Time (ms)", ylabel="Neuron", dpi=500,
         left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin,
         xtickfont=font(10), ytickfont=font(10), legendfontsize=10, titlefontsize=10, xlabelfontsize=10, ylabelfontsize=10)

for ci in 1:4000
    vals = times[ci, 1:ns[ci]]
    y = fill(ci, length(vals))
    scatter!(vals, y, markersize=1, color=:black, marker=:circle, legend=false)
end

# Define the directory to save the plot
directory_for_plot = "/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper"

# Save the raster plot as SVG, PDF, and PNG
savefig(p, joinpath(directory_for_plot, "raster_flex.svg"))
savefig(p, joinpath(directory_for_plot, "raster_flex.pdf"))
savefig(p, joinpath(directory_for_plot, "raster_flex.png"))

println("Raster plot saved.")