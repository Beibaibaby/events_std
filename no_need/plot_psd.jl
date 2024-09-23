using JLD2, FileIO
using Plots
ENV["GKSwstype"] = "100"  # Use the off-screen GKS terminal

# Paths to your files
main_dir_1 = "/gpfs/data/doiron-lab/draco/results_nn/exp_8905776"
main_dir_2 = "/gpfs/data/doiron-lab/draco/results_nn/exp_8905827"

# Load the data
@load joinpath(main_dir_1, "all_psds.jld2") all_psds
@load joinpath(main_dir_1, "mean_psd.jld2") mean_psd
@load joinpath(main_dir_1, "std_psd.jld2") std_psd
println("Loaded first set of data")
all_psds_1= all_psds
mean_psd_1 = mean_psd
std_psd_1 = std_psd
#load second set but name them differently
@load joinpath(main_dir_2, "mean_psd.jld2") mean_psd
@load joinpath(main_dir_2, "std_psd.jld2") std_psd

mean_psd_2 = mean_psd
std_psd_2 = std_psd

# Assuming the frequency vector is the same for both and stored in `all_psds_1.freq`
frequency = first(all_psds_1).freq

# Plotting
p = plot(size=(800, 600), xlabel="Frequency (Hz)", ylabel="Normalized Power (0-1)", title="Comparison of Normalized Power Spectral Density")

# Plot first set
plot!(p, frequency, mean_psd_1, ribbon=std_psd_1, label="exp_8905776 P", fillalpha=0.2, color=:red)

# Plot second set
plot!(p, frequency, mean_psd_2, ribbon=std_psd_2, label="exp_8905827 No P", fillalpha=0.2, color=:blue)

# Save the figure
output_file = "/gpfs/data/doiron-lab/draco/results_nn/comparison_psd_plot.png" # or .svg, .pdf, etc, depending on your preference
savefig(p, output_file)
