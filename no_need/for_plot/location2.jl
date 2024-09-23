using Plots
using Measures

ENV["GKSwstype"] = "100"

# Define the plot size and margins
plot_size = (350, 600)
plot_margin = 5mm

# Define the coordinates for the three dots
x = [0.001, 0.24, 0.24, 0.24,0.25,0.25]
y = [0.0, 0.85, 0.37, 0.24,0.0,0.3]

# Create the plot
p = plot(x, y, seriestype = :scatter,
    size = plot_size,
    xlabel = "SRD",
    ylabel = "SRF",
    xticks = [ 0.24],  # Only show the tick for 0.24
    xlim=(0,0.5),
    ylim=(0,1),
    yticks = :auto,  # Automatically set y-axis ticks
    legend = false,  # No legend
    grid = false,
    markersize = 1,
    xtickfont = font(18),
    ytickfont = font(18),
    legendfontsize = 18,
    titlefontsize = 18,
    xlabelfontsize = 18,
    ylabelfontsize = 18,
    foreground_color_legend=nothing,
    left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin,dti=500)

# Save the plot
directory_for_plot = "/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper"

# Save the E-I plot as SVG, PDF, and PNG
savefig(p, joinpath(directory_for_plot, "location2_ie.svg"))
savefig(p, joinpath(directory_for_plot, "location2_ie.pdf"))
savefig(p, joinpath(directory_for_plot, "location2_ie.png"))

println("IE plot saved.")



##########EE starts
# Define the plot size and margins
plot_size = (600, 600)
plot_margin = 5mm

# Define the coordinates for the three dots
x = [0.98, 0.24, 0.24, 0.24,0.15,0.95]
y = [0.01, 0.85, 0.85, 0.85,0.3,0.01]

# Create the plot
p = plot(x, y, seriestype = :scatter,
    size = plot_size,
    xlabel = "SRD",
    ylabel = "SRF",
    #xticks = [ 0.15, 0.24, 1],  # Only show the tick for 0.24
    xlim=(0,1.05),
    ylim=(0,1),
    yticks = :auto,  # Automatically set y-axis ticks
    legend = false,  # No legend
    grid = false,
    markersize = 1,
    xtickfont = font(18),
    ytickfont = font(18),
    legendfontsize = 18,
    titlefontsize = 18,
    xlabelfontsize = 18,
    ylabelfontsize = 18,
    foreground_color_legend=nothing,
    left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin,dti=500)

# Save the plot
directory_for_plot = "/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper"

# Save the E-I plot as SVG, PDF, and PNG
savefig(p, joinpath(directory_for_plot, "location2_ee.svg"))
savefig(p, joinpath(directory_for_plot, "location2_ee.pdf"))
savefig(p, joinpath(directory_for_plot, "location2_ee.png"))

println("EE plot saved.")