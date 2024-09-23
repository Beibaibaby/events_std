using Plots
using Measures

ENV["GKSwstype"] = "100"

# Define the plot size and margins
plot_size = (200, 350)
plot_margin = 5mm

# Define the coordinates for the three dots
x = [0.24, 0.24, 0.24]
y = [0.85, 0.37, 0.24]

# Create the plot
p = plot(x, y, seriestype = :scatter,
    size = plot_size,
    xlabel = "SRD",
    ylabel = "SRF",
    xticks = [0.24],  # Only show the tick for 0.24
    xlim=(0.15,0.33),
    ylim=(0,1),
    yticks = :auto,  # Automatically set y-axis ticks
    legend = false,  # No legend
    grid = false,
    markersize = 1,
    xtickfont = font(10),
    ytickfont = font(10),
    legendfontsize = 10,
    titlefontsize = 10,
    xlabelfontsize = 12,
    ylabelfontsize = 12,
    foreground_color_legend=nothing,
    left_margin=plot_margin, right_margin=plot_margin, top_margin=plot_margin, bottom_margin=plot_margin,dti=500)

# Save the plot
directory_for_plot = "/gpfs/data/doiron-lab/draco/Balanced_Sipiking/plots_for_paper"

# Save the E-I plot as SVG, PDF, and PNG
savefig(p, joinpath(directory_for_plot, "location.svg"))
savefig(p, joinpath(directory_for_plot, "location.pdf"))
savefig(p, joinpath(directory_for_plot, "location.png"))

println("Small plot saved.")