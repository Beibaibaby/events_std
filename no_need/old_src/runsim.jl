#this file is part of litwin-kumar_doiron_cluster_2012
#Copyright (C) 2014 Ashok Litwin-Kumar
#see README for more information

using Statistics
#uncomment the line below and set doplot=true to plot a raster
#using PyPlot
doplot = false

include("sim.jl")

times,ns,Ne,Ncells,T = sim()

println("mean excitatory firing rate: ",mean(1000*ns[1:Ne]/T)," Hz")
println("mean inhibitory firing rate: ",mean(1000*ns[(Ne+1):Ncells]/T)," Hz")

if doplot 
	println("creating plot")
	figure(figsize=(4,4))
	for ci = 1:Ne
		vals = times[ci,1:ns[ci]]
		y = ci*ones(length(vals))
		scatter(vals,y,s=.3,c="k",marker="o",linewidths=0)
	end
	xlim(0,T)
	ylim(0,Ne)
	ylabel("Neuron")
	xlabel("Time")
	tight_layout()
	savefig("output.png",dpi=150)
end
