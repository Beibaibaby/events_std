using Statistics
using Plots  
using Dates  # for generating timestamps
using JSON3  # Use JSON3 package for JSON handling
using Random
using JLD2
using Distributed
using SharedArrays
using Measures
using Distributions
using SharedArrays

ENV["GKSwstype"] = "100"  # Use the off-screen GKS terminal to aviod running messages

include("sim_7_t150.jl") # include the network simulation code

top_n_e_neurons = load_data_from_jld2("/gpfs/data/doiron-lab/draco/results_150_old/d_ee=0.24+f_ie=0.0+d_ie=0.24+2024-03-20_12-52-39/top_n_e_neurons_noise.jld2", "top_n_e_neurons")

print(size(top_n_e_neurons))