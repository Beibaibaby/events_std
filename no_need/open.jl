using JLD2

function load_data_from_jld2(file_path, data_key)
    # Open the .jld2 file and read the data associated with the given key
    data = jldopen(file_path, "r") do file
        return file[data_key]
    end

    return data
end

# Specify the file path
file_path = "/gpfs/data/doiron-lab/draco/results_new/d_ee=0.24+f_ie=0.0+d_ie=0.24+2024-01-10_15-07-18/top_n_e_neurons_noise.jld2"

# Specify the data key
data_key = "top_n_e_neurons"

# Call the function to load data
loaded_data = load_data_from_jld2(file_path, data_key)

# Now you can process or display 'loaded_data' as needed.
# For example, if it's an array, you can print it or perform calculations.
# If it's more complex, adjust your processing accordingly.

println("Loaded data: ")
println(loaded_data)

ci = 158
if ci in loaded_data
    # Your statements here
    println("ci is in the top_n_e_neurons list!")
    # Add more statements as needed
end