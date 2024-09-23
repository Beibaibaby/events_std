using JLD2
using FilePathsBase
using Plots
include("sim_7_t150.jl") 

using JLD2
using Measures
using DSP 
using CSV
using DataFrames


function average_correlations(main_dir::String, file_name::String)
    total_corr = Dict()
    count = 0

    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            file_path = joinpath(sub_dir, file_name)

            if isfile(file_path)
  
                corr_data = load_object(file_path)

                for (key, value) in corr_data
                    total_corr[key] = get(total_corr, key, 0) + value
                end
                count += 1
            end
        end
    end

    # Average the aggregated data
    if count > 0
        for key in keys(total_corr)
            total_corr[key] /= count
        end
    end

    return total_corr
end



function get_arg(key, default)
    index = findfirst(==(key), ARGS)
    if index !== nothing && index < length(ARGS)
        return ARGS[index + 1]
    end
    return default
end

# Retrieve the main directory from command-line arguments
main_dir = get_arg("--main_dir", "/default/path/to/main_dir")
event_thre = parse(Float64, get_arg("--event_thre", "0.8"))

# Average each type of correlation
avg_cross_corr_E_E = average_correlations(main_dir, "cross_corr_E_E.jld2")
avg_cross_corr_I_I = average_correlations(main_dir, "cross_corr_I_I.jld2")
avg_cross_corr_E_I = average_correlations(main_dir, "cross_corr_E_I.jld2")
avg_cross_corr_I_E = average_correlations(main_dir, "cross_corr_I_E.jld2")
avg_cross_corr_C_C = average_correlations(main_dir, "cross_corr_C_C.jld2")

# Average each type of event-based correlation
avg_cross_corr_E_E_event = average_correlations(main_dir, "cross_corr_E_E_event.jld2")
avg_cross_corr_I_I_event = average_correlations(main_dir, "cross_corr_I_I_event.jld2")
avg_cross_corr_E_I_event = average_correlations(main_dir, "cross_corr_E_I_event.jld2")
avg_cross_corr_I_E_event = average_correlations(main_dir, "cross_corr_I_E_event.jld2")
avg_cross_corr_C_C_event = average_correlations(main_dir, "cross_corr_C_C_event.jld2")

# Average each type of nonevent-based correlation
avg_cross_corr_E_E_nonevent = average_correlations(main_dir, "cross_corr_E_E_nonevent.jld2")
avg_cross_corr_I_I_nonevent = average_correlations(main_dir, "cross_corr_I_I_nonevent.jld2")
avg_cross_corr_E_I_nonevent = average_correlations(main_dir, "cross_corr_E_I_nonevent.jld2")
avg_cross_corr_I_E_nonevent = average_correlations(main_dir, "cross_corr_I_E_nonevent.jld2")
avg_cross_corr_C_C_nonevent = average_correlations(main_dir, "cross_corr_C_C_nonevent.jld2")

# Plot the averaged correlations for overall
plot_correlations(avg_cross_corr_E_E, avg_cross_corr_I_I, avg_cross_corr_E_I, avg_cross_corr_I_E, avg_cross_corr_C_C, joinpath(main_dir, "plot_corr_overall.png"),event_thre)

# Plot the averaged correlations for event-based
plot_correlations(avg_cross_corr_E_E_event, avg_cross_corr_I_I_event, avg_cross_corr_E_I_event, avg_cross_corr_I_E_event, avg_cross_corr_C_C_event, joinpath(main_dir, "plot_corr_events.png"),event_thre)

# Plot the averaged correlations for nonevent-based
plot_correlations(avg_cross_corr_E_E_nonevent, avg_cross_corr_I_I_nonevent, avg_cross_corr_E_I_nonevent, avg_cross_corr_I_E_nonevent, avg_cross_corr_C_C_nonevent, joinpath(main_dir, "plot_corr_nonevents.png"),event_thre)



if all([haskey(avg_cross_corr_E_E_event, 0), haskey(avg_cross_corr_C_C_event, 0), haskey(avg_cross_corr_I_I_event, 0),
    haskey(avg_cross_corr_E_E_nonevent, 0), haskey(avg_cross_corr_C_C_nonevent, 0), haskey(avg_cross_corr_I_I_nonevent, 0)])
    value_E_E_event = avg_cross_corr_E_E_event[0]
    value_C_C_event = avg_cross_corr_C_C_event[0]
    value_I_I_event = avg_cross_corr_I_I_event[0]

    value_E_E_nonevent = avg_cross_corr_E_E_nonevent[0]
    value_C_C_nonevent = avg_cross_corr_C_C_nonevent[0]
    value_I_I_nonevent = avg_cross_corr_I_I_nonevent[0]
    categories = ["E_E", "C_C", "I_I"]
    event_data = [value_E_E_event, value_C_C_event, value_I_I_event]
    nonevent_data = [value_E_E_nonevent, value_C_C_nonevent, value_I_I_nonevent]

    # Create the plot
    p = plot(xticks = (1:3, categories))
    scatter!(p, 1:3, event_data, label = "Event", color = :blue)
    scatter!(p, 1:3, nonevent_data, label = "Non-event", color = :red)
    plot!(p, 1:3, event_data, label = "", color = :blue, line = :line)
    plot!(p, 1:3, nonevent_data, label = "", color = :red, line = :line)

    # Save the plot to a file
    savefig(p, joinpath(main_dir, "event_vs_nonevent.png"))

end



function load_and_merge_data(main_dir::String)
    e_rates_all = Float64[]
    i_rates_all = Float64[]

    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            # Check for the existence of cross_corr_E_E.jld2 in the subdirectory
            cross_corr_file_path = joinpath(sub_dir, "cross_corr_E_E.jld2")

            e_rate_after_peak_file_path = joinpath(sub_dir, "e_rate_after_peak.jld2")
            
            #cross_corr_file_path = joinpath(sub_dir, "directory_name.txt")
            if isfile( e_rate_after_peak_file_path)
                e_file_path = joinpath(sub_dir, "e_rate_after_peak.jld2")
                i_file_path = joinpath(sub_dir, "i_rate_after_peak.jld2")

                if isfile(e_file_path) && isfile(i_file_path)
                    e_data = load(e_file_path, "e_rate_after_peak")
                    i_data = load(i_file_path, "i_rate_after_peak")

                    #disregard the entire e_data and i_data is there is rate higher than 30 in e_data
                    #println("e_data",e_data)
                    if e_data == []
                        println("e_data is empty")
                        continue
                    end

                    if maximum(e_data) > 30
                        #println("e_data",e_data)
                        #println("skiping")
                        continue
                    end

                    
                    #print(sub_dir)
                    #println(e_data)

                    append!(e_rates_all, e_data)
                    append!(i_rates_all, i_data)
                    
                end
            end
        end
    end

    return e_rates_all, i_rates_all
end



e_rates_merged, i_rates_merged = load_and_merge_data(main_dir)
#println("e_rates_merged",e_rates_merged)
@save joinpath(main_dir, "merged_e_rates.jld2") e_rates_merged
@save joinpath(main_dir, "merged_i_rates.jld2") i_rates_merged

function plot_rates_with_stats(e_rates, i_rates, output_file)
    e_rates_filtered = filter!(x -> !isnan(x), e_rates)
    i_rates_filtered = filter!(x -> !isnan(x), i_rates)
    
    e_mean = mean(e_rates_filtered)
    i_mean = mean(i_rates_filtered)
    
    e_std = std(e_rates_filtered)
    i_std = std(i_rates_filtered)
    
    
    plot_size = (600, 400)  # Adjust the size as needed
    left_margin = 10mm  # Increase left margin to ensure y-axis labels are visible
    mean_marker_size = 12  # Increase the size of the marker for the mean
    errorbar_width = 2  # Increase error bar width for better visibility

    p = scatter(ones(length(e_rates)), e_rates, label="E data", color=:pink, alpha=0.5, size=plot_size, left_margin=left_margin, xlim=(0.5,2.5))
    scatter!(p, 2*ones(length(i_rates)), i_rates, label="I data", color=:lightblue, alpha=0.5)
    #println("E Rates Mean is")
    #print(skipmissing(e_rates))
    #println(e_mean)
    #println("I Rates Mean is")
    # print(i_mean)
    # Overlay the mean markers
    scatter!(p, [1], [e_mean], yerr=e_std, label="E mean", color=:red, markershape=:cross, markersize=mean_marker_size, markeralpha=1, linealpha=1, linewidth=errorbar_width)
    scatter!(p, [2], [i_mean], yerr=i_std, label="I mean", color=:blue, markershape=:cross, markersize=mean_marker_size, markeralpha=1, linealpha=1, linewidth=errorbar_width)
    
    xticks!(p, [1, 2], ["E", "I"])
    ylabel!(p, "Rate")
    xlabel!(p, "Neuron Type")
    savefig(p, output_file)
end



plot_rates_with_stats(e_rates_merged, i_rates_merged, joinpath(main_dir, "rates_plot.png"))



file_e = "e_rate_raw_after_peak.jld2"
file_i = "i_rate_raw_after_peak.jld2"
max_trial_length = 85 # Define buffer_before and win_buff as needed
time_step = 5 # Define your time step size in ms

# ... [previous parts of the script] ...

function should_skip_directory(sub_dir::String, file_name::String, value_threshold::Float64)
    file_path = joinpath(sub_dir, file_name)
    if isfile(file_path)
        raw_activity_data = load(file_path)
        key = file_name == "e_rate_raw_after_peak.jld2" ? "e_rate_raw_after_peak" : "i_rate_raw_after_peak"
        raw_activity = raw_activity_data[key]
        # Check if any sublist in raw_activity contains a value greater than the threshold
        return any(sublist -> any(value -> value > value_threshold, sublist), raw_activity)
    end
    return false
end




function find_skippable_directories(main_dir::String, file_name::String, value_threshold::Float64)
    skippable_directories = Set{String}()
    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir) && should_skip_directory(sub_dir, file_name, value_threshold)
            push!(skippable_directories, sub_dir)
        end
    end
    return skippable_directories
end



function average_raw_activity_and_std(main_dir::String, file_name::String, max_length::Int64)

    skippable_directories = find_skippable_directories(main_dir, "e_rate_raw_after_peak.jld2", 30.0)

    total_activity = Vector{Float64}()
    squared_activity = Vector{Float64}()
    count = 0

    
    p = plot(size=(800, 600), xlabel="Time (ms)", ylabel="Rate", title="Raw Neural Activity")
    save_list =[]
    for sub_dir in readdir(main_dir, join=true) 
        if isdir(sub_dir) && !(sub_dir in skippable_directories)
            file_path = joinpath(sub_dir, file_name)
            if isfile(file_path)
                cross_corr_file_path = joinpath(sub_dir, "cross_corr_E_E.jld2")
                
                #cross_corr_file_path = joinpath(sub_dir, "guagua.txt")
                e_rate_after_peak_file_path = joinpath(sub_dir, "e_rate_after_peak.jld2")
            
            #cross_corr_file_path = joinpath(sub_dir, "directory_name.txt")
            if true #isfile(e_rate_after_peak_file_path)

                raw_activity_data = load(file_path)
                key = file_name == "e_rate_raw_after_peak.jld2" ? "e_rate_raw_after_peak" : "i_rate_raw_after_peak"
                raw_activity = raw_activity_data[key]
                
                # Filter the activities in raw_activity, discard if length > max_length
                filtered_activity = filter(a -> 12001 <= length(a) <= 12001, raw_activity)
                 
  
                #skip the current iteration of the loop iterating through subdirectories if any sublist in filtered_activity contains a value greater than 30. 
                if any(sublist -> any(value -> value > 30, sublist), filtered_activity)
                    #println("skiping")
                    continue # Skip this iteration if the condition is met
                end

                #plot all the activities on a same plot and save it
                
                for activity in filtered_activity



 
                    if length(total_activity) < length(activity)
                        
                        println("!!!!")
                        print(length(total_activity))
                        print(length(activity))

                        resize!(total_activity, length(activity))
                        
                        total_activity[(end - length(activity) + 1):end] .= 0.0
                        resize!(squared_activity, length(activity))
                        squared_activity[(end - length(activity) + 1):end] .= 0.0
                    end

                    #plot to p
                    plot!(p, 0:5:(length(activity) - 1) * 5, activity, label="Trial $count", legend=false)
                    push!(save_list,activity)

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
        
    end

    if file_name == "e_rate_raw_after_peak.jld2"
        println(save_list)
        savefig(p, joinpath(main_dir, "e_raw_activity_plot.png"))

        
        csv_file_path=joinpath(main_dir, "all_raw_activities.csv")

        open(csv_file_path, "w") do file
            for list in save_list
                # Join the elements of the list into a comma-separated string and write to the file
                write(file, join(list, ",") * "\n")
            end
        end
        

        
    end

    
    

    # Calculate mean and standard deviation
    mean_activity = count > 0 ? total_activity ./ count : total_activity
    variance_activity = count > 0 ? (squared_activity ./ count) .- mean_activity.^2 : squared_activity
    std_activity = sqrt.(variance_activity)

    return mean_activity, std_activity
end


# ... [rest of the script] ...



# Function to plot average raw activities with error ribbons
function plot_avg_raw_activity(time_step::Int, avg_e_activity, std_e_activity, avg_i_activity, std_i_activity, output_file::String)
    println("plotting")
    #println("avg_e_activity",avg_e_activity)
    time_axis = 0:time_step:(length(avg_e_activity) - 1) * time_step
    
    left_margin = 20mm
    p = plot(size=(800, 400), left_margin=left_margin, ylim=(0, 15), bottom_margin=5mm)  # Adjust bottom margin as needed

    ribbon_alpha = 0.1
    # Plot I rates with error ribbon
    plot!(p, time_axis, avg_i_activity, ribbon = std_i_activity, label="I Rates", color=:blue, fillalpha=ribbon_alpha)
    plot!(p, time_axis, avg_e_activity, ribbon = std_e_activity, label="E Rates", color=:red, fillalpha=ribbon_alpha)
    # Increase the size of x and y tick labels
    plot!(p, xtickfontsize=20, ytickfontsize=20)
    
    # Set the labels with increased size and add the title
    xlabel!(p, "Time (ms)", xguidefontsize=20)  # Increase label size for x-axis
    ylabel!(p, "Average Rate", yguidefontsize=20)  # Increase label size for y-axis
    title!(p, "Average Neural Activity Over Time")

    
    savefig(p, output_file)
end

# Main script execution



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
#zzz=zeros(size(std_i_rate_raw_after_peak))
#plot_avg_raw_activity(time_step, avg_e_rate_raw_after_peak, zzz, avg_i_rate_raw_after_peak, zzz, plot_file)



function compute_and_plot_psd(main_dir::String, file_name::String, max_length::Int64, output_file::String, output_file_2::String)
    all_psds = []  # To store the PSD of each trial
    count = 0

    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            file_path = joinpath(sub_dir, file_name)
            if isfile(file_path)
                cross_corr_file_path = joinpath(sub_dir, "cross_corr_E_E.jld2")
                
                #cross_corr_file_path = joinpath(sub_dir, "guagua.txt")
                e_rate_after_peak_file_path = joinpath(sub_dir, "e_rate_after_peak.jld2")
            
            #cross_corr_file_path = joinpath(sub_dir, "directory_name.txt")
            if isfile(e_rate_after_peak_file_path)
                raw_activity_data = load(file_path)
                key = file_name == "e_rate_raw_after_peak.jld2" ? "e_rate_raw_after_peak" : "i_rate_raw_after_peak"
                raw_activity = raw_activity_data[key]
                for i in 1:length(raw_activity)
                    #println(length(raw_activity[i]))
                end



                filtered_activity = filter(a -> 12001 <= length(a) <= 12001, raw_activity)

                if any(sublist -> any(value -> value > 30, sublist), filtered_activity)
                    #println("skiping")
                    continue # Skip this iteration if the condition is met
                end

                

                for activity in filtered_activity
                    # Calculate Power Spectral Density (PSD)
                    println(activity)
                    activity_r = activity[1:end] #remove the first 20 values
                    psd = periodogram(activity_r,fs=200)
                    #print(psd)
                    

                    psd_power_db = 10*log10.(psd.power[1:20])
                    psd_power_db_normalized = (psd_power_db .- minimum(psd_power_db)) ./ (maximum(psd_power_db) - minimum(psd_power_db))
                    push!(all_psds, (freq=psd.freq[1:20], power_db=psd_power_db_normalized))
                    count += 1
                end
            end
            end
        end
    end

    # Plotting the PSDs
    # Plotting the normalized PSDs
    p = plot(size=(800, 600), xlabel="Frequency (Hz)", ylabel="Normalized Power (0-1)", title="Normalized Power Spectral Density")
    for psd in all_psds
        plot!(p, psd.freq, psd.power_db, label="Trial $count", legend=false)
    end
    savefig(p, output_file)


    num_freqs = length(first(all_psds).freq)
    mean_psd = zeros(num_freqs)
    std_psd = zeros(num_freqs)

    for i in 1:num_freqs
        freq_psds = [psd.power_db[i] for psd in all_psds]
        mean_psd[i] = mean(freq_psds)
        std_psd[i] = std(freq_psds)
    end

    # Plotting the mean PSD with standard deviation as a ribbon
    p2 = plot(size=(800, 600), xlabel="Frequency (Hz)", ylabel="Mean Normalized Power (0-1)", title="Mean Normalized Power Spectral Density with STD",fillalpha=0.2,color="red")
    plot!(p2, first(all_psds).freq, mean_psd, ribbon=std_psd, label="Mean PSD", legend=false)
    savefig(p2, output_file_2)
    #save the data for plotting for the furture use
    @save joinpath(main_dir, "all_psds.jld2") all_psds
    @save joinpath(main_dir, "mean_psd.jld2") mean_psd
    @save joinpath(main_dir, "std_psd.jld2") std_psd




end


#compute_and_plot_psd(main_dir, file_e, max_trial_length,joinpath(main_dir, "power_spectra_e.png")) 
#compute_and_plot_psd(main_dir, file_i, max_trial_length,joinpath(main_dir, "power_spectra_i.png"))
println("computing psd for e")
#println("file_e",file_e)
compute_and_plot_psd(main_dir, file_e, max_trial_length,joinpath(main_dir, "power_spectra_e.png"),joinpath(main_dir, "power_spectra_e_mean.png")) 

compute_and_plot_psd(main_dir, file_i, max_trial_length,joinpath(main_dir, "power_spectra_i.png"),joinpath(main_dir, "power_spectra_i_mean.png"))



function convert_to_2d_matrix(sorted_psd_matrix)
    if isempty(sorted_psd_matrix)
        return Matrix{Float64}(undef, 0, 0)  # Return an empty matrix if the input is empty
    end

    # Determine the length of individual PSD vectors
    psd_length = length(first(sorted_psd_matrix))

    # Initialize a matrix to store the concatenated PSD data
    psd_2d_matrix = Matrix{Float64}(undef, length(sorted_psd_matrix), psd_length)

    for (i, psd_vector) in enumerate(sorted_psd_matrix)
        psd_2d_matrix[i, :] = psd_vector
    end

    return psd_2d_matrix
end



function compute_and_plot_psd_heatmap(main_dir::String, file_name::String, output_file::String, fs::Float64)
    psd_matrix = []  # To store the PSD (in dB) of each trial
    total_energy_per_trial = []  # To store the total energy of each trial for sorting

    for sub_dir in readdir(main_dir, join=true)
        if isdir(sub_dir)
            file_path = joinpath(sub_dir, file_name)
            if isfile(file_path)
                cross_corr_file_path = joinpath(sub_dir, "cross_corr_E_E.jld2")
                
                #cross_corr_file_path = joinpath(sub_dir, "guagua.txt")
            
            #cross_corr_file_path = joinpath(sub_dir, "directory_name.txt")
            if true #isfile(cross_corr_file_path)
                raw_activity_data = load(file_path)
                key = file_name == "e_rate_raw_after_peak.jld2" ? "e_rate_raw_after_peak" : "i_rate_raw_after_peak"
                raw_activity = raw_activity_data[key]

                filtered_activity = filter(a -> 12001 <= length(a) <= 12001, raw_activity)

                if any(sublist -> any(value -> value > 30, sublist), filtered_activity)
                    #println("skiping")
                    continue # Skip this iteration if the condition is met
                end
                
                for activity in filtered_activity
                    # Calculate Power Spectral Density (PSD)
                    #println(activity)
                    activity_r = activity[1:end] #remove the first 20 values
                    psd = periodogram(activity_r, fs=200)
                    #psd_db = 10 * log10.(psd.power)
                    #push!(psd_matrix, psd_db)
                    
                    psd_power_db = 10*log10.(psd.power[1:20])
                    psd_power_db_normalized = (psd_power_db .- minimum(psd_power_db)) ./ (maximum(psd_power_db) - minimum(psd_power_db))
                    push!(psd_matrix, psd_power_db_normalized)



                    # Calculate total energy (sum of PSD values)
                    total_energy = sum(psd_power_db)
                    push!(total_energy_per_trial, total_energy)
                end
             end

            end

        end
    end

    # Sort trials by total energy (higher energy at the bottom)
    sorted_indices = sortperm(total_energy_per_trial, rev=true)
    sorted_psd_matrix = psd_matrix[sorted_indices]
    #println(size(sorted_psd_matrix))
    #println(sorted_psd_matrix)
    psd_2d_matrix = convert_to_2d_matrix(sorted_psd_matrix)
    
    # Create the heatmap
    # Is there a way to change the x axis?
    p = heatmap(psd_2d_matrix, color=:viridis, xlabel="Frequency (Hz)", ylabel="Trials", title="PSD Heatmap")
    savefig(p, output_file)
end


compute_and_plot_psd_heatmap(main_dir, file_e,joinpath(main_dir, "power_spectra_heat_e.png"),200.0) 
compute_and_plot_psd_heatmap(main_dir, file_i,joinpath(main_dir, "power_spectra_heat_i.png"),200.0)
