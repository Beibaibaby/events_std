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

include("sim_7_t20.jl") # include the network simulation code

struct NetworkParameters #struct for the network parameters
    Ncells::Int
    Ne::Int
    Ni::Int
    T::Int
    taue::Int
    taui::Int
    pei::Float64
    pie::Float64
    pii::Float64
    pee::Float64
    K::Int
    jie_para::Float64
    jei_para::Float64
    jii_para::Float64
    jee_para::Float64
    Nstim::Int
    stimstr_para::Float64
    stim_duration::Int
    stim_start_time::Int
    ie_sign::Bool
    ee_sign::Bool
    corr_flag::Bool
    add_noise::Bool
    sigma_noise::Float64
    scale_noise::Float64
    d_ee::Float64
    f_ee::Float64
    d_ie::Float64
    f_ie::Float64
    stim_duration_2::Int
    stim_start_2::Int
    stimstr_2::Float64
    c_noise::Float64
    peak_ratio::Float64
    large_peak_mean::Float64
end

# Define a function to retrieve a value from ARGS or return a default value if not present.
function get_arg(key, default)
    index = findfirst(==(key), ARGS)
    if index !== nothing && index < length(ARGS)
        return ARGS[index + 1]
    end
    return default
end


function run_experiment(;
    Ncells,
    Ne,
    Ni,
    T,
    taue,
    taui,
    pei,
    pie,
    pii,
    pee,
    K,
    jie,
    jei,
    jii,
    jee,
    Nstim,
    stimstr,
    stim_duration,
    stim_start_time,
    ie_sign,
    ee_sign,
    corr_flag,
    low_plot,
    add_noise,
    sigma_noise,
    scale_noise,
    env,
    d_ee,
    f_ee,
    d_ie,
    f_ie,
    stim_duration_2,
    stim_start_2,
    stimstr_2,
    c_noise,
    dir_name_in,
    corr_sign,
    event_thre,
    peak_ratio,
    large_peak_mean,
    use_init_weights, 
    init_weights
)
        
        doplot = true
        do_v_corr = false
        do_repeat_v_corr=false

        #Setting the parameters
        ############################
        # Now, use the provided values to create an instance of the struct:
        
        params = NetworkParameters(Ncells, Ne, Ni, T, taue, taui, pei, pie, pii, pee, K, jie, jei, jii, jee, Nstim, stimstr,stim_duration, stim_start_time,ie_sign,ee_sign,corr_flag, add_noise,
        sigma_noise,scale_noise,d_ee,f_ee,d_ie,f_ie,stim_duration_2,stim_start_2,stimstr_2,c_noise,peak_ratio,large_peak_mean)
         
        timestamp = Dates.now()
        timestamp_str = Dates.format(timestamp, "yyyy-mm-dd_HH-MM-SS")

        # These codes tell us which machine we are using 
        if env == 1 
            dir_name ="/root/autodl-tmp/d_ee=$d_ee+f_ie=$f_ie+d_ie=$d_ie+$timestamp_str"
            if !isdir(dir_name_in)
                mkpath(dir_name_in)
            end
        elseif env == 2
           #dir_name = "/gpfs/data/doiron-lab/draco/results/d_ee=$d_ee+f_ie=$f_ie+d_ie=$d_ie+$timestamp_str"
           if !isdir(dir_name_in)
            mkpath(dir_name_in)
        end
        else
            dir_name=dir_name_in
            println("nice")
            if !isdir(dir_name_in)
                mkpath(dir_name_in)
            end
        end

        file_path_test = joinpath(dir_name, "directory_name.txt")


        # Open the file in write mode, write the dir_name, and close the file
        open(file_path_test, "w") do file
            write(file, dir_name)
        end

        #store it
        #run the stimulus
        #times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights = sim_old()
        global times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, weights_D_ee_track, weights_F_ee_track,weights_IE_mean_history,weights_EE_mean_history
        times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, weights_D_ee_track, weights_F_ee_track,weights_IE_mean_history,weights_EE_mean_history,weights_D_ie_track, weights_F_ie_track, record_kick, weights_initial=sim_dynamic(
            params.Ne,params.Ni,params.T,params.taue,params.taui,params.pei,params.pie,params.pii,params.pee,params.K,
            params.stimstr_para,params.Nstim,params.jie_para,params.jei_para,params.jii_para,params.jee_para,
            params.stim_duration,params.stim_start_time,params.ie_sign,params.ee_sign,params.corr_flag,
            params.add_noise, params.sigma_noise, params.scale_noise, params.d_ee,params.f_ee,params.d_ie,params.f_ie,
            params.stim_duration_2,params.stim_start_2,params.stimstr_2,params.c_noise,peak_ratio,large_peak_mean,
            use_init_weights, init_weights)
            
        # Here, we check whether we are using existing initial weights. If not, we save the initial weights.
        if !use_init_weights  # Save only if new weights are generated
            print("no using initials")
            @save "/gpfs/data/doiron-lab/draco/weights_$timestamp_str.jld2" weights_initial
        end

        overall_mean_rate_E =  mean(1000 * ns[1:params.Ne] / params.T)
        overall_mean_rate_I = mean(1000 * ns[(params.Ne+1):Ncells] / params.T)

        println("mean excitatory firing rate: ", overall_mean_rate_E, " Hz")
        println("mean inhibitory firing rate: ", overall_mean_rate_I, " Hz")
        
        product_weights_ee = weights_D_ee_track .* weights_F_ee_track


        if doplot 

            # Generate a timestamp for unique filenames
                println("creating plot")
           
                # Define plot size: (width, height)
                plot_size = (2500, 1100) 
                plot_margin = 50mm

                # Parameters for sliding window
                window_size = 25  # in ms
                step_size = 5     # in ms

                if add_noise
                    
                    neighborhood_size = 2
                    
                    times_modified = remove_spikes_near_kicks(times, record_kick, neighborhood_size)
                    println("kick time is")
                    println(record_kick)
                    println("times")
                    println(times[1:50,:])
                    println("times motified")
                    println(times_modified[1:50,:])
                    
                    e_rate_cons = compute_sliding_rate(times_modified[1:params.Ne, :], window_size, step_size, params.T)
                    i_rate_cons = compute_sliding_rate(times_modified[(params.Ne+1):Ncells, :], window_size, step_size, params.T)

                    
                    win_buff = 400
                    buffer_before=10
                    #e_rate_after_peak = compute_average_activity_post_kick(times_modified[1:params.Ne, :], record_kick, win_buff,T)
                    #i_rate_after_peak = compute_average_activity_post_kick(times_modified[(params.Ne+1):Ncells, :], record_kick, win_buff,T)
                    
                    # if kick_record is non-empty, then we compute the average activity post-kick
                    # condition check for non-empty kick_record
                    if record_kick != []
                        e_rate_after_peak = compute_average_activity_post_kick_2(e_rate_cons, record_kick, win_buff, step_size, params.T)
                        i_rate_after_peak = compute_average_activity_post_kick_2(i_rate_cons, record_kick, win_buff, step_size, params.T)
                        
                        e_rate_raw_after_peak=collect_raw_activity_clips(e_rate_cons, record_kick, buffer_before, win_buff, step_size, T)
                        i_rate_raw_after_peak=collect_raw_activity_clips(i_rate_cons, record_kick, buffer_before, win_buff, step_size, T)

                        println("events average for E")
                        println(e_rate_after_peak)
                        println("events average for I")
                        println(i_rate_after_peak)

                        @save joinpath(dir_name, "e_rate_after_peak.jld2") e_rate_after_peak
                        @save joinpath(dir_name, "i_rate_after_peak.jld2") i_rate_after_peak
                        @save joinpath(dir_name, "e_rate_raw_after_peak.jld2") e_rate_raw_after_peak
                        @save joinpath(dir_name, "i_rate_raw_after_peak.jld2") i_rate_raw_after_peak
                

                        plot_ef = plot(e_rate_after_peak, title = "Excitatory Rate After Peak", xlabel = "Index", ylabel = "Rate", legend = false)
                        plot_if = plot(i_rate_after_peak, title = "Inhibitory Rate After Peak", xlabel = "Index", ylabel = "Rate", legend = false)
                        savefig(plot_ef, joinpath(dir_name, "e_rate_after_peak_plot.png"))
                        savefig(plot_if, joinpath(dir_name, "i_rate_after_peak_plot.png"))
                    end
                end 

                e_rate = compute_sliding_rate(times[1:params.Ne, :], window_size, step_size, params.T)
                i_rate = compute_sliding_rate(times[(params.Ne+1):Ncells, :], window_size, step_size, params.T)

                
                @save joinpath(dir_name, "e_rate.jld2") e_rate
                @save joinpath(dir_name, "i_rate.jld2") i_rate

                #compute avg_E_rate, which is the average rate of each e neuron
                T_seconds = T / 1000

                avg_E_rate = ns[1:Ne] / T_seconds

                #compute avg_I_rate, which is the average rate of each i neuron
                avg_I_rate = ns[(Ne+1):Ncells] / T_seconds



                N_boss = 50 # The boss in the population

                top_n_e_neurons = sortperm(avg_E_rate, rev=true)[1:N_boss]
            

                # Get indices of top N inhibitory neurons
                top_n_i_neurons = sortperm(avg_I_rate, rev=true)[1:N_boss]

                @save joinpath(dir_name, "top_n_e_neurons_noise.jld2") top_n_e_neurons
                @save joinpath(dir_name, "top_n_i_neurons_noise.jld2") top_n_i_neurons

                # Sort the average rates in descending order
                sorted_avg_E_rate = sort(avg_E_rate, rev=true)
                sorted_avg_I_rate = sort(avg_I_rate, rev=true)

                
            
                bar(1:length(sorted_avg_E_rate), sorted_avg_E_rate, 
                title="Sorted Average Rate of Excitatory Neurons", 
                xlabel="Neuron Index (Sorted)", ylabel="Average Rate (Hz)", 
                legend=false, color=:red,lw=0)
                savefig(joinpath(dir_name, "sorted_average_rate_excitatory_neurons_bar.png"))

                
                bar(1:length(sorted_avg_I_rate), sorted_avg_I_rate, 
                title="Sorted Average Rate of Inhibitory Neurons", 
                xlabel="Neuron Index (Sorted)", ylabel="Average Rate (Hz)", 
                legend=false, color=:blue,lw=0)

                savefig(joinpath(dir_name, "sorted_average_rate_inhibitory_neurons_bar.png"))
                if corr_sign
                    if any(x -> x > 150, e_rate)
                        println("Element in e_rate greater than 150 found. Exiting script.")

                        open(file_path_test, "w") do file
                            #write "Element in e_rate greater than 150 found. Exiting script."
                            write(file, "Element in e_rate greater than 150 found. Exiting script.")
                        end

                        exit()
                    end
                end

                # Compute the time values based on window_size and step_size
                n_steps = length(e_rate)  # or length(i_rate), assuming they have the same length
                #time_values = [i * step_size + (window_size / 2) for i in 1:n_steps]
                time_values = [i * step_size + window_size  for i in 1:n_steps]

                # Add a code to detect low rate or not 
                if add_noise
                    if low_plot ##this para control whether focus on the zoom in low activity
                        p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)", label="Excitatory", lw=2, linecolor=:red, size=plot_size, title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie Noise=$scale_noise sigma=$sigma_noise c_noise=$c_noise event_thre=$event_thre p_ratio=$peak_ratio large_p_mean=$large_peak_mean" , ylim=(0,5))
                        plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2,left_margin=plot_margin)
                    else
                        p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)", label="Excitatory", lw=5, linecolor=:red, size=plot_size, title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie Noise=$scale_noise sigma=$sigma_noise c_noise=$c_noise event_thre=$event_thre p_ratio=$peak_ratio large_p_mean=$large_peak_mean",xtickfont=font(12), ytickfont=font(12),legendfontsize=12,
                        titlefontsize=14, xlabelfontsize=12, ylabelfontsize=12)
                        plot!(time_values, i_rate, label="Inhibitory", lw=5, linecolor=:deepskyblue2,left_margin=plot_margin,bottom_margin=plot_margin)
                    end 

                    if low_plot # this parameter controls whether to focus on the zoom in low activity
                        p2_cons = plot(time_values, e_rate_cons, xlabel="Time (ms)", ylabel="Firing rate (Hz)", label="Excitatory", lw=2, linecolor=:red, size=plot_size, title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie Noise=$scale_noise sigma=$sigma_noise c_noise=$c_noise event_thre=$event_thre p_ratio=$peak_ratio large_p_mean=$large_peak_mean", ylim=(0,5))
                        plot!(time_values, i_rate_cons, label="Inhibitory", lw=2, linecolor=:deepskyblue2, left_margin=plot_margin)
                    else
                        p2_cons = plot(time_values, e_rate_cons, xlabel="Time (ms)", ylabel="Firing rate (Hz)", label="Excitatory", lw=2, linecolor=:red, size=plot_size, title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie Noise=$scale_noise sigma=$sigma_noise c_noise=$c_noise event_thre=$event_thre p_ratio=$peak_ratio large_p_mean=$large_peak_mean")
                        plot!(time_values, i_rate_cons, label="Inhibitory", lw=2, linecolor=:deepskyblue2, left_margin=plot_margin)
                    end
                    
                else
                    if low_plot ##this para control whether focus on the zoom in low activity
                        p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)", label="Excitatory", lw=2, linecolor=:red, size=plot_size, title="Firing rate(d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie stim=$stimstr event_thre=$event_thre)", ylim=(0,5))
                        plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2,left_margin=plot_margin)
                    else
                        p2 = plot(time_values, e_rate, xlabel="Time (ms)", ylabel="Firing rate (Hz)", label="Excitatory", lw=2, linecolor=:red, size=plot_size, title="Firing rate(d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie stim=$stimstr event_thre=$event_thre)")
                        plot!(time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2,left_margin=plot_margin)
                    end 

                end
        
                          

                plot_curve(d_ee, f_ee, d_ie, f_ie, "$dir_name/plot_curve_$timestamp_str.png")
                
                #@save "$dir_name/times.jld2" times
                #@save "$dir_name/ns.jld2" ns
                #@save "$dir_name/Ne.jld2" Ne
                #@save "$dir_name/Ncells.jld2" Ncells
                #@save "$dir_name/T.jld2" T
                #@save "$dir_name/v_history.jld2" v_history
                #@save "$dir_name/E_input.jld2" E_input
                #@save "$dir_name/I_input.jld2" I_input
                #@save "$dir_name/weights.jld2" weights
                #@save "$dir_name/weights_D_ee_track.jld2" weights_D_ee_track
                #@save "$dir_name/weights_F_ee_track.jld2" weights_F_ee_track
                #@save "$dir_name/weights_IE_mean_history.jld2" weights_IE_mean_history
                #@save "$dir_name/weights_EE_mean_history.jld2" weights_EE_mean_history
                


                fig_filename = "$dir_name/plot_FR_$timestamp_str.png"
                savefig(p2, fig_filename)
                
                if add_noise 
                   fig_filename_cons = "$dir_name/plot_FR_cons_$timestamp_str.png"
                   savefig(p2_cons, fig_filename_cons)
                end


                # Generate scaled x-values
                x_values = 0:0.1:(length(product_weights_ee) - 1) * 0.1

                # Create individual subplots
                # Adjust title font size and margins for the plots
                title_fontsize = 12  # You can adjust this value as needed
                left_margin = 10   # Adjust this value as needed

                # Create individual subplots with the modifications
                fig_width = 1800  # Width in pixels
                fig_height = 600  # Height in pixels
                title_fontsize = 14  # Title font size
                

                p3 = plot(
                    x_values, weights_EE_mean_history, 
                    label="EE Mean History", 
                    xlabel="Time (ms)", 
                    ylabel="Value", 
                    title="E->E strength d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # And also here in pixels
                )
                
                p4 = plot(
                    x_values, weights_IE_mean_history, 
                    label="IE Mean History", 
                    xlabel="Time (ms)", 
                    ylabel="Value", 
                    title="E->I strength d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # And also here in pixels
                )
                
                # Now combine the subplots
                p5 = plot(p3, p4, layout=(2,1), size=(fig_width, fig_height*2))

                # Save the plot as an image
                savefig(p5,"$dir_name/plot_syn_$timestamp_str.png") 

                println("Combined figure saved as $dir_name/plot_SYN_$timestamp_str.png") 

                p6 = plot(
                    x_values, product_weights_ee, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="Product Value", 
                    title="Product of D and F (EE) d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                p7 = plot(
                    x_values, weights_D_ee_track, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="D Value", 
                    title="D factor (EE) d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie ", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                p8 = plot(
                    x_values, weights_F_ee_track, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="F Value", 
                    title="F factor (EE) d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                combined_plot = plot(p6, p7, p8, layout = (3, 1), size = (fig_width*1.2, fig_height*3.2))

                savefig(combined_plot,"$dir_name/plot_DF_EE_$timestamp_str.png") 
                

                product_weights_ie = weights_D_ie_track .* weights_F_ie_track

                p9 = plot(
                    x_values, product_weights_ie, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="Product Value", 
                    title="Product of D and F (IE) d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                p10 = plot(
                    x_values, weights_D_ie_track, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="D Value", 
                    title="D factor (IE) d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie ", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                p11 = plot(
                    x_values, weights_F_ie_track, 
                    label="Product over time", 
                    xlabel="Time (ms)", 
                    ylabel="F Value", 
                    title="F factor (IE) d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie", 
                    titlefontsize=title_fontsize,
                    size=(fig_width, fig_height),
                    left_margin=plot_margin # Apply the margin adjustment here in pixels
                )

                combined_plot_2 = plot(p9, p10, p11, layout = (3, 1), size = (fig_width*1.2, fig_height*3.2))

                savefig(combined_plot_2,"$dir_name/plot_DF_IE_$timestamp_str.png") 
                

                println("Figure saved as $fig_filename")  


                json_filename = "$dir_name/$timestamp_str.json"



                param_dict = Dict(
                    "Ncells" => params.Ncells,
                    "Ne" => params.Ne,
                    "Ni" => params.Ni,
                    "T" => params.T,
                    "taue" => params.taue,
                    "taui" => params.taui,
                    "pei" => params.pei,
                    "pie" => params.pie,
                    "pii" => params.pii,
                    "pee" => params.pee,
                    "K" => params.K,
                    "jie_para" => params.jie_para,
                    "jee_para" => params.jee_para,
                    "jei_para" => params.jei_para,
                    "jii_para" => params.jii_para,
                    "Nstim" => params.Nstim,
                    "stimstr_para" => params.stimstr_para,
                    "stim_duration" => params.stim_duration,
                    "stim_start_time" => params.stim_start_time,
                    "ie_sign"=> params.ie_sign,
                    "ee_sign"=> params.ee_sign,
                    "corr_flag" =>  params.corr_flag,
                    "sigma_noise" => params.sigma_noise,
                    "add_noise" => params.add_noise,
                    "scale_noise" => params.scale_noise,
                    "d_ee" => params.d_ee,
                    "f_ee" => params.f_ee,
                    "d_ie" => params.d_ie,
                    "f_ie" => params.f_ie,
                    "stim_duration_2"=>params.stim_duration_2,
                    "stim_start_2"=>params.stim_start_2,
                    "stimstr_2"=>params.stimstr_2,
                    "c_noise"=>params.c_noise,
                    "peak_ratio"=>peak_ratio,
                    "large_peak_mean"=>large_peak_mean
                )

                # Now, you can access any of these values using the dictionary's keys, e.g., param_dict["Ne"] or param_dict["jie"].


                JSON3.open(json_filename, "w") do io
                JSON3.print(io, param_dict)
                end

                println("Parameters saved as $json_filename")
                println("done")
                #println(cross_EPSP_EPSP)
                #println("cross E-E: ", cross_corr_E_E)



               if corr_sign
                
                        #### Codes for searching events
                        # Assuming step_size in ms and time_values array is in place
                        buffer_before = round(Int, 150 / step_size)  # 100 ms before in steps, rounded to nearest integer
                        buffer_after = round(Int,  350/ step_size)   # 200 ms after in steps, rounded to nearest integer
 
                        # Calculate overall average for excitatory rates
                        overall_avg_e_rate = mean(e_rate)

                        # Initialize array to store excursion periods
                        excursions_e = []

                        # Identify excursions for excitatory rates
                        start_index_e = -1
                        for (index, rate) in enumerate(e_rate)
                            if rate > overall_avg_e_rate + event_thre
                                if start_index_e == -1
                                    start_index_e = max(1, index - buffer_before)
                                end
                            else
                                if start_index_e != -1
                                    end_index_e = min(length(e_rate), index + buffer_after)
                                    push!(excursions_e, (time_values[start_index_e], time_values[end_index_e]))
                                    start_index_e = -1  # Reset start index
                                end
                            end
                        end

                        # Check for an excursion at the end of the time series
                        if start_index_e != -1
                            end_index_e = min(length(e_rate), length(time_values))
                            push!(excursions_e, (time_values[start_index_e], time_values[end_index_e]))
                        end

                        sorted_excursions_e = sort(excursions_e, by = x -> x[1])  # Sort by start time
                        merged_excursions_e = []

                        for excursion in sorted_excursions_e
                            if isempty(merged_excursions_e) || merged_excursions_e[end][2] < excursion[1]
                                push!(merged_excursions_e, excursion)
                            else
                                # Merge overlapping intervals
                                merged_excursions_e[end] = (merged_excursions_e[end][1], max(merged_excursions_e[end][2], excursion[2]))
                            end
                        end
                        excursions_e=merged_excursions_e
                        # Excursions including buffers are now stored in excursions_e

                        # Assuming excursions_e contains the start and end times of event periods

                        ########################Plot firing rates#############################
 
                        p_event = plot(size=plot_size, left_margin=plot_margin)
                     
                        # Plot grey rectangles for event periods
                        for (start_time, end_time) in excursions_e
                            vspan!(p_event, [start_time, end_time], color=:grey, alpha=0.3,label=false)
                            open(file_path_test, "w") do file
                                write(file, "start time= $start_time\n")
                                write(file, "end time = $end_time\n")
                            end
                        end



                        # Add the excitatory and inhibitory rates to the plot
                        plot!(p_event, time_values, e_rate, label="Excitatory", lw=2, linecolor=:red, title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f=$f_ie sigma=$sigma_noise c_noise=$c_noise event_thre=$event_thre p_ratio=$peak_ratio large_p_mean=$large_peak_mean")
                        plot!(p_event, time_values, i_rate, label="Inhibitory", lw=2, linecolor=:deepskyblue2)

                        # Set plot labels and title
                        xlabel!(p_event, "Time (ms)")
                        ylabel!(p_event, "Firing rate (Hz)")

                        fig_filename = "$dir_name/plot_FR_event_$timestamp_str.png"

                        # Save the figure
                        savefig(p_event, fig_filename)




                        #########################ploting end###########################################
                        # Initialize arrays for event and nonevent segments
                        E_input_event = []
                        E_input_nonevent = []

                        E_input_event_test = []
                        E_input_nonevent_test = []
                       
                        start_nonevent_index = 1
                         
                        index_test = 0 
                        # Loop through each excursion period
                        for (start_time, end_time) in excursions_e
                            
                            start_index = start_time*10
                            end_index = min(end_time * 10, T * 10 - 1)

                            # Extract event segments
                            event_segment = E_input[:, start_index:end_index]
                            push!(E_input_event, event_segment)

                            open(file_path_test, "w") do file
                                write(file, "event start:= $start_index\n")
                                write(file, "event end:= $end_index\n")
                            end

                            # Extract nonevent segments
                            nonevent_segment = E_input[:, start_nonevent_index:start_index-1]
                            push!(E_input_nonevent, nonevent_segment)

                            open(file_path_test, "w") do file
                                write(file, "non_event start:= $start_nonevent_index\n")
                                write(file, "non_event end:= $start_index\n")
                            end

                            if index_test == 0
                                push!(E_input_event_test, event_segment)
                                push!(E_input_nonevent_test, nonevent_segment)
                            end

                            # Update the start index for the next nonevent segment
                            start_nonevent_index = end_index + 1
                            index_test = index_test + 1
                        end

                        # Handle the final nonevent segment after the last excursion
                        final_nonevent_segment = E_input[:, start_nonevent_index:end]
                        push!(E_input_nonevent, final_nonevent_segment)
                        

                        # Concatenate the segments
            ##################just for testing########################
                        #E_input_event = E_input_event_test
                        #E_input_nonevent = E_input_nonevent_test
           ########comment out when not testing##########################

                        E_input_event = hcat(E_input_event...)
                        E_input_nonevent = hcat(E_input_nonevent...)

                        # Initialize arrays for event and nonevent segments
                        I_input_event = []
                        I_input_nonevent = []

                        I_input_event_test = []
                        I_input_nonevent_test = []

                        # Initialize the starting index for nonevent segments
                        start_nonevent_index = 1

                        index_test_2 = 0

                        # Loop through each excursion period
                        for (start_time, end_time) in excursions_e
                            # Convert times to indices (assuming time_values is in sync with I_input)
                            start_index = start_time*10
                            end_index = min(end_time * 10, T * 10 - 1)

                            # Extract event segments
                            event_segment = I_input[:, start_index:end_index]
                            push!(I_input_event, event_segment)

                            open(file_path_test, "w") do file
                                write(file, "event start:= $start_index\n")
                                write(file, "event end:= $end_index\n")
                            end

                            # Extract nonevent segments
                            nonevent_segment = I_input[:, start_nonevent_index:start_index-1]
                            push!(I_input_nonevent, nonevent_segment)

                            open(file_path_test, "w") do file
                                write(file, "non_event start:= $start_nonevent_index\n")
                                write(file, "non_event end:= $start_index\n")
                            end

                            if index_test_2 == 0
                                push!(I_input_event_test, event_segment)
                                push!(I_input_nonevent_test, nonevent_segment)
                            end

                            # Update the start index for the next nonevent segment
                            start_nonevent_index = end_index + 1
                            index_test_2=index_test_2+1
                        end

                        # Handle the final nonevent segment after the last excursion
                        final_nonevent_segment = I_input[:, start_nonevent_index:end]
                        push!(I_input_nonevent, final_nonevent_segment)

            ##################just for testing########################
                        #I_input_event = I_input_event_test
                        #I_input_nonevent = I_input_nonevent_test
           ########comment out when not testing##########################

                        # Concatenate the segments
                        I_input_event = hcat(I_input_event...)
                        I_input_nonevent = hcat(I_input_nonevent...)
                        
                        C_input = E_input .+ I_input

                        event_length = size(E_input_event)
                        nonevent_length = size(E_input_nonevent)
                        overall_length = size(E_input)

                        open(file_path_test, "w") do file
                            write(file, "event_size = $event_length\n")
                            write(file, "nonevent_size = $nonevent_length\n")
                            write(file, "overall_size = $overall_length\n")
                        end

                        avg_correlation_E_I=compute_correlation(E_input, I_input)
                        avg_correlation_E_E=compute_correlation(E_input, E_input)
                        avg_correlation_I_I=compute_correlation(I_input, I_input)
                        avg_correlation_C_C=compute_correlation(C_input, C_input)

                        println("avg correlation(E-I): ", avg_correlation_E_I)
                        println("avg correlation(E-E): ", avg_correlation_E_E)
                        println("avg correlation(I-I): ", avg_correlation_I_I)
                        println("avg correlation(C-C): ", avg_correlation_C_C)

                        cross_corr_E_E=compute_cross_correlation(E_input, E_input)
                        cross_corr_I_I=compute_cross_correlation(I_input, I_input)
                        cross_corr_E_I=compute_cross_correlation(E_input, I_input)
                        cross_corr_I_E=compute_cross_correlation(I_input, E_input)
                        cross_corr_C_C=compute_cross_correlation(C_input, C_input)
                       
                        @save joinpath(dir_name_in, "cross_corr_E_E.jld2") cross_corr_E_E
                        @save joinpath(dir_name_in, "cross_corr_I_I.jld2") cross_corr_I_I
                        @save joinpath(dir_name_in, "cross_corr_E_I.jld2") cross_corr_E_I
                        @save joinpath(dir_name_in, "cross_corr_I_E.jld2") cross_corr_I_E
                        @save joinpath(dir_name_in, "cross_corr_C_C.jld2") cross_corr_C_C

                        plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I, cross_corr_I_E, cross_corr_C_C, joinpath(dir_name, "plot_corr_overall.png"),event_thre)
                            

                        # Create combined input for events
                        if isempty(E_input_event) || isempty(I_input_event)
                           println("skip Event")
                        else
                            C_input_event = E_input_event .+ I_input_event
                        
                            # Compute average correlations for event-based inputs
                            avg_correlation_E_I_event = compute_correlation(E_input_event, I_input_event)
                            avg_correlation_E_E_event = compute_correlation(E_input_event, E_input_event)
                            avg_correlation_I_I_event = compute_correlation(I_input_event, I_input_event)
                            avg_correlation_C_C_event = compute_correlation(C_input_event, C_input_event)

                            println("Event-based avg correlation (E-I): ", avg_correlation_E_I_event)
                            println("Event-based avg correlation (E-E): ", avg_correlation_E_E_event)
                            println("Event-based avg correlation (I-I): ", avg_correlation_I_I_event)
                            println("Event-based avg correlation (C-C): ", avg_correlation_C_C_event)

                            # Compute cross-correlations for event-based inputs
                            cross_corr_E_E_event = compute_cross_correlation(E_input_event, E_input_event)
                            cross_corr_I_I_event = compute_cross_correlation(I_input_event, I_input_event)
                            cross_corr_E_I_event = compute_cross_correlation(E_input_event, I_input_event)
                            cross_corr_I_E_event = compute_cross_correlation(I_input_event, E_input_event)
                            cross_corr_C_C_event = compute_cross_correlation(C_input_event, C_input_event)

                                                    # Saving event-based cross-correlations
                            @save joinpath(dir_name_in, "cross_corr_E_E_event.jld2") cross_corr_E_E_event
                            @save joinpath(dir_name_in, "cross_corr_I_I_event.jld2") cross_corr_I_I_event
                            @save joinpath(dir_name_in, "cross_corr_E_I_event.jld2") cross_corr_E_I_event
                            @save joinpath(dir_name_in, "cross_corr_I_E_event.jld2") cross_corr_I_E_event
                            @save joinpath(dir_name_in, "cross_corr_C_C_event.jld2") cross_corr_C_C_event

                            plot_correlations(cross_corr_E_E_event, cross_corr_I_I_event, cross_corr_E_I_event, cross_corr_I_E_event, cross_corr_C_C_event, joinpath(dir_name, "plot_corr_events.png"),event_thre)
                        end 
                        # Plot correlations for event-based inputs
                        
                        if isempty(E_input_nonevent) || isempty(I_input_nonevent)
                            println("skip non-Event")
                         else
                            # Create combined input for nonevents
                            C_input_nonevent = E_input_nonevent .+ I_input_nonevent

                            # Compute average correlations for nonevent-based inputs
                            avg_correlation_E_I_nonevent = compute_correlation(E_input_nonevent, I_input_nonevent)
                            avg_correlation_E_E_nonevent = compute_correlation(E_input_nonevent, E_input_nonevent)
                            avg_correlation_I_I_nonevent = compute_correlation(I_input_nonevent, I_input_nonevent)
                            avg_correlation_C_C_nonevent = compute_correlation(C_input_nonevent, C_input_nonevent)

                            println("Nonevent-based avg correlation (E-I): ", avg_correlation_E_I_nonevent)
                            println("Nonevent-based avg correlation (E-E): ", avg_correlation_E_E_nonevent)
                            println("Nonevent-based avg correlation (I-I): ", avg_correlation_I_I_nonevent)
                            println("Nonevent-based avg correlation (C-C): ", avg_correlation_C_C_nonevent)

                            # Compute cross-correlations for nonevent-based inputs
                            cross_corr_E_E_nonevent = compute_cross_correlation(E_input_nonevent, E_input_nonevent)
                            cross_corr_I_I_nonevent = compute_cross_correlation(I_input_nonevent, I_input_nonevent)
                            cross_corr_E_I_nonevent = compute_cross_correlation(E_input_nonevent, I_input_nonevent)
                            cross_corr_I_E_nonevent = compute_cross_correlation(I_input_nonevent, E_input_nonevent)
                            cross_corr_C_C_nonevent = compute_cross_correlation(C_input_nonevent, C_input_nonevent)
                
                            # Plot correlations for nonevent-based inputs


                            plot_correlations(cross_corr_E_E_nonevent, cross_corr_I_I_nonevent, cross_corr_E_I_nonevent, cross_corr_I_E_nonevent, cross_corr_C_C_nonevent, joinpath(dir_name, "plot_corr_nonevents.png"),event_thre)
                            

                            # Saving nonevent-based cross-correlations
                            @save joinpath(dir_name_in, "cross_corr_E_E_nonevent.jld2") cross_corr_E_E_nonevent
                            @save joinpath(dir_name_in, "cross_corr_I_I_nonevent.jld2") cross_corr_I_I_nonevent
                            @save joinpath(dir_name_in, "cross_corr_E_I_nonevent.jld2") cross_corr_E_I_nonevent
                            @save joinpath(dir_name_in, "cross_corr_I_E_nonevent.jld2") cross_corr_I_E_nonevent
                            @save joinpath(dir_name_in, "cross_corr_C_C_nonevent.jld2") cross_corr_C_C_nonevent

                         end
                end

        end



end



# Extract parameters from ARGS or use default values
Ncells = parse(Int, get_arg("--Ncells", "5000"))
Ne = parse(Int, get_arg("--Ne", "4000"))
Ni = parse(Int, get_arg("--Ni", "1000"))
T = parse(Int, get_arg("--T", "1000"))
taue = parse(Int, get_arg("--taue", "15"))
taui = parse(Int, get_arg("--taui", "10"))
pei = parse(Float64, get_arg("--pei", "0.5"))
pie = parse(Float64, get_arg("--pie", "0.5"))
pii = parse(Float64, get_arg("--pii", "0.5"))
pee = parse(Float64, get_arg("--pee", "0.2"))
K = parse(Int, get_arg("--K", "800"))
jie = parse(Float64, get_arg("--jie", "4.0"))
jei = parse(Float64, get_arg("--jei", string(-18.0 * 1.2)))
jii = parse(Float64, get_arg("--jii", "-16.0"))
jee = parse(Float64, get_arg("--jee", "10.0"))
Nstim = parse(Int, get_arg("--Nstim", "4000"))


ie_sign = parse(Bool, get_arg("--ie_sign", "true")) #controal E->I is dynamic or not 
ee_sign = parse(Bool, get_arg("--ee_sign", "true")) #controal E->E is dynamic or not 
corr_flag = parse(Bool, get_arg("--corr_flag", "false")) ##wether compute and plot EPSP and IPSP
low_plot = parse(Bool, get_arg("--low_plot", "false")) #contronl whether manully plot a low ativity regime

sigma_noise = parse(Float64, get_arg("--sigma_noise", "1.0"))
add_noise = parse(Bool, get_arg("--add_noise", "false"))
scale_noise = parse(Float64, get_arg("--scale_noise", "1"))
c_noise = parse(Float64, get_arg("--c_noise", "0.01"))

env = parse(Int, get_arg("--env", "3"))

d_ee = parse(Float64, get_arg("--d_ee", "0.24"))
f_ee = parse(Float64, get_arg("--f_ee", "0.85"))
d_ie = parse(Float64, get_arg("--d_ie", "0.24"))
f_ie = parse(Float64, get_arg("--f_ie", "0.0"))

stimstr = parse(Float64, get_arg("--stimstr", "0.0"))
stim_duration= parse(Int, get_arg("--stim_duration", "5"))
stim_start_time= parse(Int, get_arg("--stim_start_time", "200"))

stimstr_2 = parse(Float64, get_arg("--stimstr_2", "0.0"))
stimstr_2 = parse(Float64, get_arg("--stimstr_2", "0.0"))
stim_duration_2 = parse(Int, get_arg("--stim_duration_2 ", "200"))
stim_start_2 = parse(Int, get_arg("--stim_start_2", "400"))

event_thre = parse(Float64, get_arg("--event_thre", "0.8"))
peak_ratio = parse(Float64, get_arg("--peak_ratio", "10000"))
large_peak_mean = parse(Float64, get_arg("--large_peak_mean", "300"))

timestamp = Dates.now()
timestamp_str = Dates.format(timestamp, "yyyy-mm-dd_HH-MM-SS")

dir_name_in = get_arg("--dir_name_in", "/gpfs/data/doiron-lab/draco/results_20/d_ee=$d_ee+f_ie=$f_ie+d_ie=$d_ie+$timestamp_str")

corr_sign = parse(Bool, get_arg("--corr_sign", "true")) ##New sign for correlation

use_init_weights = parse(Bool, get_arg("--use_init_weights", "true"))
init_weights_name = get_arg("--weights_file", "/gpfs/data/doiron-lab/draco/weights_2024-03-20_16-03-16.jld2")  # Adjust the default as needed

# Initialize variable for initial weights
#init_weights = nothing
#init_weights = 1
jldopen(init_weights_name, "r") do file
    println(keys(file))  # This will print all the variables (keys) in the file
end

if use_init_weights
    # Load the weights from the specified file
    @load init_weights_name weights_initial
    init_weights = weights_initial
else
    init_weights = nothing
end


run_experiment(;Ncells,
    Ne,
    Ni,
    T,
    taue,
    taui,
    pei,
    pie,
    pii,
    pee,
    K,
    jie,
    jei,
    jii,
    jee,
    Nstim,
    stimstr,
    stim_duration,
    stim_start_time,
    ie_sign,
    ee_sign,
    corr_flag,
    low_plot,
    add_noise,
    sigma_noise,
    scale_noise,
    env,
    d_ee,
    f_ee,
    d_ie,
    f_ie,
    stim_duration_2,
    stim_start_2,
    stimstr_2,
    c_noise,
    dir_name_in,
    corr_sign,
    event_thre,
    peak_ratio,
    large_peak_mean,
    use_init_weights, 
    init_weights
)