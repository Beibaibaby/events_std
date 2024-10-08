#Copyright (C) 2024 Draco(Yunlong) Xu
#see README for more information
using Statistics
using ProgressMeter
using Plots

ENV["GKSwstype"] = "100"  # Use the off-screen GKS terminal

include("./LIF_1.0.0_compare.jl")

function sim_dynamic(Ne,Ni,T,taue,taui,pei,pie,pii,pee,K,stimstr_para,Nstim,jie_para,
    jei_para,jii_para,jee_para, stim_duration,stim_start_time,ie_sign,ee_sign,corr_flag,
    add_noise,sigma_noise,scale_noise,d_ee,f_ee,d_ie,f_ie,
    stim_duration_2,stim_start_2,stimstr_2,c_noise,peak_ratio,large_peak_mean,
    use_init_weights=false, init_weights=nothing)
    println("Setting up parameters")
    #corr_flag=false
    # Network parameters

    Ncells = Ne+Ni
    sqrtK = sqrt(K)
    if add_noise

       print("noise level is ")
       println(sigma_noise)

       print("noise lamba is ")
       println(scale_noise)
    end

    # Synaptic weights
    jie = jie_para / (taui * sqrtK)
    jei = jei_para / (taue * sqrtK)
    jii = jii_para / (taui * sqrtK)
    jee = jee_para / (taue * sqrtK)
   

    # Stimulation
    stimstr = stimstr_para/taue
    stimstr_2 = stimstr_2/taue
    print("stim is ")
    println(stimstr)
    print("stim 2 is ")
    println(stimstr_2)

    stimstart = stim_start_time
    stimend = stim_start_time+stim_duration

    stim_end_2 = stim_duration_2 + stim_start_2
    
    print("stim starting time is ")
    println(stim_start_time)

    print("stim ending time is ")
    println(stimend)

    print("stim 2 starting time is ")
    println(stim_start_2)

    print("stim 2 ending time is ")
    println(stim_end_2)


    #d = 0.15       # depression fraction upon a spike
    #f = 0.92  (orignal)      # facilitation increment upon a spike
    #f = 0.5
    #tau_d = 103.0     # time constant for D to recover to 1 (ms)
    #tau_f = 96.0     # time constant for F to recover to 1 (ms)
    
    tau_d = 1030 # time constant for D to recover to 1 (step;not ms)
    tau_f = 960 # time constant for D to recover to 1 (step;not ms)

    #todo: might need adjust the time scale

    # Constants and thresholds
    muemin = 1.1
    muemax = 1.2
    muimin = 1
    muimax = 1.05


    vre = 0.0
    threshe = 1
    threshi = 1

    dt = 0.1
    refrac = 5

    tauerise = 1
    taurise = 1
    tauedecay = 3
    tauidecay = 2

    maxrate = 500  # Maximum average firing rate (Hz) was 100 orignally

    # Initialize parameters
    mu = zeros(Ncells)
    thresh = zeros(Ncells)
    tau = zeros(Ncells)

    if use_init_weights && init_weights !== nothing

        mu=init_weights[2]


    else

        mu[1:Ne] .= (muemax - muemin) .* rand(Ne) .+ muemin
        mu[(Ne + 1):Ncells] .= (muimax - muimin) .* rand(Ni) .+ muimin

    end
    ##########
    #mu = zeros(Ncells)
    ##########
    thresh[1:Ne] .= threshe
    thresh[(Ne + 1):Ncells] .= threshi

    tau[1:Ne] .= taue
    tau[(Ne + 1):Ncells] .= taui

    weights = zeros(Ncells, Ncells)
    weights_D_ee = ones(Ne) #the sending D (E->E)
    weights_F_ee = ones(Ne) #the sending F (E->E)

    weights_D_ie = ones(Ne) #the sending F (E->I)
    weights_F_ie = ones(Ne) #the sending F (E->I)

    # Here we only need one decay/facilitation factor for one given neuron i, the factors from i to j are all the same

    # Initialize or load weights
    if use_init_weights && init_weights !== nothing
        weights = init_weights[1]
    else
        weights = zeros(Ncells, Ncells)
        weights[1:Ne, 1:Ne] .= jee .* (rand(Ne, Ne) .< pee)
        weights[1:Ne, (1 + Ne):Ncells] .= jei .* (rand(Ne, Ni) .< pei)
        weights[(1 + Ne):Ncells, 1:Ne] .= jie .* (rand(Ni, Ne) .< pie)
        weights[(1 + Ne):Ncells, (1 + Ne):Ncells] .= jii .* (rand(Ni, Ni) .< pii)
        
        for ci = 1:Ncells
            weights[ci, ci] = 0
        end
    end

    
    weights_initial = [copy(weights),copy(mu)] #save the weights_initial

    maxTimes = round(Int, maxrate * T / 1000)
    times = zeros(Ncells, maxTimes)
    ns = zeros(Int, Ncells)
    Nsteps = round(Int, T / dt)


    weights_D_ee_track = zeros(Nsteps)  
    weights_F_ee_track = zeros(Nsteps)

    weights_D_ie_track = zeros(Nsteps)  
    weights_F_ie_track = zeros(Nsteps)

    weights_IE_mean_history = zeros(Nsteps)
    weights_EE_mean_history = zeros(Nsteps)

    v_history = zeros(Ncells, Nsteps)  # Nsteps because we're saving at each time step, not just spikes
    E_input=zeros(Ncells, Nsteps) 
    I_input=zeros(Ncells, Nsteps)
    E_input_raw=zeros(Ncells, Nsteps) 
    I_input_raw=zeros(Ncells, Nsteps)
    
    forwardInputsE = zeros(Ncells)
    forwardInputsI = zeros(Ncells)

    forwardInputsEPrev = zeros(Ncells)
    forwardInputsIPrev = zeros(Ncells)

    xerise = zeros(Ncells)
    xedecay = zeros(Ncells)
    xirise = zeros(Ncells)
    xidecay = zeros(Ncells)

    v = rand(Ncells)
    #v[1:3] .=0
    #v = vre*ones(Ncells)
    lastSpike = -100 * ones(Ncells)
    
    
    println("Starting simulation")
    pl = Progress(Nsteps, 5)
    
    corr_pairs=100
    weights_copy = copy(weights)
    record_kick = Int[]
    
    gaussian_noise_global=0
    large_noise_flag=false

    top_n_e_neurons = load_data_from_jld2("./hyperparameter/top_n_e_neurons.jld2", "top_n_e_neurons")
    println(top_n_e_neurons)
    large_noise_last = 0
    for ti = 1:Nsteps
        t = dt * ti
        forwardInputsE[:] .= 0
        forwardInputsI[:] .= 0
        
        weights_D_ee .+= (1 .- weights_D_ee) ./ tau_d
        weights_F_ee .+= (1 .- weights_F_ee) ./ tau_f

        weights_D_ie .+= (1 .- weights_D_ie) ./ tau_d
        weights_F_ie .+= (1 .- weights_F_ie) ./ tau_f

        weights_D_ee_track[ti] = mean(weights_D_ee)
        weights_F_ee_track[ti] = mean(weights_F_ee)
        weights_D_ie_track[ti] = mean(weights_D_ie)
        weights_F_ie_track[ti] = mean(weights_F_ie)
        
        if ie_sign == true
           weights[(1 + Ne):Ncells, 1:Ne] = weights_copy[(1 + Ne):Ncells, 1:Ne] .* (weights_D_ie .* weights_F_ie)'
        end 

        if ee_sign == true
            weights[1:Ne, 1:Ne] = weights_copy[1:Ne, 1:Ne] .* (weights_D_ee .* weights_F_ee)'
        end 

        W_sub_view_ie = @view weights[(1 + Ne):Ncells, 1:Ne]
        weights_IE_mean_history[ti] = mean(W_sub_view_ie)
        W_sub_view_ee = @view weights[1:Ne, 1:Ne]
        weights_EE_mean_history[ti] = mean(W_sub_view_ee)

        if add_noise
            
           #gaussian_noise_global = sqrt(c_noise) * scale_noise * randn() * sigma_noise * dt

           large_peak_mean = 5.4/taue
           small_variance = sigma_noise*0.01
           large_variance = sigma_noise*0.005
           #peak_ratio = 2000
           
           #longer noise or not
           
           gaussian_noise_global, large_noise_flag = bimodal_gaussian_noise(c_noise, scale_noise, sigma_noise, dt, large_peak_mean, small_variance, large_variance, peak_ratio)
           
            
            if large_noise_flag
              push!(record_kick, ti)
              large_noise_last = gaussian_noise_global
              println("noise")
              println("large_noise_last is",large_noise_last)
            else
                if record_kick != []
                    if abs(ti - record_kick[end]) < 800 #800 
                                #println("ti is",ti)
                                #println("record_kick is",record_kick[end])
                                
                                gaussian_noise_global = large_noise_last
                                #println("large_noise_last is",large_noise_last)
                        
                    else
                       gaussian_noise_global=0
                    end

                end
            end
        end

        for ci = 1:Ncells
            
            xerise[ci] += -dt * xerise[ci] / tauerise + forwardInputsEPrev[ci]
            xedecay[ci] += -dt * xedecay[ci] / tauedecay + forwardInputsEPrev[ci]
            xirise[ci] += -dt * xirise[ci] / taurise + forwardInputsIPrev[ci]
            xidecay[ci] += -dt * xidecay[ci] / tauidecay + forwardInputsIPrev[ci]
            
            

            #synInput = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            synInput_E = (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise)
            synInput_I = (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            E_input[ci, ti] = synInput_E #* 0.1
            I_input[ci, ti] = synInput_I #* 0.1
            E_input_raw[ci, ti] = forwardInputsEPrev[ci]
            I_input_raw[ci, ti] = forwardInputsIPrev[ci]
           

            synInput = synInput_E+synInput_I # (xedecay[ci] - xerise[ci]) / (tauedecay - tauerise) + (xidecay[ci] - xirise[ci]) / (tauidecay - taurise)
            

            #if (ci < Nstim) && (t > stimstart) && (t < stimend)
            #    synInput += stimstr
            #end

            #if (ci in top_n_e_neurons) && (t > stimstart) && (t < stimend)
            #   synInput += stimstr
            #end 

            #Here we stimulate all the E neurons

            if (ci < Nstim) && (t > stimstart) && (t < stimend)
                synInput += stimstr
             end 
            

            if (ci < Nstim) && (t > stim_start_2) && (t < stim_end_2)
                synInput += stimstr_2
            end

            gaussian_noise_local = sqrt(1-c_noise) * randn() * sigma_noise * dt
            
            ##########
            synInput += gaussian_noise_local
            ###########

            if add_noise



               if ci in top_n_e_neurons
                synInput += gaussian_noise_global
               end 

            end
           
            if t > (lastSpike[ci] + refrac)
                v[ci] += dt * ((1 / tau[ci]) * (mu[ci] - v[ci]) + synInput)     
                v_history[ci, ti] = v[ci]
                if v[ci] > thresh[ci]
                    if corr_flag == true
                        if (ci>=corr_pairs && ci<500 )||(ci>=600 && ci<=1000)  || ci>=1100
                           v[ci] = vre
                        end
                    else
                        v[ci] = vre
                    end

                    lastSpike[ci] = t
                    ns[ci] += 1
                    if ns[ci] <= maxTimes
                        times[ci, ns[ci]] = t
                    end
                    
                    #test if remove
                   if ci < Ne
                    weights_D_ee[ci] = d_ee * weights_D_ee[ci]
                    weights_F_ee[ci] = f_ee + weights_F_ee[ci]
                    weights_D_ie[ci] = d_ie * weights_D_ie[ci]
                    weights_F_ie[ci] = f_ie + weights_F_ie[ci]
                   end                             

                    for j = 1:Ncells
                        if weights[j, ci] > 0  # E synapse
                            forwardInputsE[j] += weights[j, ci]
                        elseif weights[j, ci] < 0  # I synapse
                            forwardInputsI[j] += weights[j, ci]
                        end
                    end
                end
            else
                v_history[ci, ti] = v[ci]#
            end
        end

        forwardInputsEPrev .= forwardInputsE
        forwardInputsIPrev .= forwardInputsI
        next!(pl)
    end

    print("\r")
    println(maximum(ns))
    println(maximum(ns))
    if maximum(ns)<=maxTimes
         times = times[:, 1:maximum(ns)]
    else
       println("no over max")
    end
    
    return times, ns, Ne, Ncells, T, v_history, E_input, I_input, weights, weights_D_ee_track, weights_F_ee_track , weights_IE_mean_history, weights_EE_mean_history, weights_D_ie_track, weights_F_ie_track, record_kick, weights_initial, E_input_raw, I_input_raw
end

function load_data_from_jld2(file_path, data_key)
    # Open the .jld2 file and read the data associated with the given key
    data = jldopen(file_path, "r") do file
        return file[data_key]
    end

    return data
end

function bimodal_gaussian_noise(c_noise, scale_noise, sigma_noise, dt, large_peak_mean, small_variance, large_variance, peak_ratio)
    # Define two normal distributions
    small_peak = Normal(0, small_variance)
    large_peak = Normal(large_peak_mean, large_variance)

    # Determine which distribution to sample from
    if rand() < peak_ratio / (peak_ratio + 1)
        noise = rand(small_peak) * 0.001
        noise = sqrt(c_noise) * scale_noise * noise * dt
        noise = 0
        large_noise_flag = false  # Flag for large noise not generated
        
    else
        noise = rand(large_peak)
        large_noise_flag = true  # Flag for large noise generated
    end



    return noise, large_noise_flag
end


function compute_correlation(E_input::Matrix, I_input::Matrix, num_pairs::Int=100)
    # Ensure input matrices have the same dimensions
    if size(E_input) != size(I_input)
        error("E_input and I_input must have the same dimensions.")
    end

    Ncells, Nsteps = size(E_input)

    # Ensure the number of pairs is less than or equal to Ncells
    if num_pairs > Ncells
        error("The number of pairs cannot exceed the number of cells.")
    end

    # Randomly select indices for both E_input and I_input
    random_e_indices = rand(1:Ncells, num_pairs)
    random_i_indices = rand(1:Ncells, num_pairs)

    total_correlation = 0.0

    for i = 1:num_pairs
        e_data = E_input[random_e_indices[i], :]
        i_data = I_input[random_i_indices[i], :]

        # Compute Pearson correlation and add to the total
        total_correlation += cor(e_data, i_data)
    end

    # Compute average correlation
    avg_correlation = total_correlation / num_pairs

    return avg_correlation
end

function compute_cross_correlation_old(E_input::Matrix, I_input::Matrix, tau_range::UnitRange{Int}=-30:30, num_pairs::Int=50)
    # Ensure input matrices have the same dimensions
    if size(E_input) != size(I_input)
        error("E_input and I_input must have the same dimensions.")
    end

    Ncells, Nsteps = size(E_input)
    
    num_pairs=round(Int,Ncells/1.1)


    # Ensure the number of pairs is less than or equal to Ncells
    if num_pairs > Ncells * Ncells
        error("The number of pairs cannot exceed the square of the number of cells.")
    end
    

    # Randomly select indices for both E_input and I_input
    random_e_indices = rand(1:Ncells, num_pairs)
    random_i_indices = rand(1:Ncells, num_pairs)

    avg_correlations = Dict{Int, Float64}()

    for tau in tau_range
        if abs(tau) >= Nsteps
            avg_correlations[tau] = 0
        else
            total_correlation = 0.0
            for i = 1:num_pairs
                
                e_data = E_input[random_e_indices[i], :]
                i_data = I_input[random_i_indices[i], :]
                
                # Shift data according to tau
                if tau > 0
                    e_data = e_data[1:end-tau]
                    i_data = i_data[tau+1:end]
                elseif tau < 0
                    e_data = e_data[-tau+1:end]
                    i_data = i_data[1:end+tau]
                end
                
                # Compute Pearson correlation for this tau and add to the total
                total_correlation += cor(e_data, i_data)
            end

            # Compute average correlation for this tau
            avg_correlations[tau] = total_correlation / num_pairs
        end
    end

    return avg_correlations
end


function compute_cross_correlation(E_input::Matrix, I_input::Matrix, tau_range::UnitRange{Int}=-30:30, num_pairs::Int=50)
    # Ensure input matrices have the same dimensions
    if size(E_input) != size(I_input)
        error("E_input and I_input must have the same dimensions.")
    end

    Ncells, Nsteps = size(E_input)
    
    num_pairs = round(Int, Ncells / 1.5)

    # Ensure the number of pairs is less than or equal to Ncells
    if num_pairs > Ncells
        error("The number of pairs cannot exceed the number of cells.")
    end
    
    avg_correlations = Dict{Int, Float64}()

    for tau in tau_range
        if abs(tau) >= Nsteps
            avg_correlations[tau] = 0
        else
            total_correlation = 0.0
            for i = 1:num_pairs
                # Randomly select distinct indices for E_input and I_input
                e_index = rand(1:Ncells)
                i_index = rand(setdiff(1:Ncells, e_index))

                e_data = E_input[e_index, :]
                i_data = I_input[i_index, :]
                
                # Shift data according to tau
                if tau > 0
                    e_data = e_data[1:end-tau]
                    i_data = i_data[tau+1:end]
                elseif tau < 0
                    e_data = e_data[-tau+1:end]
                    i_data = i_data[1:end+tau]
                end
                
                # Compute Pearson correlation for this tau and add to the total
                total_correlation += cor(e_data, i_data)
            end

            # Compute average correlation for this tau
            avg_correlations[tau] = total_correlation / num_pairs
        end
    end

    return avg_correlations
end

    
#tune_parameters()



function compute_sliding_rate(spiketimes, window_size, step_size, T)
    println("starting computing firing rate")
    n_neurons, _ = size(spiketimes)
    n_steps = floor(Int, (T - window_size) / step_size) + 1
    rates = zeros(n_steps)

    p = Progress(n_steps, 2) # Initialize the progress bar

    for i = 1:n_steps
        next!(p) # Update the progress bar

        t_start = (i-1) * step_size + 1  # Ensure the start time is non-zero
        t_end = t_start + window_size - 1  # Adjust end time based on start
        
        # Check for out-of-bounds
        if t_end > T
            break
        end

        n_spikes = sum([sum((t_start .<= spiketimes[j, :]) .& (spiketimes[j, :] .< t_end)) for j=1:n_neurons])
        rates[i] = n_spikes / (n_neurons * window_size) * 1000  # rate in Hz
    end
    #println(rates)
    return rates
end

function remove_spikes_near_kicks(times, record_kick, neighborhood_size)
    # Deep copy to preserve the original times array
    times_copy = deepcopy(times)

    # First loop: Iterate over each neuron and each kick time
    for neuron_index in 1:size(times, 1)
        for kick_time in record_kick
            adjusted_kick_time = kick_time / 10

            # Find indices of spikes to be removed
            spikes_indices = findall(x -> abs(x - adjusted_kick_time) <= neighborhood_size && x != 0, times[neuron_index, :])

            # Set the identified spikes to zero
            for index in spikes_indices
                times_copy[neuron_index, index] = 0
            end
        end
    end

    # Second loop: Shift non-zero elements to the left for each neuron
    for neuron_index in 1:size(times, 1)
        non_zero_elements = filter(x -> x != 0, times_copy[neuron_index, :])
        num_non_zero = length(non_zero_elements)
        times_copy[neuron_index, 1:num_non_zero] = non_zero_elements
        times_copy[neuron_index, (num_non_zero+1):end] .= 0
    end

    # Return the modified times array
    return times_copy
end




function plot_correlations(cross_corr_E_E, cross_corr_I_I, cross_corr_E_I, cross_corr_I_E, cross_corr_C_C, output_path, event_thre )
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]
    sorted_values_C_C = [cross_corr_C_C[k] for k in sorted_keys]
    plot_margin = 25mm
    p_size = (1200, 800)

    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright,size=p_size,left_margin=plot_margin)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    plot!(sorted_keys, sorted_values_C_C, label="C-C", linewidth=2, marker=:circle, color=:green)

    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation | Event threshold = $event_thre")
    savefig(output_path)
end


function plot_correlations_mem(cross_corr_E_E, cross_corr_I_I, cross_corr_T_T, cross_corr_E_I, cross_corr_I_E, output_file)
    # Assuming the dictionaries have the same keys
    sorted_keys = sort(collect(keys(cross_corr_E_E)))
    
    sorted_values_E_E = [cross_corr_E_E[k] for k in sorted_keys]
    sorted_values_I_I = [cross_corr_I_I[k] for k in sorted_keys]
    sorted_values_T_T = [cross_corr_T_T[k] for k in sorted_keys]
    sorted_values_E_I = [cross_corr_E_I[k] for k in sorted_keys]
    sorted_values_I_E = [cross_corr_I_E[k] for k in sorted_keys]

    # Create the plot
    plot(sorted_keys, sorted_values_E_E, label="E-E", linewidth=2, marker=:circle, color=:blue, legend=:topright)
    plot!(sorted_keys, sorted_values_I_I, label="I-I", linewidth=2, marker=:circle, color=:red)
    plot!(sorted_keys, sorted_values_E_I, label="E-I", linewidth=2, marker=:circle, color=:yellow)
    plot!(sorted_keys, sorted_values_I_E, label="I-E", linewidth=2, marker=:circle, color=:orange)
    plot!(sorted_keys, sorted_values_T_T, label="T-T", linewidth=2, marker=:circle, color=:green)

    # Set labels and title
    xlabel!("Tau")
    ylabel!("Correlation")
    title!("Cross-correlation")

    # Save the figure to the specified output file
    savefig(output_file)
end


function plot_cells(v_history, cells)
    p = plot(layout=(3,1), size=(600,800))

    for (index, cell) in enumerate(cells)
        plot!(p[index], v_history[cell, :], label="Cell $cell", xlabel="Time", ylabel="Voltage", legend=:topright)
    end
    savefig("v_history.png")
end

function plot_curve(d_ee, f_ee, d_ie, f_ie, output_filename)
    A = 20.0
    tau_d = 103.0
    tau_f = 96.0
    dt = 1.0
    T = 3000.0
    first_spike_time = 50.0
    temporal_frequencies = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 20, 25, 30, 40, 50, 60, 75, 80, 100, 110, 130, 150]

    depression_ratios = []
    depression_ratios_fixed_F = []

    for tf in temporal_frequencies
        S_input = generate_spike_train(T, dt, first_spike_time, tf)

        # Simulate with variable F
        _, _, _, Hs = simulate_LIF_neuron(A, d_ee, f_ee, tau_d, tau_f, dt, T, S_input)
        if length(Hs) >= 2
            push!(depression_ratios, Hs[2] / Hs[1])
        end

        # Simulate with fixed F
        _, _, _, Hs_fixed_F = simulate_LIF_neuron_fixed_F(A, d_ie, f_ie, tau_d, tau_f, dt, T, S_input)
        if length(Hs_fixed_F) >= 2
            push!(depression_ratios_fixed_F, Hs_fixed_F[2] / Hs_fixed_F[1])
        end
    end

    # Plot and save results
    p = plot(temporal_frequencies, depression_ratios, xlabel="Temporal Frequency (Hz)", title="d_ee=$d_ee f_ee=$f_ee d_ie=$d_ie f_ie=$f_ie", ylabel="Depression Ratio", label="E->E", marker=:circle, markercolor=:red, linewidth=2, linecolor=:red)
    plot!(p, temporal_frequencies, depression_ratios_fixed_F, label="E->I", marker=:circle, linewidth=2, linecolor=:deepskyblue2, markercolor=:deepskyblue2)
    savefig(p, output_filename)
end

function compute_average_activity_post_kick(times, record_kick, buffer_time,T)
    # Sort kick times to handle them in order
    sort!(record_kick)
    n_neurons = size(times, 1)
    # Create and merge intervals
    intervals = []
    current_interval = (record_kick[1]+300, min(record_kick[1] + buffer_time*10, T*10))
    for kick_time in record_kick[2:end]
        if kick_time <= current_interval[2]
            # Extend the current interval if there is an overlap
            current_interval = (current_interval[1]+300, min(kick_time + buffer_time*10, T*10))
        else
            # Save the current interval and start a new one
            push!(intervals, current_interval)
            current_interval = (kick_time, min(kick_time + buffer_time*10, T*10))
        end
    end
    push!(intervals, current_interval)  # Add the last interval
    println(intervals)

    # Compute average activity for each interval
    average_activities = []
    for interval in intervals
        total_spikes = 0
        interval_duration = interval[2]/10 - interval[1]/10

        n_spikes = sum([sum(((interval[1]/10) .<= times[j, :]) .& (times[j, :] .< (interval[2]/10))) for j=1:n_neurons])


        #for neuron_index in 1:size(times, 1)
            # Count spikes in the current interval
        #    total_spikes += count(x -> (interval[1]/10) <= x <= (interval[2]/10), times[neuron_index, :])
        #end
        # Calculate the average activity
        # Converting interval duration from ms to seconds for Hz calculation
        average_activity = n_spikes / (n_neurons * interval_duration) * 1000
        #(total_spikes / interval_duration) * 1000
        push!(average_activities, average_activity)
        print(average_activity)
    end

    return average_activities
end


function compute_average_activity_post_kick_2(rate_cons, record_kick, buffer_time, step_size, T)
    sort!(record_kick)
    println(record_kick)
    println("step_size",step_size)

    # Convert kick times to indices in the rate arrays, adjusting for the time scale
    kick_indices = [round(Int, (kick_time / 10 ) / step_size) for kick_time in record_kick]
    buffer_steps = round(Int, buffer_time / step_size)
    println("kick_indices",kick_indices)
    println("buffer_steps",buffer_steps)

    # Create and merge intervals
    intervals = []
    current_interval = (kick_indices[1], min(kick_indices[1] + buffer_steps, round(Int, T / step_size)))
    println("current_interval",current_interval)
    for kick_index in kick_indices[2:end]
        print("kick_index",kick_index)
        if kick_index <= current_interval[2]
            # Extend the current interval if there is an overlap
            current_interval = (current_interval[1], min(kick_index + buffer_steps, round(Int, T / step_size)))
        else
            # Save the current interval and start a new one
            current_interval = (kick_index, min(kick_index + buffer_steps, round(Int, T / step_size)))
            push!(intervals, current_interval)
            
        end
        println("current_interval",current_interval)

    end
    push!(intervals, current_interval)  # Add the last interval

    # Compute average activity for each interval
    average_activities = []
    #exculde the interval less than 5 steps
    intervals = filter(x -> x[2]-x[1]>49, intervals)
    println("intervals",intervals)

    for interval in intervals
        println("interval",interval)
        interval_rates = rate_cons[interval[1]:min(interval[2], length(rate_cons))]
        println("interval_rates",interval_rates)
        average_activity = mean(interval_rates)
        push!(average_activities, average_activity)
    end

    return average_activities
end


function collect_raw_activity_clips(rate_cons, record_kick, buffer_before, buffer_after, step_size, T)
    sort!(record_kick)

    # Convert kick times to indices in the rate arrays
    kick_indices = [round(Int, (kick_time / 10 ) / step_size) for kick_time in record_kick]
    buffer_before_steps = round(Int, buffer_before / step_size)
    buffer_after_steps = round(Int, buffer_after / step_size)

    # Create and merge intervals
    clips = []
    current_clip = (max(1, kick_indices[1] - buffer_before_steps), min(kick_indices[1] + buffer_after_steps, round(Int, T / step_size)))
    for kick_index in kick_indices[2:end]
        if kick_index - buffer_before_steps <= current_clip[2]
            # Extend the current interval if there is an overlap
            current_clip = (current_clip[1], min(kick_index + buffer_after_steps, round(Int, T / step_size)))
        else
            # Save the current interval and start a new one
            clip = rate_cons[current_clip[1]:min(current_clip[2], length(rate_cons))]
            push!(clips, clip)
            current_clip = (max(1, kick_index - buffer_before_steps), min(kick_index + buffer_after_steps, round(Int, T / step_size)))
        end
    end
    # Add the last interval
    clip = rate_cons[current_clip[1]:min(current_clip[2], length(rate_cons))]
    push!(clips, clip)

    return clips
end














