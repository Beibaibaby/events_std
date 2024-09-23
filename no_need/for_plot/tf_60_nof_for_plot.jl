ENV["GKSwstype"] = "100"  # Use the off-screen GKS terminal to aviod running messages
using Plots
using Measures

function simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input)
    τ_m = 10.0       
    V_thresh = -50.0  
    V_reset = -75.0   
    V_rest = -75.0    
    R_m = 10.0        

    time = 0:dt:T      
     
    V = V_rest                
    Vs = Float64[]           
    spike_times = []         
    D = 1.0                  
    F = 1.0                  
    Hs = Float64[] 

    Ds = Float64[]
    Fs = Float64[] 

    prev_V = V 
    spike_count = 0  

    for (idx, t) in enumerate(time)
        # Store the D and F at each timestep
        push!(Ds, D)
        push!(Fs, F)

        W = A * D * F
        dV = (-(V - V_rest) + R_m * W * S_input[idx]) / τ_m
        V += dV * dt

        if S_input[idx] == 1
            spike_count += 1
            if spike_count == 1 || spike_count == 5
                # Record the change in potential due to the input spike
                push!(Hs, V - prev_V)
            end
            D *= d  
            F += f  
        end

        dD = (1 - D) / tau_d
        D += dD * dt
    
        dF = (1 - F) / tau_f
        F += dF * dt
        
        prev_V = V
        push!(Vs, V)
    end
    
    return time, Vs, spike_times, Hs, Ds, Fs  # Added Ds and Fs to the return values
end




using Plots

# Define the RGB color using Plots' RGB function
#orangex = RGB(242/255, 163/255, 58/255)
purplex = RGB(134/255, 99/255, 152/255)

function plot_LIF_simulation(time, Vs, spike_times, S_input, Ds, Fs, filename)
    # Plot the membrane potential
    p1 = plot(
        time[1:200], Vs[1:200], 
        lw=2, label="Membrane Potential", 
        ylabel="Voltage (mV)", ylims=(-80, -50), 
        legendfontsize=12, linecolor=purplex, 
        grid=false, dpi=400, left_margin=15mm, 
        linewidth=12, legend=:none, 
        tickfont=font(15), guidefont=font(15), 
        bottom_margin=5mm
    )

    # Plot the depression and facilitation factors
    p3 = plot(
        time[1:200], Ds[1:200], 
        lw=2, label="Depression Factor D", 
        xlabel="Time (ms)", ylabel="Value", 
        color=:black, legend=:topright, 
        legendfontsize=6, linewidth=12, 
        grid=false, dpi=400, left_margin=15mm, 
        ylim=(-0.05,5), legend_background_color=:transparent, 
        foreground_color_legend=nothing, 
        tickfont=font(15), guidefont=font(15), 
        bottom_margin=5mm
    )

    plot!(
        p3, time[1:200], Fs[1:200], 
        lw=2, label="Facilitation Factor F", 
        color=:darkgreen, legend=:topleft, 
        legendfontsize=12, linewidth=12, 
        legend_background_color=:transparent, 
        foreground_color_legend=nothing, 
        tickfont=font(15), guidefont=font(15), 
        bottom_margin=5mm
    )

    # Combine the plots
    p = plot(
        p1, p3, 
        layout=(2,1), link=:x, 
        size=(800, 400), dpi=400
    )

    # Save the figure
    savefig(p, filename)
end



function generate_spike_train(T, dt, initial_spike_time, tf)
    # T: total simulation time
    # dt: time step
    # initial_spike_time: time for the first spike
    # tf: temporal frequency, indicating how often a spike appears in the train
    
    # Calculate the number of time steps
    num_steps = convert(Int, T/dt) + 1
    
    # Initialize the spike train with all zeros
    S_input = zeros(Int, num_steps)
    
    # Set the initial spike
    S_input[convert(Int, initial_spike_time/dt)] = 1
    
    # Calculate the time interval between spikes based on the temporal frequency
    spike_interval = round(Int, 1000/tf)
    
    # Generate following spikes at evenly spaced intervals
    for i in 2:5  # since the first spike is already set, we generate the next 4 spikes
        next_spike_time = initial_spike_time + (i-1)*spike_interval
        if next_spike_time <= T  # only set spike if it is within the total time T
            S_input[convert(Int, next_spike_time/dt)] = 1
        end
    end
    
    return S_input
end

# Example usage
A = 20.0        # fixed parameter A
d = 0.24       # depression fraction upon a spike
f = 0.00         # facilitation increment upon a spike
tau_d = 103.0     # time constant for D to recover to 1 (ms)
tau_f = 96.0     # time constant for F to recover to 1 (ms)
dt = 1.0        # time step (ms)
T = 1000.0       # total time to simulate (ms)

# Generate spike train
first_spike_time = 50.0  # ms
temporal_frequency = 60.0  # Hz
S_input = generate_spike_train(T, dt, first_spike_time, temporal_frequency)

# Run simulation
time, Vs, spike_times, Hs, Ds, Fs = simulate_LIF_neuron(A, d, f, tau_d, tau_f, dt, T, S_input) 

# Calculate depression ratio H_5/H_1
depression_ratio = Hs[2] / Hs[1]
println("Depression ratio H_5/H_1: ", depression_ratio)


# Plot results and save to file
output_filename = "../plots_for_paper/tf=60_noF.png"
plot_LIF_simulation(time, Vs, spike_times, S_input, Ds, Fs, output_filename)  # Added Ds and Fs
output_filename = "../plots_for_paper/tf=60_noF.svg"
plot_LIF_simulation(time, Vs, spike_times, S_input, Ds, Fs, output_filename)  # Added Ds and Fs