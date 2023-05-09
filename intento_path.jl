using JSON
using Plots
# DEFINE THE STRUCTURE OF THE NEURAL NETWORK:
# TAKING INTO ACCOUNT THE NUMBER OF NEURONS, THEIR CONNECTIVITY,
# AND THEIR WEIGHTS (VALUE).

# Define connectivity between neurons
# read in the data from the file
syn_data = JSON.parsefile("chem.json")

nodes = syn_data["nodes"]
links = syn_data["links"]

connectivity = zeros(302, 302)
weights = zeros(302, 302)
for link in links
    from_index = link["source"] + 1
    to_index = link["target"] + 1
    value = link["value"]
    connectivity[from_index, to_index] = 1
    weights[from_index, to_index] = value
    
end

spy(connectivity)

# SET THE INITIAL STATE OF THE NEURONS:
# INITIAL MEMBRANE POTENTIALS AND CONDUCTANCES OF ION CHANNELS


# Set the initial membrane potentials of the neurons to -74 mV
# and the initial conductances *******
initial_membrane_potentials = -74 * connectivity
initial_conductances = 0 * connectivity   #PREGUNTAR





function kunert_method(t::Float64, dt::Float64, f::Function)
    t_mid = t + dt / 2
    k1 = f(t)
    k2 = f(t_mid)
    return dt * k2 + (dt^2 / 4) * (k2 - k1)
end

# Define the time step and simulation duration
dt = 0.01
tmax = 100
num_neurons = 302

# Create an array to store the membrane potentials of all neurons
V = zeros(Float64, num_neurons, Int(tmax / dt) + 1)

# Initialize the membrane potentials to their resting values
for i in 1:num_neurons
    V[i,1] = V_rest[i]
end

# Define a function that computes the time derivative of the membrane potential
function dV_dt(V::Vector{Float64}, t::Float64)
    # Compute the total synaptic input to the neuron
    I_syn = sum(I_synapses(V, gap_junctions, i) for i in 1:n_neurons if i != n)

    # Compute the time derivative of the membrane potential
    dVdt = (-g_leak * (V - V_rest) .- I_syn) / Cm
    return dVdt
end

# Use the Kunert method to update the membrane potentials of all neurons over time
for j in 2:Int(tmax / dt) + 1
    t = (j - 1) * dt
    for i in 1:num_neurons
        V[i,j] = V[i,j-1] + kunert_method(t, dt, (V_i, t) -> dV_dt(V_i, t))
    end
end



# Create an array of colors for each neuron
colors = [:blue, :green, :red, :orange, :purple, :brown, :pink, :gray, :black, :cyan]

# Create a plot with subplots for each neuron
plot_array = plot(layout=(num_neurons,1), size=(800,800))
for i in 1:num_neurons
    plot!(plot_array[i], V[i,:], color=:blue, xlabel="Time (ms)", ylabel="Membrane Potential (mV)")
end

# Display the plot
display(plot_array)


"""


# Define the dynamics of the ion channels and the neurons
function ion_channel_dynamics(conductance, V, E)
    return conductance * (E - V)
end

function neuron_dynamics(V, I_syn)
    C = 1 # membrane capacitance (nF)
    g_L = 0.1 # leak conductance (nS)
    E_L = -70 # reversal potential for the leak current (mV)
    I_L = g_L * (E_L - V) # leak current (nA)
    I_m = I_L + I_syn # total membrane current (nA)
    dV_dt = (1 / C) * I_m # membrane potential dynamics (mV/ms)
    return V + dt * dV_dt # update the membrane potential using the Euler method
end

function synaptic_current(g, Vpre, Vpost, E_syn)
    g_syn = sum(g)
    I_syn = g_syn * (E_syn - Vpre)
    return I_syn
end


# Define the parameters of the simulation
dt = 0.01 # time step (ms)
t_end = 100 # end time (ms)

# Define the time of the external impulse
t_impulse = 0.04 # time of the external impulse (ms)

# Define the simulation loop
num_steps = Int(t_end / dt)
V = initial_membrane_potentials
g = initial_conductances
num_neurons = 302
num_channels = 302
E_syn = 0 
E_k = -90 * connectivity
membrane_potential_history = Vector{Matrix{Float64}}() 
#conductance_history = zeros(num_neurons, num_channels, num_steps)
for i in 1:num_steps
    t = i * dt # current time (ms)

    I_syn = zeros(num_neurons)
    # Inject an external current into neuron 1 at time t_impulse
    if t == t_impulse
        I_ext = 0 # amplitude of the external current (nA)
        I_syn[1] += I_ext # add the external current to the synaptic current of neuron 1
    end
    
    # Compute the synaptic currents for each neuron
    for j in 1:num_neurons
        for k in 1:num_neurons
            if connectivity[k,j] == 1
                I_syn[j] += synaptic_current(g[k,j], V[k], V[j], E_syn)
            end
        end
    end
    
    # Compute the ion channel currents for each neuron
    I_ion = zeros(num_neurons, num_channels)
    for j in 1:num_neurons
        for k in 1:num_channels
            I_ion[j,k] = ion_channel_dynamics(g[j,k], V[j], E_k[k])
        end
    end
    
    # Compute the new membrane potentials for each neuron
    for j in 1:num_neurons
        V[j] = neuron_dynamics(V[j], sum(I_syn) + sum(I_ion[j,:]))
    end

    
    # Record the membrane potentials and conductances for each neuron
    push!(membrane_potential_history, V)
    #conductance_history[:,:,i] = g
end

# Determine which path is taken
threshold = -40  # in mV
println(maximum(membrane_potential_history[1,end]))
if maximum(membrane_potential_history[1,end]) > threshold
    println("The impulse has taken the path that includes neuron 1.")
else
    println("The impulse has taken a different path.")
end



# Plot with subplots for each neuron
plot_array = plot(layout=(5,1))
for i in 1:5
    plot!(plot_array[i], V[i,:], color=:purple, xlabel="Time (ms)", ylabel="Membrane Potential (mV)")
end
display(plot_array)


"""