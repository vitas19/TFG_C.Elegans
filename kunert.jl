using DifferentialEquations
using Plots
using JSON

# Define the Neuron type
struct Neuron
    id::Int
    Cm::Float64
    gL::Float64
    EL::Float64
    gNa::Float64
    ENa::Float64
    gK::Float64
    EK::Float64
    gCa::Float64
    ECa::Float64
    gm::Float64
    Em::Float64
    Ih::Float64
    Eh::Float64
    V::Float64
    m::Float64
    h::Float64
    n::Float64
    a::Float64
    b::Float64
    p::Float64
    q::Float64
end

# Define the Network type
struct Network
    neurons::Vector{Neuron}
    weights::Matrix{Float64}
end


function kunert_model(du, u, p, t)
    # Extract variables from input vector
    V, n, m, h = u
    
    # Extract parameters from input vector
    gna, gk, gl, ena, ek, el, Cm, Iapp = p
    
    # Compute membrane currents
    ina = gna * m^3 * h * (V - ena)
    ik = gk * n^4 * (V - ek)
    il = gl * (V - el)
    I = Iapp
    
    # Compute derivatives of gating variables
    ninf = 1 / (1 + exp((-V - 34) / 6.5))
    taun = 0.5 / (0.0197 * (1 + exp(-0.0337 * (V + 14.56))))
    minf = 1 / (1 + exp((-V - 30.6) / 6.25))
    tauh = 0.25 + 4.35 / (1 + exp(-0.1 * (V + 21)))
    mtau = 1 / (1.32 * exp((V + 119.0) / (-40.0)) + 0.23 * exp((V - 31.8) / 20.5))
    hinf = 1 / (1 + exp((V + 50) / 5.5))
    
    # Compute derivatives of membrane potential
    dVdt = (I - ina - ik - il) / Cm
    dndt = (ninf - n) / taun
    dmdt = (minf - m) / mtau
    dhdt = (hinf - h) / tauh
    
    # Assign output vector
    du[1] = dVdt
    du[2] = dndt
    du[3] = dmdt
    du[4] = dhdt
end



# Define network architecture and connectivity
network = Network()
neurons = []
for i in 1:302
    neuron = Neuron()
    add_node!(network, neuron)
    push!(neurons, neuron)
end

# Define connectivity between neurons
# read in the data from the file
syn_data = JSON.parsefile("chem.json")

nodes = syn_data["nodes"]
links = syn_data["links"]

strings=[]

num_neurons = length(nodes)
link_matrix = zeros(Float32, num_neurons, num_neurons)
connectivity = zeros(302, 302)
for link in links
    from_index = link["source"] +1
    to_index = link["target"] +1
    value = link["value"]
    connectivity[from_index, to_index] = value
    push!(strings, value)
    
end

spy(connectivity)
savefig("spy_plot.png")

weights = connectivity

# Define initial conditions for neurons
u0 = [0.0, 0.0, 0.0, 0.0]  # initial membrane potential, gating variables

# Define time span of simulation
tspan = (0.0, 10.0)  # start and end times

# Define parameters for simulation
gna = 120.0  # maximum sodium conductance (mS/cm^2)
gk = 36.0  # maximum potassium conductance (mS/cm^2)
gl = 0.3  # leak conductance (mS/cm^2)
ena = 55.0  # sodium reversal potential (mV)
ek = -90.0  # potassium reversal potential (mV)
el = -54.4  # leak reversal potential (mV)
Cm = 1.0  # membrane capacitance (uF/cm^2)
Iapp = 0.0  # applied current (uA/cm^2)
p = [gna, gk, gl, ena, ek, el, Cm, Iapp]  # parameters as a vector




solver = Tsit5()


# Solve differential equations for each neuron
for i in eachindex(neurons)
    neuron = neurons[i]
    p[8] = neuron.iapp  # set applied current for neuron
    prob = ODEProblem(kunert_model, u0, tspan, p)
    sol = solve(prob, solver)
    neuron.v = sol[1:end, 1]  # store membrane potential for neuron
    neuron.n = sol[1:end, 2]  # store gating variables for neuron
    neuron.m = sol[1:end, 3]
    neuron.h = sol[1:end, 4]
end


# Plot membrane potential for each neuron
plot()
for i in eachindex(neurons)
    neuron = neurons[i]
    plot!(sol.t, neuron.v, label="Neuron $i")
end
xlabel!("Time (ms)")
ylabel!("Membrane potential (mV)")
title!("Membrane potential of 302 neurons in C. elegans")
