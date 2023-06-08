using DifferentialEquations
using Plots
include("data_extraction.jl")


# Define the parameters of the model
num_neurons = 302


# Initialize the state of the model
state = zeros(num_neurons)

function model(du, u, t)
    # Calculate the rate of change of each neuron
    for n in 1:num_neurons
      # Calculate the sum of the incoming connections
      sum_inputs = 0
      for j in 1:num_neurons
        if u[j] > 0 && synaptic_number[n,j] > 0
          sum_inputs += u[j]
        end
      end
  
      # Calculate the rate of change of the neuron
      du[n] = sum_inputs - u[n]
    end
  
    return du
  end

# Create the ODE system
ode = ODESystem(model, state)

# Solve the ODE system
tspan = (0.0, 3.0)
sol = solve(ode, tspan)

# Plot the results
plot(sol)
