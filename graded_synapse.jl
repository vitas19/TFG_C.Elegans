using DifferentialEquations
using Plots

# Define the ODE function
function synapse_dynamics!(du, u, p, t)
    # Unpack parameters
    gs, Vpost, Es, Vpre, Vth, Δ, k = p
    
    # Unpack state variables
    s = u[1]
    
    ss = 1 / (1 + exp((Vth - Vpre) / Δ))

    # Compute ts and s_dot
    ts = (1 - ss) / k
    s_dot = (ss - s) / ts
    
    # Compute IS
    IS = gs * s * (Vpost - Es)
    
    # Update state variables
    du[1] = s_dot
    
    return IS
end

# Set initial conditions and parameters
s0 = [0.1]
p = (0.47, -90.0, 0.0, -74.0, -37.0, 20.0, 2*log(0.9/0.1), 1.0)  # gs, Vpost, Es, Vpre, Vth, Δ, k, ts

#Vth, Δ, ts NO SE DONDE ESTA SU VALOR

# Set time span and solve the system of equations
tspan = (0.0, 100.0)
prob = ODEProblem(synapse_dynamics!, s0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# Plot the results
plot(sol.t, sol[1,:],xlimits=(0,3), xlabel="Time (ms)", ylabel="s", legend=false)
