using DifferentialEquations

function kunert_equations(v, t, params)
    C, g_L, E_L, g_I, E_I, I_app = params
    dvdt = -g_L (v - E_L) - g_I (v - E_I) + I_app
    return dvdt
end

function simulate_kunert_equations(params, tspan, v0)
    prob = ODEProblem(kunert_equations, v0, tspan, params)
    sol = solve(prob)
    return sol
end

params = (C = 1e-12, g_L = 1e-3, E_L = -60 mV, g_I = 1e-2, E_I = 0 mV, I_app = 0 uA)
tspan = (0 s, 100 ms)
v0 = 0 mV

# Initialize the state vector
state = zeros(302, 1)

# Simulate the activity of all 302 neurons
for i in 1:302
    state[i] = simulate_kunert_equations(params, tspan, v0)[i]
end

# Plot the membrane potentials of all 302 neurons
plot(state)
