# KUNERT MODEL 
include("data_preprocessing.jl") 
include("plots.jl")
using Plots
# Constants 

" Number of neurons in the Celegans connectome."
const nNeurons = 302

"""
mV, Reversal potential of the synapse (-65 mV)

But Wicks et al, 1996, reports a private communication with R.E.Davis:
- excitatory synapses: 0mV
- inhibitory synapses: -45mV

"""
const E_exc = 0#mV
const E_inh = -48#mv 
const E_rev = -65.0 #mV, Reversal potential of the synapse (-65 mV)

const spikeThresh = 0 #mV Spiking threshold?????
const specific_capacitance = 1 #uF/cm2
const intracellular_resistivity = 0.03 #kΩ*cm

"""
100pS, Conductance gap and synaptic (Varshney et al., 2011)

But Wicks et al, 1996 suggests 5nS, 10x for sensibility analysis.
"""
const gᵞ = 100 #pS#FVA OK

const Gc = 10 #pS, Cell membrane conductance

"""
Membrane capacitance 1(pF) (Varshney et al., 2011)

"""
const C = 1 #pF, Membrane capacitance 

"""
Leakage potential (-35mV, Wicks et al, 1996, from David, pers.comm.)
"""
const Ecell = -35.0 #mV, Leakage potential

"""
Beta constant(mV−1) in sigmoid function.
"""
const beta = 0.125  #mV−1, const

"""
Synaptic activity rise time ( Kunert et al, 2014)
"""
const ar = 1 

"""
Synaptic activity decay time ( Kunert et al, 2014)
"""
const ad = 5 


# THRESHOLD POTENTIAL COMPUTATION
"""
Threshold potential for each neuron is computed by imposing dVi/dt=0 (Equation 2 for C. elegans) and solving for Vi. This is equivalent to Solving the following system of linear equations
    Ax = b ===> A = M1+M2+M3;
where the solution x is N × 1 vector with each entry being the threshold potential Vthreshold for the ith neuron.
M1 is a matrix of size N × N where N is the number of neurons (279 for C. elegans) with its diagonal terms populated with −Gc (cell membrane capacitance).
M2 is a diagonal matrix where diagonal term in ith row corresponds to −∑jGgij i.e., the sum of total conductivity of gap junctions for the ith neuron.
M3 is a diagonal matrix where its ith diagonal term corresponds to −∑jseqGsij, where seq=arar+2ad and Gsij is maximum total conductivity of synapses to i from j. Note that seq is obtained by imposing dsidt=0 and synaptic activation Φ = 1/2 in Equation 5.
b1=Gc∗Ec where Ec is a 1D vector of size N × 1 in which all of its elements are Ec (leakage potential).
b3=Gs⋅(s∗eqEj) where Ej is a 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −48 mV if inhibitory).
Iext is the input stimuli vector where its ith element determines the input current amplitude for the ith neuron.
"""

# M1 computation:
M1 = Diagonal(fill(-Gc, nNeurons))

# M2 computation:
M2 = Diagonal(-gᵞ * vec(sum((gᵞ * gap_number), dims=2)))
#@assert isless(vec(sum(gapJunctionConnectome, dims=2)) - gapJunctionConnectome * ones(nNeurons), fill(eps(Float64), nNeurons))

# M3 computation:
s_eq = ar / (ar + 2 * ad)       # Formula to compute the synaptic activity in eq# uilibrium (is obtained by imposing dsi/dt=0 and synaptic activation Φ = 1/2 in Kunert equations)
M3 = Diagonal(- gᵞ * s_eq * vec(sum(synaptic_number, dims=2)));

# b1 computation
b1 = Gc * fill(Ecell, nNeurons) # 1D vector of size N × 1 in which all of its elements are Ec (leakage potential)

# b3 computation
#Ej -> 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −48 mV if inhibitory)
neurotransmitter1 = neurotransmitter[!, "Neurotransmitter1"]
Ej = [neurotrans == "GABA" || neurotrans == "Dopamine" ? -48 : 0 for neurotrans in neurotransmitter1]

b3 = gᵞ * (s_eq * Ej)   # 1D vector of size N × 1 

"""
Excitatory current for sensory experiment on CO2 sensing on BAG in C.elegans. In mV(?)
"""
Iext = zeros(nNeurons);#zero-input condition

# Compute A and b
A = M1 + M2 + M3;
b = - b1 - b3 - Iext;

# Solve Ax=b for x
A_inv = inv(A);      # Inverse of a
"""
Vth is a N × 1 vector that contains the threshold potential value for each neuron
"""
Vth = A_inv * b

# Plot the vector Vth
gr()        
Plots.plot(Vth)           
histogram(Vth,
          bins=:scott, 
          xticks = unique(Vth), legend = false, normalize=:probability,
          xlabel = "Vth of sigmoid",
          ylabel = "Frequency",
          title = "Threshold voltage histogram")
histogram!(size=(1500,500))
###

#---------------------------------------------------------


###############################################################################
# Model 1 using naked Dif Equations
#using DifferentialEquations
using Pkg;Pkg.instantiate()
using OrdinaryDiffEq
##############################################################################

# Define the derivatives of the states with respect to time for each neuron.
# - u is the neuron voltage
# - i is the neuron activation state/level
# all constants imported from Celegans.Constants
Gᵞ = gᵞ * gap_number #The matrix of GJ conductances

Sᵞ = gᵞ * synaptic_number#Matrix of Synaptic conductances
EJ = repeat(Ej',nNeurons,1) #Same ROW is the same
spy(EJ)
@assert all(vec(EJ[rand(1:nNeurons),:] .== EJ[rand(1:nNeurons),:]))
#typeof(EJ[rand(1:nNeurons),:] .== EJ[rand(1:nNeurons),:])
function derivatives(nStates, p, t)
    #dv = (-Gc * (nStates[:, 1] - Ecell) - eq_IiGap[1] - eq_IiSyn[1] +  ) / C
    allV = repeat(nStates[:,1],1,nNeurons)#Makes the space for it. 
    dv = (-Gc * (nStates[:, 1] .- Ecell) #- eq_IiGap[1]
          -  (sum( Gᵞ .* (allV - allV'), dims=2))
          -  ( Sᵞ .* (allV - EJ)) * nStates[:, 2]) ./ C  #- eq_IiSyn[1] +  ) / C
    di = ar * (1 ./ (1 .+ exp.(-beta * (nStates[:, 1] .- Vth)))) .* (1 .- nStates[:, 2]) .- ad * nStates[:, 2]
  return hcat(dv,di)
end

"""
In place derivative function for the easier solvers.

- Note: excitations Iₑ have are encoded here as a current vector for each neuron
@example
S = zeros(nNeurons,2)
dS = zeros(nNeurons,2)
Iₑ = Ie #see below
#Iₑ = Array{Any}(undef,nNeurons);
#null = t -> 0.0
#Iₑ[:] .= null
#Iₑ = t -> sin(t)
t = 0.010
"""
function celegansChangeState!(dS, S, Iₑ, t)
    allV = repeat(S[:,1],1,nNeurons);#Makes the space for it. 
    dS[:,1] .= (-Gc * (S[:, 1] .- Ecell) 
                -  (sum( Gᵞ .* (allV - allV'), dims=2))
                -  (Sᵞ .* (allV - EJ)) * S[:, 2]
                + Iₑ(t)
                ) ./ C
    dS[:,2] .= ar * (1 ./ (1 .+ exp.(-beta * (S[:, 1] .- Vth)))) .* (1 .- S[:, 2]) .- ad * S[:, 2]
  #return hcat(dv,di)
end

"""
A curried function for two parameters takes:
- fVector:a list of functions of time (WARNING!: Unchecked!)
- t:a time instant
and evaluates all the functions, returning a vector of the same length
"""
vectorInput = fVector ->
    function(t)
        return([f(v) for (f, v) in zip(fVector,fill(t,size(fVector)))])
    end
#vectorInput("hello")(1.0)

"""
        null(t)

The zero function on a real domain
"""
null = (t::Real) -> 0.0

"""
        pulse(A,D,τ)

A function to generate a symmetric pulse of height [A], and duration [D] with
positive delay [τ]. The pulse is only symmetric if τ=0.

"""
pulse = function(A,D,τ)
    function(t::Real)
        return(abs(t - τ) <= D/2 ? A : 0.0)
    end
end

pulse(20.0, 1.0, 0.5)(-0.01)
pulse(20.0, 1.0, 0.5)(0.0)
pulse(20.0, 1.0, 0.5)(0.5)
pulse(20.0, 1.0, 0.5)(1.5)


# 5.1. Put together a dataframe with the timeseries of the affected neurons
# Experiment from Kunert, Shlizerman and Kutz, 2014
# 1. Get the neuron names of some given neuron name patterns.
inputNeuronPatterns = ["PLM"]
outputNeuronPatterns = ["AVA", "AVB", "AVD", "AVE", "PVC"]
focusNeuronPatterns =  vcat(inputNeuronPatterns, outputNeuronPatterns)


function findNeuronByIndices(NeuronPatterns)
    results = []
    for pattern in NeuronPatterns
        indices = findall(contains(pattern), data_type_sort_soma.Neuron)
        push!(results, indices)
    end
    return vcat(results...)
end

neuronIndices = findNeuronByIndices(focusNeuronPatterns)



"""
Example vector of excitation functions for an experiment.
"""
If = Array{Function}(undef,nNeurons);
If[:] .= null
A = 100.0#pA, input level for the pulse
If[findNeuronByIndices(inputNeuronPatterns)] .= pulse(A, 1.0, 0.5)
Ie = vectorInput(If)#This is the matrix of excitations
# FVA: some assertions should follow.
matrixInputs = mapreduce(Ie, hcat, 0.0:0.05:1.0);
[findNeuronByIndices(inputNeuronPatterns), 1:10]

# dState = ODEFunction{true}(
#     derivatives;
#     syms=repeat(neuron_list.Neuron, 1,2)
# )
dStateIIP = ODEFunction{true}(
    celegansChangeState!;
    syms=repeat(data_type_sort_soma.Neuron, 1,2)    
)
#######################
TEST = true
#TEST = false
#######################
# Create a vector of nNeurons neurons, each with two states (u and i).
# First in pair is neuron voltage, second is activation
"""
nState - Neuron voltage in column 1 and activation state of synapses in 2.
"""
nStates = zeros(nNeurons, 2)
# if !TEST#set the states already in their equilibrium
#     nStates[:,2] = Vth#Let's place the synapses in their original activation
# end

"""
For the worm we use:
- time steps of 10ms
- time spans of 10 m = 600s
For tests, timespans of 0.5m = 30s
 """
if !TEST
    tspan = (0.0, 100.0)
else
    tspan = (0.0, 30.0)#only for half a minute.
end
tspan = (0.0, 2.0)#for inputs lasting like a minute.
#tspan = (0.0, 0.5)#for inputs lasting like a minute.
# gusano1 = ODEProblem(dState, nStates, tspan)#Not working
# sol1 = solve(gusano1)
gusano2 =
    ODEProblem(celegansChangeState!,#dStateIIP,
               nStates,
#               tspan
               tspan,
               Ie
               )
fieldnames(ODEProblem)

"""
        algo - the algorithm for solving the ODE
In DifferentialEquations.jl, some good “go-to” choices for ODEs are:

* AutoTsit5(Rosenbrock23()) handles both stiff and non-stiff equations. This is a good algorithm to use if you know nothing about the equation.
* AutoVern7(Rodas5()) handles both stiff and non-stiff equations in a way that's efficient for high accuracy.
* Tsit5() for standard non-stiff. This is the first algorithm to try in most cases.
* BS3() for fast low accuracy non-stiff.
* Vern7() for high accuracy non-stiff.
* Rodas4() or Rodas5() for small stiff equations with Julia-defined types, events, etc.
* KenCarp4() or TRBDF2() for medium-sized (100-2000 ODEs) stiff equations
* RadauIIA5() for really high accuracy stiff equations
* QNDF() for large stiff equations
    """
timestep = 0.02# 20ms to sample at a rate of 50Hz
#algo = KenCarp4()#Is this too demandind?
algo = Tsit5()
using TimerOutputs
sol2 = @time "solve worm" solve(gusano2, algo;
             alg_hints = [:stiff],#want more accuracy in solutions
             abstol = 1e-8,#This and previous param to demand more accuracy
             reltol = 1e-8,#time tolerance
             saveat = timestep#interpolate solutions
             )
fieldnames(ODESolution)

# Plot the results of the simulation.
# TODO: tspan must be transformed to the time scale
########################################################################
# using Plots
########################################################################
# gr()
# focusNeuronPattern = "BAG"
# focusNeuronIndices =
#     findall(map(x -> str_detect(x, focusNeuronPattern), newNameByIndex))
# focusNeurons = newNameByIndex[focusNeuronIndices]
# #sol[2, 1, :], is the timeseries for the component, which is the 2nd row and 1 column.
# sol2[focusNeuronIndices[1],1,:]
# display(plot(sol2.t, sol2[focusNeuronIndices,1,:]'))

########################################################################
using Gadfly# For the plotting
########################################################################
#display(Gadfly.plot(x=sol2.t,  y=sol2[focusNeuronIndices[1],1,:]', Geom.line))
# 2. Get the timeseries for these neurons and create a dataframe with the ordered focussed Neuron set
df = 
    DataFrame(
        Dict(#"Time" => sol2.t,
            zip(
                data_type_sort_soma.Neuron[neuronIndices],
                map(id -> sol2[id,1,:], neuronIndices)
            )
        )
    );
df[!, "Time"] = sol2.t;#Adding the time as an independent criterion
describe(df)
using Chain
#dfs = stack(df, newNameByIndex[neuronIndices], "Time");
dfs = @chain stack(df, Not("Time")) begin
    rename(:variable => :Neuron, :value => :Voltage)
end

show(dfs, allrows = true)


# 3. Plot the timeseries according to whatever criterion
p = Gadfly.plot(dfs, x=:Time, y=:Voltage, color=:Neuron, Geom.line)
img = SVG("100mV_PLM_input.svg")
import Cairo, Fontconfig
img = PDF("100mV_PLM_input.pdf", 10cm,dpi=300)
draw(img, p)
display(p)

"""
println("4. Environment description")
using Pkg;Pkg.status()
exit(0)
"""