using DataFrames
using CSV
using XLSX
using LinearAlgebra


data_connect_phar = CSV.read("RawData/CONEXIONESPHARINGEAL.csv", DataFrame)
"""
Sending: Name of sending neuron
Receiving: Name of receiving neuron
Number: Number of synapses between the given neuron pair.
Type: Type of synapse: S: synaptic; G: gap.
"""


data_connect_neuron = DataFrame(XLSX.readtable("RawData/NeuronConnect.xlsx", "Sheet1"))
"""
N1: Neuron 1 name
N2: Neuron 2 name
Type: Type of synapse: 
    S: Send or output (Neuron 1 pre-synaptic to Neuron 2); 
    Sp: Send-poly (Neuron 1 is pre-synaptic to more than one postsynaptic partner. Neuron 2 is just one of these post-synaptic neurons); 
    R: Receive or input (Neuron 1 is post-synaptic to Neuron 2); 
    Rp: Receive-poly (Neuron 1 is one of several post-synaptic partners of Neuron 2.); 
    EJ: Electric junction; NMJ: Neuromuscular junction (only reconstructed NMJ's are represented).
It must be noted that at polyadic synaptic sites, not all “send-poly” were faithfully labeled 
    as such in White et al, 1986. Some pre-synaptic connections were labeled simply as “sends”. 
    Reconciliation of chemical synapses did not previously distinguish between send from send-poly 
    and receive from receive-poly. In this new reconciliation, the total number of send and send-poly 
    is equal to the total number of receive and receive-poly (S+Sp=R+Rp). Every documented synapse is 
    now listed in this Table, both with respect to the sending neuron and with respect to the receiving neuron(s).
Nbr: Number of synapses between the given neuron pair.
"""

data_type_neuron = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet1"))
"""
Neuron: Name of neuron
Soma Position: Position of cell body along the AP axis of worm body. 0=tip of nose; 1=tail tip.
Soma region: Cell body position by head, mid-body, or tail region.
Span: Length of neuron span. Neurons spanning <25% of worm body (e.g., motor neurons in the ventral cord, 
    neurons with processes confined to the nerve ring and neurons confined in the mid-body) are defined to 
    have short spans (S). All other neurons are defined to have long spans (L).
Ambiguity: If applicable, code for the type of ambiguity. Codes beginning with M denote ambiguity 
    citedin White et al, 1986. Codes beginning with R denote ambiguity found in reconstructions during 
    update of wiring diagram (MB=cell body position ambiguous, MTS=tail synapses ambiguous and/or sparse 
    connections in the tail; MAS=anterior body ambiguous and/or sparse connections in the anterior; MD=dorsal 
    side ambiguous; MAD=anterior and dorsal side ambiguous; MS=neurons with sublateral processes not covered 
    by reconstructions. RDI=dorsal reconstruction incomplete; RDM=dorsal reconstruction completely missing; 
    RVI=ventral reconstruction incomplete.)
TotHead: Total number of synapses in the head including EJ and NMJ.
TotTail: Total number of synapses in the tail including EJ and NMJ.
TotMid: Total number of synapses in the mid-body including EJ and NMJ.
S_Head: Number of “sends” or output synapses in the head, includes send polyadic synapses (see Figure 1) and NMJ.
R_Head: Number of “receives” or input synapses (includes polyadic synapses) in the head.
S_Mid: Number of “sends” or output synapses in the mid-body, includes polyadic synapses and NMJ.
R_Mid: Number of “receives” or input synapses (includes polyadic synapses) in the mid-body.
S_Tail: Number of “sends” or output synapses in the tail, includes polyadic synapses and NMJ.
R_Tail: Number of “receives” or input synapses (includes polyadic synapses) in the tail.
AY NeuronType: Letter codes denoting ganglion group as defined by Achacoso and Yamamoto W.S., 1991, where 
    A=anterior ganglion, B=dorsal ganglion, C=lateral ganglion, D=ventral ganglion, E=retrovesicular ganglion, 
    F=posterolateral ganglion, G=ventral cord neuron group, H=pre-anal ganglion, J=dorsorectal ganglion, 
    K=lumbar ganglion.
AYNbr: Numeric identifier given by AY for each neuron.
Note:  Sum of S_Head and R_Head does not include electrical junctions (EJ), therefore, does not equal TotHead.  Similar is true for mid-body and tail.
"""


data_type_phar = DataFrame(XLSX.readtable("RawData/NeuronType.xlsx", "Sheet2"))
"""
Neuron: Name of neuron
Soma Position: Position of cell body along the AP axis of worm body. 0=tip of nose; 1=tail tip.
Soma region: Cell body position by head, mid-body, or tail region.
"""

data_connect_monoamine = CSV.read("RawData/MonoaminesConnect.csv", DataFrame)
"""
Neuron1: Name of sending neuron
Neuron2: Name of receiving neuron
Type: Name of the monoamine
Monoamine: MA stands for monoamine
Specific: The specific type of monoamine
"""




""" IMPORTANT VARIABLES:
data_connect_phar
data_connect_neuron
data_type_neuron
data_type_phar
data_connect_monoamine"""




# DATAFRAMES CLEANING
# GAP AND SYNAPTIC
# From all the dataframe of data_type_neuron only select the first three columns
data_type_neuron = data_type_neuron[:, 1:3]
# Concatenate the two dataframes about neuron type
data_type = vcat(data_type_phar, data_type_neuron)
# Add the two neurons that do not have any connections
push!(data_type,["CANL", 0.61, "M"])
push!(data_type,["CANR", 0.61, "M"])
# Sort the dataframe by the position of the neuron in the soma
data_type_sort_soma = sort!(data_type, [:"Soma Position"])
# Add a column to know for the future the number of the neuron
data_type_sort_soma.Place = 1:302

" IMPORTANT VARIABLES:
data_type_sort_soma
"

# Interchange the "Number" and "Type" columns
select!(data_connect_phar, [:Sending,:Receiving,:Type, :Number])
# Rename columns so they match the names of "data_connect_phar"
rename!(data_connect_neuron,[:Sending,:Receiving,:Type, :Number])
# Concatenate the two dataframes about neuron connectivity
data_connect = vcat(data_connect_phar, data_connect_neuron)

# Eliminate the rows in which the "Type" is either receive, receive-poly or NMJ.
data_connect = data_connect[data_connect.Type .!= "R", :]
data_connect = data_connect[data_connect.Type .!= "Rp", :]
data_connect = data_connect[data_connect.Type .!= "NMJ", :]


# TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR GAP AND SYNAPTIC
# Create a new dictionary to append the indexes
data_index = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# Create two vectors to store the indeces before appending
sending_value = Vector{Int64}()
receiving_value = Vector{Int64}()
i = 1
for row in eachrow(data_connect)
    for value in eachrow(data_type_sort_soma)
        if row[:Sending] == value[:Neuron]
            append!(sending_value, value[:Place])
        end
        if row[:Receiving] == value[:Neuron]
            append!(receiving_value, value[:Place])
        end    
    end
    push!(data_index, (sending_value[i], receiving_value[i]))
    i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indeces
data_connect = hcat(data_connect, data_index)

# Select only the synaptic connections and join them
data_connect_S = data_connect[data_connect.Type .== "S", :]
data_connect_Sp = data_connect[data_connect.Type .== "Sp", :]
data_connect_synaptic = vcat(data_connect_S, data_connect_Sp)

# Select the gap juntions and join them
data_connect_G = data_connect[data_connect.Type .== "G", :]
data_connect_EJ = data_connect[data_connect.Type .== "EJ", :]
data_connect_gap = vcat(data_connect_G, data_connect_EJ)



# DATAFRAMES CLEANING
# MONOAMINES AND NUEROPEPTIDES
# From all the dataframe of data_type_neuron only select the first three columns
data_connect_monoamine = data_connect_monoamine[:, 1:3]

# TRANSFORM THE CONNECTIONS FROM NEURON NAMES TO NUMBERS FOR MONOAMINES
# Create a new dictionary to append the indexes
data_index_mono = DataFrame(IndexSending=Any[], IndexReceiving=Any[])
# Create two vectors to store the indeces before appending
sending_value_mono = Vector{Int64}()
receiving_value_mono = Vector{Int64}()
i = 1
for row in eachrow(data_connect_monoamine)
    for value in eachrow(data_type_sort_soma)
        if row[:Neuron1] == value[:Neuron]
            append!(sending_value_mono, value[:Place])
        end
        if row[:Neuron2] == value[:Neuron]
            append!(receiving_value_mono, value[:Place])
        end    
    end
    push!(data_index_mono, (sending_value_mono[i], receiving_value_mono[i]))
    i = i+1
end
# Concatenate horizontally the dictionary of data_connect and indeces
data_connect_monoamine = hcat(data_connect_monoamine, data_index_mono)


" IMPORTANT VARIABLES:
data_connect_synaptic
data_connect_gap
"


# Plot the number of connections of the synaptic and gap
synaptic_number = zeros(302, 302)
synaptic_connections = zeros(302, 302)
for link in eachrow(data_connect_synaptic)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    value = link[:Number]
    synaptic_number[from_index, to_index] = value   
    synaptic_connections[from_index, to_index] = 1
end
spy(synaptic_number, plot_title= "Number of synaptic connections among neurons", xlabel = "Neuron index", ylabel = "Neuron index")
spy(synaptic_connections, plot_title= "Synaptic connections among neurons", xlabel = "Neuron index", ylabel = "Neuron index")

gap_number = zeros(302, 302)
gap_connections = zeros(302,302)
for link in eachrow(data_connect_gap)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    value = link[:Number]
    gap_number[from_index, to_index] = value    
    gap_connections[from_index, to_index] = 1
end
spy(gap_number, plot_title= "Number of gap connections among neurons", xlabel = "Neuron index", ylabel = "Neuron index")
spy(gap_connections, plot_title= "Gap connections among neurons", xlabel = "Neuron index", ylabel = "Neuron index")
xticks!(1:50:302)
yticks!(1:50:302)

# Monoamines
mono_connections_tyr = zeros(302,302)
mono_connections_oct = zeros(302,302)
mono_connections_dop = zeros(302,302)
mono_connections_ser = zeros(302,302)
for link in eachrow(data_connect_monoamine)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    if link[:Type] == "tyramine"
        mono_connections_tyr[from_index, to_index] = 1
    elseif link[:Type] == "octopamine"
        mono_connections_oct[from_index, to_index] = 1
    elseif link[:Type] == "dopamine"
        mono_connections_dop[from_index, to_index] = 1
    elseif link[:Type] == "serotonin"
        mono_connections_ser[from_index, to_index] = 1
    end
end
spy(mono_connections_tyr, color = :orange, plot_title= "Monoamine connections among neurons", xlabel = "Neuron index", ylabel = "Neuron index")
spy!(mono_connections_oct, color = :red)
spy!(mono_connections_dop, color = :green)
spy!(mono_connections_ser, color = :blue)
xticks!(1:50:302)
yticks!(1:50:302)



# Set variables
N = 302 # number of neurons
E_rev = -65.0 #mV Reversal potential of the synapse
spikeThresh = 0 #mV Spiking threshold
specific_capacitance = 1 #uF/cm2
intracellular_resistivity = 0.03 #kΩ*cm
g = 100 #pS Conductance (Varshney et al., 2011)
Gc = 10 #pS Cell membrane conductance 
C = 1.5 #pF Membrane capacitance (Varshney et al., 2011)
Ecell = -35.0 #mV Leakage potential 
#while reversal potential E_j = 0mV for excitatory synapses and −48 mV for inhibitory synapses (Wicks et al., 1996). For the synaptic activity variable, we take ar=11.5, ad=51.5 and width of the sigmoid β = 0.125mV−1 (Wicks et al., 1996). Also for the initial condition of the membrane voltages V and synaptic activity variable s, we sample the normal distribution of μ = 0 and σ = 0.94 with size 279 * 2 (for both V and s) and multiply by 10−4
C = 0.015 # Cell Membrane Capacitance
ar = 1/1.5
ad = 5/1.5

"""
Threshold potential for each neuron is computed by imposing dVi/dt=0 (Equation 2 for C. elegans) and solving for Vi. This is equivalent to Solving the following system of linear equations
Ax=b
(11)
A=M1+M2+M3; b=−b1−b3−Iext,
(12)
where the solution x is N × 1 vector with each entry being the threshold potential Vthreshold for the ith neuron.
M1 is a matrix of size N × N where N is the number of neurons (279 for C. elegans) with its diagonal terms populated with −Gc (cell membrane capacitance).
M2 is a diagonal matrix where diagonal term in ith row corresponds to −∑jGgij i.e., the sum of total conductivity of gap junctions for the ith neuron.
M3 is a diagonal matrix where its ith diagonal term corresponds to −∑jseqGsij, where seq=arar+2ad and Gsij is maximum total conductivity of synapses to i from j. Note that seq is obtained by imposing dsidt=0 and synaptic activation Φ = 1/2 in Equation 5.
b1=Gc∗Ec where Ec is a 1D vector of size N × 1 in which all of its elements are Ec (leakage potential).
b3=Gs⋅(s∗eqEj) where Ej is a 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −48 mV if inhibitory).
Iext is the input stimuli vector where its ith element determines the input current amplitude for the ith neuron.
"""

println("----------------------------")

show(data_type_sort_soma[!,:Neuron])


"""
# Ax = b ===> A = M1+M2+M3; b = -b1-b3-Iext
function threshold_potential_computation()
    # M1 computation
    Gc_diag = zeros(302,302)  # diagonal matrix with membrane conductance
    for i in 1:N
        Gc_diag[i, i] = Gc
    end
    M1 = -Gc_diag

    # M2 computation
    g_gap = 1.0 # condctance of each gap channel
    gap_number_cond = g_gap * gap_number  # multiply the conductance of each channel by the number of connections
    gap_diag_array = sum(gap_number_cond, dims = 2)   # sum of total conductivity for each neuron in its row
    gap_diag = Diagonal(vec(gap_diag_array))  # diagonal matrix with the values of above
    M2 = -gap_diag  # minus the above

    # M3 computation
    g_syn = 1.0 # conductance of each synapsis
    s_eq = ar / (ar + 2 * ad)
    synaptic_number_cond = g_syn * synaptic_number
    mult = s_eq * synaptic_number_cond
    mult_diag_array = sum(mult, dims = 2)
    mult_diag = Diagonal(vec(mult_diag_array))
    M3 = -mult_diag

    # b1 computation
    Ec = fill(Ecell, N) # 1D vector of size N × 1 in which all of its elements are Ec (leakage potential)
    b1 = Gc *Ec  # cell membrane conductance * Ec

    # b3 computation
    #Ej # 1D vector of size N × 1 that enlists the directionality of each neuron (0 if excitatory or −48 mV if inhibitory)
    b3 = synaptic_number_cond .* (s_eq * Ej)


    M = M1 + M2 + M3
    b = - b1 - b3 # - Iext

end

"""





