# Define connectivity between neurons
# read in the data from the file
syn_data = JSON.parsefile("gap.json")

nodes = syn_data["nodes"]
links = syn_data["links"]

strings=[]

num_neurons = length(nodes)
gap_connectivity = zeros(302, 302)
for link in links
    from_index = link["source"] +1
    to_index = link["target"] +1
    value = link["value"]
    gap_connectivity[from_index, to_index] = value
    push!(strings, value)
    
end

spy(gap_connectivity)


# Define the gap junction formula function
function igap(ggap, vpost, vpre)
    return ggap * (vpost - vpre)
end

# Set the neuron voltages
voltages = zeros(302) # Initialize all voltages to 0
voltages[1] = -50.0 # Set the voltage of the first neuron to -50 mV


# Calculate the total gap junction current
total_igap = 0
println(gap_connectivity, 1)
for i in eachindex(gap_connectivity, 1)
    pre = gap_connectivity[i, 1]
    post = gap_connectivity[i, 2]
    ggap = 0.1 # Set the gap junction conductance to a constant value
    total_igap += igap(ggap, voltages[post], voltages[pre])
end

println("Total gap junction current = ", total_igap)