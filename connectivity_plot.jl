# Define connectivity between neurons
# read in the data from the file
syn_data = JSON.parsefile("chem.json")

nodes = syn_data["nodes"]
links = syn_data["links"]

strings=[]

num_neurons = length(nodes)
connectivity = zeros(302, 302)
for link in links
    from_index = link["source"] +1
    to_index = link["target"] +1
    value = link["value"]
    connectivity[from_index, to_index] = value
    push!(strings, value)
    
end

spy(connectivity)
#savefig("spy_plot.png")