using CSV
using CSVFiles
using DataFrames
using Plots


include("data_preprocessing.jl")


#%% Data exploration 
"Plot the number of connections of the synaptic and gap"

# Create two matrices for storing the connectivity information
synaptic_number = zeros(302, 302)           # To store the number of connections
synaptic_connections = zeros(302, 302)      # To store if there is a connection
# Iterate through the datatframe and add the information to the corresponding point
for link in eachrow(data_connect_synaptic)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    value = link[:Number]
    synaptic_number[from_index, to_index] = value
    synaptic_connections[from_index, to_index] = 1
end
# Plot both matrices
spy(synaptic_number, plot_title= "Number of synaptic connections", xlabel = "Sending neuron index", ylabel = "Rceiving neuron index", color=:berlin)
savefig(plotsdir("num_syn_conn.png"))
spy(synaptic_connections, plot_title= "Synaptic connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")
savefig(plotsdir("syn_conn.png"))

#%% For viewing the information of the frequency of the number of connections
data_s = vec(synaptic_number)           # Flatten the matrix into a 1D array
data_s = filter(x -> x != 0, data_s)    # Take out the values that are 0 for a true result
# Generate histogram
histogram(data_s, bins=:scott, 
    xticks = unique(data_s), legend = false, normalize=:probability,
    xlabel = "Number of connections between two neurons", ylabel = "Frequency", title = "Synaptic frequency histogram")
histogram!(size=(760,500))
savefig(plotsdir("syn_freq_histogram.png"))


# Repeat the previous process for GAP junctions
gap_number = zeros(302, 302)
gap_connections = zeros(302,302)
for link in eachrow(data_connect_gap)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    value = link[:Number]
    gap_number[from_index, to_index] = value
    gap_connections[from_index, to_index] = 1
end

spy(gap_number, plot_title= "Number of gap connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index", color=:berlin)
savefig(plotsdir("num_gap_conn.png"))
spy(gap_connections, plot_title= "Gap connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")
savefig(plotsdir("gap_conn.png"))


data = vec(gap_number)
data = filter(x -> x != 0, data)
# Generate histogram
histogram(data, bins=:scott, 
    xticks = unique(data), legend = false, normalize=:probability,
    xlabel = "Number of connections between two neurons", ylabel = "Frequency", title = "Gap frecuency histogram")
histogram!(size=(800,500))
savefig(plotsdir("gap_freq_histogram.png"))



"Plot the connections by Monoamines"
# Create a matrix for each type
mono_connections_tyr = zeros(302,302)
mono_connections_oct = zeros(302,302)
mono_connections_dop = zeros(302,302)
mono_connections_ser = zeros(302,302)
# Add each connection to the corresponding matrix
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
# Plot all in the same graph
spy(mono_connections_tyr, color = :orange, plot_title= "Monoamine connections among neurons", xlabel = "Sending neuron index", ylabel = "Receiving neuron index")
spy!(mono_connections_oct, color = :red)
spy!(mono_connections_dop, color = :green)
spy!(mono_connections_ser, color = :blue,legend=:outertopright)
savefig(plotsdir("mono_conn.png"))




"Plot the connections by Neuropeptides"
# Create one matrix for each column of the dataframe
neuropep_connections_Type1 = zeros(302,302)
neuropep_connections_Type2 = zeros(302,302)
# Append the value to the matrix
for link in eachrow(data_connect_neuropep)
    from_index = link[:IndexSending]
    to_index = link[:IndexReceiving]
    if link[:Type1] !== nothing
        neuropep_connections_Type1[from_index, to_index] = 1
    end
    if link[:Type2] !== nothing
        neuropep_connections_Type2[from_index, to_index] = 1
    end
end
spy(neuropep_connections_Type1, color = :blue, plot_title= "Neuropeptides connections among neurons", xlabel = "Neuron index", ylabel = "Neuron index")
spy!(neuropep_connections_Type2, color = :red)
savefig(plotsdir("neuropep_conn.png"))
