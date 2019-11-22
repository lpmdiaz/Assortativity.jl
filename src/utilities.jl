using CSV: read
using DelimitedFiles: readdlm
using JSON

# make an InferredNetwork from a network text file
function load_network(filename::String; ignoreotherline = true)

    # read only every other line in case of symmetric networks
    factor = ignoreotherline == true ? 2 : 1

    # read data
    network = readdlm(filename)
    labels = string.(unique(network[:, 1:2]))
    nodes_dict = Dict{String,Node}()
    for label in labels
        nodes_dict[label] = Node(label, Float64[], 0, Float64[])
    end

    # construct a NetworkInference.Edge
    function make_edge(net, pos)
    	NetworkInference.Edge(
    		[nodes_dict[ string(net[pos,:][1]) ],
    		 nodes_dict[ string(net[pos,:][2]) ]],
    		net[pos, :][3])
    end

    # collect nodes and edges
	nodes = collect(values(nodes_dict))
    edges = [make_edge(network, i) for i in 1 : factor : size(network)[1]]

    InferredNetwork(nodes, edges)

end

# get a dictionary of nodes labels names to group labels
function get_labels_to_groups(nodes::Array{Node}, groups_filename::String)

    # read in node labels as sorted by group
    groups = read(groups_filename, delim='\t')

    # declare dictionary
    labels_to_groups = Dict{String,Symbol}()

    # fill in dictionary with labels in each group
    for group in names(groups)
        labels = collect(skipmissing(groups[:,group])) # ignore missing values
        for label in labels
            labels_to_groups[label] = group
        end
    end

    # assign any labels not in the groups file to the :Other group
    for node in nodes
        if !(node.label in keys(labels_to_groups))
            labels_to_groups[node.label] = :Other
        end
    end

    labels_to_groups

end

# build a dictionary of group labels to indices
function get_groups_to_indices(groups)

    # declare dictionary
    groups_to_indices = Dict{Symbol,Int}()

    # assign an index to each group
	[groups_to_indices[group] = i for (i, group) in enumerate(sort(unique(values(groups))))]

    groups_to_indices

end

# make an InferredNetwork at a given edge threshold
function set_threshold(network::InferredNetwork, threshold::Int)

    threshold > length(network.edges) ? error("edge treshold is greater than edge number (in this case: max. $(length(network.edges)))") : nothing

    new_network_edges = network.edges[1:threshold]
    nodes = unique(vcat([edge.nodes[1] for edge in new_network_edges], [edge.nodes[2] for edge in new_network_edges]))

    InferredNetwork(nodes, new_network_edges)

end

# make a LightGraph from an InferredNetwork
function InferredNetwork_to_LightGraph(network::InferredNetwork)

    # make dictionaries to keep track of the IDs assigned to nodes
    labels_to_ids = Dict{String,Int}()
    for (i, node) in enumerate(network.nodes)
        labels_to_ids[node.label] = i
    end
    ids_to_labels = Dict(value => key for (key, value) in labels_to_ids) # reverse dictionary

    # make an adjacency matrix
    adjacency_matrix = zeros(Int, (length(network.nodes), length(network.nodes)))
    for edge in network.edges
        node1_id = labels_to_ids[edge.nodes[1].label]
        node2_id = labels_to_ids[edge.nodes[2].label]
        adjacency_matrix[node1_id, node2_id] += 1
        adjacency_matrix[node2_id, node1_id] += 1
    end

    # return a SimpleGraph made from the adjacency matrix, and the ids dictionary
    LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix), ids_to_labels

end

# convert InferredNetwork to a JSON object that can then be exported
function InferredNetwork_to_JSON(network::InferredNetwork, labels_to_groups, groups_to_indices)

   # collect relevant data
   nodes_labels = [node.label for node in network.nodes]
   nodes_groups = [labels_to_groups[node_label] for node_label in nodes_labels]
   nodes_groups_indices = [groups_to_indices[node_group] for node_group in nodes_groups]
   edges_sources = [node.label for node in [edge.nodes[1] for edge in network.edges]]
   edges_destinations = [node.label for node in [edge.nodes[2] for edge in network.edges]]
   edges_weight = [edge.weight for edge in network.edges]

   # make nodes and edges
   JSON_nodes = [JSON_node(nodes_labels[i], nodes_groups[i], nodes_groups_indices[i]) for i in 1:length(network.nodes)]
   JSON_edges = [JSON_edge(edges_sources[i], edges_destinations[i], edges_weight[i]) for i in 1:length(network.edges)]

   # return dictionary of nodes and edges
   Dict("nodes" => JSON_nodes, "edges" => JSON_edges)

end

# export JSON network. use formatted to have tabs
function write_JSON_network(JSON_net::Dict, out_path::String; formatted = true)
    formatted ? data = JSON_net : data = JSON.json(JSON_net)
    open(out_path,"w") do f
        JSON.print(f, data, 4)
    end
end

# filter connectivity of an AssortativityObject to remove empty columns and rows
function filter_connectivity(obj::AssortativityObject)
	if typeof(obj.groups) <: Dict # ignore degree assortativity
		sorted_groups = sort(obj.groups)
		indices = collect(values(sorted_groups))
		new_connectivity = obj.connectivity[indices,indices]
		new_groups = Dict(keys(sorted_groups) .=> 1:length(indices))
		AssortativityObject(obj.value, new_connectivity, new_groups)
	else
		obj
	end
end
