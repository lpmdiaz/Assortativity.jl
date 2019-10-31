using CSV
using DelimitedFiles
using JSON

# make an InferredNetwork from a network text file
function load_network(filename::String; ignoreotherline = true)

    # read only every other line in case of symmetric networks
    factor = ignoreotherline == true ? 2 : 1

    # read data
    network = readdlm(filename)
    gene_names = uppercase.(string.(unique(network[:, 1:2])))
    genes_dict = Dict{String,Node}()
    for gene_name in gene_names
        genes_dict[gene_name] = Node(gene_name, Float64[], 0, Float64[])
    end

    # collect genes
    genes = collect(values(genes_dict))

    # construct a NetworkInference.Edge
    function make_edge(net, pos)
    	NetworkInference.Edge(
    		[genes_dict[ uppercase(string(net[pos,:][1])) ],
    		 genes_dict[ uppercase(string(net[pos,:][2])) ]],
    		net[pos, :][3])
    end

    # collect edges
    edges = [make_edge(network, i) for i in 1 : factor : size(network)[1]]

    InferredNetwork(genes, edges)

end

# get a dictionary of gene names to group labels
function get_genes_to_groups(genes::Array{Node}, groups_filename::String)

    # read in genes as sorted by group
    groups = CSV.read(groups_filename, delim='\t')

    # declare dictionary
    genes_to_groups = Dict{String,Symbol}()

    # fill in dictionary with genes in each gene group
    for group_label in names(groups)
        gene_names = collect(skipmissing(groups[:,group_label])) # ignore missing values
        for gene_name in gene_names
            genes_to_groups[uppercase(gene_name)] = group_label
        end
    end

    # assign any genes not in the groups file to the :Other group
    for gene in genes
        if !(uppercase(gene.label) in keys(genes_to_groups))
            genes_to_groups[uppercase(gene.label)] = :Other
        end
    end

    genes_to_groups

end

# build a dictionary of group labels (= gene symbols or attribute) to indices
function get_groups_to_indices(groups)

    # declare dictionary
    groups_to_indices = Dict{Symbol,Int}()

    # assign an index to each group
    for (i, group) in enumerate(groups)
        groups_to_indices[group] = i
    end

    groups_to_indices

end

# make a LightGraph from an InferredNetwork
function InferredNetwork_to_LightGraph(network::InferredNetwork)

    # make dictionaries to keep track of the IDs assigned to nodes
    labels_to_ids = Dict{String,Int}()
    for (i, node) in enumerate(network.nodes)
        labels_to_ids[uppercase(node.label)] = i
    end
    ids_to_labels = Dict(value => key for (key,value) in labels_to_ids) # reverse dictionary

    # make an adjacency matrix
    adjacency_matrix = zeros(Int, (length(network.nodes), length(network.nodes)))
    for edge in network.edges
        node1_id = labels_to_ids[uppercase(edge.nodes[1].label)]
        node2_id = labels_to_ids[uppercase(edge.nodes[2].label)]
        adjacency_matrix[node1_id, node2_id] += 1
        adjacency_matrix[node2_id, node1_id] += 1
    end

    # return a SimpleGraph made from the adjacency matrix, and the ids dictionary
    LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix), ids_to_labels

end

# make an InferredNetwork at a given edge threshold
function set_threshold(network::InferredNetwork, threshold::Int)

    threshold > length(network.edges) ? error("edge treshold is greater than edge number (in this case: max. $(length(network.edges)))") : nothing

    new_network_edges = network.edges[1:threshold]
    genes = unique(vcat([edge.nodes[1] for edge in new_network_edges], [edge.nodes[2] for edge in new_network_edges]))

    InferredNetwork(genes, new_network_edges)

end

# convert InferredNetwork to a JSON object that can then be exported
function InferredNetwork_to_JSON(network::InferredNetwork, genes_to_groups, groups_to_indices)

   # collect relevant data
   nodes_labels = [node.label for node in network.nodes]
   nodes_groups = [genes_to_groups[node_label] for node_label in nodes_labels]
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
function write_json_network(json_net::Dict, out_path::String; formatted=true)
    formatted ? data = json_net : data = JSON.json(json_net)
    open(out_path,"w") do f
        JSON.print(f, data, 4)
    end
end
