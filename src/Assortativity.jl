module Assortativity

using LightGraphs
using NetworkInference

export
	# common
	JSON_node,
	JSON_edge,
	# measures
	assortativity,
	second_neighbour_assortativity,
	excess_degree_assortativity,
	excess_degree_second_neighbour_assortativity,
	get_communities_number,
	get_modularity,
	# noise
	random_edge_rewiring,
	random_network,
	random_node_deletion,
	randomise_annotations,
	# utilities
	load_network,
	get_genes_to_groups,
	get_groups_to_indices,
	set_threshold,
	InferredNetwork_to_LightGraph,
	InferredNetwork_to_JSON,
	write_json_network

include("common.jl")
include("measures.jl")
include("noise.jl")
include("utilities.jl")

end # module
