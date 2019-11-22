module Assortativity

using LightGraphs, NetworkInference

export
	# common
	AssortativityObject,
	JSON_node,
	JSON_edge,
	# measures
	assortativity,
	second_neighbour_assortativity,
	get_communities_number,
	get_modularity,
	# noise
	random_edge_rewiring,
	random_node_deletion,
	random_network,
	randomise_annotations,
	# utilities
	load_network,
	get_labels_to_groups,
	get_groups_to_indices,
	set_threshold,
	InferredNetwork_to_LightGraph,
	InferredNetwork_to_JSON,
	write_JSON_network,
	filter_connectivity

include("common.jl")
include("measures.jl")
include("noise.jl")
include("utilities.jl")

end # module
