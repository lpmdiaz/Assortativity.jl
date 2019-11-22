using Combinatorics: combinations

# introduce noise in networks for a specified number of edges
function random_edge_rewiring(network::InferredNetwork, number_of_random::Int)

	# retrieve positions of all unique edges in an adjacency matrix
	function get_edges_indices(adjacency_matrix::Array{Int64,2})

		indices = Array{Array{Int64,2},1}()

		for x in 1:size(adjacency_matrix)[1]
			for y in 1:x-1 # -1 to ignore the diagonal
				if adjacency_matrix[x,y] == 1
					push!(indices, [x y])
				end
			end
		end

		indices

	end

	# update adjacency matrix with given coordinates
	function update_adjacency_matrix!(adjacency_matrix::Array{Int,2}, direction::Int, node1x::Int, node1y::Int, node2x::Int, node2y::Int)

		# decrement adjacency matrix to remove edges being swapped
		adjacency_matrix[node1x, node1y] -= 1
		adjacency_matrix[node1y, node1x] -= 1
		adjacency_matrix[node2x, node2y] -= 1
		adjacency_matrix[node2y, node2x] -= 1

		# increment adjacency matrix to add the new edges
		if direction == 1 # cross swap
			adjacency_matrix[node1x, node2y] += 1
			adjacency_matrix[node2y, node1x] += 1
			adjacency_matrix[node2x, node1y] += 1
			adjacency_matrix[node1y, node2x] += 1
		elseif direction == 2 # vertical swap
			adjacency_matrix[node1x, node2x] += 1
			adjacency_matrix[node2x, node1x] += 1
			adjacency_matrix[node1y, node2y] += 1
			adjacency_matrix[node2y, node1y] += 1
		end

	end

	# initialise var to count total number of rewired edges (should end up being number_of_random)
	rewired_edges = 0

	# start as for converting an InferredNetwork to a LightGraph --
	# make dictionaries to keep track of the IDs assigned to nodes
	labels_to_ids = Dict{String,Int}()
	for (i, node) in enumerate(network.nodes)
		labels_to_ids[node.label] = i
	end
	ids_to_labels = Dict(value => key for (key,value) in labels_to_ids) # reverse dictionary

	# make an adjacency matrix
	adjacency_matrix = zeros(Int, (length(network.nodes), length(network.nodes)))
	for edge in network.edges
		label1id = labels_to_ids[edge.nodes[1].label]
		label2id = labels_to_ids[edge.nodes[2].label]
		adjacency_matrix[label1id, label2id] += 1
		adjacency_matrix[label2id, label1id] += 1
	end

	# define function return behaviour
	function final_output()
		LightGraphs.SimpleGraphs.SimpleGraph(adjacency_matrix), ids_to_labels, rewired_edges
	end

	# collect the position of edges in the adjacency matrix
	indices = get_edges_indices(adjacency_matrix)

	# check that edges can possibly be rewired; if true, mo two edges can be rewired in the network e.g. complete network
	if length(indices) > ((((size(adjacency_matrix)[1] * size(adjacency_matrix)[2]) - size(adjacency_matrix)[1]) / 2) - 2)
		return final_output()
	end

	# retrieve the number of edges in the network
	length(indices) != length(network.edges) ? error("non symmetric network") : number_of_edges = length(network.edges)

	# select two edges at random
	function select_two(number_of_edges::Int64)
		edge1 = rand(1:number_of_edges)
		edge2 = rand(filter(x -> x ≠ edge1, collect(1:number_of_edges))) # avoid selecting the same node twice
		edge1, edge2
	end

	# retrieve the degree sequence of the original adjacency matrix
	orig_degree_seq = sum(adjacency_matrix,dims=1)

	# iterate over number of edges that need to be rewired
	for iteration in 1:number_of_random

		# initialise variables (direction of edge rewiring, node coordinates, exit loop control)
		direction, node1x, node1y, node2x, node2y, control = 0, 0, 0, 0, 0, 0

		# pick edges and check if suitable for a swap
		while direction == 0

			# select edges
			edge1, edge2 = select_two(number_of_edges)

			# retrieve corresponding nodes
			node1x = indices[edge1][1]
			node1y = indices[edge1][2]
			node2x = indices[edge2][1]
			node2y = indices[edge2][2]

			# check if first direction is suitable (cross swap)
			if node1x != node2y && node2x != node1y # avoid the diagonal
				if adjacency_matrix[node1x, node2y] == 0 && adjacency_matrix[node2x, node1y] == 0
					direction = 1
				end
			end

			# check if second direction is suitable (vertical swap)
			if node1x != node2x && node1y != node2y # avoid the diagonal
				if adjacency_matrix[node1x, node2x] == 0 && adjacency_matrix[node1y, node2y] == 0
					direction = 2
				end
			end

			# if can't exit loop
			if control > 1000
				println("loop stuck -- exiting")
				return final_output()
			end

			# increment loop control
			control += 1

		end

		# swap the selected edges
		update_adjacency_matrix!(adjacency_matrix, direction, node1x, node1y, node2x, node2y)

		# increment number of rewired edges
		rewired_edges += 2

		# check that swap kept the same degree sequence
		if sum(adjacency_matrix,dims=1) != orig_degree_seq
			error("incorrect degree sequence")
		end

		# collect the positions of edges in the new adjacency matrix
		indices = get_edges_indices(adjacency_matrix)

	end

	final_output()

end

# introduce noise in networks by deleting a specified number of nodes and their edges
function random_node_deletion(network::InferredNetwork, number_of_delete::Int64)

	if number_of_delete >= length(network.nodes)
		error("cannot delete more nodes than present in the network")
	end

	# copy nodes and edges from input network
	new_nodes = deepcopy(network.nodes)
	new_edges = deepcopy(network.edges)

	# remove the desired number of nodes
	for i in 1:number_of_delete

		number_of_nodes = length(new_nodes)
		del_index = rand(1:number_of_nodes)
		del_node = new_nodes[del_index]
		deleteat!(new_nodes, del_index)

		# collect all node labels
		all_node_labels = [node.label for node in new_nodes]

		# collect remaning edges after node deletion
		src_labels = [node[1].label for node in [edge.nodes for edge in new_edges]]
		dst_labels = [node[2].label for node in [edge.nodes for edge in new_edges]]
		srclabelsin = collect(label in all_node_labels for label in (src_labels))
		dstlabelsin = collect(label in all_node_labels for label in (dst_labels))
		bothlabelsin = collect(srclabelsin .== dstlabelsin .== true)
		new_edges = new_edges[bothlabelsin]

	end

	# make a LightGraph from an InferredNetwork
	InferredNetwork_to_LightGraph(InferredNetwork(new_nodes, new_edges))

end

# creates a random graph from the given input network
function random_network(network::InferredNetwork; strict_rand = true)

	# number of edges in input network
	number_of_edges = length(network.edges)

	# declare array to store edges
	new_edges = Array{NetworkInference.Edge}(undef,number_of_edges)

	if strict_rand

		# make a dictionary of all possible edges
		possible_edges = collect(combinations(deepcopy(network.nodes), 2))

		# collect edges
		for i in 1:number_of_edges

			# select one possible edge at random
			new_edge_index = rand(1:length(possible_edges))

			# construct edge from the selected two nodes
			new_edges[i] = NetworkInference.Edge([possible_edges[new_edge_index][1], possible_edges[new_edge_index][2]], 1.0)

			# remove the selected edge from the list of possible edges
			deleteat!(possible_edges, new_edge_index)

		end

	else

		# construct edge from two nodes at random
		make_random_edge(nodes::Array{NetworkInference.Node,1}) = NetworkInference.Edge([nodes[rand(1:length(nodes))], nodes[rand(1:length(nodes))]], 1.0)

		# make edges
		for i in 1:number_of_edges
			new_edges[i] = make_random_edge(network.nodes)
		end

	end

	# make a LightGraph from an InferredNetwork
	InferredNetwork_to_LightGraph(InferredNetwork(network.nodes, new_edges))

end

# introduce noise in group annotations dictionary
function randomise_annotations(in_dict::Dict, number_of_random::Int; behaviour = 5)

	# count how many keys / values are different between two dictionaries
	function get_number_differences(old_dict, new_dict)
		a = convert(Vector{Bool}, collect(keys(old_dict)) .!= collect(keys(new_dict)))
		b = convert(Vector{Bool}, collect(values(old_dict)) .!= collect(values(new_dict)))
		union(collect(1:length(old_dict))[a], collect(1:length(old_dict))[b])
	end

	# pick a random key and attribute it a random value
	function relaxed_randomising()
		for i in 1:number_of_random
			dict[rand(collect(keys(dict)))] = rand(collect(values(dict)))
		end
		dict
	end

	# pick two keys at random and swap their values
	function swapped_randomising()
		for i in 1:number_of_random
			key1 = rand(collect(keys(dict)))
			key2 = rand(collect(keys(dict)))
			value1 = dict[key1]
			value2 = dict[key2]
			dict[key1] = value2
			dict[key2] = value1
		end
		dict
	end

	# pick a key at random and randomly assign it a new value that is not the one it already had
	function strict_randomising()
		for i in 1:number_of_random
			key = rand(collect(keys(dict)))
			new_value = rand(filter(x -> x ≠ dict[key], collect(values(dict))))
			dict[key] = new_value
		end
		dict
	end

	# pick a key at random and randomly assign it a new value that is not the one it originally had
	function stricter_randomising()
		for i in 1:number_of_random
			key = rand(collect(keys(dict)))
			new_value = rand(filter(x -> x ≠ in_dict[key], collect(values(dict))))
			dict[key] = new_value
		end
		dict
	end

	# make sure dictionary is changes as many times as number_of_random
	function strictest_randomising()
		while length(get_number_differences(in_dict, dict)) != number_of_random
			key = rand(collect(keys(dict)))
			new_value = rand(filter(x -> x ≠ in_dict[key], collect(values(dict))))
			dict[key] = new_value
		end
		dict
	end

	# prevent modification of the input dictionary
	dict = deepcopy(in_dict)

	# choose algorithm behaviour
	if behaviour == 1
		relaxed_randomising()
	elseif behaviour == 2
		swapped_randomising()
	elseif behaviour == 3
		strict_randomising()
	elseif behaviour == 4
		stricter_randomising()
	elseif behaviour == 5
		number_of_random <= length(dict) ? strictest_randomising() : error("number of random exceeds dictionary length")
	end

	dict, length(get_number_differences(in_dict, dict))

end
