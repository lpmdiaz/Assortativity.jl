using LinearAlgebra

function calculate_assortativity(connectivity_matrix)

	# normalise matrix
	connectivity_matrix = connectivity_matrix / sum(connectivity_matrix)

	# calculate assortativity
	ai = sum(connectivity_matrix,dims=1)
	aibi_summed = sum(ai.^2) # e matrix is symmetric thus ai = bi
	r = (tr(connectivity_matrix) - aibi_summed) / (1 - aibi_summed)

	isnan(r) ? .0 : r, connectivity_matrix

end

function assortativity(graph::SimpleGraph)

    # make a matrix containing the strength of connection between each group
    e_matrix = zeros(Δ(graph), Δ(graph))

    for edge in edges(graph)

        # retrieve the degree of the two nodes
        degree1 = degree(graph,edge.src)
        degree2 = degree(graph,edge.dst)

        # increment connectivity between the two degrees
        e_matrix[degree1, degree2] += 1
        e_matrix[degree2, degree1] += 1

    end

	calculate_assortativity(e_matrix)

end

function assortativity(graph::SimpleGraph, genes_to_groups, groups_to_indices, ids_to_labels)

	# make a matrix containing the strength of connection between each group
    number_of_indices = length(unique(values(groups_to_indices)))
    e_matrix = zeros(number_of_indices, number_of_indices)

    for edge in edges(graph)

        # recover node labels stored as ids by LightGraphs
        node1 = ids_to_labels[edge.src]
        node2 = ids_to_labels[edge.dst]

        # retrieve the corresponding groups
        group1 = genes_to_groups[node1]
        group2 = genes_to_groups[node2]

        # increment connectivity between the two groups
        e_matrix[groups_to_indices[group1], groups_to_indices[group2]] += 1
        e_matrix[groups_to_indices[group2], groups_to_indices[group1]] += 1

    end

	calculate_assortativity(e_matrix)

end

# degree assortativity from InferredNetwork (at threshold)
function assortativity(network::InferredNetwork)
end
function assortativity(network::InferredNetwork, threshold::Int)
end
function assortativity(network::InferredNetwork, thresholds::AbstractRange)
end

# label assortativity from InferredNetwork (at threshold)
function assortativity(network::InferredNetwork, nodes_to_groups)
end
function assortativity(network::InferredNetwork, threshold::Int, nodes_to_groups)
end
function assortativity(network::InferredNetwork, thresholds::AbstractRange, nodes_to_groups)
end

function excess_degree_assortativity(graph::SimpleGraph)

    # make a matrix containing the strength of connection between each group
    e_matrix = zeros(Δ(graph) - 1, Δ(graph) - 1)

    for edge in edges(graph)

        # retrieve the degree of the two nodes. uses the
		# excess degree, i.e. disregards the current edge
        degree1 = degree(graph,edge.src) - 1
        degree2 = degree(graph,edge.dst) - 1

		# check that excess degree isn't null
		if degree1 != 0 && degree2 != 0

	        # increment connectivity between the two degrees
	        e_matrix[degree1, degree2] += 1
	        e_matrix[degree2, degree1] += 1

		end
    end

	# omitted edges
	(omitted_edges = Int(ne(graph) - sum(e_matrix)/2)) == 0 ? nothing : @warn "$(omitted_edges) edges were omitted out of $(ne(graph)) considered (at least one node with excess degree of 0)"

	calculate_assortativity(e_matrix)

end

function second_neighbour_assortativity(graph::SimpleGraph)

	# make a matrix containing the strength of connection between each group
	e_matrix = zeros(Δ(graph), Δ(graph))

	# iterate over nodes
	for node in vertices(graph)

		# retrieve the degree of that node
		degree1 = degree(graph,node)

		# iterate over the neighbours of that first node
		for one_walk in outneighbors(graph,node)

			# iterate over the neighbours of that second node (excluding the starting node)
			for two_walk in filter(x -> x ≠ node, outneighbors(graph,one_walk))

				# retrieve the degree of the second neighbour
				degree2 = degree(graph,two_walk)

				# increment connectivity between the two degrees
				e_matrix[degree1, degree2] += 1
				e_matrix[degree2, degree1] += 1

			end
		end
	end

	calculate_assortativity(e_matrix)

end

function second_neighbour_assortativity(graph::SimpleGraph, genes_to_groups, groups_to_indices, ids_to_labels)

	    # make a matrix containing the strength of connection between each group
	    number_of_indices = length(unique(values(groups_to_indices)))
	    e_matrix = zeros(number_of_indices, number_of_indices)

	    # iterate over nodes
	    for node in vertices(graph)

	        # retrieve the group of current node
	        node1_id = ids_to_labels[node]
	        group1 = genes_to_groups[node1_id]

	        # iterate over the neighbours of that first node
	        for one_walk in outneighbors(graph,node)

	            # iterate over the neighbours of that second node (excluding the starting node)
	            for two_walk in filter(x -> x ≠ node, outneighbors(graph,one_walk))

	                # retrieve the group of the second neighbour
	                node2_id = ids_to_labels[two_walk]
	                group2 = genes_to_groups[node2_id]

	                # increment connectivity between the two groups
	                e_matrix[groups_to_indices[group1], groups_to_indices[group2]] += 1
	                e_matrix[groups_to_indices[group2], groups_to_indices[group1]] += 1

	            end
	        end
	    end

	calculate_assortativity(e_matrix)

end

function excess_degree_second_neighbour_assortativity(graph::SimpleGraph)

	# make a matrix containing the strength of connection between each group
	e_matrix = zeros(Δ(graph) - 1, Δ(graph) - 1)

	# keep track of omitted edges
	omitted_edges = 0

	# count edges
	total_edges = 0

	# iterate over nodes
	for node in vertices(graph)

		# retrieve the degree of that node
		degree1 = degree(graph,node) - 1

		# iterate over the neighbours of that first node
		for one_walk in outneighbors(graph,node)

			# iterate over the neighbours of that second node (excluding the starting node)
			for two_walk in filter(x -> x ≠ node, outneighbors(graph,one_walk))

				# retrieve the degree of the second neighbour
				degree2 = degree(graph,two_walk) - 1

				# check that excess degree isn't null
				if degree1 != 0 && degree2 != 0

					# increment connectivity between the two degrees
					e_matrix[degree1, degree2] += 1
					e_matrix[degree2, degree1] += 1

				else
					omitted_edges += 1
				end
				total_edges += 1
			end
		end
	end

	# omitted edges
	(omitted_edges)  == 0 ? nothing : @warn "$(omitted_edges) edges were omitted out of $(total_edges) considered (at least one node with excess degree of 0)"

	calculate_assortativity(e_matrix)

end

# second neighbour degree assortativity from InferredNetwork (at threshold)
function second_neighbour_assortativity(network::InferredNetwork)
end
function second_neighbour_assortativity(network::InferredNetwork, threshold::Int)
end
function second_neighbour_assortativity(network::InferredNetwork, thresholds::AbstractRange)
end

# second neighbour label assortativity from InferredNetwork (at threshold)
function second_neighbour_assortativity(network::InferredNetwork, nodes_to_groups)
end
function second_neighbour_assortativity(network::InferredNetwork, threshold::Int, nodes_to_groups)
end
function second_neighbour_assortativity(network::InferredNetwork, thresholds::AbstractRange, nodes_to_groups)
end

# retrieve community number and communities using the label propagation algorithm
function get_communities(graph::SimpleGraph)
	communities, convergence = label_propagation(graph)
	communities_number = 0
	for i in 1:maximum(communities)
		if length(findall(x -> x == i, communities)) >= 1
			communities_number += 1
		end
	end
	communities_number, communities
end

# LightGraphs syntax helpers
get_communities_number(graph::SimpleGraph) = get_communities(graph)[1]
get_modularity(graph::SimpleGraph) = modularity(graph, get_communities(graph)[2])
