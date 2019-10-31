using LinearAlgebra

# """
# 	assortativity(network::InferredNetwork, threshold<:Real, groups_filename::String)
#
# description.
#
# arguments:
# - 'network': network.
# - 'threshold': the proportion of network edges to include in the analysis.
# - 'annotations_filename': path to annotations file.
# """

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

function excess_degree_assortativity(graph::SimpleGraph)

    # make a matrix containing the strength of connection between each group
    e_matrix = zeros(Δ(graph) - 1, Δ(graph) - 1)

    for edge in edges(graph)

        # retrieve the degree of the two nodes. uses the
		# excess degree, i.e. disregards the current edge
        degree1 = degree(graph,edge.src) - 1
        degree2 = degree(graph,edge.dst) - 1

		if degree1 != 0 && degree2 != 0

	        # increment connectivity between the two degrees
	        e_matrix[degree1, degree2] += 1
	        e_matrix[degree2, degree1] += 1

		end
    end

	# omitted edges
	@warn "$(Int(ne(graph) - sum(e_matrix)/2)) edges were omitted (at least one node with excess degree of 0)"

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
