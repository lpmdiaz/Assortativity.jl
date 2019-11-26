using LinearAlgebra: tr

# helper to increment a matrix given coordinates
increment_connectivity_matrix!(mat::Matrix, i::Int, j::Int) = (mat[i,j] += 1; mat[j,i] += 1; mat)

# return the assortativity coefficient of a connectivity matrix
function calculate_assortativity(e_matrix)

	# normalise matrix
	e_matrix = e_matrix / sum(e_matrix)

	# calculate assortativity
	ai = sum(e_matrix,dims=1)
	aibi_summed = sum(ai.^2) # e matrix is symmetric thus ai = bi
	r = ((ndims(e_matrix) > 1 ? tr(e_matrix) : e_matrix[1]) - aibi_summed) / (1 - aibi_summed)

	# return a coefficient of 1 if r is NaN -- happens because the e_matrix has only one entry which is 1,
	# meaning that there is only one group and all edges are between nodes of that group
	isnan(r) ? 1.0 : r

end

# format output into an AssortativityObject
make_AssortativityObject(e_matrix) = AssortativityObject(calculate_assortativity(e_matrix),e_matrix,nothing)
function make_AssortativityObject(e_matrix, labels_to_groups, groups_to_indices, ids_to_labels)

	# filter dictionaries for groups in current graph
	curr_groups = unique([labels_to_groups[value] for value in collect(values(ids_to_labels))])
	curr_groups_to_indices = collect(Iterators.flatten([filter(p -> p.first == curr_groups[i], groups_to_indices) for i in 1:length(curr_groups)])) |> Dict{Symbol, Int}

	AssortativityObject(calculate_assortativity(e_matrix),e_matrix,curr_groups_to_indices)

end

# calculate the degree assortativity of a graph
function assortativity(graph::SimpleGraph; excess_degree = false)

	# set variable to offset e_matrix size and node degrees to return excess degree assortativity
	excess_degree ? offset = 1 : offset = 0

    # make a matrix containing the strength of connection between each group
    e_matrix = zeros(Int, Δ(graph) - offset, Δ(graph) - offset)

    for edge in edges(graph)

        # retrieve the degree of the two nodes
        degree1 = degree(graph,edge.src) - offset
        degree2 = degree(graph,edge.dst) - offset

		# if calculating excess degree, check that resulting excess degree isn't null
		if degree1 != 0 && degree2 != 0

	        # increment connectivity between the two degrees
			e_matrix = increment_connectivity_matrix!(e_matrix, degree1, degree2)

		end
    end

	# send a warning with omitted edges if returning the excess degree assortativity
	if excess_degree
		(omitted_edges = Int(ne(graph) - sum(e_matrix)/2)) == 0 ? nothing : @warn "$(omitted_edges) edges were omitted out of $(ne(graph)) considered (at least one node with excess degree of 0)"
	end

	make_AssortativityObject(e_matrix)

end

# calculate the label assortativity of a graph
function assortativity(graph::SimpleGraph, labels_to_groups, groups_to_indices, ids_to_labels)

	# make a matrix containing the strength of connection between each group
    number_of_indices = length(unique(values(groups_to_indices)))
    e_matrix = zeros(Int, number_of_indices, number_of_indices)

    for edge in edges(graph)

        # recover node labels stored as ids by LightGraphs
        node1 = ids_to_labels[edge.src]
        node2 = ids_to_labels[edge.dst]

        # retrieve the corresponding groups
        group1 = labels_to_groups[node1]
        group2 = labels_to_groups[node2]

        # increment connectivity between the two groups
		e_matrix = increment_connectivity_matrix!(e_matrix,groups_to_indices[group1],groups_to_indices[group2])

    end

	make_AssortativityObject(e_matrix, labels_to_groups, groups_to_indices, ids_to_labels)

end

### dispatch assortativity on InferredNetwork type
# calculate the degree assortativity from an InferredNetwork (at threshold)
assortativity(network::InferredNetwork; excess_degree = false) = assortativity(InferredNetwork_to_LightGraph(set_threshold(network, length(network.edges)))[1], excess_degree = excess_degree)
assortativity(network::InferredNetwork, threshold::Int; excess_degree = false) = assortativity(InferredNetwork_to_LightGraph(set_threshold(network, threshold))[1], excess_degree = excess_degree)
function assortativity(network::InferredNetwork, thresholds::AbstractRange; excess_degree = false)
	assortativities = Vector{AssortativityObject}()
	for threshold in thresholds
		push!(assortativities, assortativity(network, threshold, excess_degree = excess_degree))
	end
	assortativities
end
# calculate the label assortativity from an InferredNetwork (at thresholds)
function assortativity(network::InferredNetwork, nodes_to_groups, groups_to_indices)
	graph, ids_to_labels = InferredNetwork_to_LightGraph(set_threshold(network, length(network.edges)))
	assortativity(graph, nodes_to_groups, groups_to_indices, ids_to_labels)
end
function assortativity(network::InferredNetwork, nodes_to_groups, groups_to_indices, threshold::Int)
	graph, ids_to_labels = InferredNetwork_to_LightGraph(set_threshold(network, threshold))
	assortativity(graph, nodes_to_groups, groups_to_indices, ids_to_labels)
end
function assortativity(network::InferredNetwork, nodes_to_groups, groups_to_indices, thresholds::AbstractRange)
	assortativities = Vector{AssortativityObject}()
	for threshold in thresholds
		push!(assortativities, assortativity(network, nodes_to_groups, groups_to_indices, threshold))
	end
	assortativities
end

function second_neighbour_assortativity(graph::SimpleGraph; excess_degree = false)

	# if calculating excess degree:
	# - set variable to offset e_matrix size and node degrees to return excess degree assortativity;
	# - keep track of omitted edges;
	# - count all edges traversed (both from the first or the second neighbour).
	excess_degree ? (offset = 1; omitted_edges = 0; total_edges = 0) : offset = 0

	# make a matrix containing the strength of connection between each group
	e_matrix = zeros(Int, Δ(graph) - offset, Δ(graph) - offset)

	# iterate over nodes
	for node in vertices(graph)

		# retrieve the degree of that node
		degree1 = degree(graph,node) - offset

		# iterate over the neighbours of that first node
		for neighbour_one in outneighbors(graph,node)

			# iterate over the neighbours of that second node (excluding the starting node)
			for neighbour_two in filter(x -> x ≠ node, outneighbors(graph,neighbour_one))

				# retrieve the degree of the second neighbour
				degree2 = degree(graph,neighbour_two) - offset

				# if calculating excess degree, check that resulting excess degree isn't null
				if degree1 != 0 && degree2 != 0

					# increment connectivity between the two degrees
					e_matrix = increment_connectivity_matrix!(e_matrix,degree1,degree2)

				elseif excess_degree
					omitted_edges += 1
				end
				excess_degree ? total_edges += 1 : nothing
			end
		end
	end

	# send a warning with omitted edges if returning the excess degree assortativity
	if excess_degree
		(omitted_edges)  == 0 ? nothing : @warn "$(omitted_edges) edges were omitted out of $(total_edges) considered (at least one node with excess degree of 0)"
	end

	make_AssortativityObject(e_matrix)

end

# calculate the second neighbour label assortativity
function second_neighbour_assortativity(graph::SimpleGraph, labels_to_groups, groups_to_indices, ids_to_labels)

	    # make a matrix containing the strength of connection between each group
	    number_of_indices = length(unique(values(groups_to_indices)))
	    e_matrix = zeros(Int, number_of_indices, number_of_indices)

	    # iterate over nodes
	    for node in vertices(graph)

	        # retrieve the group of current node
	        node1_id = ids_to_labels[node]
	        group1 = labels_to_groups[node1_id]

	        # iterate over the neighbours of that first node
	        for neighbour_one in outneighbors(graph,node)

	            # iterate over the neighbours of that second node (excluding the starting node)
	            for neighbour_two in filter(x -> x ≠ node, outneighbors(graph,neighbour_one))

	                # retrieve the group of the second neighbour
	                node2_id = ids_to_labels[neighbour_two]
	                group2 = labels_to_groups[node2_id]

	                # increment connectivity between the two groups
					e_matrix = increment_connectivity_matrix!(e_matrix,groups_to_indices[group1],groups_to_indices[group2])

	            end
	        end
	    end

	make_AssortativityObject(e_matrix, labels_to_groups, groups_to_indices, ids_to_labels)

end

### dispatch second_neighbour_assortativity on InferredNetwork type
# calculate the second neighbour degree assortativity from an InferredNetwork (at threshold)
second_neighbour_assortativity(network::InferredNetwork; excess_degree = false) = second_neighbour_assortativity(InferredNetwork_to_LightGraph(set_threshold(network, length(network.edges)))[1], excess_degree = excess_degree)
second_neighbour_assortativity(network::InferredNetwork, threshold::Int; excess_degree = false) = second_neighbour_assortativity(InferredNetwork_to_LightGraph(set_threshold(network, threshold))[1], excess_degree = excess_degree)
function second_neighbour_assortativity(network::InferredNetwork, thresholds::AbstractRange; excess_degree = false)
	assortativities = Vector{AssortativityObject}()
	for threshold in thresholds
		push!(assortativities, second_neighbour_assortativity(network, threshold, excess_degree = excess_degree))
	end
	assortativities
end
# calculate the second neighbour label assortativity from an InferredNetwork (at thresholds)
function second_neighbour_assortativity(network::InferredNetwork, nodes_to_groups, groups_to_indices)
	graph, ids_to_labels = InferredNetwork_to_LightGraph(set_threshold(network, length(network.edges)))
	second_neighbour_assortativity(graph, nodes_to_groups, groups_to_indices, ids_to_labels)
end
function second_neighbour_assortativity(network::InferredNetwork, nodes_to_groups, groups_to_indices, threshold::Int)
	graph, ids_to_labels = InferredNetwork_to_LightGraph(set_threshold(network, threshold))
	second_neighbour_assortativity(graph, nodes_to_groups, groups_to_indices, ids_to_labels)
end
function second_neighbour_assortativity(network::InferredNetwork, nodes_to_groups, groups_to_indices, thresholds::AbstractRange)
	assortativities = Vector{AssortativityObject}()
	for threshold in thresholds
		push!(assortativities, second_neighbour_assortativity(network, nodes_to_groups, groups_to_indices, threshold))
	end
	assortativities
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
