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

# function assortativity(network::InferredNetwork, threshold::Int, annotations_filename::String)

function assortativity(graph::SimpleGraph, genes_to_groups, groups_to_indices, ids_to_labels)

	# make a matrix containing the strength of connection between each group
    number_of_indices = length(unique(values(groups_to_indices)))
    e_matrix = zeros(number_of_indices, number_of_indices)

    for edge in edges(graph)

        # recover node labels stored as ids by LightGraphs
        gene1 = ids_to_labels[edge.src]
        gene2 = ids_to_labels[edge.dst]

        # get the corresponding group of each label
        group1 = genes_to_groups[gene1]
        group2 = genes_to_groups[gene2]

        # increment connectivity between these two groups
        e_matrix[groups_to_indices[group1], groups_to_indices[group2]] += 1
        e_matrix[groups_to_indices[group2], groups_to_indices[group1]] += 1

    end

	# normalise matrix
    e_matrix = e_matrix / sum(e_matrix)

    # calculate and return assortativity
    ai = sum(e_matrix,dims=1)
    aibi_summed = sum(ai.^2) # e matrix is symmetric, thus ai = bi
    r = (tr(e_matrix) - aibi_summed) / (1 - aibi_summed)
    isnan(r) ? .0 : r

end
