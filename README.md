# Assortativity.jl

## Description

This package implements the calculation of the assortativity coefficient for inferred networks as used in a forthcoming paper (*Gaining confidence in inferred networks*). The assortativity coefficient was introduced in [[1]](#references) and represents the tendency of nodes with similar properties to be connected by an edge. These properties are freely defined, with the most simple being the node degree (this then yields degree assortativity, answering the question: do nodes that have the same degree tend to be more often connected to each other than to nodes with different degrees?).

This package is based on code originally written by [Thalia Chan](https://github.com/Tchanders). It supports **undirected** networks and primarily relies on [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [NetworkInference.jl](https://github.com/Tchanders/NetworkInference.jl) for inferring and storing networks.

The approach taken to look at the assortativity coefficient is to gradually increase the number of edges in a network, i.e. the edge threshold, starting from 1, and to calculate the assortativity coefficient at each step. More details can be found in the supporting publication. The package also includes additional functionalities not discussed in the publication (e.g. second neighbour assortativity coefficient).

This package was originally written for gene regulatory networks (node labels then refer to gene names and groups to e.g. biological function) but can work with any type of network.

### Prerequisites

The following packages need to be installed: `NetworkInference`, `LightGraphs`, `Combinatorics`, `CSV`, `JSON`. If using correlation for inferring networks, please use [this fork](https://github.com/lpmdiaz/NetworkInference.jl) of `NetworkInference.jl` (i.e. run `using Pkg; Pkg.add("https://github.com/lpmdiaz/NetworkInference.jl")`).

Using this package requires the source data (e.g. gene expression values) to be in the following format (with the first line being disregarded):

| nodes_labels | condition1 | condition2 | ... |
|:------------:|:----------:| :---------:|:---:|
| node_label1  |   value1   |   value2   | ... |
| node_label2  |   value1   |   value2   | ... |
|      ...     |     ...    |     ...    | ... |

If calculating label assortativity i.e. assortativity based on custom labelling of the nodes, a file that links each node label (e.g. gene name) to a property of that node (e.g. biological function) is required. For clarity, these properties are sorted in 'groups', and label assortativity thus refers to the assortativity coefficient calculated from the attachment patterns of nodes belonging to these groups with custom labels. This file should be tab delimited and in the following format:

|       group1       |       group2       |       group3       | ... |
|:------------------:|:------------------:|:------------------:|:---:|
| group1_node_label1 | group2_node_label1 | group3_node_label1 | ... |
| group1_node_label2 |        n/a         | group3_node_label2 | ... |
|        ...         |        n/a         |        ...         | ... |

This means all node labels associated with a given group should be in the same column as that group and that columns can have different lengths depending on the number of node labels associated with each group. Note that the node labels (e.g. the gene names) in both files must be consistent (e.g. code is case sensitive).

## Installation

`using Pkg; Pkg.add("https://github.com/lpmdiaz/Assortativity.jl")`

## Functionalities

### Network inference

Network inference is supported via integration with `NetworkInference.jl`. As described in this package, the inference algorithms currently implemented are MI, CLR, PUC, PIDC, as well as two types of correlation (Pearson and Spearman, both signed and unsigned) if using [this fork](https://github.com/lpmdiaz/NetworkInference.jl)⁠—see usage [details](https://github.com/lpmdiaz/NetworkInference.jl#inference-algorithms-currently-implemented) for how to use them and the [reference publication](https://github.com/lpmdiaz/NetworkInference.jl#references) for details on the first four algorithms.

### Assortativity and other  graph measures

The main aim of this package is to calculate the assortativity coefficient; two flavours of this measure are currently implemented: degree assortativity and label assortativity (the latter only if provided a file containing node annotations as described in the [Prerequisites](#prerequisites) section).

Degree assortativity is also available as excess degree assortativity as originally used by Newman [[1]](#references)—the excess degree of a given node is its degree minus one as this measure only considers edges other than the one used to reach that node.

Both flavours of assortativity are also available in a second neighbour implementation⁠—for a given node, all second neighbours of this node i.e. all nodes that can be reached by traversing two edges (excluding the walk back to the original node) are considered to calculate the assortativity coefficient. This approach iterate over nodes and thus has a different behaviour from the regular assortativity that iterates over edges.

Other measures are also available via `LightGraphs.jl` to inspect relevant graph properties and compare them to the assortativity coefficient: the clustering coefficient, the number of communities, the graph modularity, and the centrality coefficient.

### Network and graph utilities

Previously inferred networks can be *imported*, and all networks can also be *exported* to file as undirected, symmetric networks.

This package allows to *convert* an `InferredNetwork` from the `NetworkInference` package to a `SimpleGraph` from the `LightGraphs` package. This is done for simple calculation of the assortativity coefficient. Note that this conversion loses track of the edges order i.e. edges in the `SimpleGraph` have no associated score or weight (setting an edge threshold must be done on an `InferredNetwork` since it does rank edges based on the score given by the inference algorithm used).

Future additions include plotting utilities, as well as graph manipulation utilities (e.g. the number of edges shared between networks, an implementation of the Borda count election approach to merge networks, etc.).

### Noise \& randomness

Four ways to introduce noise in networks are implemented in order to study the behaviour of the assortativity coefficient and to know what is expected in the random case.

`random_edge_rewiring` allows to introduce noise in edges; given a number `n` of edges to be rewired, will pick `2n` edges and swap the nodes they are connected to.

`random_node_deletion` will randomly delete the given number of nodes from the network.

`random_network` will return a graph with as many nodes and edges than the network given as an input but with edges populated at random.

`randomise_annotations` introduces noise in the group dictionary; it has five different behaviours:
- (1) relaxed — fully random and thus does not prevent attributing a key the same value it already had after a few iterations.
- (2) swapped ⁠— swaps two random keys' values with the possibility that the swapped values are the same than before.
- (3) strict ⁠— changes a random key's value to one that must be different from the one it previously had. However, several runs could revert back to the original value.
- (4) stricter ⁠— changes a random key's value to one that must be different from the one it originally had. This behaviour prevents reverting back to the original value of the keys. Will result in a number of changes in the dictionary that is close to the number of iterations.
- (5) strictest (**default**) ⁠— same as (4) but will always result in as many changes in the dictionary as iterations.

## Example usage and data

All basic usage is described in the `basic_usage.ipynb` Jupyter notebook within the `examples` directory.

The dataset used in the supporting publication is from [[2]](#references) and is included in the `examples/data` subdirectory. It has been filtered for genes expressed in less than 80% cells and is also available as two smaller files where the data has been divided according to the two cell lines used in the experiment.

## References

[1] Newman, M. E. J. (2003). Mixing patterns in networks. *Physical Review E*, 67(2), 13. ([doi](https://doi.org/10.1103/PhysRevE.67.026126))

[2] Stumpf, P. S., Smith, R. C. G., Lenz, M., Schuppert, A., Müller, F.-J., Babtie, A. C., Thalia, E. C., Stumpf, M. P. H., Please, C. P., Howison, S. D., Arai, F. and MacArthur, B. D. (2017). Stem cell differentiation as a non-Markov stochastic process. *Cell Systems*, 5(3), 268–282.  ([doi](https://doi.org/10.1016/j.cels.2017.08.009))
