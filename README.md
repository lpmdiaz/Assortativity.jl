# Assortativity.jl

(code coverage etc.)

## Description

Implementation of the assortativity coefficient for inferred networks. Assortativity is explained in [[1]](#references).
Based on code / work by [Thalia Chan](https://github.com/Tchanders). Supports paper. Based on edge threshold approach. Gaining confidence in inferred networks. Also includes additional stuff not discussed in paper (e.g. second neighbour measures). Relies primarily on [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [NetworkInference.jl](https://github.com/Tchanders/NetworkInference.jl).

If using correlation for inferring networks, use [this fork](https://github.com/lpmdiaz/NetworkInference.jl) of NetworkInference.jl (i.e. `using Pkg; Pkg.clone("https://github.com/lpmdiaz/NetworkInference.jl.git")`).

### Requirements

Requires data + file describing gene annotations (or annotations for any type of data not just gene as long as described by graph).

## Installation

`using Pkg; Pkg.clone("https://github.com/lpmdiaz/Assortativity.jl.git")`

## Usage

First include the package: `using Assortativity`.

### Infer network

Inference algorithms: MI, CLR, PUC, PIDC and two types of correlation coefficients (Pearson and Spearman) both signed and unsigned -- see usage [details](https://github.com/lpmdiaz/NetworkInference.jl#inference-algorithms-currently-implemented) and [reference](https://github.com/lpmdiaz/NetworkInference.jl#references) publication for the first four algorithms.

inference algorithms + export net (symmetric)
`nodes_number = length(network.nodes) # number of nodes in the full InferredNetwork`
`edges_number = length(network.edges) # number of edges in the full InferredNetwork`

### Import inferred network

...

### Conversions

Lose edge order when convert InferredNetwork (at threshold or not) to LightGraph. (So always must start from InferredNetwork when setting a threshold.)

### Calculate assortativity and other measures

Degree assortativity: also have excess degree.
Second neighbour: does not consider walk back to the starting node. Different behaviour as iterate over nodes and not edges.

### Noise \& randomness

`randomise_annotations` behaviours:
- 1: relaxed -- Fully random and thus allows attributing the key the same value it already had.
- 2: swapped -- Swaps two random keys' values with the possibility that the swapped values are the same than before.
- 3: strict -- Changes a random key's value to one that must be different from the one it had. However, several runs could revert back to the original value.
- 4: stricter -- Changes a random key's value to one that must be different from the one it originally had. This behaviour prevents reverting back to the original value of the keys. Will often result in as many changes in the dictionary as iterations.
- 5: strictest -- Will always result in as many changes in the dictionary as iterations.

### Graph manipulation utilities

Future addition.

### Plotting utilities

Future addition.

## References

[1] Newman, M. E. J. (2003). Mixing patterns in networks. _Physical Review E_, 67(2), 13. ([doi])(https://doi.org/10.1103/PhysRevE.67.026126)
