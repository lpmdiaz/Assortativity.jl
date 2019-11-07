# Assortativity.jl

(code coverage etc.)

## Description

Based on code / work by [Thalia Chan](https://github.com/Tchanders). Supports paper. Based on edge threshold approach. Gaining confidence in inferred networks. Also includes additional stuff not discussed in paper (e.g. second neighbour measures). Relies primarily on [LightGraphs.jl](https://github.com/JuliaGraphs/LightGraphs.jl) and [NetworkInference.jl](https://github.com/Tchanders/NetworkInference.jl).

(use [fork](https://github.com/lpmdiaz/NetworkInference.jl) of NetworkInference.jl for correlation (for now))

### Requirements

Requires data + file describing gene annotations (or annotations for any type of data not just gene as long as described by graph).

## Installation

`]
clone
`

## Usage

### Infer network

inference algorithms + export net (symmetric)
`nodes_number = length(network.nodes) # number of nodes in the full InferredNetwork`
`edges_number = length(network.edges) # number of edges in the full InferredNetwork`

### Import inferred network

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

Later / WIP.

### Plotting utilities

Later / WIP.
