# Assortativity.jl

(code coverage etc.)

## Description

Based on code / work by [Thalia Chan](https://github.com/Tchanders). Supports paper. Based on edge threshold approach. Gaining confidence in inferred networks. Also includes additional stuff not discussed in paper (e.g. two-walk measures). Relies primarily on [LightGraph](https://github.com/JuliaGraphs/LightGraphs.jl) and [NetworkInference](https://github.com/Tchanders/NetworkInference.jl).

### Requirements

Requires data + file describing gene annotations (or annotations for any type of data not just gene as long as described by graph).

## Installation

'
]
clone?
'

## Usage

### Infer network

inference algorithms + export net (symmetric)
'nodes_number = length(network.nodes) # number of nodes in the full InferredNetwork'
'edges_number = length(network.edges) # number of edges in the full InferredNetwork'

### Import inferred network

### Conversions

Lose edge order when convert InferredNetwork (at threshold or not) to LightGraph. (So always must start from InferredNetwork when setting a threshold.)

### Calculate assortativity and other measures

### Noise \& randomness

### Graph manipulation utilities

### Plotting utilities
