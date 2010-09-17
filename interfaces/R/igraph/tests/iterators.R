
library(igraph)

## Create a small ring graph, assign attributes
ring <- graph.formula( A-B-C-D-E-F-G-A )
E(ring)$weight <- seq_len(ecount(ring))

## Selection based on attributes
E(ring)[ weight < 4 ]
V(ring)[ c("A", "C") ]

## TODO: %--%, %->%, other special functions

