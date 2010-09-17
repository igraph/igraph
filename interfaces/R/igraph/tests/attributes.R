
library(igraph)

## Create a small ring graph, assign attributes
ring <- graph.formula( A-B-C-D-E-F-G-A )
E(ring)$weight <- seq_len(ecount(ring))

## Query attributes
V(ring)$name
E(ring)$weight

## TODO: subsetting


