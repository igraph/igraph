
library(igraph)
library(Matrix, quietly=TRUE, warn.conflicts=FALSE)

## Are these vertices connected?
g <- graph.tree(20)
g[1,2]
g[c(1,1,7), c(2,3,14)]
g[c(1,1,7), c(5,3,12)]
g[c(1,1,1,1), c(2,3,2,2)]
g[c(8,17), c(17,8)]

## The same with symbolic names
V(g)$name <- letters[1:vcount(g)]
g['a','b']
g[c('a','a','g'), c('b','c','n')]
g[c('a','a','g'), c('e','c','l')]
g[c('a','a','a','a'), c('b','c','b','b')]
g[c('h','q'), c('q','h')]

## Logical vectors
g[degree(g,mode="in")==0,2]
g[2:3,TRUE]

## Negative indices
g[2:3,-1]

## Weighted graphs
el <- get.edgelist(g, names=FALSE)
E(g)$weight <- el[,1] * el[,2]
g[1,2]
g[c(1,1,7), c(2,3,14)]
g[c(1,1,7), c(5,3,12)]
g[c(1,1,1,1), c(2,3,2,2)]
g[c(8,17), c(17,8)]

## Weighted graph, with symbolic names
g['a','b']
g[c('a','a','g'), c('b','c','n')]
g[c('a','a','g'), c('e','c','l')]
g[c('a','a','a','a'), c('b','c','b','b')]
g[c('h','q'), c('q','h')]

################################################################

## Adjacent vertices
g[[1, ]]
g[[, 2]]
g[[, 2, directed=FALSE]]
g[[2, directed=FALSE]]

g[[1:3, ]]
g[[, 1:3]]

# Same with vertex names
g[['a', ]]
g[[, 'b']]
g[[, 'b', directed=FALSE]]
g[['b', directed=FALSE]]

g[[letters[1:3],]]
g[[, letters[1:3]]]

# Logical vectors
g[[degree(g,mode="in")==0,]]

# Filtering on both ends
g[[1:10, 1:10]]

################################################################

## Query edge ids
g[1,2, edges=TRUE]
g[c(1,1,7), c(2,3,14), edges=TRUE]
g[c(1,1,7), c(5,3,12), edges=TRUE]
g[c(1,1,1,1), c(2,3,2,2), edges=TRUE]
g[c(8,17), c(17,8), edges=TRUE]

## The same with symbolic names
g['a','b', edges=TRUE]
g[c('a','a','g'), c('b','c','n'), edges=TRUE]
g[c('a','a','g'), c('e','c','l'), edges=TRUE]
g[c('a','a','a','a'), c('b','c','b','b'), edges=TRUE]
g[c('h','q'), c('q','h'), edges=TRUE]

################################################################

## Incident edges of vertices
g[[1, , edges=TRUE]]
g[[, 2, edges=TRUE]]
g[[, 2, directed=FALSE, edges=TRUE]]
g[[2, directed=FALSE, edges=TRUE]]

g[[1:3, , edges=TRUE]]
g[[, 1:3, edges=TRUE]]

# Same with vertex names
g[['a', , edges=TRUE]]
g[[, 'b', edges=TRUE]]
g[[, 'b', directed=FALSE, edges=TRUE]]
g[['b', directed=FALSE, edges=TRUE]]

g[[letters[1:3],, edges=TRUE]]
g[[, letters[1:3], edges=TRUE]]

# Filtering on both ends
g[[1:10, 1:10, edges=TRUE]]

#################################################################

## from & to
g <- graph.tree(20)
g[from=c(1,2,2,3), to=c(3,4,8,7)]

V(g)$name <- letters[1:20]
g[from=c("a","b","b","c"), to=c("c","d","h","g")]

E(g)$weight <- (1:ecount(g)) ^ 2 
g[from=c("a","b","b","c"), to=c("c","d","h","g")]

g[from=c("a","b","b","c"), to=c("c","d","h","g"), edges=TRUE]
