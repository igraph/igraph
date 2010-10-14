
library(igraph)

## Are these vertices connected?
g <- graph.tree(20)
g[1,2]
g[c(1,1,7), c(2,3,14)]
g[c(1,1,7), c(5,3,12)]
g[c(1,1,1,1), c(2,3,2,2)]
g[c(1,1,1,1), c(2,3,2,2), multi=TRUE]
g[c(8,17), c(17,8)]
g[c(8,17), c(17,8), directed=FALSE]
g[c(8,17), c(17,8), directed=FALSE, multi=TRUE]

## The same with symbolic names
V(g)$name <- letters[1:vcount(g)]
g['a','b']
g[c('a','a','g'), c('b','c','n')]
g[c('a','a','g'), c('e','c','l')]
g[c('a','a','a','a'), c('b','c','b','b')]
g[c('a','a','a','a'), c('b','c','b','b'), multi=TRUE]
g[c('h','q'), c('q','h')]
g[c('h','q'), c('q','h'), directed=FALSE]
g[c('h','q'), c('q','h'), directed=FALSE, multi=TRUE]

## Weighted graphs
el <- get.edgelist(g, names=FALSE)
E(g)$weight <- el[,1] * el[,2]
g[1,2]
g[c(1,1,7), c(2,3,14)]
g[c(1,1,7), c(5,3,12)]
g[c(1,1,1,1), c(2,3,2,2)]
g[c(1,1,1,1), c(2,3,2,2), multi=TRUE]
g[c(8,17), c(17,8)]
g[c(8,17), c(17,8), directed=FALSE]
g[c(8,17), c(17,8), directed=FALSE, multi=TRUE]

## Weighted graph, with symbolic names
g['a','b']
g[c('a','a','g'), c('b','c','n')]
g[c('a','a','g'), c('e','c','l')]
g[c('a','a','a','a'), c('b','c','b','b')]
g[c('a','a','a','a'), c('b','c','b','b'), multi=TRUE]
g[c('h','q'), c('q','h')]
g[c('h','q'), c('q','h'), directed=FALSE]
g[c('h','q'), c('q','h'), directed=FALSE, multi=TRUE]

################################################################

## Adjacent vertices
g[1, ]
g[, 2]
g[, 2, directed=FALSE]
g[2, directed=FALSE]

g[1:3, ]
g[, 1:3]
g[, 1:3, unique=FALSE]
g[1:3, simplify=FALSE]
g[, 1:3, simplify=FALSE]

# Same with vertex names
g['a', ]
g[, 'b']
g[, 'b', directed=FALSE]
g['b', directed=FALSE]

g[letters[1:3],]
g[, letters[1:3]]
g[, letters[1:3], unique=FALSE]
g[letters[1:3], simplify=FALSE]
g[, letters[1:3], simplify=FALSE]

################################################################

# Incident vertices of some edges
get.edges(g, 1:3)
g[edges=1:3]
g[edges=1:3, simplify=FALSE]
g[edges=1:3, unique=FALSE]

# Same with edge names
E(g)$name <- LETTERS[1:ecount(g)]
g[edges=LETTERS[1:3]]
g[edges=LETTERS[1:3], simplify=FALSE]
g[edges=LETTERS[1:3], unique=FALSE]

# Edge names based on vertex names
g[edges=c("a|b", "a|c", "b|d")]
g[edges=c("a|b", "a|c", "b|d"), simplify=FALSE]
g[edges=c("a|b", "a|c", "b|d"), unique=FALSE]

################################################################

## Query edge ids
g[[1,2]]
g[[c(1,1,7), c(2,3,14)]]
g[[c(1,1,7), c(5,3,12)]]
g[[c(1,1,1,1), c(2,3,2,2)]]
g[[c(1,1,1,1), c(2,3,2,2), multi=TRUE]]
g[[c(8,17), c(17,8)]]
g[[c(8,17), c(17,8), directed=FALSE]]
g[[c(8,17), c(17,8), directed=FALSE, multi=TRUE]]

## The same with symbolic names
g[['a','b']]
g[[c('a','a','g'), c('b','c','n')]]
g[[c('a','a','g'), c('e','c','l')]]
g[[c('a','a','a','a'), c('b','c','b','b')]]
g[[c('a','a','a','a'), c('b','c','b','b'), multi=TRUE]]
g[[c('h','q'), c('q','h')]]
g[[c('h','q'), c('q','h'), directed=FALSE]]
g[[c('h','q'), c('q','h'), directed=FALSE, multi=TRUE]]

################################################################

## Incident edges of vertices
g[[1, ]]
g[[, 2]]
g[[, 2, directed=FALSE]]
g[[2, directed=FALSE]]

g[[1:3, ]]
g[[, 1:3]]
g[[, 1:3, unique=FALSE]]
g[[1:3, simplify=FALSE]]
g[[, 1:3, simplify=FALSE]]

# Same with vertex names
g[['a', ]]
g[[, 'b']]
g[[, 'b', directed=FALSE]]
g[['b', directed=FALSE]]

g[[letters[1:3],]]
g[[, letters[1:3]]]
g[[, letters[1:3], unique=FALSE]]
g[[letters[1:3], simplify=FALSE]]
g[[, letters[1:3], simplify=FALSE]]

