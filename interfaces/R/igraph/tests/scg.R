
library(igraph)

tree <- graph.tree(10, 3, "undirected")
treeM <- get.adjacency(tree)
treeM2 <- get.adjacency(tree, sparse=FALSE)

g <- scg(tree, ev=1, nt=3, mtype="laplacian")
M <- scg(treeM, ev=1, nt=3, mtype="laplacian")
M2 <- scg(treeM2, ev=1, nt=3, mtype="laplacian")
graph.laplacian(g$Xt)
M$Xt
M2$Xt

g <- scg(tree, ev=1, nt=3, mtype="symmetric")
M <- scg(treeM, ev=1, nt=3, mtype="symmetric")
M2 <- scg(treeM2, ev=1, nt=3, mtype="symmetric")
get.adjacency(g$Xt, attr="weight")
M$Xt
M2$Xt

g <- scg(tree, ev=1, nt=3, mtype="stochastic")
M <- scg(treeM, ev=1, nt=3, mtype="stochastic")
M2 <- scg(treeM2, ev=1, nt=3, mtype="stochastic")
get.adjacency(g$Xt, attr="weight")
M$Xt
M2$Xt


############

library(igraph)

dirtree <- graph.tree(10, 3, "out")
dirtreeM <- get.adjacency(dirtree)
dirtreeM2 <- get.adjacency(dirtree, sparse=FALSE)

g <- scg(dirtree, ev=3, nt=2, mtype="laplacian")
M <- scg(dirtreeM, ev=3, nt=2, mtype="laplacian")
M2 <- scg(dirtreeM2, ev=3, nt=2, mtype="laplacian")
graph.laplacian(g$Xt)
M$Xt
M2$Xt

g <- scg(dirtree, ev=3, nt=2, mtype="symmetric")
M <- scg(dirtreeM, ev=3, nt=2, mtype="symmetric")
M2 <- scg(dirtreeM2, ev=3, nt=2, mtype="symmetric")
get.adjacency(g$Xt, attr="weight")
M$Xt
M2$Xt

dirtree2 <- graph.tree(10, 3, "in") %u% graph.tree(10, 3, "out")
dirtree2 <- delete.edges(dirtree2, E(dirtree2)[1 %->% 2])
dirtree2M <- get.adjacency(dirtree2)
dirtree2M2 <- get.adjacency(dirtree2, sparse=FALSE)

g <- scg(dirtree2, ev=3, nt=2, mtype="stochastic")
M <- scg(dirtree2M, ev=3, nt=2, mtype="stochastic")
M2 <- scg(dirtree2M2, ev=3, nt=2, mtype="stochastic")
get.adjacency(g$Xt, attr="weight")
M$Xt
M2$Xt

