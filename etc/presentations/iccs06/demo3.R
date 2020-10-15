library(igraph)
g <- \emph{erdos.renyi.game}(5000, 0.8/5000)

cl <- \emph{clusters}(g)
large <- which(cl$csize > 3)-1
g2 <- \emph{subgraph}(g, which(cl$membership %in% large)-1)

graphs <- \emph{decompose.graph}(g2)
layouts <- lapply(graphs, \emph{layout.kamada.kawai})
coords <- \emph{layout.merge}(graphs, layouts)
g3 <- \emph{graph.disjoint.union}(graphs)

cl3 <- \emph{clusters}(g3)
cl.no <- length(cl3$csize)
colorbar <- heat.colors(cl.no)
\emph{V(g3)$color} <- colorbar[cl3$membership+1]

\emph{plot}(g3, layout=coords, vertex.size=1, labels=NA, vertex.color="a:color")
