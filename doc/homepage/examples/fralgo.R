library(igraph)

source("frgraphs.R")
graphs <- frgraphs()
lay <- lapply(graphs, layout.fruchterman.reingold, niter=3000)

# png(file="frplots.png")
par(mai=c(0,0,0,0))
layout(matrix(1:16, nr=4, byrow=TRUE))
layout.show(16)
for (i in seq(along=graphs)) {
  plot(graphs[[i]], layout=lay[[i]],
       vertex.label=NA, vertex.size=13, edge.color="black",
       vertex.color="red")
}
# dev.off()

