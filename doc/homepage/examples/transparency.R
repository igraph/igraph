library(igraph)

g <- barabasi.game(100, power=0.6, m=10)
V(g)$name <- seq(vcount(g))
g <- simplify(g)
l <- layout.fruchterman.reingold(g)
l <- layout.norm(l, -1,1, -1,1) 

fcs <- fastgreedy.community(simplify(as.undirected(g)))
Q <- round(modularity(fcs), 3)

# png(file="fastgreedy.png")

plot(g, layout=l, vertex.size=3, vertex.label=NA, vertex.color="#ff000033",
     vertex.frame.color="#ff000033", edge.color="#55555533", edge.arrow.size=0.3,
     rescale=FALSE, xlim=range(l[,1]), ylim=range(l[,2]),
     main=paste(sep="", "Fast greedy community detection,\nQ=", Q))

g2 <- induced.subgraph(g, communities(fcs)[[1]])
l2 <- l[ communities(fcs)[[1]], ]

plot(g2, layout=l2, vertex.size=3, vertex.label=V(g2)$name,
     vertex.color="#ff0000", vertex.frame.color="#ff0000", edge.color="#555555",
     vertex.label.dist=0.5, vertex.label.cex=0.8, vertex.label.font=2,
     edge.arrow.size=0.3, add=TRUE, rescale=FALSE)

nodes <- communities(fcs)[[1]]
nodes <- paste(collapse=", ", nodes)
text(0,-1.3, nodes, col="blue")

# dev.off()
