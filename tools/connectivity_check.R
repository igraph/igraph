
R

library(igraph)

# Create the graph

g <- simplify(barabasi.game(10, m=5, directed=F))
g <- simplify(erdos.renyi.game(20, p=6/20, directed=F))
g <- simplify(graph.union(graph.star(20, mode="undirected"),
                          erdos.renyi.game(20, p=1/20)))
g <- graph( c(0,1, 0,1, 0,2, 0,2, 1,2, 1,2, 2,3, 2,3, 2,4, 2,4, 3,4),
           direc=FALSE)
el <- get.edgelist(g)+1
write(t(el), "edgelist.txt", ncol=2)

math.command <- ' <<DiscreteMath`Combinatorica`;
edges = ReadList["edgelist.txt", {Number, Number}];
g = FromOrderedPairs[edges, Type -> Undirected];
VertexConnectivity[g]
EdgeConnectivity[g]'

command <- paste("echo ", "'", math.command, "'", "|MathKernel")
output <- system(command, int=T)
output <- output[ grep("^Out", output) ]
output <- sub("^.* ", "", output)
output <- as.numeric(output)

igraph.output <- c(vertex.connectivity(g), edge.connectivity(g))

print(output)
print(igraph.output)

#######################################

g <- erdos.renyi.game(50, p=15/50)
date()
vertex.connectivity(g)
edge.connectivity(g)
date()

el <- get.edgelist(g)+1
write(t(el), "edgelist.txt", ncol=2)

math.command <- ' <<DiscreteMath`Combinatorica`;
edges = ReadList["edgelist.txt", {Number, Number}];
g = FromOrderedPairs[edges, Type -> Undirected];
VertexConnectivity[g]
EdgeConnectivity[g]'

command <- paste("echo ", "'", math.command, "'", "|MathKernel")
date()
output <- system(command, int=T)
date()
output <- output[ grep("^Out", output) ]
output <- sub("^.* ", "", output)
output <- as.numeric(output)

