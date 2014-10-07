
R

library(igraph)

# Create the graph

g <- simplify(barabasi.game(10, m=5, directed=F))
g <- simplify(g_np(20, p=6/20, directed=F))
g <- simplify(union(g_star(20, mode="undirected"),
                          g_np(20, p=1/20)))
g <- graph( c(0,1, 0,1, 0,2, 0,2, 1,2, 1,2, 2,3, 2,3, 2,4, 2,4, 3,4),
           direc=FALSE)
el <- edgelist(g)+1
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

igraph.output <- c(vertex_conn(g), edge_conn(g))

print(output)
print(igraph.output)

#######################################

g <- g_np(50, p=15/50)
date()
vertex_conn(g)
edge_conn(g)
date()

el <- edgelist(g)+1
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

