# Read the files first
traits <- read.csv("traits.csv", head=FALSE)
rel <- read.csv("relations.csv", head=FALSE)

# Create the graph, add the vertices
library(igraph)
g <- graph.empty()
g <- add.vertices(g, nrow(traits), 
                     name=as.character(traits[,1]), age=traits[,2],
                     gender=as.character(traits[,3]))

# Extract first names from the full names
names <- sapply(strsplit(V(g)$name, " "), "[",1)
ids <- 1:length(names)-1
names(ids) <- names

# Create the edges
from <- as.character(rel[,1])
to <- as.character(rel[,2])
edges <- matrix(c(ids[from], ids[to]), nc=2)

# Add the edges
g <- add.edges(g, t(edges), 
               room=as.character(rel[,3]),
	       friend=rel[,4], advice=rel[,5])
