
url <- "http://www-personal.umich.edu/~mejn/netdata/adjnoun.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "adjnoun.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip adjnoun.zip"))

txt <- readLines(paste(sep="", tmp, "/adjnoun.txt"))
gml <- paste(sep="", tmp, "/adjnoun.gml")

library(igraph)
g <- read.graph(gml, format="gml")

g <- remove.vertex.attribute(g, "id")

V(g)$Adjnoun <- V(g)$value
g <- remove.vertex.attribute(g, "value")

V(g)$name <- V(g)$label
g <- remove.vertex.attribute(g, "label")

g$name <- "Word adjacencies from David Copperfield"
g$Author <- "Mark E. J. Newman"
g$Citation <- "M. E. J. Newman: Finding community structure in networks using the eigenvectors of matrices. Phys. Rev. E 74, 036104 (2006)"
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- "The network of common adjective and noun adjacencies for the novel 'David Copperfield' by Charles Dickens, as described by M. Newman.  Nodes represent the most commonly occurring adjectives and nouns in the book.  Node values are 0 for adjectives and 1 for nouns.  Edges connect any pair of words that occur in adjacent position in the text of the book.  Please cite M. E. J. Newman, Finding community structure in networks using the eigenvectors of matrices, Phys. Rev. E 74, 036104 (2006)"

adjnoun <- g
save(adjnoun, file="/tmp/adjnoun.Rdata")

