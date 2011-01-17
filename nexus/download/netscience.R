
url <- "http://www-personal.umich.edu/~mejn/netdata/netscience.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "netscience.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip netscience.zip"))

txt <- readLines(paste(sep="", tmp, "/netscience.txt"))
gml <- paste(sep="", tmp, "/netscience.gml")

library(igraph)
g <- read.graph(gml, format="gml")

g <- remove.vertex.attribute(g, "id")

V(g)$name <- V(g)$label
g <- remove.vertex.attribute(g, "label")

E(g)$weight <- E(g)$value
g <- remove.edge.attribute(g, "value")

g$name <- "Coauthorships in network science"
g$Author <- "Mark E. J. Newman"
g$Citation <- "M. E. J. Newman, Finding community structure in networks using the eigenvectors of matrices, Phys. Rev. E 74, 036104 (2006)"
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

netscience <- g
save(netscience, file="/tmp/netscience.Rdata")

