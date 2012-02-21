
url <- "http://www-personal.umich.edu/~mejn/netdata/celegansneural.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "celegansneural.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip celegansneural.zip"))

txt <- readLines(paste(sep="", tmp, "/celegansneural.txt"))
gml <- paste(sep="", tmp, "/celegansneural.gml")

library(igraph)
g <- read.graph(gml, format="gml")

g <- remove.vertex.attribute(g, "id")

V(g)$name <- V(g)$label
g <- remove.vertex.attribute(g, "label")

E(g)$weight <- E(g)$value
g <- remove.edge.attribute(g, "value")

g$name <- "C. Elegans neural network"
g$Author <- "Duncan Watts and Steven Strogatz from original experimental
data by White et al."
g$Citation <- "J. G. White, E. Southgate, J. N. Thompson, and S. Brenner, The structure of the nervous system of the nematode C. Elegans, Phil. Trans. R. Soc. London 314, 1-340 (1986).\n\nD. J. Watts and S. H. Strogatz, Collective dynamics of `small-world' networks, Nature 393, 440-442 (1998)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

celegansneural <- g
save(celegansneural, file="/tmp/celegansneural.Rdata")

