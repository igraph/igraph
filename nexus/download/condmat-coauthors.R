
url <- "http://www-personal.umich.edu/~mejn/netdata/cond-mat.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "cond-mat.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip cond-mat.zip"))

txt <- readLines(paste(sep="", tmp, "/cond-mat.txt"))
gml <- paste(sep="", tmp, "/cond-mat.gml")

library(igraph)
g <- read.graph(gml, format="gml")

g <- remove.vertex.attribute(g, "id")

V(g)$name <- V(g)$label
g <- remove.vertex.attribute(g, "label")

E(g)$weight <- E(g)$value
g <- remove.edge.attribute(g, "value")

g$name <- "Condensed matter collaborations"
g$Author <- "Mark E. J. Newman"
g$Citation <- "M. E. J. Newman, The structure of scientific collaboration networks, Proc. Natl. Acad. Sci. USA 98, 404-409 (2001).\n\nM. E. J. Newman, Scientific collaboration networks: I. Network construction and fundamental results, Phys. Rev. E 64, 016131 (2001).\n\nM. E. J. Newman, Scientific collaboration networks: II. Shortest paths, weighted networks, and centrality, Phys. Rev. E 64, 016132 (2001)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

condmatcollab <- g
save(condmatcollab, file="/tmp/condmatcollab.Rdata")
