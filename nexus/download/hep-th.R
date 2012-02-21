
url <- "http://www-personal.umich.edu/~mejn/netdata/hep-th.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "hep-th.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip hep-th.zip"))

txt <- readLines(paste(sep="", tmp, "/hep-th.txt"))
gml <- paste(sep="", tmp, "/hep-th.gml")

library(igraph)
g <- read.graph(gml, format="gml")

g <- remove.vertex.attribute(g, "id")

V(g)$name <- V(g)$label
g <- remove.vertex.attribute(g, "label")

E(g)$weight <- E(g)$value
g <- remove.edge.attribute(g, "value")

g$name <- "High energy physics collaborations"
g$Author <- "Mark E. J. Newman"
g$Citation <- "M. E. J. Newman, The structure of scientific collaboration networks, Proc. Natl. Acad. Sci. USA 98, 404-409 (2001)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

hepcollab <- g
save(hepcollab, file="/tmp/hepcollab.Rdata")
