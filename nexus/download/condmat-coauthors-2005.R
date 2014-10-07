
url <- "http://www-personal.umich.edu/~mejn/netdata/cond-mat-2005.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "cond-mat-2005.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip cond-mat-2005.zip"))

txt <- readLines(paste(sep="", tmp, "/cond-mat-2005.txt"))
gml <- paste(sep="", tmp, "/cond-mat-2005.gml")

library(igraph)
g <- read_graph(gml, format="gml")

g <- delete_vertex_attr(g, "id")

V(g)$name <- V(g)$label
g <- delete_vertex_attr(g, "label")

E(g)$weight <- E(g)$value
g <- delete_edge_attr(g, "value")

g$name <- "Condensed matter collaborations, 2005"
g$Author <- "Mark E. J. Newman"
g$Citation <- "M. E. J. Newman, The structure of scientific collaboration networks, Proc. Natl. Acad. Sci. USA 98, 404-409 (2001).\n\nM. E. J. Newman, Fast algorithm for detecting community structure in networks, Phys. Rev. E 69, 066133 (2004)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

condmatcollab2005 <- g
save(condmatcollab2005, file="/tmp/condmatcollab2005.Rdata")
