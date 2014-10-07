
url <- "http://www-personal.umich.edu/~mejn/netdata/astro-ph.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "astro-ph.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip astro-ph.zip"))

txt <- readLines(paste(sep="", tmp, "/astro-ph.txt"))
gml <- paste(sep="", tmp, "/astro-ph.gml")

library(igraph)
g <- read_graph(gml, format="gml")

g <- delete_vertex_attr(g, "id")

V(g)$name <- V(g)$label
g <- delete_vertex_attr(g, "label")

E(g)$weight <- E(g)$value
g <- delete_edge_attr(g, "value")

g$name <- "Astrophysics collaborations"
g$Author <- "Mark E. J. Newman"
g$Citation <- "M. E. J. Newman, The structure of scientific collaboration networks, Proc. Natl. Acad. Sci. USA 98, 404-409 (2001)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

astrocollab <- g
save(astrocollab, file="/tmp/astrocollab.Rdata")
