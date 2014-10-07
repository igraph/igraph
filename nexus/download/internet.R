
url <- "http://www-personal.umich.edu/~mejn/netdata/as-22july06.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "as-22july06.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip as-22july06.zip"))

txt <- readLines(paste(sep="", tmp, "/as-22july06.txt"))
gml <- paste(sep="", tmp, "/as-22july06.gml")

library(igraph)
g <- read_graph(gml, format="gml")

g <- delete_vertex_attr(g, "id")
g <- delete_vertex_attr(g, "label")

g$name <- "Internet at the AS level, July 22, 2006"
g$Author <- "Mark E. J. Newman"
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

as2006 <- g
save(as2006, file="/tmp/as2006.Rdata")

