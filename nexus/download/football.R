
url <- "http://www-personal.umich.edu/~mejn/netdata/football.zip"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "football.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip football.zip"))

txt <- readLines(paste(sep="", tmp, "/football.txt"))
gml <- paste(sep="", tmp, "/football.gml")

library(igraph)
g <- read.graph(gml, format="gml")
g <- remove.vertex.attribute(g, "id")

V(g)$name <- V(g)$label
g <- remove.vertex.attribute(g, "label")

V(g)$Conference <- V(g)$value
g <- remove.vertex.attribute(g, "value")

g$name <- "Network of American college football games"
g$Author <- "Michelle Girvan and Mark EJ Newman"
g$Citation <- "M. Girvan and M. E. J. Newman, Community structure in social and biological networks, Proc. Natl. Acad. Sci. USA 99, 7821-7826 (2002)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(txt, collapse="\n"), "\n", sep="")

football <- g
save(football, file="/tmp/football.Rdata")

