
url <- "http://www-personal.umich.edu/~mejn/netdata/dolphins.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "dolphins.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip dolphins.zip"))

txt <- readLines(paste(sep="", tmp, "/dolphins.txt"))
gml <- paste(sep="", tmp, "/dolphins.gml")

library(igraph)
g <- read.graph(gml, format="gml")

g <- remove.vertex.attribute(g, "id")

V(g)$name <- V(g)$label
g <- remove.vertex.attribute(g, "label")

g$name <- "Dophin social network"
g$Author <- "D. Lusseau et al."
g$Citation <- "D. Lusseau, K. Schneider, O. J. Boisseau, P. Haase, E. Slooten, and S. M. Dawson, The bottlenose dolphin community of Doubtful Sound features a large proportion of long-lasting associations, Behavioral Ecology and Sociobiology 54, 396-405 (2003)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

dolphins <- g
save(dolphins, file="/tmp/dolphins.Rdata")

