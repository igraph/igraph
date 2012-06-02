
url <- "http://www-personal.umich.edu/~mejn/netdata/power.zip"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "power.zip")

download.file(url, dest)
system(paste("cd ", tmp, "; unzip power.zip"))

txt <- readLines(paste(sep="", tmp, "/power.txt"))
gml <- paste(sep="", tmp, "/power.gml")

library(igraph)
g <- read.graph(gml, format="gml")

g <- remove.vertex.attribute(g, "id")

g$name <- "Western states power grid"
g$Author <- "Duncan Watts and Steven Strogatz"
g$Citation <- "D. J. Watts and S. H. Strogatz, Collective dynamics of `small-world' networks, Nature 393, 440-442 (1998)."
g$URL <- "http://www-personal.umich.edu/~mejn/netdata/"
g$Description <- paste(paste(collapse="\n", txt), sep="", "\n")

powergrid <- g
save(powergrid, file="/tmp/powergrid.Rdata")


