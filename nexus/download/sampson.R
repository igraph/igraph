
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/sampson.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "sampson.dat")

download.file(url, dest)

l <- readLines(paste(sep="", tmp, "/sampson.dat"))
rowlab <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
collab <- l[ (grep("^COLUMN LABELS:", l)+1):(grep("^LEVEL LABELS:", l)-1) ]
llab <- l[ (grep("^LEVEL LABELS:", l)+1):(grep("^DATA:$", l)-1)]
llab[10] <- "SAMPNPR"

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- matrix(scan(tc), nr=length(llab), byrow=TRUE)
close(tc)

mat <- apply(mat, 1, list)
mat <- lapply(mat, "[[", 1)
mat <- lapply(mat, matrix, nc=18, byrow=TRUE)

for (i in seq_along(mat)) {
  rownames(mat[[i]]) <- colnames(mat[[i]]) <- rowlab
}

library(igraph)
mat <- lapply(mat, graph.adjacency, weighted=TRUE)
names(mat) <- llab

for (i in seq_along(mat)) {
  mat[[i]]$name <- paste("Sampson monastery,", llab[i])
  mat[[i]]$Author <- "S. Sampson"
  mat[[i]]$Citation <- "Breiger R., Boorman S. and Arabie P. (1975). An algorithm for clustering relational data with applications to social network analysis and comparison with multidimensional scaling. Journal of Mathematical Psychology, 12, 328-383.\n\nSampson, S. (1969). Crisis in a cloister. Unpublished doctoral dissertation, Cornell University."
  mat[[i]]$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"
}

sampson <- mat

save(sampson, file="/tmp/sampson.Rdata.gz")

