
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wiring.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "wiring.dat")

download.file(url, dest)

l <- readLines(paste(sep="", tmp, "/wiring.dat"))
rowlab <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
collab <- l[ (grep("^COLUMN LABELS:", l)+1):(grep("^LEVEL LABELS:", l)-1) ]
llab <- l[ (grep("^LEVEL LABELS:", l)+1):(grep("^DATA:", l)-1) ]

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- matrix(scan(tc), nc=14, byrow=TRUE)
close(tc)

wiring <- list(mat[1:14,], mat[15:28,], mat[29:42,], mat[43:56,],
               mat[57:70,], mat[71:84,])

for (i in seq_along(wiring)) {
  rownames(wiring[[i]]) <- colnames(wiring[[i]]) <- rowlab
}

library(igraph)
wiring <- list(graph.adjacency(wiring[[1]], mode="undirected"),
               graph.adjacency(wiring[[2]], mode="undirected"),
               graph.adjacency(wiring[[3]], mode="undirected"),
               graph.adjacency(wiring[[4]], mode="undirected"),
               graph.adjacency(wiring[[5]], mode="directed"),
               graph.adjacency(wiring[[3]], mode="directed", weighted=TRUE))
names(wiring) <- llab

for (i in seq_along(wiring)) {
  wiring[[i]]$name <- paste("Roethlisberger & Dickson bank wiring room,",
                            llab[i])
  wiring[[i]]$Author <- "Roethlisberger F. and Dickson W."
  wiring[[i]]$Citation <- "Breiger R., Boorman S. and Arabie P. (1975). An algorithm for clustering relational data with applications to social network analysis and comparison with multidimensional scaling. Journal of Mathematical Psychology, 12, 328-383.\n\nHomans G. (1950). The human group. New York: Harcourt-Brace.\n\nRoethlisberger F. and Dickson W. (1939). Management and the worker. Cambridge: Cambridge University Press."
  wiring[[i]]$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"
}

save(wiring, file="/tmp/wiring.Rdata.gz")

