
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/kaptail.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "kaptail.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/kaptail.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

labs <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
netlabs <- l[ (grep("^LEVEL LABELS:", l)+1):(grep("^DATA:", l)-1) ]
  
mat1 <- mat[1:(length(mat)/4)]
mat2 <- mat[(length(mat)/4+1):(length(mat)/4*2)]
mat3 <- mat[(length(mat)/4*2+1):(length(mat)/4*3)]
mat4 <- mat[(length(mat)/4*3+1):length(mat)]

mat1 <- matrix(mat1, sqrt(length(mat1)), byrow=TRUE)
mat2 <- matrix(mat2, sqrt(length(mat2)), byrow=TRUE)
mat3 <- matrix(mat3, sqrt(length(mat3)), byrow=TRUE)
mat4 <- matrix(mat4, sqrt(length(mat4)), byrow=TRUE)

library(igraph)
kaptail <- list(graph.adjacency(mat1, mode="undirected"),
                graph.adjacency(mat2, mode="undirected"),
                graph.adjacency(mat3, mode="directed"),
                graph.adjacency(mat4, mode="directed"))
names(kaptail) <- netlabs
for (i in seq_along(kaptail)) V(kaptail[[i]])$name <- labs
for (i in seq_along(kaptail)) {
  kaptail[[i]]$name <- paste("Kapferer tailor shop,", netlabs[[i]])
  kaptail[[i]]$Author <- "Bruce Kapferer"
  kaptail[[i]]$Citation <- "Kapferer B. (1972). Strategy and transaction in an African factory. Manchester: Manchester University Press."
  kaptail[[i]]$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"
}

save(kaptail, file="/tmp/kaptail.Rdata.gz")
