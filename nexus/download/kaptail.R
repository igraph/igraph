
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/kaptail.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "kaptail.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/kaptail.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

mat1 <- mat[1:(length(mat)/4)]
mat2 <- mat[(length(mat)/4+1):(length(mat)/4*2)]
mat3 <- mat[(length(mat)/4*2+1):(length(mat)/4*3)]
mat4 <- mat[(length(mat)/4*3+1):length(mat)]

mat1 <- matrix(mat1, sqrt(length(mat1)), byrow=TRUE)
mat2 <- matrix(mat2, sqrt(length(mat2)), byrow=TRUE)
mat3 <- matrix(mat3, sqrt(length(mat3)), byrow=TRUE)
mat4 <- matrix(mat4, sqrt(length(mat4)), byrow=TRUE)

library(igraph)
g1 <- graph.adjacency(mat1)
g2 <- graph.adjacency(mat2)
g3 <- graph.adjacency(mat3)
g4 <- graph.adjacency(mat4)
el <- data.frame(rbind(get.edgelist(g1), get.edgelist(g2),
                       get.edgelist(g3), get.edgelist(g4)))
el <- cbind(el, Type=rep(c("sociational1", "sociational2",
                  "instrumental1", "instrumental2"),
                  c(ecount(g1), ecount(g2), ecount(g3), ecount(g4))))
g <- graph.data.frame(el, directed=TRUE)
g <- remove.vertex.attribute(g, "name")

g$name <- "Kapferer tailor shop"
g$Author <- "Bruce Kapferer"
g$Citation <- "Kapferer B. (1972). Strategy and transaction in an African factory. Manchester: Manchester University Press."
g$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

kaptail <- g
save(kaptail, file="/tmp/kaptail.Rdata.gz")
