
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/kapmine.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "kapmine.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/kapmine.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

mat1 <- mat[1:(length(mat)/2)]
mat2 <- mat[(length(mat)/2+1):length(mat)]

mat1 <- matrix(mat1, sqrt(length(mat1)), byrow=TRUE)
mat2 <- matrix(mat2, sqrt(length(mat2)), byrow=TRUE)

library(igraph)
g1 <- graph.adjacency(mat1, mode="undirected")
g2 <- graph.adjacency(mat2, mode="undirected")
el <- data.frame(rbind(get.edgelist(g1), get.edgelist(g2)))
el <- cbind(el, Type=rep(c("multi", "uni"), c(ecount(g1), ecount(g2))))
g <- graph.data.frame(el, directed=FALSE)
g <- remove.vertex.attribute(g, "name")

g$name <- "Kapferer mine"
g$Author <- "Bruce Kapferer"
g$Citation <- "Kapferer B. (1969). Norms and the manipulation of relationships in a work context. In J Mitchell (ed), Social networks in urban situations. Manchester: Manchester University Press.\n\nDoreian P. (1974). On the connectivity of social networks. Journal of Mathematical Sociology, 3, 245-258."
g$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

kapmine <- g
save(kapmine, file="/tmp/kapmine.Rdata.gz")
