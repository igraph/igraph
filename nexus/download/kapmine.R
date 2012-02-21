
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/kapmine.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "kapmine.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/kapmine.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]
adata <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

mat1 <- mat[1:(length(mat)/2)]
mat2 <- mat[(length(mat)/2+1):length(mat)]

mat1 <- matrix(mat1, sqrt(length(mat1)), byrow=TRUE)
mat2 <- matrix(mat2, sqrt(length(mat2)), byrow=TRUE)
colnames(mat1) <- rownames(mat1) <- colnames(mat2) <- rownames(mat2) <- adata

library(igraph)
g1 <- graph.adjacency(mat1, mode="undirected")
g2 <- graph.adjacency(mat2, mode="undirected")

g1$name <- "Kapferer mine, KAPFMM"
g2$name <- "Kapferer mine, KAPFMU"
g1$Author <- g2$Author <- "Bruce Kapferer"
g1$Citation <- g2$Citation <- "Kapferer B. (1969). Norms and the manipulation of relationships in a work context. In J Mitchell (ed), Social networks in urban situations. Manchester: Manchester University Press.\n\nDoreian P. (1974). On the connectivity of social networks. Journal of Mathematical Sociology, 3, 245-258."
g1$URL <- g2$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

kapmine <- list(KAPFMM=g1, KAPFMU=g2)
save(kapmine, file="/tmp/kapmine.Rdata.gz")
