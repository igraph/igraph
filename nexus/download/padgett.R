
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/padgett.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "padgett.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/padgett.dat"))

url2 <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/padgw.dat"
dest <- paste(sep="", tmp, "/", "padgw.dat")
download.file(url2, dest)
l2 <- readLines(paste(sep="", tmp, "/padgw.dat"))

rowlab <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
collab <- l[ (grep("^COLUMN LABELS:", l)+1):(grep("^LEVEL LABELS:", l)-1) ]

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- as.matrix(read.table(tc))
close(tc)

mat1 <- mat[1:16,]
mat2 <- mat[17:32,]
colnames(mat1) <- rownames(mat1) <- colnames(mat2) <- rownames(mat2) <- rowlab

library(igraph)
g1 <- graph.adjacency(mat1, mode="undirected")
g2 <- graph.adjacency(mat2, mode="undirected")

mlabs <- l2[ (grep("^COLUMN LABELS:", l2)+1):(grep("^DATA", l2)-1) ]
alabs <- l2[ (grep("^ROW LABELS:", l2)+1):(grep("^COLUMN LABELS", l2)-1) ]
data2 <- l2[ (grep("^DATA:", l2)+1):length(l2)]
data2 <- paste(data2, collapse="\n")

tc <- textConnection(data2)
meta <- as.matrix(read.table(tc))
close(tc)

rownames(meta) <- alabs
colnames(meta) <- mlabs
meta <- meta[V(g1)$name,]

g1$name <- "Padgett Frorentine families, business"
g2$name <- "Padgett Frorentine families, marriage"
g1$Author <- g2$Author <- "John Padgett"
g1$Citation <- g2$Citation <- "Kent D. (1978). The rise of the Medici: Faction in Florence, 1426-1434. Oxford: Oxford University Press."
g1$URL <- g2$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"
V(g1)$Wealth <- V(g2)$Wealth <- meta[,1]
V(g1)$Priors <- V(g2)$Priors <- meta[,2]
V(g1)$Ties <- V(g2)$Ties <- meta[,3]

padgett <- list(PADGB=g1, PADGM=g2)
save(padgett, file="/tmp/padgett.Rdata.gz")

