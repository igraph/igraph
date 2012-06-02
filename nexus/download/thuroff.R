
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/thuroff.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "thuroff.dat")

download.file(url, dest)

l <- readLines(paste(sep="", tmp, "/thuroff.dat"))
rowlab <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
collab <- l[ (grep("^COLUMN LABELS:", l)+1):(grep("^LEVEL LABELS:", l)-1) ]

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- matrix(scan(tc), nc=length(rowlab), byrow=TRUE)
close(tc)

mat1 <- mat[1:15,]
mat2 <- mat[16:30,]
colnames(mat1) <- rownames(mat1) <- colnames(mat2) <- rownames(mat2) <- rowlab

library(igraph)
g1 <- graph.adjacency(mat1, mode="directed")
g2 <- graph.adjacency(mat2, mode="undirected")

g1$name <- "Thurman office, organizational chart"
g2$name <- "Thurman office, multiplex ties"

g1$Author <- g2$Author <- "B. Thurman"
g1$Citation <- g2$Citation <- "Thurman B. (1979). In the office: Networks and coalitions. Social Networks, 2, 47-63."
g1$URL <- g2$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm" 

thuroff <- list(THURA=g1, THURM=g2)
save(thuroff, file="/tmp/thuroff.Rdata.gz")

