
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wolf.dat"
url2 <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/wolfi.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "wolf.dat")
download.file(url, dest)
l <- readLines(dest)

dest <- paste(sep="", tmp, "/", "wolfi.dat")
download.file(url2, dest)
l2 <- readLines(dest)

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- matrix(scan(tc), nc=20, byrow=TRUE)
close(tc)

mat1 <- mat[1:20,]
mat2 <- mat[21:40,]

library(igraph)
g1 <- graph.adjacency(mat1, mode="directed")
g2 <- graph.adjacency(mat2, mode="undirected")

g1$name <- "Wolfe primates, putative kin relationships"
g2$name <- "Wolfe primates, interactions"
g1$Author <- g2$Author <- "Linda Wolfe"
g1$URL <- g2$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

sex <- l2[ (grep("^ROW LABELS:", l2)+1):(grep("^COLUMN LABELS:", l2)-1) ]
sex <- sub("^[0-9]+ ", "", sex)
age <- l2[ (grep("^DATA:", l2)+1):length(l2) ]
age <- as.numeric(sub("[ ][0-9]+$", "", sub("^[ ]*", "", age)))
rank <- l2[ (grep("^DATA:", l2)+1):length(l2) ]
rank <- as.numeric(sub("^[ ]*[0-9]+[ ]*", "", rank))

V(g1)$Sex <- V(g2)$Sex <- sex
V(g1)$Age <- V(g2)$Age <- age
V(g1)$Rank <- V(g2)$Rank <- rank

wolfe <- list(WOLFK=g1, WOLFN=g2)
save(wolfe, file="/tmp/wolfe.Rdata.gz")

