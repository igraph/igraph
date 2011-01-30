
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/bktec.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "bktec.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/bktec.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

mat1 <- mat[1:(length(mat)/2)]
mat2 <- mat[(length(mat)/2+1):length(mat)]

mat1 <- matrix(mat1, sqrt(length(mat1)), byrow=TRUE)
mat2 <- matrix(mat2, sqrt(length(mat2)), byrow=TRUE)
diag(mat2) <- 0

library(igraph)
g1 <- graph.adjacency(mat1, mode="undirected", weighted=TRUE)
g2 <- graph.adjacency(mat2, mode="directed", weighted=TRUE)

g1$name <- "Bernard & Killworth technical research group, BKTECB"
g2$name <- "Bernard & Killworth technical research group, BKTECC"
g1$Author <- g2$Author <- "H. Bernard H and P. Killworth"
g1$Citation <- g2$Citation <- "Killworth B and Bernard H. (1976). Informant accuracy in social network data. Human Organization, 35, 269-286.\n\nBernard H and Killworth P. (1977). Informant accuracy in social network data II. Human Communication Research, 4, 3-18.\n\nKillworth P and Bernard H. (1979). Informant accuracy in social network data III. Social Networks, 2, 19-46."
g1$URL <- g2$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

bktec <- list(BKTECB=g1, BKTECC=g2)
save(bktec, file="/tmp/bktec.Rdata.gz")

