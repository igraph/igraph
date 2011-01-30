
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/bkfrat.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "bkfrat.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/bkfrat.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

mat1 <- mat[1:(length(mat)/2)]
mat2 <- mat[(length(mat)/2+1):length(mat)]

mat1 <- matrix(mat1, sqrt(length(mat1)), byrow=TRUE)
mat2 <- matrix(mat2, sqrt(length(mat2)), byrow=TRUE)

library(igraph)
g1 <- graph.adjacency(mat1, weighted=TRUE, mode="undirected")
g2 <- graph.adjacency(mat1, weighted=TRUE, mode="directed")

g1$name <- "Bernard & Killworth fraternity, BKFRAB"
g2$name <- "Bernard & Killworth fraternity, BKFRAC"
g1$Author <- g2$Author <- "H. Bernard H and P. Killworth"
g1$Citation <- g2$Citation <- "Bernard H, Killworth P and Sailer L. (1980). Informant accuracy in social network data IV. Social Networks, 2, 191-218.\n\nBernard H, Killworth P and Sailer L. (1982). Informant accuracy in social network data V. Social Science Research, 11, 30-66.\n\nRomney K and Weller S. (1984). Predicting informant accuracy from patterns of recall among individuals. Social Networks, 6, 59-78."
g1$URL <- g2$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

bkfrat <- list(BKFRAB=g1, BKFRAC=g2)
save(bkfrat, file="/tmp/bkfrat.Rdata.gz")


