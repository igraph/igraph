
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/knokbur.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "knokbur.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/knokbur.dat"))

names <- strsplit(l[ grep("^information", l)+2 ], ",")[[1]]
data1 <- l[ (grep("^information", l)+3):(grep("\032", l)[1]-1) ]
data2 <- l[ (grep("^money", l)+3):(grep("\032", l)[2]-1) ]
data3 <- l[ (grep("ID,", l)+1):(grep("\032", l)[3]-1) ]
data3 <- gsub("[ ]", "", data3)

tc <- textConnection(data1)
mat1 <- scan(tc)
close(tc)
tc <- textConnection(data2)
mat2 <- scan(tc)
close(tc)
tc <- textConnection(data3)
met <- read.delim(tc, header=FALSE)
close(tc)
colnames(met) <- c("Id", "Sector", "Rating")

mat1 <- matrix(mat1, sqrt(length(mat1)), byrow=TRUE)
mat2 <- matrix(mat2, sqrt(length(mat2)), byrow=TRUE)
colnames(mat1) <- rownames(mat1) <- colnames(mat2) <- rownames(mat2) <- names

library(igraph)
g1 <- graph.adjacency(mat1)
g2 <- graph.adjacency(mat2)
V(g1)$Sector <- V(g2)$Sector <- as.character(met$Sector)
V(g1)$Rating <- V(g2)$Rating <- met$Rating

g1$name <- "Knoke bureaucracies, KNOKM"
g2$name <- "Knoke bureaucracies, KNOKI"
g1$Author <- g2$Author <- "D. Knoke and J. Wood"
g1$Citation <- g2$Citation <- "Knoke D. and Wood J. (1981). Organized for action: Commitment in voluntary associations. New Brunswick, NJ: Rutgers University Press.\n\n
Knoke D. and Kuklinski J. (1982). Network analysis, Beverly Hills, CA: Sage."
g1$URL <- g2$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

knokbur <- list(KNOKM=g1, KNOKI=g2)
save(knokbur, file="/tmp/knokbur.Rdata.gz")
