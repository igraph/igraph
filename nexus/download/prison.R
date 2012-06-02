
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/prison.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "prison.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/prison.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

mat <- matrix(mat, sqrt(length(mat)), byrow=TRUE)

library(igraph)
g <- graph.adjacency(mat)

g$name <- "Gagnon & MacRae prison"
g$Author <- "John Gagnon"
g$Citation <- "MacRae J. (1960). Direct factor analysis of sociometric data. Sociometry, 23, 360-371."
g$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

prison <- g
save(prison, file="/tmp/prison.Rdata.gz")


