
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/taro.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "taro.dat")

download.file(url, dest)

l <- readLines(paste(sep="", tmp, "/taro.dat"))

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- matrix(scan(tc), nc=22, byrow=TRUE)
close(tc)

library(igraph)
taro <- graph.adjacency(mat, mode="undirected")

taro$name <- "Schwimmer taro exchange"
taro$Author <- "E. Schwimmer"
taro$Citation <- "Hage P. and Harary F. (1983). Structural models in anthropology. Cambridge: Cambridge University Press.\n\nSchwimmer E. (1973). Exchange in the social structure of the Orokaiva. New York: St Martins."
taro$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

save(taro, file="/tmp/taro.Rdata.gz")

