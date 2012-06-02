
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/gama.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "gama.dat")

download.file(url, dest)

l <- readLines(paste(sep="", tmp, "/gama.dat"))
rowlab <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
collab <- l[ (grep("^COLUMN LABELS:", l)+1):(grep("^LEVEL LABELS:", l)-1) ]

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- matrix(scan(tc), nc=16, byrow=TRUE)
close(tc)

mat <- mat[1:16,] - mat[17:32,]
rownames(mat) <- colnames(mat) <- rowlab

library(igraph)
gama <- graph.adjacency(mat, weighted=TRUE, mode="undirected")

gama$name <- "Read highland tribes"
gama$Author <- "K. Read"
gama$Citation <- "Hage P. and Harary F. (1983). Structural models in anthropology. Cambridge: Cambridge University Press. (See p 56-60).\n\nRead K. (1954). Cultures of the central highlands, New Guinea. Southwestern Journal of Anthropology, 10, 1-43."
gama$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

save(gama, file="/tmp/gama.Rdata.gz")

