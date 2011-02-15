
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/davis.dat"
tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "davis.dat")

download.file(url, dest)

l <- readLines(paste(sep="", tmp, "/davis.dat"))
rowlab <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
collab <- l[ (grep("^COLUMN LABELS:", l)+1):(grep("^DATA:", l)-1) ]

data <- l[ (grep("^DATA:", l)+1):length(l) ]
data <- paste(data, collapse="\n")

tc <- textConnection(data)
mat <- as.matrix(read.table(tc))
close(tc)

rownames(mat) <- rowlab
colnames(mat) <- collab

library(igraph)
g <- graph.incidence(mat)

g$name <- "Davis's Southern club women"
g$Author <- "Allison Davis"
g$Citation <- "Allison Davis, Burleigh Bradford Gardner, Mary R. Gardner, Deep south: a social anthropological study of caste and class, University of Chicago Press (1941)."
g$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"
g$Description <- "These data were collected by Davis et al in the 1930s. They represent observed attendance at 14 social events by 18 Southern women."

V(g)$type <- as.numeric(grepl("^E[0-9][0-9]?$", V(g)$name))

Davis <- g
save(Davis, file="/tmp/Davis.Rdata.gz")

