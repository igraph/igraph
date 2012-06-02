
urls <- c("http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/szcid.dat",
          "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/szcig.dat")

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "szcid.dat")
download.file(urls[1], dest)
l1 <- readLines(dest)
dest <- paste(sep="", tmp, "/", "szcig.dat")
download.file(urls[2], dest)
l2 <- readLines(dest)

getnet <- function(l) {
  n <- as.numeric(sub("^N=", "", l[grep("^N=", l)]))
  rl <- l[ (grep("^ROW LABELS:", l)+1):(grep("^COLUMN LABELS:", l)-1) ]
  cl <- l[ (grep("^COLUMN LABELS:", l)+1):(grep("^DATA:", l)-1) ]
  if (length(rl) != length(cl) || any(rl != cl)) {
    stop("Row and column names no not match")
  }
  data <- l[ (grep("^DATA:", l)+1):length(l) ]
  tc <- textConnection(data)
  mat <- matrix(scan(tc), nc=n, byrow=TRUE)
  colnames(mat) <- rownames(mat) <- rl
  graph.adjacency(mat, mode="undirected", weighted=TRUE)
}

library(igraph)
szcid <- getnet(l1)
szcig <- getnet(l2)

szcid$name <- "Stokman-Ziegler corporate interlocks, Holland"
szcig$name <- "Stokman-Ziegler corporate interlocks, Germany"
szcid$Author <- szcig$Author <- "F. Stokman and R. Ziegler"
szcid$Citation <- szcig$Citation <- "Ziegler R., Bender R. and Biehler H. (1985). Industry and banking in the German corporate network. In F. Stokman, R. Ziegler & J. Scott (eds), Networks of corporate power. Cambridge: Polity Press, 1985.\n\nStokman F., Wasseur F. and Elsas D. (1985). The Dutch network: Types of interlocks and network structure. In F. Stokman, R. Ziegler & J. Scott (eds), Networks of corporate power. Cambridge: Polity Press, 1985."
szcid$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"

szci <- list(SZCID=szcid, SZCIG=szcig)
save(szci, file="/tmp/szci.Rdata.gz")



