
url1 <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/krackad.dat"
url2 <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/krackfr.dat"

tmp <- tempdir()

dest <- paste(sep="", tmp, "/", "krackad.dat")
download.file(url1, dest)
l1 <- readLines(paste(sep="", tmp, "/krackad.dat"))

dest <- paste(sep="", tmp, "/", "krackfr.dat")
download.file(url2, dest)
l2 <- readLines(paste(sep="", tmp, "/krackfr.dat"))

data1 <- l1[ (grep("^DATA:", l1)+1):length(l1) ]
data2 <- l2[ (grep("^DATA:", l2)+1):length(l2) ]

tc <- textConnection(data1)
mat1 <- scan(tc)
close(tc)

tc <- textConnection(data2)
mat2 <- scan(tc)
close(tc)

mat1 <- apply(matrix(mat1, nrow=21, byrow=TRUE), 1, list)
mat2 <- apply(matrix(mat2, nrow=21, byrow=TRUE), 1, list)

mat1 <- lapply(mat1, "[[", 1)
mat2 <- lapply(mat2, "[[", 1)

mat1 <- lapply(mat1, matrix, nc=21, byrow=TRUE)
mat2 <- lapply(mat2, matrix, nc=21, byrow=TRUE)

library(igraph)
KRACKAD <- lapply(mat1, graph.adjacency, mode="directed")
KRACKFR <- lapply(mat2, graph.adjacency, mode="directed")

for (i in seq_along(KRACKAD)) {
  KRACKAD[[i]]$name <- paste(sep="", "Krackhardt office CSS, advice, #", i)
  KRACKFR[[i]]$name <-
    paste(sep="", "Krackhardt office CSS, friendship, #", i)
  KRACKAD[[i]]$Author <- KRACKFR[[i]]$Author <- "David Krackhardt"
  KRACKAD[[i]]$Citation <- KRACKFR[[i]]$Citation <- "Krackhardt D. (1987). Cognitive social structures. Social Networks, 9, 104-134."
  KRACKAD[[i]]$URL <- KRACKFR[[i]]$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"
}

names(KRACKAD) <- paste(sep="", "KRACKAD-", seq_along(KRACKAD))
names(KRACKFR) <- paste(sep="", "KRACKFR-", seq_along(KRACKFR))

krackoffice <- c(KRACKAD, KRACKFR)

save(krackoffice, file="/tmp/krackoffice.Rdata.gz")

stop()

library(RSQLite)
con <- dbConnect("SQLite", dbname="../db/test.db")
for (i in seq_along(KRACKAD)) {
  
  sid <- paste(sep="", "'KRACKAD-", i, "'")
  desc <- paste(sep="", "'Advice, as viewed by actor #", i, "'")
  values <- paste(sep=",", 36, i, sid, desc, vcount(KRACKAD[[i]]),
                  ecount(KRACKAD[[i]]), sid, 1, 0, 0)
  dbSendQuery(con, paste(sep="",
                    'INSERT INTO network (dataset, id, sid, description,
                                    vertices, edges, filename, directed,
                                    bipartite, weighted)
                     VALUES (', values, ')'))
}

for (i in seq_along(KRACKFR)) {
  
  sid <- paste(sep="", "'KRACKFR-", i, "'")
  desc <- paste(sep="", "'Friendship, as viewed by actor #", i, "'")
  values <- paste(sep=",", 36, i+21, sid, desc, vcount(KRACKFR[[i]]),
                  ecount(KRACKFR[[i]]), sid, 1, 0, 0)
  dbSendQuery(con, paste(sep="",
                    'INSERT INTO network (dataset, id, sid, description,
                                    vertices, edges, filename, directed,
                                    bipartite, weighted)
                     VALUES (', values, ')'))
}

