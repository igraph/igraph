
url <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/newfrat.dat"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "newfrat.dat")
download.file(url, dest)
l <- readLines(paste(sep="", tmp, "/newfrat.dat"))

lab <- l[ (grep("^LEVEL LABELS:", l)+1):(grep("^DATA", l)-1) ]
data <- l[ (grep("^DATA:", l)+1):length(l) ]

tc <- textConnection(data)
mat <- scan(tc)
close(tc)

mat <- apply(matrix(mat, nrow=15, byrow=TRUE), 1, list)
mat <- lapply(mat, "[[", 1)
mat <- lapply(mat, matrix, nc=17, byrow=TRUE)

library(igraph)
newfrat <- lapply(mat, graph.adjacency, mode="directed", weighted=TRUE)
names(newfrat) <- lab

for (i in seq_along(newfrat)) {
  newfrat[[i]]$name <- paste(sep="", "Newcomb fraternity, ", lab[i])
  newfrat[[i]]$Author <- "Theodore Newcomb"
  newfrat[[i]]$Citation <- "Newcomb T. (1961). The acquaintance process. New York: Holt, Reinhard & Winston"
  newfrat[[i]]$URL <- "http://vlado.fmf.uni-lj.si/pub/networks/data/ucinet/ucidata.htm"
}

save(newfrat, file="/tmp/newfrat.Rdata.gz")

stop()

library(RSQLite)
con <- dbConnect("SQLite", dbname="../db/test.db")

for (i in seq_along(newfrat)) {
  
  sid <- paste(sep="", "'", lab[i], "'")
  values <- paste(sep=",", 37, i, sid, vcount(newfrat[[i]]),
                  ecount(newfrat[[i]]), sid, 1, 0, 1)
  dbSendQuery(con, paste(sep="",
                    'INSERT INTO network (dataset, id, sid, 
                                    vertices, edges, filename, directed,
                                    bipartite, weighted)
                     VALUES (', values, ')'))
}

