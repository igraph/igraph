
url <- "http://www.itp.uni-bremen.de/complex/ebel_bornholdt_email.net.gz"

tmp <- tempdir()
dest <- paste(sep="", tmp, "/", "kielemail.txt.gz")

download.file(url, dest)
system(paste("cd ", tmp, "; zcat kielemail.txt.gz | sed -e 's/  */ /g' >kielemail.txt"))

df <- read.table(paste(sep="", tmp, "/kielemail.txt"), header=FALSE,
                 sep=" ", row.names=NULL, colClasses="character",
                 as.is=TRUE, stringsAsFactors=FALSE)

edges <- data.frame(source=I(paste(sep="", df$V5, df$V6, df$V7)),
                    target=I(paste(sep="", df$V9, df$V10, df$V11)),
                    Date=I(paste(df$V1, df$V2, df$V3)))
vertices <- data.frame(id=I(c(edges$source, edges$target)),
                       Internal=c(df$V4, df$V8))
vertices <- aggregate(vertices$Internal,
    by=list(id=vertices$id), FUN=function(xs) { paste(unique(xs), collapse=" ") })
names(vertices) <- c("id", "Internal")
vertices$Internal <- vertices$Internal == "s"

library(igraph)
g <- graph.data.frame(edges, vertices=vertices)
E(g)$weight <- 1
g <- simplify(g, edge.attr.comb=list(weight="sum"), remove.loops=FALSE)

g$Name <- "University of Kiel email dataset"
g$Author <- "Holger Ebel, Lutz-Ingo Mielsch and Stefan Bornholdt"
g$Citation <- "H. Ebel, L-I. Mielsch and S. Bornholdt: Scale-free topology of e-mail networks. Phys. Rev. E 66, 035103(R) (2002)."
g$URL <- "http://www.itp.uni-bremen.de/complex/email_net.html"

kielemail <- g
save(kielemail, file="/tmp/kielemail.Rdata.gz")
