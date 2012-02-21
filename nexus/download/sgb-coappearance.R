
getnet <- function(url) {
  tmp <- tempdir()
  dest <- paste(sep="", tmp, "/", "anna.dat")

  download.file(url, dest)
  l <- readLines(dest)
  l <- l[!grepl("^\\*", l)]
  
  emptyline <- grep("^$", l)
  name <- l[1:(emptyline-1)]
  secs <- l[(emptyline+1):length(l)]
  
  mon <- sub(" .*$", "", name)
  name <- sub("^.. ", "", name)
  
  meet <- strsplit(sub("^[0-9\\.]+:?", "", secs), ";")
  meet2 <- lapply(meet, strsplit, ",")
  
  edges <- character()
  for (i in seq_along(meet)) {
    for (j in seq_along(meet[[i]])) {
      ev <- meet2[[i]][[j]]
      if (length(ev) >= 2) {
        edges <- c(edges, combn(ev, 2))
      }
    }
  }
  
  edges <- matrix(edges, ncol=2, byrow=TRUE)
  
  names(meet2) <- sub("([0-9\\.]+).*$", "\\1", secs)
  
  library(igraph)
  g <- graph.edgelist(edges, directed=FALSE)
  E(g)$weight <- 1
  g <- simplify(g)

  V(g)$Description <- name[match(V(g)$name, mon)]
  g$Citation <- "Donald E. Knuth, The Stanford GraphBase: A Platform for Combinatorial Computing, Addison-Wesley, Reading, MA (1993)"
  g$Author <- "Donald E. Knuth"
  g$URL <- "http://www.ctan.org/tex-archive/support/graphbase/"

  g$Coappearance <- meet2

  g
}
