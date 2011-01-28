
library(igraph)

r.to.graphml <- function(input_file, output_file) {
  E <- new.env()
  load(input_file, envir=E)
  g <- get(ls(E), envir=E)
  if (is.igraph(g)) {
    g <- list(g)
  }
  f=file(output_file, open="wb")
  for (i in seq_along(g)) {
    write.graph(g[[i]], file=f, format="GraphML")
  }
  close(f)
  NULL
}
 
r.to.pajek <- function(input_file, output_file, dataset_name) {
  E <- new.env()
  load(input_file, envir=E)
  g <- get(ls(E), envir=E)
  if (is.igraph(g)) {
    g <- list(g)
  }
  f=file(output_file, open="wb")
  for (i in seq_along(g)) {
    cat(sep="", "*Network ", dataset_name, ".", names(g)[i], "\n", file=f)
    write.graph(g[[i]], file=f, format="Pajek")
  }
  close(f)
  NULL
}
