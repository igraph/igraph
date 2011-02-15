
library(igraph)

get.graph.list <- function(input_file) {
  E <- new.env()
  load(input_file, envir=E)
  g <- get(ls(E), envir=E)
  if (is.igraph(g)) {
    g <- list(g)
    names(g) <- strsplit(basename(input_file), ".", fixed=TRUE)[[1]][1]
  }
  g
}

zip.output <- function(tmp, output_file) {
  system(paste(sep="", "cd ", tmp, ";", "zip ", basename(output_file), " *"))
  file.copy(paste(sep="", tmp, "/", basename(output_file)), output_file)
  NULL
}

create.tmp <- function() {
  tmp <- tempdir()
  ntmp <- paste(sep="", tmp, "/nexus")
  unlink(ntmp, recursive=TRUE)
  dir.create(ntmp)
  ntmp
}

r.to.graphml <- function(input_file, output_file) {

  g <- get.graph.list(input_file)

  tmp <- create.tmp()
  sapply(names(g), function(n) {
    f <- paste(sep="", tmp, "/", n, ".GraphML")
    write.graph(g[[n]], file=f, format="GraphML")
  })

  zip.output(tmp, output_file)
  
  NULL
}
 
r.to.pajek <- function(input_file, output_file) {
  g <- get.graph.list(input_file)

  tmp <- create.tmp()
  sapply(names(g), function(n) {
    f <- paste(sep="", tmp, "/", n, ".net")
    write.graph(g[[n]], file=f, format="Pajek")
  })

  zip.output(tmp, output_file)  
  NULL
}

print_check <- function(res) {
  sapply(names(res), function(x) {
    cat(sep="", x, ":", paste(res[[x]], collapse=";"), "\n")
  })
  invisible(NULL)
}

r.check <- function(input_file, tags, netnames, netno, meta) {

  res <- list()

  ## Parse input
  tags=strsplit(tags, ";")[[1]]
  
  ## File exists
  if (! (res$`Data file exists` <- file.exists(input_file)) ) {
    return(print_check(res))
  }

  ## Can be loaded
  env <- new.env()
  r <- try( { load(input_file, envir=env) }, silent=TRUE )
  if (! (res$`Data file can be loaded` <- ! inherits(r, "try-error")) ) {
    return(print_check(res))
  }
  
  ## Contains a single object, an igraph graph or a list of igraph graphs
  r <- length(ls(env)) == 1
  if (r) {
    G <- get(ls(env), envir=env)
    r <- r && (is.igraph(G) || is.list(G))
    if (r && !is.igraph(G)) { r <- r && all(sapply(G, is.igraph)) }
  }
  if (! (res$`File contains proper data` <- r) ) {
    return(print_check(res))
  }

  ## Number of vertices and edges
  if (is.igraph(G)) {
    G <- list(G)
  } else {
    G <- G[strsplit(netnames, ";")[[1]]]
  }
  res$`Vertices` <- sapply(G, vcount)
  res$`Edges` <- sapply(G, ecount)
  
  ## Tags
  if ('directed' %in% tags) {
    res$`Tags, directed` <- any(sapply(G, is.directed))
  }
  if ('undirected' %in% tags) {
    res$`Tags, undirected` <- any(!sapply(G, is.directed))
  }
  if ('weighted' %in% tags) {
    res$`Tags, weighted` <- any(sapply(G, function(x)
                                       'weight' %in% list.edge.attributes(x)))
  }
  if ('bipartite' %in% tags) {
    res$`Tags, bipartite` <-any(sapply(G, function(x)
                                       'type' %in% list.vertex.attributes(x)))
  }

  ## Network properties, directed, weighted and bipartite, TODO
  
  ## Metadata
  meta <- lapply(strsplit(meta, ";;")[[1]],
                 function(x) strsplit(x, ";")[[1]])
  va <- lapply(G, list.vertex.attributes)
  ea <- lapply(G, list.edge.attributes)
  ga <- lapply(G, list.graph.attributes)
  netids <- strsplit(netno, ";")[[1]]
  for (m in meta) {
    aa <- switch(m[[1]], "vertex"=va, "edge"=ea, "graph"=ga)
    res[[paste(sep="", "Metadata, ", m[[1]], ", ", m[[2]])]] <-
      if (m[[3]] == "NULL") {
        all(sapply(aa, function(x) m[[2]] %in% x))
      } else {
        m[[2]] %in% aa[[ match(m[[3]], netids) ]]
      }
    
  }
  
  ## Finished
  print_check(res)
}
