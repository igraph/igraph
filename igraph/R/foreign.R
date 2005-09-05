
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

###################################################################
# Reading foreign file formats
###################################################################

read.graph <- function(file, format="edgelist", ...) {

  res <- switch(format,
#                "pajek"=read.graph.pajek(file, format, ...),
                "edgelist"=read.graph.edgelist(file, format, ...),
                stop(paste("Unknown file format:",format))
                )
  res
}

write.graph <- function(graph, file, format="edgelist", ...) {
  
  res <- switch(format,
#                "pajek"=write.graph.pajek(graph, file, format, ...),
                "edgelist"=write.graph.edgelist(graph, file, format, ...),
                "lgl"=write.graph.lgl(graph, file, format, ...),
                stop(paste("Unknown file format:",format))
                )
  res
}

################################################################
# Pajek format
# see http://vlado.fmf.uni-lj.si/pub/networks/pajek/
################################################################

read.graph.pajek <- function(filename, format="pajek", attributes=TRUE, ...) {

##   lines <- readLines(filename)
##   res <- .Call("REST_import_pajek", igraph.c.interface, lines, list(...), attributes)

  res
}

################################################################
# Plain edge list format, not sorted
################################################################

read.graph.edgelist <- function(filename, format="edgelist", ...) {

  res <- graph.empty(...)

  if (is.character(filename)) {
    filename <- file(filename)
  }
  if (!isOpen(filename)) {
    open(filename)
  }
  
  edges <- scan(filename, nmax=20000)
  while(length(edges)>0) {
    m <- max(edges)+1
    v <- vcount(res)
    if (m>v) {
      res <- add.vertices(res, m-v)
    }
    res <- add.edges(res, edges)
    edges <- scan(filename, nmax=20000)
  }

  res
}

################################################################
# Write an edgelist format, sorted
################################################################

write.graph.edgelist <- function(graph, file, format="edgelist",
                                 ...) {

  if (is.character(file)) {
    file <- file(file, open="w+")
  }
  
  closeit <- FALSE
  if (!isOpen(file)) {
    file <- open(file)
    closeit <- TRUE
  }
  
  vc <- vcount(graph)
  for (i in 0:(vc-1)) {
    neis <- neighbors(graph, i, "out")
    if (!is.directed(graph)) {
      no.loops <- sum(neis==i)
      neis <- c(neis[neis > i], rep(i, no.loops/2))
    }
    neis <- sort(neis)
    for (n in neis) {
      cat(i,    file=file)
      cat(" ",  file=file)
      cat(n,    file=file)
      cat("\n", file=file)
    }
  }

  if (closeit) {
    close(file)
  }

  invisible(NULL)
}

################################################################
# Write an LGL format, quite simple
################################################################

write.graph.lgl <- function(graph, file, format="lgl", ...) {
  
  if (is.character(file)) {
    file <- file(file, open="w+")
  }
  
  closeit <- FALSE
  if (!isOpen(file)) {
    file <- open(file)
    closeit <- TRUE
  }

  vc <- vcount(graph)
  for (i in 0:(vc-1)) {
    neis <- neighbors(graph, i, "out")
    if (!is.directed(graph)) {
      neis <- neis [ neis < i ]
    }
    if (length(neis) > 0) {
      neis <- sort(neis)
      cat("# ", file=file)
      cat(i+1,    file=file)
      cat("\n", file=file)
      for (n in neis) {
        cat(n+1,    file=file)
        cat("\n", file=file)
      }
    }
  }

  invisible(NULL)
}
      
