
#   IGraph R package
#   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

plot.igraph <- function(x, layout=layout.random, layout.par=list(),
                       labels=NULL, label.color="darkblue",
                       label.font=NULL, label.degree=-pi/4, label.dist=0,
                       vertex.color="SkyBlue2", vertex.size=15,
                       edge.color="darkgrey", edge.width=1,
                       edge.labels=NA, 
                       vertex.frame.color="black", loop.cex=1.0,
                       # SPECIFIC: #####################################
                       axes=FALSE, xlab="", ylab="",
                       xlim=c(-1,1), ylim=c(-1,1),
                       ...) {

  graph <- x
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  # Interpret parameters
  layout <- i.get.layout(graph, layout, layout.par)
  vertex.color <- i.get.vertex.color(graph, vertex.color)
  vertex.frame.color <- i.get.vertex.frame.color(graph, vertex.frame.color)
  vertex.size <- (1/200) * i.get.vertex.size(graph, vertex.size)
  edge.color <- i.get.edge.color(graph, edge.color)
  edge.width <- i.get.edge.width(graph, edge.width)
  label.degree <- i.get.label.degree(graph, label.degree)
  labels <- i.get.labels(graph, labels)
  edge.labels <- i.get.edge.labels(graph, edge.labels)

  # create the plot
  plot(0, 0, type="n", xlab=xlab, ylab=ylab, asp=1, xlim=xlim, ylim=ylim,
       axes=axes, ...)

  # norm layout to (-1, 1)
  layout <- i.layout.norm(layout, -1, 1, -1, 1)
  
  # add the edges
  el <- get.edgelist(graph)
  loops.v <- el[,1] [ el[,1] == el[,2] ] + 1
  loops.e <- which(el[,1] == el[,2])
  el <- el[el[,1] != el[,2],]
  dim(el) <- c(length(el)/2, 2)
  x0 <- layout[,1][el[,1]+1]
  y0 <- layout[,2][el[,1]+1]
  x1 <- layout[,1][el[,2]+1]
  y1 <- layout[,2][el[,2]+1]

  # we do this for undirected graphs also because some
  # graphics drivers do not handle 'depth' properly (or at all)
  if (length(vertex.size)!=1) {
    vsize.from <- vertex.size[el[,1]+1]
    vsize.to   <- vertex.size[el[,2]+1]
  } else {
    vsize.from <- vsize.to <- vertex.size
  }
  rm (el)
  phi <- atan2(y1-y0, x1-x0)
  r <- sqrt( (x1-x0)^2 + (y1-y0)^2 )
  x1 <- x0 + (r-vsize.to)*cos(phi)
  y1 <- y0 + (r-vsize.to)*sin(phi)
  x0 <- x0 + vsize.from*cos(phi)
  y0 <- y0 + vsize.from*sin(phi)
  
  # add the loop edges
  if (length(loops.e) > 0) {
    ec <- edge.color
    if (length(ec)>1) { ec <- ec[loops.e] }
    symbols(x=layout[loops.v,1]+.08*loop.cex,
            y=layout[loops.v,2], fg=ec, bg="white",
            circles=rep(.08*loop.cex, length(loops.e)), add=TRUE, inches=FALSE)
    if (is.directed(graph)) {
      # arrows
      lr0 <- vertex.size
      if (length(vertex.size)>1) { lr0 <- vertex.size[loops.v] }
      lx0 <- layout[loops.v,1]
      ly0 <- layout[loops.v,2]
      lr1 <- .08*loop.cex
      lcx <- (lr0^2+2*lx0*lr1)/(2*lr1)
      lcy <- ly0+sqrt(ly0^2-((lcx-lx0)^2+ly0^2-lr0^2))
      arrows(lcx+.001, lcy+.001, lcx, lcy, angle=20, length=0.2,
             col=edge.color, lwd=edge.width)
    }
  }
  
  arrow.code <- ifelse(is.directed(graph), 2, 0)
  if (length(x0) != 0) {    
    arrows(x0, y0, x1, y1, angle=20, length=0.2, code=arrow.code,
           col=edge.color, lwd=edge.width)
    x <- (x0+x1)/2
    y <- (y0+y1)/2
    text(x, y, labels=edge.labels, col=label.color)
  }
  
  rm(x0, y0, x1, y1)
  
  # add the vertices
  if (length(vertex.size)==1) { vertex.size <- rep(vertex.size, nrow(layout)) }
  symbols(x=layout[,1], y=layout[,2], bg=vertex.color, fg=vertex.frame.color,
          circles=vertex.size, add=TRUE, inches=FALSE)

  # add the labels
  if (!is.null(label.font)) par(family=label.font)
  par(xpd=TRUE)
  x <- layout[,1]+label.dist*cos(-label.degree)* 
    (vertex.size+6*8*log10(nchar(labels)+1))/200
  y <- layout[,2]+label.dist*sin(-label.degree)*
    (vertex.size+6*8*log10(nchar(labels)+1))/200
  text(x, y, labels=labels, col=label.color)
  rm(x, y)
  invisible(NULL)
}

rglplot        <- function(x, layout=layout.random, layout.par=list(),
                           labels=NULL, label.color="darkblue",
                           label.font=NULL, label.degree=-pi/4, label.dist=0,
                           vertex.color="SkyBlue2", vertex.size=15,
                           edge.color="darkgrey", edge.width=1,
                           edge.labels=NA, 
                           # SPECIFIC: #####################################
                           ...)
  UseMethod("rglplot", x)


rglplot.igraph <- function(x, layout=layout.random, layout.par=list(),
                           labels=NULL, label.color="darkblue",
                           label.font=NULL, label.degree=-pi/4, label.dist=0,
                           vertex.color="SkyBlue2", vertex.size=15,
                           edge.color="darkgrey", edge.width=1,
                           edge.labels=NA, 
                           # SPECIFIC: #####################################
                           ...) {

  require(rgl)
  
  graph <- x
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  # Interpret parameters
  layout <- i.get.layout(graph, layout, layout.par)
  vertex.color <- i.get.vertex.color(graph, vertex.color)
  vertex.size <- (1/200) * i.get.vertex.size(graph, vertex.size)
  edge.color <- i.get.edge.color(graph, edge.color)
  edge.width <- i.get.edge.width(graph, edge.width)
  label.degree <- i.get.label.degree(graph, label.degree)
  labels <- i.get.labels(graph, labels)
  edge.labels <- i.get.edge.labels(graph, edge.labels)

  # norm layout to (-1, 1)
  layout <- i.layout.norm(layout, -1, 1, -1, 1, -1, 1)
  
  # add the edges
  # TODO: loops
  x0 <- layout[,1][get.edgelist(graph)[,1]+1]
  y0 <- layout[,2][get.edgelist(graph)[,1]+1]
  z0 <- layout[,3][get.edgelist(graph)[,1]+1]
  x1 <- layout[,1][get.edgelist(graph)[,2]+1]
  y1 <- layout[,2][get.edgelist(graph)[,2]+1]
  z1 <- layout[,3][get.edgelist(graph)[,2]+1]

  # we do this for undirected graphs also because some
  # graphics drivers do not handle 'depth' properly (or at all)
  if (length(vertex.size)!=1) {
    vsize.from <- vertex.size[get.edgelist(graph)[,1]+1]
    vsize.to   <- vertex.size[get.edgelist(graph)[,2]+1]
  } else {
    vsize.from <- vsize.to <- vertex.size
  }

  rgl.lines(as.numeric(t(matrix( c(x0,x1), nc=2))),
            as.numeric(t(matrix( c(y0,y1), nc=2))),
            as.numeric(t(matrix( c(z0,z1), nc=2))),
            col=edge.color, size=edge.width)

  # add the vertices
  if (length(vertex.size)==1) { vertex.size <- rep(vertex.size, nrow(layout)) }
  rgl.spheres(layout[,1], layout[,2], layout[,3], radius=vertex.size,
              col=vertex.color)

  # add the labels
  if (!is.na(labels)) {
    x <- layout[,1]+label.dist*cos(-label.degree)* 
      (vertex.size+6*10*log10(nchar(labels)+1))/200
    y <- layout[,2]+label.dist*sin(-label.degree)*
      (vertex.size+6*10*log10(nchar(labels)+1))/200
    z <- layout[,3]
    rgl.texts(x,y,z, labels, col=label.color, justify="left")
  }
  
  if (!is.na(edge.labels)) {
    rgl.texts((x0+x1)/2, (y0+y1)/2, (z0+z1)/2, edge.labels,
              col=label.color)
  }
  
  invisible(NULL)
}

