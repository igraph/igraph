
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

plot.graph <- function(graph, layout=layout.random, layout.par=list(),
                       labels=NULL, label.color="darkblue",
                       label.font=NULL, label.degree=-pi/4, label.dist=0,
                       vertex.color="SkyBlue2", vertex.size=15,
                       edge.color="darkgrey", edge.width=1,
                       # SPECIFIC: #####################################
                       axes=FALSE, xlab="", ylab="",
                       ...) {

  # Interpret parameters
  layout <- i.get.layout(graph, layout, layout.par)
  vertex.color <- i.get.vertex.color(graph, vertex.color)
  vertex.size <- (1/200) * i.get.vertex.size(graph, vertex.size)
  edge.color <- i.get.edge.color(graph, edge.color)
  edge.width <- i.get.edge.width(graph, edge.width)
  label.degree <- i.get.label.degree(graph, label.degree)
  labels <- i.get.labels(graph, labels)

  # create the plot
  plot(0, 0, type="n", xlab=xlab, ylab=ylab, asp=1, xlim=c(-1,1), ylim=c(-1,1),
       axes=axes, ...)

  # norm layout to (-1, 1)
  layout <- i.layout.norm(layout, -1, 1, -1, 1)
  
  # add the edges
  # TODO: loops
  x0 <- layout[,1][get.edgelist(graph)[,1]]
  y0 <- layout[,2][get.edgelist(graph)[,1]]
  x1 <- layout[,1][get.edgelist(graph)[,2]]
  y1 <- layout[,2][get.edgelist(graph)[,2]]

  # we do this for undirected graphs also because some
  # graphics drivers do not handle 'depth' properly (or at all)
  if (length(vertex.size)!=1) {
    vsize.from <- vertex.size[get.edgelist(graph)[,1]]
    vsize.to   <- vertex.size[get.edgelist(graph)[,2]]
  } else {
    vsize.from <- vsize.to <- vertex.size
  }
  phi <- atan2(y1-y0, x1-x0)
  r <- sqrt( (x1-x0)^2 + (y1-y0)^2 )
  x1 <- x0 + (r-vsize.to)*cos(phi)
  y1 <- y0 + (r-vsize.to)*sin(phi)
  x0 <- x0 + vsize.from*cos(phi)
  y0 <- y0 + vsize.from*sin(phi)
  
  arrow.code <- ifelse(is.directed(graph), 2, 0)
  arrows(x0, y0, x1, y1, angle=20, length=0.2, code=arrow.code,
         col=edge.color, lwd=edge.width)
  rm(x0, y0, x1, y1)
  
  # add the vertices
  if (length(vertex.size)==1) { vertex.size <- rep(vertex.size, nrow(layout)) }
  symbols(x=layout[,1], y=layout[,2], bg=vertex.color,
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
}
