
#   IGraph R package
#   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
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

###################################################################
## API design
##
## A vertex shape is defined by two functions: the clipping function and
## the plotting function.
##
## The clipping function is called to determine where to put the
## arrowhead of a potential (incoming) incident edge. Its signature is
##    function(coords, el, params, end=c("both", "from", "to"))
## where the arguments are:
##    coords    A matrix with one row for each edge, and four columns.
##              It contains the coordinates of the end points of all
##              edges. The first two columns are the coordinates of the
##              first end points (sources, if the graph is directed),
##              the last two columns are for the other end points
##              (targets if the graph is directed).
##    el        The edge list itself, with vertex ids.
##    params    A function object to query plotting parameters.
##    end       Which end points to calculate. "both" means both,
##              "from" means the first end point, "to" the second.
## The clipping function must return the new version of "coords",
## modified according to the vertex sizes/shapes, with proper positions
## for the potential arrow heads. The positions are for the tips of the
## arrows.
##
## The plotting function plots the vertex. Its signature is
##    function(coords, v=NULL, params)
## where the arguments are
##    coords    Two column matrix, the coordinates for the vertices to draw.
##    v         The vertex ids of the vertices to draw. If NULL, then all
##              vertices are drawn.
##    params    A function object to query plotting parameters.
## 
## vertex.shapes()         - lists all vertex shapes
## vertex.shapes(shape)    - returns the clipping and plotting functions
##                           for a given vertex shape
## add.vertex.shape()      - adds a new vertex shape, the clipping and
##                           plotting functions must be given, and
##                           optionally the newly introduced plotting
##                           parameters. This function can also be used
##                           to overwrite a given vertex shape.
##
## Examples:
## add.vertex.shapes("image", clip=image.clip, plot=image.plot,
##                   parameters=c("filename"))
##
## add.vertex.shapes("triangle", clip=vertex.shapes("circle")$clip,
##                   plot=triangle.plot)
##
## add.vertex.shapes("polygon", clip=vertex.shapes("circle")$clip,
##                   plot=polygon.plot)
##
###################################################################

vertex.shapes <- function(shape=NULL) {
  if (is.null(shape)) {
    ls(.igraph.shapes)
  } else {
    ## checkScalarString(shape)
    .igraph.shapes[[shape]]
  }
}

igraph.shape.noclip <- function(coords, el, params,
                                end=c("both", "from", "to")) {
  end <- igraph.match.arg(end)

  if (end=="both") {
    coords
  } else if (end=="from") {
    coords[,1:2,drop=FALSE]
  } else {
    coords[,3:4,drop=FALSE]
  }
}

igraph.shape.noplot <- function(coords, v=NULL, params) {
  invisible(NULL)
}

add.vertex.shape <- function(shape, clip=igraph.shape.noclip,
                             plot=igraph.shape.noplot,
                             parameters=character()) {

  ## TODO
  ## checkScalarString(shape)
  ## checkFunction(clip)
  ## checkFunction(plot)
  ## checkCharacter(parameters)

  assign(shape, value=list(clip=clip, plot=plot), envir=.igraph.shapes)

  invisible(TRUE)
}

## These are the predefined shapes

.igraph.shape.circle.clip <- function(coords, el, params,
                                      end=c("both", "from", "to")) {
  
  end <- match.arg(end)

  if (length(coords)==0) { return (coords) }     

  vertex.size <- 1/200 * params("vertex", "size")

  if (end=="from") {
    phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
    vsize.from <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    res <- cbind(coords[,1] + vsize.from*cos(phi),
                 coords[,2] + vsize.from*sin(phi) )
  } else if (end=="to") {
    phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
    r <- sqrt( (coords[,3] - coords[,1])^2 + (coords[,4] - coords[,2])^2 )
    vsize.to <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    res <- cbind(coords[,1] + (r-vsize.to)*cos(phi),
                 coords[,2] + (r-vsize.to)*sin(phi) )
  } else if (end=="both") {
    phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
    r <- sqrt( (coords[,3] - coords[,1])^2 + (coords[,4] - coords[,2])^2 )
    vsize.from <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    vsize.to <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    res <- cbind(coords[,1] + vsize.from*cos(phi),
                 coords[,2] + vsize.from*sin(phi),
                 coords[,1] + (r-vsize.to)*cos(phi),
                 coords[,2] + (r-vsize.to)*sin(phi) )
  }
  res
}

.igraph.shape.circle.plot <- function(coords, v=NULL, params) {
  
  vertex.color       <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size        <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.size <- rep(vertex.size, length=nrow(coords))
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          circles=vertex.size, add=TRUE, inches=FALSE)
}

.igraph.shape.square.clip <- function(coords, el, params,
                                      end=c("both", "from", "to")) {
  end <- match.arg(end)

  if (length(coords)==0) { return (coords) }     

  vertex.size <- 1/200 * params("vertex", "size")

  square.shift <- function(x0, y0, x1, y1, vsize) {
    m <- (y0-y1)/(x0-x1)
    l <- cbind(x1-vsize/m , y1-vsize,
               x1-vsize , y1-vsize*m,
               x1+vsize/m, y1+vsize,
               x1+vsize , y1+vsize*m )
    
    v <- cbind(x1-vsize <= l[,1] & l[,1] <= x1+vsize &
               y1-vsize <= l[,2] & l[,2] <= y1+vsize,
               x1-vsize <= l[,3] & l[,3] <= x1+vsize &
               y1-vsize <= l[,4] & l[,4] <= y1+vsize,
               x1-vsize <= l[,5] & l[,5] <= x1+vsize &
               y1-vsize <= l[,6] & l[,6] <= y1+vsize,
               x1-vsize <= l[,7] & l[,7] <= x1+vsize &
               y1-vsize <= l[,8] & l[,8] <= y1+vsize)
    
    d <- cbind((l[,1]-x0)^2 + (l[,2]-y0)^2,
               (l[,3]-x0)^2 + (l[,4]-y0)^2,
               (l[,5]-x0)^2 + (l[,6]-y0)^2,
               (l[,7]-x0)^2 + (l[,8]-y0)^2)
    
    t(sapply(seq(length=nrow(l)), function(x) {
      d[x,][!v[x,]] <- Inf
      m <- which.min(d[x,])
      l[x, c(m*2-1, m*2)]
    }))
  }
    
  if (end %in% c("from", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    res <- res1 <- square.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                                vsize)
  }
  if (end %in% c("to", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    res <- res2 <- square.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                                vsize)
  }
  if (end=="both") {
    res <- cbind(res1, res2)
  }
  
  res
}

.igraph.shape.square.plot <- function(coords, v=NULL, params) {

  vertex.color       <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size        <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.size <- rep(vertex.size, length=nrow(coords))
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          squares=2*vertex.size, add=TRUE, inches=FALSE)
}  

.igraph.shape.csquare.clip <- function(coords, el, params,
                                       end=c("both", "from", "to")) {

  end <- match.arg(end)

  if (length(coords)==0) { return (coords) }     

  vertex.size <- 1/200 * params("vertex", "size")
  
  square.shift <- function(x0, y0, x1, y1, vsize) {

    l <- cbind(x1,       y1-vsize,
               x1-vsize, y1,
               x1,       y1+vsize,
               x1+vsize, y1)
    
    d <- cbind((l[,1]-x0)^2 + (l[,2]-y0)^2,
               (l[,3]-x0)^2 + (l[,4]-y0)^2,
               (l[,5]-x0)^2 + (l[,6]-y0)^2,
               (l[,7]-x0)^2 + (l[,8]-y0)^2)
    
    t(sapply(seq(length=nrow(l)), function(x) {
      m <- which.min(d[x,])
      l[x, c(m*2-1, m*2)]
    }))
  }
  
  if (end %in% c("from", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    res <- res1 <- square.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                                vsize)
  }
  if (end %in% c("to", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    res <- res2 <- square.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                                vsize)
  }
  if (end=="both") {
    res <- cbind(res1, res2)
  }
  
  res
}

.igraph.shape.csquare.plot <- .igraph.shape.square.plot

.igraph.shape.rectangle.clip <- function(coords, el, params,
                                         end=c("both", "from", "to")) {

  end <- match.arg(end)

  if (length(coords)==0) { return (coords) }     

  vertex.size <- 1/200 * params("vertex", "size")
  vertex.size2 <- 1/200 * params("vertex", "size2")
  
  rec.shift <- function(x0, y0, x1, y1, vsize, vsize2) {
    m <- (y0-y1)/(x0-x1)
    l <- cbind(x1-vsize/m,  y1-vsize2,
               x1-vsize,    y1-vsize*m,
               x1+vsize2/m, y1+vsize2,
               x1+vsize,    y1+vsize*m )
    
    v <- cbind(x1-vsize <= l[,1] & l[,1] <= x1+vsize &
               y1-vsize2 <= l[,2] & l[,2] <= y1+vsize2,
               x1-vsize <= l[,3] & l[,3] <= x1+vsize &
               y1-vsize2 <= l[,4] & l[,4] <= y1+vsize2,
               x1-vsize <= l[,5] & l[,5] <= x1+vsize &
               y1-vsize2 <= l[,6] & l[,6] <= y1+vsize2,
               x1-vsize <= l[,7] & l[,7] <= x1+vsize &
               y1-vsize2 <= l[,8] & l[,8] <= y1+vsize2)
    
    d <- cbind((l[,1]-x0)^2 + (l[,2]-y0)^2,
               (l[,3]-x0)^2 + (l[,4]-y0)^2,
               (l[,5]-x0)^2 + (l[,6]-y0)^2,
               (l[,7]-x0)^2 + (l[,8]-y0)^2)
    
    t(sapply(seq(length=nrow(l)), function(x) {
      d[x,][!v[x,]] <- Inf
      m <- which.min(d[x,])
      l[x, c(m*2-1, m*2)]
    }))
  }
  
  if (end %in% c("from", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    vsize2 <- if (length(vertex.size2)==1) {
      vertex.size2
    } else {
      vertex.size2[ el[,1] ]
    }
    res <- res1 <- rec.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                             vsize, vsize2)
  }
  if (end %in% c("to", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    vsize2 <- if (length(vertex.size2)==1) {
      vertex.size2
    } else {
      vertex.size2[ el[,1] ]
    }
    res <- res2 <- rec.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                             vsize, vsize2)
  }
  if (end=="both") {
    res <- cbind(res1, res2)
  }
  
  res
}

.igraph.shape.rectangle.plot <- function(coords, v=NULL, params) {

  vertex.color       <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size        <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.size <- rep(vertex.size, length=nrow(coords))   
  vertex.size2       <- 1/200 * params("vertex", "size2")
  if (length(vertex.size2) != 1 && !is.null(v)) {
    vertex.size2 <- vertex.size2[v]
  }
  vertex.size <- cbind(vertex.size, vertex.size2)
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          rectangles=2*vertex.size, add=TRUE, inches=FALSE)
}

.igraph.shape.crectangle.clip <- function(coords, el, params,
                                          end=c("both", "from", "to")) {

  end <- match.arg(end)

  if (length(coords)==0) { return (coords) }     

  vertex.size <- 1/200 * params("vertex", "size")
  vertex.size2 <- 1/200 * params("vertex", "size2")
    
  rec.shift <- function(x0, y0, x1, y1, vsize, vsize2) {
    
    l <- cbind(x1,       y1-vsize2,
               x1-vsize, y1,
               x1,       y1+vsize2,
               x1+vsize, y1)
    
    d <- cbind((l[,1]-x0)^2 + (l[,2]-y0)^2,
               (l[,3]-x0)^2 + (l[,4]-y0)^2,
               (l[,5]-x0)^2 + (l[,6]-y0)^2,
               (l[,7]-x0)^2 + (l[,8]-y0)^2)
      
    t(sapply(seq(length=nrow(l)), function(x) {
      m <- which.min(d[x,])
      l[x, c(m*2-1, m*2)]
    }))
  }
  
  if (end %in% c("from", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    vsize2 <- if (length(vertex.size2)==1) {
      vertex.size2
    } else {
      vertex.size2[ el[,1] ]
    }
    res <- res1 <- rec.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                             vsize, vsize2)
  }
  if (end %in% c("to", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    vsize2 <- if (length(vertex.size2)==1) {
      vertex.size2
    } else {
      vertex.size2[ el[,1] ]
    }
    res <- res2 <- rec.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                             vsize, vsize2)
  }
  if (end=="both") {
    res <- cbind(res1, res2)
  }
  
  res
}

.igraph.shape.crectangle.plot <- .igraph.shape.rectangle.plot

.igraph.shape.vrectangle.clip <- function(coords, el, params,
                                          end=c("both", "from", "to")) {

  end <- match.arg(end)

  if (length(coords)==0) { return (coords) }     

  vertex.size <- 1/200 * params("vertex", "size")
  vertex.size2 <- 1/200 * params("vertex", "size2")
  
  rec.shift <- function(x0, y0, x1, y1, vsize, vsize2) {
    
    l <- cbind(x1-vsize, y1, x1+vsize, y1)
    
    d <- cbind((l[,1]-x0)^2 + (l[,2]-y0)^2,
               (l[,3]-x0)^2 + (l[,4]-y0)^2)
    
    t(sapply(seq(length=nrow(l)), function(x) {
      m <- which.min(d[x,])
      l[x, c(m*2-1, m*2)]
    }))
  }

  if (end %in% c("from", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    vsize2 <- if (length(vertex.size2)==1) {
      vertex.size2
    } else {
      vertex.size2[ el[,1] ]
    }
    res <- res1 <- rec.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                             vsize, vsize2)
  }
  if (end %in% c("to", "both")) {
    vsize <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    vsize2 <- if (length(vertex.size2)==1) {
      vertex.size2
    } else {
      vertex.size2[ el[,1] ]
    }
    res <- res2 <- rec.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                             vsize, vsize2)
  }
  if (end=="both") {
    res <- cbind(res1, res2)
  }

  res
}    

.igraph.shape.vrectangle.plot <- .igraph.shape.rectangle.plot

.igraph.shape.none.clip <- .igraph.shape.circle.clip

.igraph.shape.none.plot <- function(coords, v=NULL, params) {
  ## does not plot anything at all
  invisible(NULL)
}

mypie <- function(x, y, values, radius, edges=200, col=NULL, angle=45,
                  density=NULL, border=NULL, lty=NULL, init.angle=90, ...) {
  values <- c(0, cumsum(values)/sum(values))
  dx <- diff(values)
  nx <- length(dx)
  twopi <- 2 * pi
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
    else par("fg")
  col <- rep(col, length.out = nx)
  border <- rep(border, length.out = nx)
  lty <- rep(lty, length.out = nx)
  angle <- rep(angle, length.out = nx)
  density <- rep(density, length.out = nx)
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p))
  }
  for (i in 1:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(values[i], values[i + 1], length.out = n))
    polygon(x+c(P$x, 0), y+c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i], ...)
  }
}

.igraph.shape.pie.clip <- function(coords, el, params,
                                   end=c("both", "from", "to")) {

  end <- match.arg(end)

  if (length(coords)==0) { return (coords) }     

  vertex.size <- 1/200 * params("vertex", "size")

  if (end=="from") {
    phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
    vsize.from <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    res <- cbind(coords[,1] + vsize.from*cos(phi),
                 coords[,2] + vsize.from*sin(phi) )
  } else if (end=="to") {
    phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
    r <- sqrt( (coords[,3] - coords[,1])^2 + (coords[,4] - coords[,2])^2 )
    vsize.to <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    res <- cbind(coords[,1] + (r-vsize.to)*cos(phi),
                 coords[,2] + (r-vsize.to)*sin(phi) )
  } else if (end=="both") {
    phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
    r <- sqrt( (coords[,3] - coords[,1])^2 + (coords[,4] - coords[,2])^2 )
    vsize.from <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,1] ]
    }
    vsize.to <- if (length(vertex.size)==1) {
      vertex.size
    } else {
      vertex.size[ el[,2] ]
    }
    res <- cbind(coords[,1] + vsize.from*cos(phi),
                 coords[,2] + vsize.from*sin(phi),
                 coords[,1] + (r-vsize.to)*cos(phi),
                 coords[,2] + (r-vsize.to)*sin(phi) )
  }

  res
}

.igraph.shape.pie.plot <- function(coords, v=NULL, params) {

  getparam <- function(pname) {
    p <- params("vertex", pname)
    if (length(p) != 1 && !is.null(v)) {
      p <- p[v]
    }
    p
  }
  vertex.color       <- getparam("color")
  vertex.frame.color <- getparam("frame.color")
  vertex.size        <- rep(1/200 * getparam("size"), length=nrow(coords))
  vertex.pie         <- getparam("pie")
  vertex.pie.color   <- getparam("pie.color")
  vertex.pie.angle   <- getparam("pie.angle")
  vertex.pie.density <- getparam("pie.density")
  vertex.pie.lty     <- getparam("pie.lty")

  for (i in seq_len(nrow(coords))) {
    pie <- if(length(vertex.pie)==1) {
      vertex.pie[[1]]
    } else {
      vertex.pie[[i]]
    }
    col <- if (length(vertex.pie.color)==1) {
      vertex.pie.color[[1]]
    } else {
      vertex.pie.color[[i]]
    }
    mypie(x=coords[i,1], y=coords[i,2], pie,
          radius=vertex.size[i], edges=200, col=col,
          angle=na.omit(vertex.pie.angle[c(i,1)])[1],
          density=na.omit(vertex.pie.density[c(i,1)])[1],
          border=na.omit(vertex.frame.color[c(i,1)])[1],
          lty=na.omit(vertex.pie.lty[c(i,1)])[1])
  }
}

.igraph.shapes <- new.env()
.igraph.shapes[["circle"]] <- list(clip=.igraph.shape.circle.clip,
                                   plot=.igraph.shape.circle.plot)
.igraph.shapes[["square"]] <- list(clip=.igraph.shape.square.clip,
                                   plot=.igraph.shape.square.plot)
.igraph.shapes[["csquare"]] <- list(clip=.igraph.shape.csquare.clip,
                                    plot=.igraph.shape.csquare.plot)
.igraph.shapes[["rectangle"]] <- list(clip=.igraph.shape.rectangle.clip,
                                      plot=.igraph.shape.rectangle.plot)
.igraph.shapes[["crectangle"]] <- list(clip=.igraph.shape.crectangle.clip,
                                       plot=.igraph.shape.crectangle.plot)
.igraph.shapes[["vrectangle"]] <- list(clip=.igraph.shape.vrectangle.clip,
                                       plot=.igraph.shape.vrectangle.plot)
.igraph.shapes[["none"]] <- list(clip=.igraph.shape.none.clip,
                                 plot=.igraph.shape.none.plot)
.igraph.shapes[["pie"]] <- list(clip=.igraph.shape.pie.clip,
                                plot=.igraph.shape.pie.plot)

