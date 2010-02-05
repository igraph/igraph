
#   IGraph R package
#   Copyright (C) 2003-2008  Gabor Csardi <csardi@rmki.kfki.hu>
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


.igraph.shape.circle <- function(coords, el=NULL, v=NULL, mode=c("clip", "plot"),
                                 params, end=c("both", "from", "to")) {

  mode=match.arg(mode)
  end =match.arg(end)

  #####################################################################
  ## clipping mode
  
  if (mode=="clip") {
    if (length(coords)==0) { return (coords) }     

    vertex.size <- 1/200 * params("vertex", "size")

    if (end=="from") {
      phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
      vsize.from <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,1]+1 ]
      }
      res <- cbind(coords[,1] + vsize.from*cos(phi),
                   coords[,2] + vsize.from*sin(phi) )
    } else if (end=="to") {
      phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
      r <- sqrt( (coords[,3] - coords[,1])^2 + (coords[,4] - coords[,2])^2 )
      vsize.to <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,2]+1 ]
      }
      res <- cbind(coords[,1] + (r-vsize.to)*cos(phi),
                   coords[,2] + (r-vsize.to)*sin(phi) )
    } else if (end=="both") {
      phi <- atan2(coords[,4] - coords[,2], coords[,3] - coords[,1])
      r <- sqrt( (coords[,3] - coords[,1])^2 + (coords[,4] - coords[,2])^2 )
      vsize.from <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,1]+1 ]
      }
      vsize.to <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,2]+1 ]
      }
      res <- cbind(coords[,1] + vsize.from*cos(phi),
                   coords[,2] + vsize.from*sin(phi),
                   coords[,1] + (r-vsize.to)*cos(phi),
                   coords[,2] + (r-vsize.to)*sin(phi) )
    }

    res

  #####################################################################
  ## plotting mode
    
  } else if (mode=="plot") {
    vertex.color       <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v+1]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v+1]
    }
    vertex.size        <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v+1]
    }
    vertex.size <- rep(vertex.size, length=nrow(coords))
    
    symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
             circles=vertex.size, add=TRUE, inches=FALSE)
  }
  
}

.igraph.shape.square <- function(coords, el, e, v=NULL, mode=c("clip", "plot"),
                                 params, end=c("both", "from", "to")) {
  mode=match.arg(mode)
  end =match.arg(end)

  #####################################################################
  ## clipping mode

  if (mode=="clip") {
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
        vertex.size[ el[,1]+1 ]
      }
      res <- res1 <- square.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                                  vsize)
    }
    if (end %in% c("to", "both")) {
      vsize <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,2]+1 ]
      }
      res <- res2 <- square.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                                  vsize)
    }
    if (end=="both") {
      res <- cbind(res1, res2)
    }

    res

  #####################################################################
  ## plotting mode
    
  } else if (mode=="plot") {
    vertex.color       <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v+1]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v+1]
    }
    vertex.size        <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v+1]
    }
    vertex.size <- rep(vertex.size, length=nrow(coords))

    symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
            squares=2*vertex.size, add=TRUE, inches=FALSE)
  }  
  
}

.igraph.shape.csquare <- function(coords, el=NULL, v=NULL, mode=c("clip", "plot"),
                                  params, end=c("both", "from", "to")) {

  mode=match.arg(mode)
  end=match.arg(end)

  #####################################################################
  ## clipping mode

  if (mode=="clip") {
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
        vertex.size[ el[,1]+1 ]
      }
      res <- res1 <- square.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                                  vsize)
    }
    if (end %in% c("to", "both")) {
      vsize <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,2]+1 ]
      }
      res <- res2 <- square.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                                  vsize)
    }
    if (end=="both") {
      res <- cbind(res1, res2)
    }

    res
    
  #####################################################################
  ## plotting mode

  } else if (mode=="plot") {

    .igraph.shape.square(coords=coords, el=el, v=v, mode=mode, params=params,
                         end=end)
    
  }
}

.igraph.shape.rectangle <- function(coords, el=NULL, v=NULL, mode=c("clip", "plot"),
                                    params, end=c("both", "from", "to")) {

  mode=match.arg(mode)
  end =match.arg(end)

  #####################################################################
  ## clipping mode

  if (mode=="clip") {
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
        vertex.size[ el[,1]+1 ]
      }
      vsize2 <- if (length(vertex.size2)==1) {
        vertex.size2
      } else {
        vertex.size2[ el[,1]+1 ]
      }
      res <- res1 <- rec.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                               vsize, vsize2)
    }
    if (end %in% c("to", "both")) {
      vsize <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,2]+1 ]
      }
      vsize2 <- if (length(vertex.size2)==1) {
        vertex.size2
      } else {
        vertex.size2[ el[,1]+1 ]
      }
      res <- res2 <- rec.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                               vsize, vsize2)
    }
    if (end=="both") {
      res <- cbind(res1, res2)
    }
    
    res

  #####################################################################
  ## plotting mode
    
  } else if (mode=="plot") {

    vertex.color       <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v+1]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v+1]
    }
    vertex.size        <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v+1]
    }
    vertex.size <- rep(vertex.size, length=nrow(coords))   
    vertex.size2       <- 1/200 * params("vertex", "size2")
    if (length(vertex.size2) != 1 && !is.null(v)) {
      vertex.size2 <- vertex.size2[v+1]
    }
    vertex.size <- cbind(vertex.size, vertex.size2)

    symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
            rectangles=2*vertex.size, add=TRUE, inches=FALSE)
  }
  
}

.igraph.shape.crectangle <- function(coords, el=NULL, v=NULL, mode=c("clip", "plot"),
                                     params, end=c("both", "from", "to")) {

  mode=match.arg(mode)
  end=match.arg(end)

  #####################################################################
  ## clipping mode

  if (mode=="clip") {
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
        vertex.size[ el[,1]+1 ]
      }
      vsize2 <- if (length(vertex.size2)==1) {
        vertex.size2
      } else {
        vertex.size2[ el[,1]+1 ]
      }
      res <- res1 <- rec.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                               vsize, vsize2)
    }
    if (end %in% c("to", "both")) {
      vsize <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,2]+1 ]
      }
      vsize2 <- if (length(vertex.size2)==1) {
        vertex.size2
      } else {
        vertex.size2[ el[,1]+1 ]
      }
      res <- res2 <- rec.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                               vsize, vsize2)
    }
    if (end=="both") {
      res <- cbind(res1, res2)
    }

    res
    
  #####################################################################
  ## plotting mode

  } else if (mode=="plot") {

    .igraph.shape.rectangle(coords=coords, el=el, v=v, mode=mode, params=params,
                            end=end)
  }
}

.igraph.shape.vrectangle <- function(coords, el=NULL, v=NULL, mode=c("clip", "plot"),
                                     params, end=c("both", "from", "to")) {

  mode=match.arg(mode)
  end=match.arg(end)

  #####################################################################
  ## clipping mode

  if (mode=="clip") {
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
        vertex.size[ el[,1]+1 ]
      }
      vsize2 <- if (length(vertex.size2)==1) {
        vertex.size2
      } else {
        vertex.size2[ el[,1]+1 ]
      }
      res <- res1 <- rec.shift(coords[,3], coords[,4], coords[,1], coords[,2],
                               vsize, vsize2)
    }
    if (end %in% c("to", "both")) {
      vsize <- if (length(vertex.size)==1) {
        vertex.size
      } else {
        vertex.size[ el[,2]+1 ]
      }
      vsize2 <- if (length(vertex.size2)==1) {
        vertex.size2
      } else {
        vertex.size2[ el[,1]+1 ]
      }
      res <- res2 <- rec.shift(coords[,1], coords[,2], coords[,3], coords[,4],
                               vsize, vsize2)
    }
    if (end=="both") {
      res <- cbind(res1, res2)
    }

    res
    
  #####################################################################
  ## plotting mode

  } else if (mode=="plot") {

    .igraph.shape.rectangle(coords=coords, el=el, v=v, mode=mode, params=params,
                            end=end)
  }
}

.igraph.shape.none <- function(coords, el=NULL, v=NULL, mode=c("clip", "plot"),
                               params, end=c("both", "from", "to")) {

  mode=match.arg(mode)

  if (mode=="clip") {
    ## clips like a circle
    .igraph.shape.circle(coords=coords, el=el, v=v, mode=mode, params=params,
                         end=end)
    
  } else if (mode=="plot") {
    ## does not plot anything at all
    invisible(NULL)
  }
}
    
.igraph.shapes <- list("circle" = .igraph.shape.circle,
                       "square" = .igraph.shape.square,
                       "csquare" = .igraph.shape.csquare,
                       "rectangle" = .igraph.shape.rectangle,
                       "crectangle" = .igraph.shape.crectangle,
                       "vrectangle" = .igraph.shape.vrectangle,
                       "none"   = .igraph.shape.none)
