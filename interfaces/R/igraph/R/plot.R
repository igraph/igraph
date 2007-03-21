
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

plot.igraph <- function(x, 
                       # SPECIFIC: #####################################
                       axes=FALSE, xlab="", ylab="",
                       xlim=c(-1,1), ylim=c(-1,1),
                       ...) {

  graph <- x
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  # Visual parameters
  params <- i.parse.plot.params(graph, list(...))
  vertex.color       <- params("vertex", "color")
  vertex.frame.color <- params("vertex", "frame.color")
  vertex.size        <- (1/200) * params("vertex", "size")
  label.font         <- params("vertex", "label.font")
  label.degree       <- params("vertex", "label.degree")
  label.color        <- params("vertex", "label.color")
  label.dist         <- params("vertex", "label.dist")
  labels             <- params("vertex", "label")

  edge.color         <- params("edge", "color")
  edge.width         <- params("edge", "width")
  edge.lty           <- params("edge", "lty")
  arrow.mode         <- params("edge", "arrow.mode")
  edge.labels        <- params("edge", "label")
  loop.angle         <- params("edge", "loop.angle")
  
  layout             <- params("plot", "layout")
  margin             <- params("plot", "margin")

  # the new style parameters can't do this yet
  arrow.mode         <- i.get.arrow.mode(graph, arrow.mode)

  # create the plot
  maxv <- max(vertex.size)
  xlim <- c(xlim[1]-margin-maxv, xlim[2]+margin+maxv)
  ylim <- c(ylim[1]-margin-maxv, ylim[2]+margin+maxv)
  plot(0, 0, type="n", xlab=xlab, ylab=ylab, asp=1, xlim=xlim, ylim=ylim,
       axes=axes)

  # norm layout to (-1, 1)
  layout <- i.layout.norm(layout, -1, 1, -1, 1)
  
  # add the edges
  el <- get.edgelist(graph, names=FALSE)
  loops.v <- el[,1] [ el[,1] == el[,2] ] + 1
  loops.e <- which(el[,1] == el[,2])
  nonloops.e <- which(el[,1] != el[,2])
  loop.labels <- edge.labels[el[,1] == el[,2]]
  edge.labels <- edge.labels[el[,1] != el[,2]]
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

    point.on.cubic.bezier <- function(cp, t) {

      c <- 3 * (cp[2,] - cp[1,])
      b <- 3 * (cp[3,] - cp[2,]) - c
      a <- cp[4,] - cp[1,] - c - b
      
      t2 <- t*t;
      t3 <- t*t*t
      
      a*t3 + b*t2 + c*t + cp[1,]
    }

    compute.bezier <- function(cp, points) {
      dt <- seq(0, 1, by=1/(points-1))
      sapply(dt, function(t) point.on.cubic.bezier(cp, t))
    }
    
    plot.bezier <- function(cp, points, color, width, arr, lty) {
      p <- compute.bezier( cp, points )
      polygon(p[1,], p[2,], border=color, lwd=width, lty=lty)
      if (arr != 0) {
        arrows(p[1,ncol(p)-1], p[2,ncol(p)-1], p[1,ncol(p)], p[2,ncol(p)],
               length=.2, angle=20, col=color, lwd=width, code=arr)
      }
    }
    
    loop <- function(x0, y0, cx=x0, cy=y0, color, angle=0, label=NA,
                     width=1, arr=2, lty=1) {
      rad <- angle/180*pi
      center <- c(cx,cy)
      cp <- matrix( c(x0,y0, x0+.4,y0+.2, x0+.4,y0-.2, x0,y0),
                   nc=2, byrow=TRUE)
      phi <- atan2(cp[,2]-center[2], cp[,1]-center[1])
      r <- sqrt((cp[,1]-center[1])**2 + (cp[,2]-center[2])**2)
      
      phi <- phi + rad

      cp[,1] <- cx+r*cos(phi)
      cp[,2] <- cy+r*sin(phi)

      plot.bezier(cp, 50, color, width, arr=arr, lty=lty)

      if (!is.na(label)) {
        lx <- x0+.3
        ly <- y0
        phi <- atan2(ly-center[2], lx-center[1])
        r <- sqrt((lx-center[1])**2 + (ly-center[2])**2)

        phi <- phi + rad

        lx <- cx+r*cos(phi)
        ly <- cy+r*sin(phi)

        text(lx, ly, label, col=label.color)
      }
    }

    ec <- edge.color
    if (length(ec)>1) { ec <- ec[loops.e] }
    vs <- vertex.size
    if (length(vertex.size)>1) { vs <- vs[loops.e] }
    ew <- edge.width
    if (length(edge.width)>1) { ew <- ew[loops.e] }
    la <- loop.angle
    if (length(loop.angle)>1) { la <- la[loops.e] }
    lty <- edge.lty
    if (length(edge.lty)>1) { lty <- lty[loops.e] }
    arr <- arrow.mode
    if (length(arrow.mode)>1) { arrow.mode <- arrow.mode[loops.e] }
    xx0 <- layout[loops.v,1] + cos(la/180*pi) * vs
    yy0 <- layout[loops.v,2] + sin(la/180*pi) * vs
    mapply(loop, xx0, yy0,
           color=ec, angle=la, label=loop.labels, lty=lty,
           width=ew, arrows=arr)
    
  }
  
  if (length(x0) != 0) {
    if (length(edge.color)>1) { edge.color <- edge.color[nonloops.e] }
    if (length(edge.width)>1) { edge.width <- edge.width[nonloops.e] }
    if (length(edge.lty)>1) { edge.lty <- edge.lty[nonloops.e] }
    if (length(arrow.mode)>1) { arrow.mode <- arrow.mode[nonloops.e] }
    if (length(unique(arrow.mode))==1) {
      arrows(x0, y0, x1, y1, angle=20, length=0.2, code=arrow.mode,
             col=edge.color, lwd=edge.width, lty=edge.lty)
      if (any(edge.lty != 1)) {
        if (arrow.mode==2 || arrow.mode==3) {
          pp <- atan2(y0-y1, x0-x1)
          xx <- x1+0.001*cos(pp)
          yy <- y1+0.001*sin(pp)
          arrows(xx, yy, x1, y1, angle=20, length=0.2, code=2,
                 col=edge.color, lwd=edge.width, lty=1)
        }
        if (arrow.mode==1 || arrow.mode==3) {          
          pp <- atan2(y1-y0, x1-x0)
          xx <- x0+0.001*cos(pp)
          yy <- y0+0.001*sin(pp)
          arrows(xx, yy, x0, y0, angle=20, length=0.2, code=2,
                 col=edge.color, lwd=edge.width, lty=1)
        }
      }
    } else {
      ## different kinds of arrow drawn separately as 'arrows' cannot
      ## handle a vector as the 'code' argument
      for (code in 0:3) {
        valid <- arrow.mode==code
        if (!any(valid)) { next }
        ec <- edge.color ; if (length(ec)>1) { ec <- ec[valid] }
        ew <- edge.width ; if (length(ew)>1) { ew <- ew[valid] }
        el <- edge.lty   ; if (length(el)>1) { el <- el[valid] }        
        arrows(x0[valid], y0[valid], x1[valid], y1[valid],
               angle=20, length=0.2, code=code,
               col=ec, lwd=ew, lty=el)
        if (any(el != 1) && code %in% c(2,3)) {
          pp <- atan2(y0-y1, x0-x1)
          xx <- x1+0.01*cos(pp)
          yy <- y1+0.01*sin(pp)
          arrows(xx[valid], yy[valid], x1[valid], y1[valid],
                 angle=20, length=0.2, code=2,
                 col=ec, lwd=ew, lty=1)
        }
        if (any(el != 1) && code %in% c(1,3)) {
          pp <- atan2(y1-y0, x1-x0)
          xx <- x0+0.01*cos(pp)
          yy <- y0+0.01*sin(pp)
          arrows(xx[valid], yy[valid], x0[valid], y0[valid],
                 angle=20, length=0.2, code=2,
                 col=ec, lwd=ew, lty=1)
        }
      }
    }
    x <- (x0+x1)/2
    y <- (y0+y1)/2
    if (!is.null(label.font)) par(family=label.font)
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
                           arrow.mode=NULL,
                           ...)
  UseMethod("rglplot", x)


rglplot.igraph <- function(x, layout=layout.random, layout.par=list(),
                           labels=NULL, label.color="darkblue",
                           label.font=NULL, label.degree=-pi/4, label.dist=0,
                           vertex.color="SkyBlue2", vertex.size=15,
                           edge.color="darkgrey", edge.width=1,
                           edge.labels=NA, 
                           # SPECIFIC: #####################################
                           arrow.mode=NULL,
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
  arrow.mode <- i.get.arrow.mode(graph, arrow.mode)

  # norm layout to (-1, 1)
  if (ncol(layout)==2) { layout <- cbind(layout, 0) }
  layout <- i.layout.norm(layout, -1, 1, -1, 1, -1, 1)
  
  # add the edges
  # TODO: loops
  el <- get.edgelist(graph, names=FALSE)
  x0 <- layout[,1][el[,1]+1]
  y0 <- layout[,2][el[,1]+1]
  z0 <- layout[,3][el[,1]+1]
  x1 <- layout[,1][el[,2]+1]
  y1 <- layout[,2][el[,2]+1]
  z1 <- layout[,3][el[,2]+1]
  
  if (length(vertex.size)!=1) {
    vsize.from <- vertex.size[get.edgelist(graph)[,1]+1]
    vsize.to   <- vertex.size[get.edgelist(graph)[,2]+1]
  } else {
    vsize.from <- vsize.to <- vertex.size
  }

  # It is faster this way
  par3d(skipRedraw=TRUE)
  
  # edges
  indices <- c(1,2,3,4,
               5,6,7,8,
               1,2,6,5,
               2,3,7,6,
               3,4,8,7,
               4,1,5,8)
  
  indices2 <- c(1,2,3,3,
                1,3,4,4,
                1,4,5,5,
                1,5,2,2,
                2,3,4,5)

  for (i in seq(x0)) {
    v1 <- c(x0[i], y0[i], z0[i])
    v2 <- c(x1[i], y1[i], z1[i])
    r <- vsize.to; if (length(r)>1) { r <- r[i] }
    alen <- sqrt(sum((v2-v1)^2))

    vertices <- c(-1,-1,alen,1,
                   1,-1,alen,1,
                   1, 1,alen,1,
                  -1, 1,alen,1,
                  -1,-1,0,   1,
                   1,-1,0,   1,
                   1, 1,0,   1,
                  -1, 1,0,   1)

    edge.prot <- qmesh3d(vertices, indices)
    ew <- edge.width; if (length(ew)>1) { ew <- ew[i] }
    edge <- transform3d(edge.prot, scaleMatrix(.01*ew, .01*ew, 1))
    phi<--atan2(v2[2]-v1[2],v1[1]-v2[1])-pi/2
    psi<-acos((v2[3]-v1[3])/alen)
    rot1 <- rbind(c(1,0,0),c(0,cos(psi),sin(psi)), c(0,-sin(psi),cos(psi)))
    rot2 <- rbind(c(cos(phi),sin(phi),0),c(-sin(phi),cos(phi),0), c(0,0,1))
    edge <- transform3d(edge, rotationMatrix(matrix=rot1))
    edge <- transform3d(edge, rotationMatrix(matrix=rot2))
    edge <- transform3d(edge, translationMatrix(v1[1], v1[2], v1[3]))
    ec <- edge.color; if (length(ec)>1) { ec <- ec[i] }
    shade3d(edge, col=ec)

    arr <- arrow.mode; if (length(arrow.mode)>1) { arr <- arr[i] }
    if (arr == 2 || arr==3) {
      ## forward
      vertices <- c( 0,  0,  0, 1,
                    -1, -1, -2, 1,
                     1, -1, -2, 1,
                     1,  1, -2, 1,
                    -1,  1, -2, 1)
      
      arrow <- qmesh3d(vertices, indices2)
      arrow <- transform3d(arrow, scaleMatrix(0.03*ew, 0.03*ew, .05*ew))
      arrow <- transform3d(arrow, translationMatrix(0,0,-r))
      arrow <- transform3d(arrow, rotationMatrix(matrix=rot1))
      arrow <- transform3d(arrow, rotationMatrix(matrix=rot2))
      arrow <- transform3d(arrow, translationMatrix(v2[1], v2[2], v2[3]))
      shade3d(arrow, col=ec)    
    }
    if (arr == 1 || arr==3) {
      ## backward
      vertices <- c( 0,  0,  0, 1,
                    -1, -1,  2, 1,
                     1, -1,  2, 1,
                     1,  1,  2, 1,
                    -1,  1,  2, 1)
      
      arrow <- qmesh3d(vertices, indices2)
      arrow <- transform3d(arrow, scaleMatrix(0.03*ew, 0.03*ew, .05*ew))
      arrow <- transform3d(arrow, translationMatrix(0,0,r))
      arrow <- transform3d(arrow, rotationMatrix(matrix=rot1))
      arrow <- transform3d(arrow, rotationMatrix(matrix=rot2))
      arrow <- transform3d(arrow, translationMatrix(v1[1], v1[2], v1[3]))
      shade3d(arrow, col=ec)      
    }
  }
    
  # add the vertices
  if (length(vertex.size)==1) { vertex.size <- rep(vertex.size, nrow(layout)) }
  rgl.spheres(layout[,1], layout[,2], layout[,3], radius=vertex.size,
              col=vertex.color)

  # add the labels, 'l1' is a stupid workaround of a mysterious rgl bug
  labels[is.na(labels)] <- ""
  x <- layout[,1]+label.dist*cos(-label.degree)* 
    (vertex.size+6*10*log10(nchar(labels)+1))/200
  y <- layout[,2]+label.dist*sin(-label.degree)*
    (vertex.size+6*10*log10(nchar(labels)+1))/200
  z <- layout[,3]
  l1 <- labels[1]
  labels[1] <- ""
  rgl.texts(x,y,z, labels, col=label.color, adj=0)
  rgl.texts(c(0,x[1]), c(0,y[1]), c(0,z[1]),
            c("",l1), col=c(label.color[1],label.color[1]), adj=0)

  edge.labels[is.na(edge.labels)] <- ""
  rgl.texts((x0+x1)/2, (y0+y1)/2, (z0+z1)/2, edge.labels,
            col=label.color)

  # draw everything
  par3d(skipRedraw=FALSE)
  
  invisible(NULL)
}

