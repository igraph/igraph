
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
                       axes=FALSE, xlab="", ylab="", add=FALSE,
                       xlim=c(-1,1), ylim=c(-1,1), main="", sub="",
                       ...) {

  graph <- x
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  ################################################################
  ## Visual parameters
  params <- i.parse.plot.params(graph, list(...))
  vertex.size        <- 1/200 * params("vertex", "size")
  label.family       <- params("vertex", "label.family")
  label.font         <- params("vertex", "label.font")
  label.cex          <- params("vertex", "label.cex")
  label.degree       <- params("vertex", "label.degree")
  label.color        <- params("vertex", "label.color")
  label.dist         <- params("vertex", "label.dist")
  labels             <- params("vertex", "label")
  shape              <- igraph.check.shapes(params("vertex", "shape"))

  edge.color         <- params("edge", "color")
  edge.width         <- params("edge", "width")
  edge.lty           <- params("edge", "lty")
  arrow.mode         <- params("edge", "arrow.mode")
  edge.labels        <- params("edge", "label")
  loop.angle         <- params("edge", "loop.angle")
  edge.label.font    <- params("edge", "label.font")
  edge.label.family  <- params("edge", "label.family")
  edge.label.cex     <- params("edge", "label.cex")
  edge.label.color   <- params("edge", "label.color")
  arrow.size         <- params("edge", "arrow.size")[1]
  arrow.width        <- params("edge", "arrow.width")[1]
  curved             <- params("edge", "curved")
  
  layout             <- params("plot", "layout")
  margin             <- params("plot", "margin")
  margin <- rep(margin, length=4)
  rescale            <- params("plot", "rescale")
  asp                <- params("plot", "asp")
  frame              <- params("plot", "frame")

  # the new style parameters can't do this yet
  arrow.mode         <- i.get.arrow.mode(graph, arrow.mode)

  ################################################################
  ## create the plot
  maxv <- max(vertex.size)
  if (rescale) {
    # norm layout to (-1, 1)
    layout <- layout.norm(layout, -1, 1, -1, 1)
    xlim <- c(xlim[1]-margin[2]-maxv, xlim[2]+margin[4]+maxv)
    ylim <- c(ylim[1]-margin[1]-maxv, ylim[2]+margin[3]+maxv)
  }
  if (!add) {
    plot(0, 0, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
         axes=axes, frame=frame, asp=asp, main=main, sub=sub)
  }
  
  ################################################################
  ## calculate position of arrow-heads
  el <- get.edgelist(graph, names=FALSE)
  loops.v <- el[,1] [ el[,1] == el[,2] ] + 1
  loops.e <- which(el[,1] == el[,2])
  nonloops.e <- which(el[,1] != el[,2])
  loop.labels <- edge.labels[el[,1] == el[,2]]
  edge.labels <- edge.labels[el[,1] != el[,2]]
  el <- el[el[,1] != el[,2],,drop=FALSE]

  edge.coords <- matrix(0, nrow=nrow(el), ncol=4)
  edge.coords[,1] <- layout[,1][ el[,1]+1 ]
  edge.coords[,2] <- layout[,2][ el[,1]+1 ]
  edge.coords[,3] <- layout[,1][ el[,2]+1 ]
  edge.coords[,4] <- layout[,2][ el[,2]+1 ]
  if ( length(unique(shape)) == 1) {
    ## same vertex shape for all vertices
    ec <- .igraph.shapes[[ shape[1] ]](edge.coords, el, mode="clip",
                                       params=params, end="both")
  } else {
    ## different vertex shapes, do it by "endpoint"
    shape <- rep(shape, length=vcount(graph))
    ec <- edge.coords
    ec[,1:2] <- t(sapply(seq(length=nrow(el)), function(x) {
      .igraph.shapes[[ shape[el[x,1]+1] ]](edge.coords[x,,drop=FALSE],
                                           el[x,,drop=FALSE],
                                           mode="clip", params=params, end="from")
    }))
    ec[,3:4] <- t(sapply(seq(length=nrow(el)), function(x) {
      .igraph.shapes[[ shape[el[x,2]+1] ]](edge.coords[x,,drop=FALSE],
                                           el[x,,drop=FALSE],
                                           mode="clip", params=params, end="to")
    }))
  }
  
  x0 <- ec[,1] ; y0 <- ec[,2] ; x1 <- ec[,3] ; y1 <- ec[,4]

  ################################################################
  ## add the loop edges
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
    
    plot.bezier <- function(cp, points, color, width, arr, lty, arrow.size, arr.w) {
      p <- compute.bezier( cp, points )
      polygon(p[1,], p[2,], border=color, lwd=width, lty=lty)
      if (arr==1 || arr==3) {
        igraph.Arrows(p[1,ncol(p)-1], p[2,ncol(p)-1], p[1,ncol(p)], p[2,ncol(p)],
                      sh.col=color, h.col=color, size=arrow.size,
                      sh.lwd=width, h.lwd=width, open=FALSE, code=2, width=arr.w)
      }
      if (arr==2 || arr==3) {
        igraph.Arrows(p[1,2], p[2,2], p[1,1], p[2,1],
                      sh.col=color, h.col=color, size=arrow.size,
                      sh.lwd=width, h.lwd=width, open=FALSE, code=2, width=arr.w)
      }
    }
    
    loop <- function(x0, y0, cx=x0, cy=y0, color, angle=0, label=NA,
                     width=1, arr=2, lty=1, arrow.size=arrow.size, arr.w=arr.w) {

      rad <- angle
      center <- c(cx,cy)
      cp <- matrix( c(x0,y0, x0+.4,y0+.2, x0+.4,y0-.2, x0,y0),
                   ncol=2, byrow=TRUE)
      phi <- atan2(cp[,2]-center[2], cp[,1]-center[1])
      r <- sqrt((cp[,1]-center[1])**2 + (cp[,2]-center[2])**2)
      
      phi <- phi + rad

      cp[,1] <- cx+r*cos(phi)
      cp[,2] <- cy+r*sin(phi)

      plot.bezier(cp, 50, color, width, arr=arr, lty=lty, arrow.size=arrow.size, arr.w=arr.w)

      if (is.language(label) || !is.na(label)) {
        lx <- x0+.3
        ly <- y0
        phi <- atan2(ly-center[2], lx-center[1])
        r <- sqrt((lx-center[1])**2 + (ly-center[2])**2)

        phi <- phi + rad

        lx <- cx+r*cos(phi)
        ly <- cy+r*sin(phi)

        text(lx, ly, label, col=edge.label.color, font=edge.label.font,
             family=edge.label.family, cex=edge.label.cex)
      }
    }

    ec <- edge.color
    if (length(ec)>1) { ec <- ec[loops.e] }
    vs <- vertex.size
    if (length(vertex.size)>1) { vs <- vs[loops.v] }
    ew <- edge.width
    if (length(edge.width)>1) { ew <- ew[loops.e] }
    la <- loop.angle
    if (length(loop.angle)>1) { la <- la[loops.e] }
    lty <- edge.lty
    if (length(edge.lty)>1) { lty <- lty[loops.e] }
    arr <- arrow.mode
    if (length(arrow.mode)>1) { arr <- arrow.mode[loops.e] }
    asize <- arrow.size
    if (length(arrow.size)>1) { asize <- arrow.size[loops.e] }
    xx0 <- layout[loops.v,1] + cos(la) * vs
    yy0 <- layout[loops.v,2] - sin(la) * vs
    mapply(loop, xx0, yy0,
           color=ec, angle=-la, label=loop.labels, lty=lty,
           width=ew, arr=arr, arrow.size=asize, arr.w=arrow.width)
  }

  ################################################################
  ## non-loop edges
  if (length(x0) != 0) {
    if (length(edge.color)>1) { edge.color <- edge.color[nonloops.e] }
    if (length(edge.width)>1) { edge.width <- edge.width[nonloops.e] }
    if (length(edge.lty)>1) { edge.lty <- edge.lty[nonloops.e] }
    if (length(arrow.mode)>1) { arrow.mode <- arrow.mode[nonloops.e] }
    if (length(arrow.size)>1) { arrow.size <- arrow.size[nonloops.e] }
    if (length(curved)>1) { curved <- curved[nonloops.e] }
    if (length(unique(arrow.mode))==1) {
      igraph.Arrows(x0, y0, x1, y1, h.col=edge.color, sh.col=edge.color,
                    sh.lwd=edge.width, h.lwd=1, open=FALSE, code=arrow.mode[1],
                    sh.lty=edge.lty, h.lty=1, size=arrow.size,
                    width=arrow.width, curved=curved)
    } else {
      ## different kinds of arrows drawn separately as 'arrows' cannot
      ## handle a vector as the 'code' argument
      curved <- rep(curved, length=ecount(graph))[nonloops.e]
      for (code in 0:3) {
        valid <- arrow.mode==code
        if (!any(valid)) { next }
        ec <- edge.color ; if (length(ec)>1) { ec <- ec[valid] }
        ew <- edge.width ; if (length(ew)>1) { ew <- ew[valid] }
        el <- edge.lty   ; if (length(el)>1) { el <- el[valid] }
        igraph.Arrows(x0[valid], y0[valid], x1[valid], y1[valid],
                      code=code, sh.col=ec, h.col=ec, sh.lwd=ew, h.lwd=1,
                      h.lty=1, sh.lty=el, open=FALSE, size=arrow.size,
                      width=arrow.width, curved=curved[valid])
      }
    }
    phi <- atan2(y1-y0, x1-x0)
    r <- sqrt( (x1-x0)^2 + (y1-y0)^2 )
    x <- x0 + 2/3*r*cos(phi)
    y <- y0 + 2/3*r*sin(phi)
    text(x, y, labels=edge.labels, col=edge.label.color, family=edge.label.family,
         font=edge.label.font, cex=edge.label.cex)
  }
  
  rm(x0, y0, x1, y1)
  
  ################################################################
  # add the vertices
  if (length(unique(shape)) == 1) {
    .igraph.shapes[[ shape[1] ]](layout, mode="plot", params=params)
  } else {
    sapply(seq(length=vcount(graph)), function(x) {
      .igraph.shapes[[ shape[x] ]](layout[x,,drop=FALSE], v=x-1,
                                   mode="plot", params=params)
    })
  }
      
  ################################################################
  # add the labels
  par(xpd=TRUE)
  x <- layout[,1]+label.dist*cos(-label.degree)* 
    (vertex.size+6*8*log10(nchar(labels)+1))/200
  y <- layout[,2]+label.dist*sin(-label.degree)*
    (vertex.size+6*8*log10(nchar(labels)+1))/200
  text(x, y, labels=labels, col=label.color, family=label.family, font=label.font,
       cex=label.cex)
  rm(x, y)
  invisible(NULL)
}

rglplot        <- function(x, ...)
  UseMethod("rglplot", x)

rglplot.igraph <- function(x, ...) {

  require(rgl)
  
  graph <- x
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  create.edge <- function(v1, v2, r1, r2, ec, ew, am, as) {
    ## these could also be parameters:
    aw <- 0.005*3*as                      # arrow width
    al <- 0.005*4*as                      # arrow length    
    
    dist <- sqrt(sum((v2-v1)^2))   # distance of the centers

    if (am==0) {
      edge <- qmesh3d(c(-ew/2,-ew/2,dist,1, ew/2,-ew/2,dist,1, ew/2,ew/2,dist,1,
                        -ew/2,ew/2,dist,1,  -ew/2,-ew/2,0,1, ew/2,-ew/2,0,1,
                        ew/2,ew/2,0,1, -ew/2,ew/2,0,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8))
    } else if (am==1) {
      edge <- qmesh3d(c(-ew/2,-ew/2,dist,1, ew/2,-ew/2,dist,1,
                        ew/2,ew/2,dist,1, -ew/2,ew/2,dist,1,
                        -ew/2,-ew/2,al+r1,1, ew/2,-ew/2,al+r1,1,
                        ew/2,ew/2,al+r1,1, -ew/2,ew/2,al+r1,1,
                        -aw/2,-aw/2,al+r1,1, aw/2,-aw/2,al+r1,1,
                        aw/2,aw/2,al+r1,1, -aw/2,aw/2,al+r1,1, 0,0,r1,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8,
                        9,10,11,12, 9,12,13,13, 9,10,13,13, 10,11,13,13,
                        11,12,13,13))
    } else if (am==2) {
      box <- dist-r2-al
      edge <- qmesh3d(c(-ew/2,-ew/2,box,1, ew/2,-ew/2,box,1, ew/2,ew/2,box,1,
                        -ew/2,ew/2,box,1,  -ew/2,-ew/2,0,1, ew/2,-ew/2,0,1,
                        ew/2,ew/2,0,1, -ew/2,ew/2,0,1,
                        -aw/2,-aw/2,box,1, aw/2,-aw/2,box,1, aw/2,aw/2,box,1,
                        -aw/2,aw/2,box,1, 0,0,box+al,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8,
                        9,10,11,12, 9,12,13,13, 9,10,13,13, 10,11,13,13,
                        11,12,13,13))
    } else {
      edge <- qmesh3d(c(-ew/2,-ew/2,dist-al-r2,1, ew/2,-ew/2,dist-al-r2,1,
                        ew/2,ew/2,dist-al-r2,1, -ew/2,ew/2,dist-al-r2,1,
                        -ew/2,-ew/2,r1+al,1, ew/2,-ew/2,r1+al,1,
                        ew/2,ew/2,r1+al,1, -ew/2,ew/2,r1+al,1,
                        -aw/2,-aw/2,dist-al-r2,1, aw/2,-aw/2,dist-al-r2,1,
                        aw/2,aw/2,dist-al-r2,1, -aw/2,aw/2,dist-al-r2,1,
                        -aw/2,-aw/2,r1+al,1, aw/2,-aw/2,r1+al,1,
                        aw/2,aw/2,r1+al,1, -aw/2,aw/2,r1+al,1,
                        0,0,dist-r2,1, 0,0,r1,1),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8,
                        9,10,11,12, 9,12,17,17, 9,10,17,17, 10,11,17,17,
                        11,12,17,17,
                        13,14,15,16, 13,16,18,18, 13,14,18,18, 14,15,18,18,
                        15,16,18,18))
    }
      

    ## rotate and shift it to its position
    phi<- -atan2(v2[2]-v1[2],v1[1]-v2[1])-pi/2
    psi<- acos((v2[3]-v1[3])/dist)    
    rot1 <- rbind(c(1,0,0),c(0,cos(psi),sin(psi)), c(0,-sin(psi),cos(psi)))
    rot2 <- rbind(c(cos(phi),sin(phi),0),c(-sin(phi),cos(phi),0), c(0,0,1))
    rot <- rot1 %*% rot2
    edge <- transform3d(edge, rotationMatrix(matrix=rot))
    edge <- transform3d(edge, translationMatrix(v1[1], v1[2], v1[3]))

    ## we are ready 
    shade3d(edge, col=ec)
  }
  
  create.loop <- function(v, r, ec, ew, am, la, la2, as) {
    aw <- 0.005*3*as
    al <- 0.005*4*as
    wi <- aw*2                          # size of the loop
    wi2 <- wi+aw-ew                     # size including the arrow heads
    hi <- al*2+ew*2
    gap <- wi-2*ew

    if (am==0) {
      edge <- qmesh3d(c(-wi/2,-ew/2,0,1, -gap/2,-ew/2,0,1,
                        -gap/2,ew/2,0,1, -wi/2,ew/2,0,1,
                        -wi/2,-ew/2,hi-ew+r,1, -gap/2,-ew/2,hi-ew+r,1,
                        -gap/2,ew/2,hi-ew+r,1, -wi/2,ew/2,hi-ew+r,1,
                        wi/2,-ew/2,0,1, gap/2,-ew/2,0,1,
                        gap/2,ew/2,0,1, wi/2,ew/2,0,1,
                        wi/2,-ew/2,hi-ew+r,1, gap/2,-ew/2,hi-ew+r,1,
                        gap/2,ew/2,hi-ew+r,1, wi/2,ew/2,hi-ew+r,1,
                        -wi/2,-ew/2,hi+r,1, -wi/2,ew/2,hi+r,1,
                        wi/2,-ew/2,hi+r,1, wi/2,ew/2,hi+r,1
                        ),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7,
                        1,4,18,17,
                        9,10,11,12, 13,14,15,16, 9,10,14,13, 10,11,15,14,
                        11,12,16,15, 9,12,20,19,
                        5,13,19,17, 17,18,20,19, 8,16,20,18, 6,7,15,14
                        ))
    } else if (am==1 || am==2) {
      edge <- qmesh3d(c(-wi/2,-ew/2,r+al,1, -gap/2,-ew/2,r+al,1,
                        -gap/2,ew/2,r+al,1, -wi/2,ew/2,r+al,1,
                        -wi/2,-ew/2,hi-ew+r,1, -gap/2,-ew/2,hi-ew+r,1,
                        -gap/2,ew/2,hi-ew+r,1, -wi/2,ew/2,hi-ew+r,1,
                        wi/2,-ew/2,0,1, gap/2,-ew/2,0,1,
                        gap/2,ew/2,0,1, wi/2,ew/2,0,1,
                        wi/2,-ew/2,hi-ew+r,1, gap/2,-ew/2,hi-ew+r,1,
                        gap/2,ew/2,hi-ew+r,1, wi/2,ew/2,hi-ew+r,1,
                        -wi/2,-ew/2,hi+r,1, -wi/2,ew/2,hi+r,1,
                        wi/2,-ew/2,hi+r,1, wi/2,ew/2,hi+r,1,
                        # the arrow
                        -wi2/2,-aw/2,r+al,1, -wi2/2+aw,-aw/2,r+al,1,
                        -wi2/2+aw,aw/2,r+al,1, -wi2/2,aw/2,r+al,1,
                        -wi2/2+aw/2,0,r,1                   
                        ),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7,
                        1,4,18,17,
                        9,10,11,12, 13,14,15,16, 9,10,14,13, 10,11,15,14,
                        11,12,16,15, 9,12,20,19,
                        5,13,19,17, 17,18,20,19, 8,16,20,18, 6,7,15,14,
                        # the arrow
                        21,22,23,24, 21,22,25,25, 22,23,25,25, 23,24,25,25,
                        21,24,25,25
                        ))
    } else if (am==3) {
      edge <- qmesh3d(c(-wi/2,-ew/2,r+al,1, -gap/2,-ew/2,r+al,1,
                        -gap/2,ew/2,r+al,1, -wi/2,ew/2,r+al,1,
                        -wi/2,-ew/2,hi-ew+r,1, -gap/2,-ew/2,hi-ew+r,1,
                        -gap/2,ew/2,hi-ew+r,1, -wi/2,ew/2,hi-ew+r,1,
                        wi/2,-ew/2,r+al,1, gap/2,-ew/2,r+al,1,
                        gap/2,ew/2,r+al,1, wi/2,ew/2,r+al,1,
                        wi/2,-ew/2,hi-ew+r,1, gap/2,-ew/2,hi-ew+r,1,
                        gap/2,ew/2,hi-ew+r,1, wi/2,ew/2,hi-ew+r,1,
                        -wi/2,-ew/2,hi+r,1, -wi/2,ew/2,hi+r,1,
                        wi/2,-ew/2,hi+r,1, wi/2,ew/2,hi+r,1,
                        # the arrows
                        -wi2/2,-aw/2,r+al,1, -wi2/2+aw,-aw/2,r+al,1,
                        -wi2/2+aw,aw/2,r+al,1, -wi2/2,aw/2,r+al,1,
                        -wi2/2+aw/2,0,r,1,
                        wi2/2,-aw/2,r+al,1, wi2/2-aw,-aw/2,r+al,1,
                        wi2/2-aw,aw/2,r+al,1, wi2/2,aw/2,r+al,1,
                        wi2/2-aw/2,0,r,1                   
                        ),
                      c(1,2,3,4, 5,6,7,8, 1,2,6,5, 2,3,7,6, 3,4,8,7,
                        1,4,18,17,
                        9,10,11,12, 13,14,15,16, 9,10,14,13, 10,11,15,14,
                        11,12,16,15, 9,12,20,19,
                        5,13,19,17, 17,18,20,19, 8,16,20,18, 6,7,15,14,
                        # the arrows
                        21,22,23,24, 21,22,25,25, 22,23,25,25, 23,24,25,25,
                        21,24,25,25,
                        26,27,28,29, 26,27,30,30, 27,28,30,30, 28,29,30,30,
                        26,29,30,30
                        ))
    }

    # rotate and shift to its position
    rot1 <- rbind(c(1,0,0),c(0,cos(la2),sin(la2)), c(0,-sin(la2),cos(la2)))
    rot2 <- rbind(c(cos(la),sin(la),0),c(-sin(la),cos(la),0), c(0,0,1))
    rot <- rot1 %*% rot2
    edge <- transform3d(edge, rotationMatrix(matrix=rot))
    edge <- transform3d(edge, translationMatrix(v[1], v[2], v[3]))

    ## we are ready
    shade3d(edge, col=ec)
  }
  
  # Visual parameters
  params <- i.parse.plot.params(graph, list(...))
  labels <- params("vertex", "label")
  label.color <- params("vertex", "label.color")
  label.font <- params("vertex", "label.font")
  label.degree <- params("vertex", "label.degree")
  label.dist <- params("vertex", "label.dist")
  vertex.color <- params("vertex", "color")
  vertex.size <- (1/200) * params("vertex", "size")
  loop.angle <- params("edge", "loop.angle")
  loop.angle2 <- params("edge", "loop.angle2")

  edge.color <- params("edge", "color")
  edge.width <- (1/200) * params("edge", "width")
  edge.labels <- params("edge","label")
  arrow.mode <- params("edge","arrow.mode")
  arrow.size <- params("edge","arrow.size")
  
  layout <- params("plot", "layout")
  rescale <- params("plot", "rescale")

  # the new style parameters can't do this yet
  arrow.mode         <- i.get.arrow.mode(graph, arrow.mode)
  
  # norm layout to (-1, 1)
  if (ncol(layout)==2) { layout <- cbind(layout, 0) }
  if (rescale) {
    layout <- layout.norm(layout, -1, 1, -1, 1, -1, 1)
  }
  
  # add the edges, the loops are handled separately
  el <- get.edgelist(graph, names=FALSE)
  
  # It is faster this way
  par3d(skipRedraw=TRUE)

  # edges first
  for (i in seq(length=nrow(el))) {
    from <- el[i,1]
    to   <- el[i,2]
    v1 <- layout[from+1,]
    v2 <- layout[to+1,]
    am <- arrow.mode; if (length(am)>1) { am <- am[i] }
    ew <- edge.width; if (length(ew)>1) { ew <- ew[i] }
    ec <- edge.color; if (length(ec)>1) { ec <- ec[i] }
    r1 <- vertex.size; if (length(r1)>1) { r1 <- r1[from+1] }
    r2 <- vertex.size; if (length(r2)>1) { r2 <- r2[to+1] }

    if (from!=to) {
      create.edge(v1,v2,r1,r2,ec,ew,am,arrow.size)
    } else {
      la <- loop.angle; if (length(la)>1) { la <- la[i] }
      la2 <- loop.angle2; if (length(la2)>1) { la2 <- la2[i] }      
      create.loop(v1,r1,ec,ew,am,la,la2,arrow.size)
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
  if (any(edge.labels != "")) {
    x0 <- layout[,1][el[,1]+1]
    x1 <- layout[,1][el[,2]+1]
    y0 <- layout[,2][el[,1]+1]
    y1 <- layout[,2][el[,2]+1]
    z0 <- layout[,3][el[,1]+1]
    z1 <- layout[,4][el[,2]+1]
    rgl.texts((x0+x1)/2, (y0+y1)/2, (z0+z1)/2, edge.labels,
              col=label.color)
  }

  # draw everything
  par3d(skipRedraw=FALSE)
  
  invisible(NULL)
}

# This is taken from the IDPmisc package,
# slightly modified: code argument added

igraph.Arrows <-
function (x1, y1, x2, y2,
                    code=2,
                    size= 1,     
                    width= 1.2/4/cin,
                    open=TRUE,
                    sh.adj=0.1, 
                    sh.lwd=1,
                    sh.col=if(is.R()) par("fg") else 1,
                    sh.lty=1,
                    h.col=sh.col,
                    h.col.bo=sh.col,
                    h.lwd=sh.lwd,
                    h.lty=sh.lty,
                    curved=FALSE)
  ## Author: Andreas Ruckstuhl, refined by Rene Locher
  ## Version: 2005-10-17
{
  cin <- size * par("cin")[2]
  width <- width * (1.2/4/cin)
  uin <- if (is.R()) 
    1/xyinch()
  else par("uin")
  x <- sqrt(seq(0, cin^2, length = floor(35 * cin) + 2))
  delta <-  sqrt(h.lwd)*par("cin")[2]*0.005      ## has been 0.05
  x.arr <- c(-rev(x), -x)
  wx2 <- width * x^2
  y.arr <- c(-rev(wx2 + delta), wx2 + delta)
  deg.arr <- c(atan2(y.arr, x.arr), NA)
  r.arr <- c(sqrt(x.arr^2 + y.arr^2), NA)

  ## backup
  bx1 <- x1 ; bx2 <- x2 ; by1 <- y1 ; by2 <- y2
  
  ## shaft
  lx <- length(x1)
  theta <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
  r.seg <- rep(cin*sh.adj, lx)
  th.seg <- theta + rep(atan2(0, -cin), lx)
  if (is.logical(curved) && all(!curved)) {
    segments(x1, y1, x2+r.seg*cos(th.seg)/uin[1], y2+r.seg*sin(th.seg)/uin[2], 
             lwd=sh.lwd, col=sh.col, lty=sh.lty)
  } else {
    if (is.numeric(curved)) {
      lambda <- curved
    } else {
      lambda <- as.logical(curved) * 0.5
    }
    c.x1 <- x1
    c.y1 <- y1
    c.x2 <- x2+r.seg*cos(th.seg)/uin[1]
    c.y2 <- y2+r.seg*sin(th.seg)/uin[2]

    midx <- (x1+x2)/2
    midy <- (y1+y2)/2  
    spx <- midx - lambda * 1/2 * (c.y2-c.y1)
    spy <- midy + lambda * 1/2 * (c.x2-c.x1)
    sh.col <- rep(sh.col, length=length(c.x1))
    sh.lty <- rep(sh.lty, length=length(c.x1))
    sh.lwd <- rep(sh.lwd, length=length(c.x1))
    for (i in seq_len(length(c.x1))) {
      spl <- xspline(x=c(c.x1[i],spx[i],c.x2[i]),
                     y=c(c.y1[i],spy[i],c.y2[i]), shape=1, draw=FALSE)
      lines(spl, lwd=sh.lwd[i], col=sh.col[i], lty=sh.lty[i])
      if (code %in% c(2,3)) {
        x1[i] <- spl$x[3*length(spl$x)/4]
        y1[i] <- spl$y[3*length(spl$y)/4]
      }
      if (code %in% c(1,3)) {
        x2[i] <- spl$x[length(spl$x)/4]
        y2[i] <- spl$y[length(spl$y)/4]
      }
    }
  }

  ## forward arrowhead
  if (code %in% c(2,3)) {    
    theta <- atan2((by2 - y1) * uin[2], (bx2 - x1) * uin[1])
    Rep <- rep(length(deg.arr), lx)
    p.x2 <- rep(bx2, Rep)
    p.y2 <- rep(by2, Rep)
    ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
    r.arr <- rep(r.arr, lx)  
    if(open) lines((p.x2 + r.arr * cos(ttheta)/uin[1]),
                   (p.y2 + r.arr*sin(ttheta)/uin[2]), 
                   lwd=h.lwd, col = h.col.bo, lty=h.lty) else
    polygon(p.x2 + r.arr * cos(ttheta)/uin[1], p.y2 + r.arr*sin(ttheta)/uin[2], 
            col = h.col, lwd=h.lwd,
            border=h.col.bo, lty=h.lty)
  }
    
  ## backward arrow head
  if (code %in% c(1,3)) {
    x1 <- bx1; y1 <- by1
    tmp <- x1 ; x1 <- x2 ; x2 <- tmp
    tmp <- y1 ; y1 <- y2 ; y2 <- tmp
    theta <- atan2((y2 - y1) * uin[2], (x2 - x1) * uin[1])
    lx <- length(x1)
    Rep <- rep(length(deg.arr), lx)
    p.x2 <- rep(x2, Rep)
    p.y2 <- rep(y2, Rep)
    ttheta <- rep(theta, Rep) + rep(deg.arr, lx)
    r.arr <- rep(r.arr, lx)

    if(open) lines((p.x2 + r.arr * cos(ttheta)/uin[1]),
                   (p.y2 + r.arr*sin(ttheta)/uin[2]), 
                   lwd=h.lwd, col = h.col.bo, lty=h.lty) else
    polygon(p.x2 + r.arr * cos(ttheta)/uin[1], p.y2 + r.arr*sin(ttheta)/uin[2], 
            col = h.col, lwd=h.lwd,
            border=h.col.bo, lty=h.lty)
  }
} # Arrows

