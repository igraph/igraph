
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
## shapes()         - lists all vertex shapes
## shapes(shape)    - returns the clipping and plotting functions
##                           for a given vertex shape
## add_shape()             - adds a new vertex shape, the clipping and
##                           plotting functions must be given, and
##                           optionally the newly introduced plotting
##                           parameters. This function can also be used
##                           to overwrite a given vertex shape.
##
## Examples:
## add_shape("image", clip=image.clip, plot=image.plot,
##                   parameters=list(filename=NA))
##
## add_shape("triangle", clip=shapes("circle")$clip,
##                   plot=triangle.plot)
##
## add_shape("polygon", clip=shapes("circle")$clip,
##                   plot=polygon.plot)
##
###################################################################

#' Various vertex shapes when plotting igraph graphs
#'
#' Starting from version 0.5.1 igraph supports different
#' vertex shapes when plotting graphs.
#'
#' @details
#' In igraph a vertex shape is defined by two functions: 1) provides
#' information about the size of the shape for clipping the edges and 2)
#' plots the shape if requested. These functions are called \dQuote{shape
#'   functions} in the rest of this manual page. The first one is the
#' clipping function and the second is the plotting function.
#'
#' The clipping function has the following arguments:
#' \describe{
#'   \item{coords}{A matrix with four columns, it contains the
#'     coordinates of the vertices for the edge list supplied in the
#'     \code{el} argument.}
#'   \item{el}{A matrix with two columns, the edges of which some end
#'     points will be clipped. It should have the same number of rows as
#'     \code{coords}.}
#'   \item{params}{This is a function object that can be called to query
#'     vertex/edge/plot graphical parameters. The first argument of the
#'     function is \dQuote{\code{vertex}}, \dQuote{\code{edge}} or
#'     \dQuote{\code{plot}} to decide the type of the parameter, the
#'     second is a character string giving the name of the
#'     parameter. E.g.
#'     \preformatted{
#'	params("vertex", "size")
#'     }
#'   }
#'   \item{end}{Character string, it gives which end points will be
#'     used. Possible values are \dQuote{\code{both}},
#'     \dQuote{\code{from}} and \dQuote{\code{to}}. If
#'     \dQuote{\code{from}} the function is expected to clip the
#'     first column in the \code{el} edge list, \dQuote{\code{to}}
#'     selects the second column, \dQuote{\code{both}} selects both.}
#' }
#'
#' The clipping function should return a matrix
#' with the same number of rows as the \code{el} arguments.
#' If \code{end} is \code{both} then the matrix must have four
#' columns, otherwise two. The matrix contains the modified coordinates,
#' with the clipping applied.
#'
#' The plotting function has the following arguments:
#' \describe{
#'   \item{coords}{The coordinates of the vertices, a matrix with two
#'     columns.}
#'   \item{v}{The ids of the vertices to plot. It should match the number
#'     of rows in the \code{coords} argument.}
#'   \item{params}{The same as for the clipping function, see above.}
#' }
#'
#' The return value of the plotting function is not used.
#'
#' \code{shapes} can be used to list the names of all installed
#' vertex shapes, by calling it without arguments, or setting the
#' \code{shape} argument to \code{NULL}. If a shape name is given, then
#' the clipping and plotting functions of that shape are returned in a
#' named list.
#'
#' \code{add_shape} can be used to add new vertex shapes to
#' igraph. For this one must give the clipping and plotting functions of
#' the new shape. It is also possible to list the plot/vertex/edge
#' parameters, in the \code{parameters} argument, that the clipping
#' and/or plotting functions can make use of. An example would be a
#' generic regular polygon shape, which can have a parameter for the
#' number of sides.
#'
#' \code{shape_noclip} is a very simple clipping function that the
#' user can use in their own shape definitions. It does no clipping, the
#' edges will be drawn exactly until the listed vertex position
#' coordinates.
#'
#' \code{shape_noplot} is a very simple (and probably not very
#' useful) plotting function, that does not plot anything.
#'
#' @aliases add.vertex.shape igraph.shape.noclip igraph.shape.noplot
#'   vertex.shapes
#'
#' @param shape Character scalar, name of a vertex shape. If it is
#'    \code{NULL} for \code{shapes}, then the names of all defined
#'    vertex shapes are returned.
#' @param clip An R function object, the clipping function.
#' @param plot An R function object, the plotting function.
#' @param parameters Named list, additional plot/vertex/edge
#'    parameters. The element named define the new parameters, and the
#'    elements themselves define their default values.
#'    Vertex parameters should have a prefix
#'    \sQuote{\code{vertex.}}, edge parameters a prefix
#'    \sQuote{\code{edge.}}. Other general plotting parameters should have
#'    a prefix \sQuote{\code{plot.}}. See Details below.
#' @param coords,el,params,end,v See parameters of the clipping/plotting
#'    functions below.
#' @return \code{shapes} returns a character vector if the
#'    \code{shape} argument is \code{NULL}. It returns a named list with
#'    entries named \sQuote{clip} and \sQuote{plot}, both of them R
#'    functions.
#'
#'    \code{add_shape} returns \code{TRUE}, invisibly.
#'
#'    \code{shape_noclip} returns the appropriate columns of its
#'    \code{coords} argument.
#' @export
#'
#' @examples
#' # all vertex shapes, minus "raster", that might not be available
#' shapes <- setdiff(shapes(), "")
#' g <- make_ring(length(shapes))
#' set.seed(42)
#' plot(g, vertex.shape=shapes, vertex.label=shapes, vertex.label.dist=1,
#'      vertex.size=15, vertex.size2=15,
#'      vertex.pie=lapply(shapes, function(x) if (x=="pie") 2:6 else 0),
#'      vertex.pie.color=list(heat.colors(5)))
#'
#' # add new vertex shape, plot nothing with no clipping
#' add_shape("nil")
#' plot(g, vertex.shape="nil")
#'
#' #################################################################
#' # triangle vertex shape
#' mytriangle <- function(coords, v=NULL, params) {
#'   vertex.color <- params("vertex", "color")
#'   if (length(vertex.color) != 1 && !is.null(v)) {
#'     vertex.color <- vertex.color[v]
#'   }
#'   vertex.size <- 1/200 * params("vertex", "size")
#'   if (length(vertex.size) != 1 && !is.null(v)) {
#'     vertex.size <- vertex.size[v]
#'   }
#'
#'   symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
#'           stars=cbind(vertex.size, vertex.size, vertex.size),
#'           add=TRUE, inches=FALSE)
#' }
#' # clips as a circle
#' add_shape("triangle", clip=shapes("circle")$clip,
#'                  plot=mytriangle)
#' plot(g, vertex.shape="triangle", vertex.color=rainbow(vcount(g)),
#'      vertex.size=seq(10,20,length=vcount(g)))
#'
#' #################################################################
#' # generic star vertex shape, with a parameter for number of rays
#' mystar <- function(coords, v=NULL, params) {
#'   vertex.color <- params("vertex", "color")
#'   if (length(vertex.color) != 1 && !is.null(v)) {
#'     vertex.color <- vertex.color[v]
#'   }
#'   vertex.size  <- 1/200 * params("vertex", "size")
#'   if (length(vertex.size) != 1 && !is.null(v)) {
#'     vertex.size <- vertex.size[v]
#'   }
#'   norays <- params("vertex", "norays")
#'   if (length(norays) != 1 && !is.null(v)) {
#'     norays <- norays[v]
#'   }
#'
#'   mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
#'          FUN=function(x, y, bg, size, nor) {
#'            symbols(x=x, y=y, bg=bg,
#'                    stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
#'                    add=TRUE, inches=FALSE)
#'          })
#' }
#' # no clipping, edges will be below the vertices anyway
#' add_shape("star", clip=shape_noclip,
#'                  plot=mystar, parameters=list(vertex.norays=5))
#' plot(g, vertex.shape="star", vertex.color=rainbow(vcount(g)),
#'      vertex.size=seq(10,20,length=vcount(g)))
#' plot(g, vertex.shape="star", vertex.color=rainbow(vcount(g)),
#'      vertex.size=seq(10,20,length=vcount(g)),
#'      vertex.norays=rep(4:8, length=vcount(g)))
#'
#' #################################################################
#' # Pictures as vertices.
#' # Similar musicians from last.fm, we start from an artist and
#' # will query two levels. We will use the XML, png and jpeg packages
#' # for this, so these must be available. Otherwise the example is
#' # skipped
#'
#' loadIfYouCan <- function(pkg) suppressWarnings(do.call(require, list(pkg)))
#'
#' if (loadIfYouCan("XML") && loadIfYouCan("png") &&
#'     loadIfYouCan("jpeg")) {
#'   url <- paste(sep="",
#'                'http://ws.audioscrobbler.com/',
#'                '2.0/?method=artist.getinfo&artist=%s',
#'                '&api_key=1784468ada3f544faf9172ee8b99fca3')
#'   getartist <- function(artist) {
#'     cat("Downloading from last.fm. ... ")
#'     txt <- readLines(sprintf(url, URLencode(artist)))
#'     xml <- xmlTreeParse(txt, useInternal=TRUE)
#'     img <- xpathSApply(xml, "/lfm/artist/image[@@size='medium'][1]",
#'                        xmlValue)
#'     if (img != "") {
#'       con <- url(img, open="rb")
#'       bin <- readBin(con, what="raw", n=10^6)
#'       close(con)
#'       if (grepl("\\\\.png$", img)) {
#'         rast <- readPNG(bin, native=TRUE)
#'       } else if (grepl("\\\\.jpe?g$", img)) {
#'         rast <- readJPEG(bin, native=TRUE)
#'       } else {
#'         rast <- as.raster(matrix())
#'       }
#'     } else {
#'       rast <- as.raster(numeric())
#'     }
#'     sim <- xpathSApply(xml, "/lfm/artist/similar/artist/name", xmlValue)
#'     cat("done.\\n")
#'     list(name=artist, image=rast, similar=sim)
#'   }
#'
#'   ego <- getartist("Placebo")
#'   similar <- lapply(ego$similar, getartist)
#'
#'   edges1 <- cbind(ego$name, ego$similar)
#'   edges2 <- lapply(similar, function(x) cbind(x$name, x$similar))
#'   edges3 <- rbind(edges1, do.call(rbind, edges2))
#'   edges <- edges3[ edges3[,1] %in% c(ego$name, ego$similar) &
#'                    edges3[,2] %in% c(ego$name, ego$similar), ]
#'
#'   musnet <- simplify(graph_from_data_frame(edges, dir=FALSE,
#'                      vertices=data.frame(name=c(ego$name, ego$similar))))
#'   str(musnet)
#'
#'   V(musnet)$raster <- c(list(ego$image), lapply(similar, "[[", "image"))
#'   plot(musnet, layout=layout_as_star, vertex.shape="raster",
#'        vertex.label=V(musnet)$name, margin=.2,
#'        vertex.size=50, vertex.size2=50,
#'        vertex.label.dist=2, vertex.label.degree=0)
#' } else {
#'   message("You need the `XML', `png' and `jpeg' packages to run this")
#' }

shapes <- function(shape=NULL) {
  if (is.null(shape)) {
    ls(.igraph.shapes)
  } else {
    ## checkScalarString(shape)
    .igraph.shapes[[shape]]
  }
}

#' @rdname shapes
#' @export

shape_noclip <- function(coords, el, params,
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

#' @rdname shapes
#' @export

shape_noplot <- function(coords, v=NULL, params) {
  invisible(NULL)
}

#' @rdname shapes
#' @export

add_shape <- function(shape, clip=shape_noclip,
                      plot=shape_noplot,
                      parameters=list()) {

  ## TODO
  ## checkScalarString(shape)
  ## checkFunction(clip)
  ## checkFunction(plot)
  ## checkList(parameters, named=TRUE)

  assign(shape, value=list(clip=clip, plot=plot), envir=.igraph.shapes)
  do.call(igraph.options, parameters)
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
      vertex.size2[ el[,2] ]
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
      vertex.size2[ el[,2] ]
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
      vertex.size2[ el[,2] ]
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

.igraph.shape.sphere.clip <- .igraph.shape.circle.clip

.igraph.shape.sphere.plot <- function(coords, v=NULL, params) {

  getparam <- function(pname) {
    p <- params("vertex", pname)
    if (length(p) != 1 && !is.null(v)) {
      p <- p[v]
    }
    p
  }
  vertex.color       <- rep(getparam("color"), length=nrow(coords))
  vertex.size        <- rep(1/200 * getparam("size"), length=nrow(coords))

  ## Need to create a separate image for every different vertex color
  allcols <- unique(vertex.color)
  images <- lapply(allcols, function(col) {
    img <- .Call("R_igraph_getsphere", pos=c(0.0,0.0,10.0), radius=7.0,
                 color=col2rgb(col)/255, bgcolor=c(0,0,0),
                 lightpos=list(c(-2,2,2)), lightcolor=list(c(1,1,1)),
                 width=100L, height=100L,
                 PACKAGE="igraph")
    as.raster(img)
  })
  whichImage <- match(vertex.color, allcols)  
  
  for (i in seq_len(nrow(coords))) {
    vsp2 <- vertex.size[i]
    rasterImage(images[[ whichImage[i] ]],
                coords[i,1]-vsp2, coords[i,2]-vsp2,
                coords[i,1]+vsp2, coords[i,2]+vsp2)
  }
}

.igraph.shape.raster.clip <- .igraph.shape.rectangle.clip

.igraph.shape.raster.plot <- function(coords, v=NULL, params) {

  getparam <- function(pname) {
    p <- params("vertex", pname)
    if (is.list(p) && length(p) != 1 && !is.null(v)) {
      p <- p[v]
    }
    p
  }

  size   <- rep(1/200 * getparam("size"), length=nrow(coords))
  size2  <- rep(1/200 * getparam("size2"), length=nrow(coords))
  raster <- getparam("raster")

  for (i in seq_len(nrow(coords))) {
    ras <- if (!is.list(raster) || length(raster)==1) raster else raster[[i]]
    rasterImage(ras, coords[i,1]-size[i], coords[i,2]-size2[i],
                coords[i,1]+size[i], coords[i,2]+size2[i])
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
.igraph.shapes[["sphere"]] <- list(clip=.igraph.shape.sphere.clip,
                                   plot=.igraph.shape.sphere.plot)
.igraph.shapes[["raster"]] <- list(clip=.igraph.shape.raster.clip,
                                   plot=.igraph.shape.raster.plot)
