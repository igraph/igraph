\name{Vertex shapes}
\alias{igraph.vertex.shapes}
\alias{vertex.shapes}
\alias{add.vertex.shape}
\alias{igraph.shape.noclip}
\alias{igraph.shape.noplot}
\concept{Vertex shapes}
\concept{Visualization}
\title{Various vertex shapes when plotting igraph graphs}
\description{Starting from version 0.5.1 igraph supports different
  vertex shapes when plotting graphs.}
\usage{
vertex.shapes (shape = NULL)
add.vertex.shape (shape, clip = igraph_shape.noclip,
       plot = igraph.shape.noplot, parameters = character())
igraph.shape.noclip (coords, el, params, end = c("both", "from", "to"))
igraph.shape.noplot (coords, v = NULL, params)
}
\arguments{
  \item{shape}{Character scalar, name of a vertex shape. If it is
    \code{NULL} for \code{vertex.shapes}, then the names of all defined
    vertex shapes are returned.} 
  \item{clip}{An R function object, the clipping function.}
  \item{plot}{An R function object, the plotting function.}
  \item{parameters}{Character vector, additional graph/vertex/edge
    parameters. Vertex parameters should have a prefix
    \sQuote{\code{vertex.}}, edge parameters a prefix
    \sQuote{\code{edge.}}.}
  \item{coords,el,params,end,v}{See parameters of the clipping/plotting
    functions below.} 
}
\details{
  In igraph a vertex shape is defined by two functions: 1) provides
  information about the size of the shape for clipping the edges and 2)
  plots the shape if requested. These functions are called \dQuote{shape
    functions} in the rest of this manual page. The first one is the
  clipping function and the second is the plotting function.
  
  The clipping function has the following arguments:
  \describe{
    \item{coords}{A matrix with four columns, it contains the
      coordinates of the vertices for the edge list supplied in the
      \code{el} argument.}
    \item{el}{A matrix with two columns, the edges of which some end
      points will be clipped. It should have the same number of rows as
      \code{coords}.}
    \item{params}{This is a function object that can be called to query
      vertex/edge/plot graphical parameters. The first argument of the
      function is \dQuote{\code{vertex}}, \dQuote{\code{edge}} or
      \dQuote{\code{plot}} to decide the type of the parameter, the
      second is a character string giving the name of the
      parameter. E.g.
      \preformatted{
	params("vertex", "size")
      }
    }
    \item{end}{Character string, it gives which end points will be
      used. Possible values are \dQuote{\code{both}},
      \dQuote{\code{from}} and \dQuote{\code{to}}. If
      \dQuote{\code{from}} the function is expected to clip the
      first column in the \code{el} edge list, \dQuote{\code{to}}
      selects the second column, \dQuote{\code{both}} selects both.}
  }

  The clipping function should return a matrix
  with the same number of rows as the \code{el} arguments.
  If \code{end} is \code{both} then the matrix must have four
  columns, otherwise two. The matrix contains the modified coordinates,
  with the clipping applied.

  The plotting function has the following arguments:
  \describe{
    \item{coords}{The coordinates of the vertices, a matrix with two
      columns.}
    \item{v}{The ids of the vertices to plot. It should match the number
      of rows in the \code{coords} argument.}
    \item{params}{The same as for the clipping function, see above.}
  }

  The return value of the plotting function is not used.

  \code{vertex.shapes} can be used to list the names of all installed
  vertex shapes, by calling it without arguments, or setting the
  \code{shape} argument to \code{NULL}. If a shape name is given, then
  the clipping and plotting functions of that shape are returned in a
  named list.

  \code{add.vertex.shape} can be used to add new vertex shapes to
  igraph. For this one must give the clipping and plotting functions of
  the new shape. It is also possible to list the plot/vertex/edge
  parameters, in the \code{parameters} argument, that the clipping
  and/or plotting functions can make use of. An example would be a
  generic regular polygon shape, which can have a parameter for the
  number of sides.

  \code{igraph.shape.noclip} is a very simple clipping function that the
  user can use in their own shape definitions. It does no clipping, the
  edges will be drawn exactly until the listed vertex position
  coordinates.

  \code{igraph.shape.noplot} is a very simple (and probably not very
  useful) plotting function, that does not plot anything.
}
\value{
  \code{vertex.shapes} returns a character vector if the \code{shape}
  argument is \code{NULL}. It returns a named list with entries named
  \sQuote{clip} and \sQuote{plot}, both of them R functions.

  \code{add.vertex.shape} returns \code{TRUE}, invisibly.

  \code{igraph.shape.noclip} returns the appropriate columns of its
  \code{coords} argument.

  \code{igraph.shape.noplot} returns \code{NULL}, invisibly.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{ \code{\link{igraph.plotting}}, \code{\link{plot.igraph}} }
\examples{
# all vertex shapes
shapes <- vertex.shapes()
g <- graph.ring(length(shapes))
set.seed(42)
plot(g, vertex.shape=shapes, vertex.label=shapes, vertex.label.dist=1,
     vertex.size=15, vertex.size2=10,
     vertex.pie=lapply(shapes, function(x) if (x=="pie") 2:6 else 0),
     vertex.pie.color=list(heat.colors(5)))

# add new vertex shape, plot nothing with no clipping
add.vertex.shape("nil")
plot(g, vertex.shape="nil")

# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size        <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }

  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add.vertex.shape("triangle", clip=vertex.shapes("circle")$clip,
                 plot=mytriangle)
plot(g, vertex.shape="triangle", vertex.color=rainbow(vcount(g)),
     vertex.size=seq(10,20,length=vcount(g)))
}
\keyword{graphs}

