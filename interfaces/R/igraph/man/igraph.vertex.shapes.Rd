\name{Vertex shapes}
\alias{igraph.vertex.shapes}
\concept{Vertex shapes}
\concept{Visualization}
\title{Various vertex shapes when plotting igraph graphs}
\description{Starting from version 0.5.1 igraph supports different
  vertex shapes when plotting graphs.}
\details{
  Note that the current vertex shape implementation is experimental and
  it might change in the future. Currently vertex shapes are implemented
  only for \code{\link{plot.igraph}}.
  
  In igraph a vertex shape is defined by a function that 1) provides
  information about the size of the shape for clipping the edges and 2)
  plots the shape if requested. These functions are called \dQuote{shape
    functions} in the rest of this manual page.

  Shape functions have a \sQuote{\code{mode}} argument that decides in
  which mode they should operate. \sQuote{\code{clip}} selects clipping
  mode and \sQuote{\code{plot}} selects plotting mode.

  In clipping mode a shape function has the following arguments:
  \itemize{
    \item{coords}{A matrix with four columns, it contains the
      coordinates of the vertices for the edge list supplied in the
      \code{el} argument.}
    \item{el}{A matrix with two columns, the edges of which some end
      points will be clipped. It should have the same number of rows as
      \code{coords}.}
    \item{mode}{\dQuote{\code{clip}} for choosing clipping mode.}
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

  In \sQuote{clipping} mode, a shape function should return a matrix
  with the same number of rows as the \code{el} arguments.
  If \code{end} is \code{both} then the matrix must have four
  columns, otherwise two. The matrix contains the modified coordinates,
  with the clipping applied.

  In plotting mode the following arguments are supplied to the shape
  function:
  \itemize{
    \item{coords}{The coordinates of the vertices, a matrix with two
      columns.}
    \item{v}{The ids of the vertices to plot. It should match the number
      of rows in the \code{coords} argument.}
    \item{mode}{\dQuote{\code{plot}} for choosing plotting mode.}
    \item{params}{The same as in clipping mode, see above.}
  }
  
  In \sQuote{plotting} mode the return value of the shape function is
  not used.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{igraph.plotting}}, \code{\link{plot.igraph}} }
\examples{
\dontrun{
g <- graph.ring(10, directed=TRUE, mut=TRUE)
plot(g, vertex.shape="rectangle", layout=layout.circle)
}
}
\keyword{graphs}

