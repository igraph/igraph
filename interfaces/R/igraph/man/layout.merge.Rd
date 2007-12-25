\name{layout.merge}
\alias{layout.merge}
\alias{piecewise.layout}
\title{Merging graph layouts}
\description{Place several graphs on the same layout}
\usage{
layout.merge(graphs, layouts, method = "dla",
       verbose = igraph.par("verbose"))
piecewise.layout(graph, layout=layout.kamada.kawai, \dots)
}
\arguments{
  \item{graphs}{A list of graph objects.}
  \item{layouts}{A list of two-column matrices.}
  \item{method}{Character constant giving the method to use. Right now
    only \code{dla} is implemented.}
  \item{verbose}{Logical constant, whether to show a progress bar while
    doing the calculation.}
  \item{graph}{The input graph.}
  \item{layout}{A function object, the layout function to use.}
  \item{\dots}{Additional arguments to pass to the \code{layout} layout
    function.}
}
\details{
  \code{layout.merge} takes a list of graphs and a list of coordinates and
  places the graphs in a common layout. The method to use is chosen via
  the \code{method} parameter, although right now only the \code{dla}
  method is implemented.

  The \code{dla} method covers the graph with circles.
  Then it sorts the graphs based on the number of vertices
  first and places the largest graph at the center of the layout. Then
  the other graphs are placed in decreasing order via a DLA (diffision
  limited aggregation) algorithm: the graph is placed randomly on a
  circle far away from the center and a random walk is conducted until
  the graph walks into the larger graphs already placed or walks too far
  from the center of the layout.

  The \code{piecewise.layout} function disassembles the graph first into
  maximal connected components and calls the supplied \code{layout}
  function for each component separately. Finally it merges the layouts
  via calling \code{layout.merge}.
}
\value{
  A matrix with two columns and as many lines as the total number of
  vertices in the graphs.
}
% \references
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{plot.igraph}}, \code{\link{tkplot}},
  \code{\link{layout}}, \code{\link{graph.disjoint.union}}
}
\examples{
# create 20 scale-free graphs and place them in a common layout
graphs <- lapply(sample(5:20, 20, replace=TRUE),
          barabasi.game, directed=FALSE)
layouts <- lapply(graphs, layout.kamada.kawai)
lay <- layout.merge(graphs, layouts)
g <- graph.disjoint.union(graphs)
\dontrun{plot(g, layout=lay, vertex.size=3, labels=NA, edge.color="black")}
}
\keyword{graphs}
