\name{conversion between igraph and graphNEL graphs}
\alias{igraph.from.graphNEL}
\alias{igraph.to.graphNEL}
\concept{Conversion}
\concept{graph package}
\concept{graphNEL object}
\title{Convert igraph graphs to graphNEL objects or back}
\description{The graphNEL class is defined in the \code{graph} package,
  it is another way to represent graphs. These functions are provided to
  convert between the igraph and the graphNEL objects.}
\usage{
igraph.from.graphNEL(graphNEL, name = TRUE, weight = TRUE,
              unlist.attrs = TRUE)
igraph.to.graphNEL(graph)
}
\arguments{
  \item{graphNEL}{The graphNEL graph.}
  \item{name}{Logical scalar, whether to add graphNEL vertex names as an
    igraph vertex attribute called \sQuote{\code{name}}.}
  \item{weight}{Logical scalar, whether to add graphNEL edge weights as
    an igraph edge attribute called \sQuote{\code{weight}}. (graphNEL
    graphs are always weighted.)}
  \item{unlist.attrs}{Logical scalar. graphNEL attribute query functions
    return the values of the attributes in R lists, if this argument is
    \code{TRUE} (the default) these will be converted to atomic vectors,
    whenever possible, before adding them to the igraph graph.}
  \item{graph}{An igraph graph object.}
}
\details{
  \code{igraph.from.graphNEL} takes a graphNEL graph and converts it to
  an igraph graph. It handles all graph/vertex/edge attributes. If the
  graphNEL graph has a vertex attribute called \sQuote{\code{name}} it
  will be used as igraph vertex attribute \sQuote{\code{name}} and the
  graphNEL vertex names will be ignored.

  Because graphNEL graphs poorly support multiple edges, the edge
  attributes of the multiple edges are lost: they are all replaced by
  the attributes of the first of the multiple edges.
  
  \code{igraph.to.graphNEL} converts and igraph graph to a graphNEL
  graph. It converts all graph/vertex/edge attributes. If the igraph
  graph has a vertex attribute \sQuote{\code{name}}, then it will be
  used to assign vertex names in the graphNEL graph. Otherwise igraph
  vertex ids will be used for this purpose.
}
\value{
  \code{igraph.from.graphNEL} returns an igraph graph object.

  \code{igraph.to.graphNEL} returns a graphNEL graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{get.adjacency}}, \code{\link{graph.adjacency}},
  \code{\link{get.adjlist}} and \code{\link{graph.adjlist}}.
}
\examples{
g <- graph.ring(10)
V(g)$name <- letters[1:10]
GNEL <- igraph.to.graphNEL(g)
g2 <- igraph.from.graphNEL(GNEL)
g2
}
\keyword{graphs}
