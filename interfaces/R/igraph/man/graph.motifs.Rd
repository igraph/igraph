\name{graph-motifs}
\alias{graph.motifs}
\alias{graph.motifs.no}
\alias{graph.motifs.est}
\concept{Graph motif}
\title{Graph motifs}
\description{Graph motifs are small subgraphs with a well-defined
  strucure. These functions search a graph for various motifs.}
\usage{
graph.motifs(graph, size = 3, cut.prob = rep(0, size)) 
graph.motifs.no(graph, size = 3, cut.prob = rep(0, size)) 
graph.motifs.est(graph, size = 3, cut.prob = rep(0, size), sample.size =
    vcount(graph)/10, sample = NULL) 
}
\arguments{
  \item{graph}{Graph object, the input graph.}
  \item{size}{The size of the motif, currently 3 and 4 are supported
    only.}
  \item{cut.prob}{Numeric vector giving the probabilities that the
    search graph is cut at a certain level. Its length should be the
    same as the size of the motif (the \code{size} argument). By default
    no cuts are made.}
  \item{sample.size}{The number of vertices to use as a starting point
    for finding motifs. Only used if the \code{sample} argument is
    \code{NULL}.}
  \item{sample}{If not \code{NULL} then it specifies the vertices to use
    as a starting point for finding motifs.}
}
\details{
  \code{graph.motifs} searches a graph for motifs of a given size and
  returns a numeric vector containing the number of different
  motifs. The order of the motifs is defined by their isomorphism class,
  see \code{\link{graph.isoclass}}.

  \code{graph.motifs.no} calculates the total number of motifs of a
  given size in graph.

  \code{graph.motifs.est} estimates the total number of motifs of a
  given size in a graph based on a sample. 
}
\value{
  \code{graph.motifs} returns a numeric vector.

  \code{graph.motifs.no} and \code{graph.motifs.est} return a numeric
  constant.
}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.isoclass}}}
\examples{
g <- barabasi.game(100)
graph.motifs(g, 3)
graph.motifs.no(g, 3)
graph.motifs.est(g, 3)
}
\keyword{graphs}

  