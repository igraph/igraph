\name{motifs}
\alias{graph.motifs}
\alias{motifs}
\alias{graph.motifs.no}
\alias{count_motifs}
\alias{graph.motifs.est}
\alias{sample_motifs}
\concept{Graph motif}
\title{Graph motifs}
\description{Graph motifs are small connected subgraphs with a
  well-defined structure. These functions search a graph for various
  motifs.}
\usage{
motifs(graph, size = 3, cut.prob = rep(0, size)) 
count_motifs(graph, size = 3, cut.prob = rep(0, size)) 
sample_motifs(graph, size = 3, cut.prob = rep(0, size), sample.size =
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
  \code{motifs} searches a graph for motifs of a given size and
  returns a numeric vector containing the number of different
  motifs. The order of the motifs is defined by their isomorphism class,
  see \code{\link{iso_class}}.

  \code{count_motifs} calculates the total number of motifs of a
  given size in graph.

  \code{sample_motifs} estimates the total number of motifs of a
  given size in a graph based on a sample. 
}
\value{
  \code{motifs} returns a numeric vector, the number of occurences
  of each motif in the graph. The motifs are ordered by their
  isomorphism classes. Note that for unconnected subgraphs, which are
  not considered to be motifs, the result will be \code{NA}.

  \code{count_motifs} and \code{sample_motifs} return a numeric
  constant.
}
% \references{}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{iso_class}}}
\examples{
g <- barabasi.game(100)
motifs(g, 3)
count_motifs(g, 3)
sample_motifs(g, 3)
}
\keyword{graphs}

  
