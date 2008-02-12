\name{reciprocity}
\alias{reciprocity}
\concept{Reciprocity}
\title{Reciprocity of graphs}
\description{Calculates the reciprocity of a directed graph.}
\usage{
reciprocity(graph, ignore.loops = TRUE) 
}
\arguments{
  \item{graph}{The graph object.}
  \item{ignore.loops}{Logical constant, whether to ignore loop edges.}
}
\details{A vertex pair (A, B) is said to be reciprocal if there are
  edges between them in both directions. The reciprocity of a directed
  graph is the proportion of all possible (A, B) pairs which are
  reciprocal, provided there is at least one edge between A and B. The
  reciprocity of an empty graph is undefined (results in an error
  code). Undirected graphs always have a reciprocity of 1.0 unless they
  are empty.}
\value{A numeric constaint between zero and one.}
% \references{}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and
  Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- random.graph.game(20, 5/20, directed=TRUE)
reciprocity(g)
}
\keyword{graphs}
