\name{reciprocity}
\alias{reciprocity}
\concept{Reciprocity}
\title{Reciprocity of graphs}
\description{Calculates the reciprocity of a directed graph.}
\usage{
reciprocity(graph, ignore.loops = TRUE, mode = c("default", "ratio")) 
}
\arguments{
  \item{graph}{The graph object.}
  \item{ignore.loops}{Logical constant, whether to ignore loop edges.}
  \item{mode}{See below.}
}
\details{
The measure of reciprocity defines the proporsion of mutual
connections, in a directed graph. It is most commonly defined as
the probability that the opposite counterpart of a directed edge is
also included in the graph. Or in adjacency matrix notation: 
\eqn{\sum_{ij} (A\cdot A')_{ij}}{sum(i, j, (A.*A')ij) / sum(i, j, Aij)},
where \eqn{A\cdot A'}{A.*A'} is the element-wise product of matrix
\eqn{A} and its transpose. This measure is calculated if the \code{mode}
argument is \code{default}.

Prior to igraph version 0.6, another measure was implemented,
defined as the probability of mutual connection between a vertex
pair, if we know that there is a (possibly non-mutual) connection
between them. In other words, (unordered) vertex pairs are
classified into three groups: (1) not-connected, (2)
non-reciprocaly connected, (3) reciprocally connected. 
The result is the size of group (3), divided by the sum of group
sizes (2)+(3). This measure is calculated if \code{mode} is
\code{ratio}.
}
\value{A numeric scalar between zero and one.}
% \references{}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and
  Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- random.graph.game(20, 5/20, directed=TRUE)
reciprocity(g)
}
\keyword{graphs}
