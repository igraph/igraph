\name{clusters}
\alias{no.clusters}
\alias{clusters}
\alias{is.connected}
\alias{cluster.distribution}
\concept{Connectedness}
\concept{Graph component}
\title{Connected components of a graph}
\description{Calculate the maximal (weakly or strongly) 
  connected components of a graph
}
\usage{
is.connected(graph, mode=c("weak", "strong"))
clusters(graph, mode=c("weak", "strong"))
no.clusters(graph, mode=c("weak", "strong"))
cluster.distribution(graph, cumulative = FALSE, mul.size = FALSE, \dots)
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{mode}{Character string, either \dQuote{weak} or \dQuote{strong}.
    For directed graphs \dQuote{weak} implies weakly, \dQuote{strong}
    strongly connected components to search. It is ignored for
    undirected graphs.}
  \item{cumulative}{Logical, if TRUE the cumulative distirubution
    (relative frequency) is calculated.}
  \item{mul.size}{Logical. If TRUE the relative frequencies will be
    multiplied by the cluster sizes.}
  \item{\dots}{Additional attributes to pass to \code{cluster}, right
    now only \code{mode} makes sense.}
}
\details{
  \code{is.connected} decides whether the graph is weakly or strongly
  connected.

  \code{clusters} finds the maximal (weakly or strongly) connected
  components of a graph.

  \code{no.clusters} does almost the same as \code{clusters} but returns
  only the number of clusters found instead of returning the actual
  clusters.

  \code{cluster.distribution} creates a histogram for the maximal
  connected component sizes.
  
  Breadth-first search is conducted from each not-yet visited
  vertex. }
\value{
  For \code{is.connected} a logical constant.
  
  For \code{clusters} a named list with three components:
  \item{membership}{numeric vector giving the cluster id to which each
    vertex belongs.}
  \item{csize}{numeric vector giving the sizes of the clusters.}
  \item{no}{numeric constant, the number of clusters.}

  For \code{no.clusters} an integer constant is returned.
  
  For \code{cluster.distribution} a numeric vector with the relative
  frequencies. The length of the vector is the size of the largest
  component plus one. Note that (for currently unknown reasons) the
  first element of the vector is the number of clusters of size zero, so
  this is always zero.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{subcomponent}}}
\examples{
g <- erdos.renyi.game(20, 1/20)
clusters(g)
}
\keyword{graphs}
