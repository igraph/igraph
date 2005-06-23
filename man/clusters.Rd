\name{clusters}
\alias{clusters}
\alias{cluster.distribution}
\title{Connected components of a graph}
\description{}
\usage{
clusters(graph, mode="weak")
cluster.distribution(graph, cumulative = FALSE, mul.size = FALSE, \dots)
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{mode}{Character string, either \dquote{weak} or \dquote{strong}.
    For directed graphs \dquote{weak} implies weakly, \dquote{strong}
    strongly connected components to search. It is ignored for
    undirected graphs.}
  \item{cumulative}{Logical, if TRUE the cumulative distirubution
    (realtive frequency) is calculated.}
  \item{mul.size}{Logical. If TRUE the relative frequencies will be
    multiplied by the cluster sizes.}
}
\details{Bread-first search is conducted from each not-yet visited
  vertex. }
\value{A named list with two components:
  \itemize{
    \item{membership}{numeric vector giving the cluster id to which each
      vertex belongs.}
    \item{csize}{numeric vector giving the sizes of the clusters.}
  }
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{subcomponent}}}
\examples{
g <- erdos.renyi.game(20, 1/20)
clusters(g)
}
\keyword{graphs}
