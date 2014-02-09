\name{local.scan}
\alias{local.scan}
\alias{local.scan1.ecount}
\alias{local.scan1.ecount.approx}
\alias{local.scan1.ecount.approx.eigen}
\concept{Scan Statistics on Graphs}
\title{Compute local scan statistics on graphs}
\description{The scan statistic is a summary of the locality statistics
  that is computed from the local neighborhood of each vertex. This 
  \code{local.scan} function
  computes the local statistics for each vertex for a given neighborhood size and
  the statistic function.} 
\usage{
local.scan (graph.us, graph.them = NULL, k = 1, FUN = NULL,
    weighted = FALSE, mode = c("out", "in", "all"),
    neighborhoods = NULL, ...) 
}
\arguments{
  \item{graph.us}{an igraph object, the graph for which the scan statistics
  will be computed}
  \item{graph.them}{an igraph object or NULL, the graph for which 
  the scan statistics will be computed for \sQuote{them} statistics}
  \item{k}{an integer scalar, The size of the local neighborhood for each vertex.
  Should be non-negative.}
  \item{FUN}{character, function name for computing the local statistic in each
  neighborhood. If \code{NULL}(the default value), \code{ecount} is used
  for unweighted graphs (if \code{weighted=FALSE}) and a function that computes
  the sum of edge weights is used for weighted graphs (if \code{weighted=TRUE}) .}
  \item{weighted}{logical scalar, TRUE if the edge weights should be
  used for  computation of the scan statistic. If TRUE, the graph should be 
  weighted.}
  \item{mode}{The mode for the graph. One of \sQuote{\code{out}}, 
  \sQuote{\code{in}}, \sQuote{\code{all}} or \sQuote{\code{total}}}
  \item{neighborhoods}{A list of neighborhoods for each vertex. 
  The default value is  \code{NULL}.  If \code{graph.them} is NULL, and
  \code{neighborhoods} is not \code{NULL}, computes the \sQuote{them} variant
  of the scan statistics using the list of vertex neighborhoods in
  \code{neighborhoods} }
  \item{\dots}{Arguments passed to \code{FUN}, the function that computes 
  the local statistics}
}
\details{
  If both \code{graph.them} and  \code{neighborhoods} are \code{NULL}, 
  \code{local.scan} computes the \sQuote{us} variant of the scan statistics. 
  Otherwise,  \code{graph.them} should be an igraph object and the  \sQuote{them}
  variant is computed using \code{graph.them} to extract   the neighborhood information.
  If \code{graph.them} is   \code{NULL}, and   \code{neighborhoods} is not 
  \code{NULL}, \code{neighborhoods} should be a list of vertex neighborhoods, and
  the \sQuote{them}  variant is computed using \code{neighborhoods}. 
  \code{neighborhoods} should a list of length \code{vcount(graph.us)}. 
}
\value{A numeric vector containing the computed local statistics for each vertex}
\references{ Priebe, C. E., Conroy, J. M., Marchette, D. J., Park, Y. (2005).
Scan Statistics on Enron Graphs. 
\emph{Computational and Mathematical Organization Theory}.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\examples{


n <- 10^3
p <- 0.1
corr.rho <- 0.8
Graph.pair <- correlated.pair.game(n=n, corr=corr.rho, p=p)
local.0.us <-local.scan(graph.us= Graph.pair$graph1, k=0)
local.1.us <-local.scan(graph.us= Graph.pair$graph1, k=1)

local.0.them <-local.scan(graph.us= Graph.pair$graph1,
                           graph.them= Graph.pair$graph2, k=0)
local.1.them <-local.scan(graph.us= Graph.pair$graph1,
                           graph.them= Graph.pair$graph2, k=1)
Neigh.1 <- neighborhood(Graph.pair$graph2, order=1)

local.1.them.nhood <- local.scan(graph.us= Graph.pair$graph1, 
                                  neighborhoods= Neigh.1)

}
\keyword{graphs,scan statistics}