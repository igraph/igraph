\name{local.scan}
\alias{local.scan}
\alias{local.scan1.ecount}
\alias{local.scan1.ecount.approx}
\alias{local.scan1.ecount.approx.eigen}
\concept{Scan Statistics on Graphs}
\title{Compute local scan statistics on graphs}
\description{The scan statistic is a summary of the locality statistics
  that is computed from the local neighborhood of each vertex. The
  \code{local.scan} function computes the local statistics for each
  vertex for a given neighborhood size and the statistic function.}
\usage{
local.scan (graph.us, graph.them = NULL, k = 1, FUN = NULL,
    weighted = FALSE, mode = c("out", "in", "all"),
    neighborhoods = NULL, ...) 
}
\arguments{
  \item{graph.us}{An igraph object, the graph for which the scan statistics
    will be computed}
  \item{graph.them}{An igraph object or \code{NULL}, if not \code{NULL},
    then the \sQuote{them} statistics is computed, i.e. the
    neighborhoods calculated from \code{graph.us} are evaluated on
    \code{graph.them}.}
  \item{k}{An integer scalar, the size of the local neighborhood for
    each vertex. Should be non-negative.}
  \item{FUN}{Character, a function name, or a function object itself,
    for computing the local statistic in each neighborhood. If
    \code{NULL}(the default value), \code{ecount} is used for unweighted
    graphs (if \code{weighted=FALSE}) and a function that computes the
    sum of edge weights is used for weighted graphs (if
    \code{weighted=TRUE}). This argument is ignored if \code{k} is
    zero.}
  \item{weighted}{Logical scalar, TRUE if the edge weights should be
    used for computation of the scan statistic. If TRUE, the graph
    should be weighted. Note that this argument is ignored if \code{FUN}
    is not \code{NULL}, \code{"ecount"} and \code{"sumweights"}.
  \item{mode}{Character scalar, the kind of neighborhoods to use for the
    calculation. One of \sQuote{\code{out}}, \sQuote{\code{in}},
    \sQuote{\code{all}} or \sQuote{\code{total}}. This argument is
    ignored for undirected graphs.}
  \item{neighborhoods}{A list of neighborhoods, one for each vertex,
    or \code{NULL}. If it is not \code{NULL}, then the function is
    evaluated on the induced subgraphs specified by these
    neighborhoods.

    In theory this could be useful if the same \code{graph.us} graph
    is used for multiple \code{graph.them} arguments. Then the
    neighborhoods can be calculated on \code{graph.us} and used with
    multiple graphs. In practice, this is currently slower than simply
    using \code{graph.them} multiple times.
  }
  \item{\dots}{Arguments passed to \code{FUN}, the function that computes 
    the local statistics.}
}
\details{
  See the given reference below for the details on the local scan
  statistics.

  If \code{graph.them} is \code{NULL}, then \code{local.scan} computes
  the \sQuote{us} variant of the scan statistics.
  Otherwise, \code{graph.them} should be an igraph object and the
  \sQuote{them} variant is computed using \code{graph.us} to extract
  the neighborhood information, and applying \code{FUN} on these
  neighborhoods in \code{graph.them}.
}
\value{Typically a numeric vector containing the computed local
  statistics for each vertex. In general a list or vector of objects, as
  returned by \code{FUN}.}
\references{Priebe, C. E., Conroy, J. M., Marchette, D. J., Park,
  Y. (2005). Scan Statistics on Enron Graphs. \emph{Computational and
  Mathematical Organization Theory}.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\examples{

Graph.pair <- correlated.pair.game(n=10^3, corr=0.8, p=0.1)
local.0.us <- local.scan(graph.us=Graph.pair$graph1, k=0)
local.1.us <- local.scan(graph.us=Graph.pair$graph1, k=1)

local.0.them <- local.scan(graph.us=Graph.pair$graph1,
                           graph.them=Graph.pair$graph2, k=0)
local.1.them <- local.scan(graph.us=Graph.pair$graph1,
                           graph.them=Graph.pair$graph2, k=1)

Neigh.1 <- neighborhood(Graph.pair$graph1, order=1)
local.1.them.nhood <- local.scan(graph.us=Graph.pair$graph1,
                                 graph.them=Graph.pair$graph2,
                                 neighborhoods=Neigh.1)
}
\keyword{graphs}
