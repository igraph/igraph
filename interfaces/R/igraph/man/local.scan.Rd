\name{local.scan}
\alias{local.scan}
\concept{Scan Statistics on Graphs}
\title{Compute local scan statistics on static and time-varying graphs}
\description{The scan statistic is a summary of the locality statistics
  that is computed from the local neighborhood of each vertex.} 
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
  \item{FUN}{character, function name for computing the local statistic in each neighborhood. If \code{NULL}(the default value), \code{ecount} is used
  for unweighted graphs (if \code{weighted=FALSE}) and a function that computes
  the sum of edge weights is used for weighted graphs (if \code{weighted=TRUE}) .}
  \item{weighted}{logical scalar, TRUE if the edge weights should be used for  computation of the scan statistic. If TRUE, the graph should be weighted  }
  \item{mode}{The mode for the graph. One of \sQuote{\code{out}}, 
  \sQuote{\code{in}},\sQuote{\code{all}} or \sQuote{\code{total}}}
  \item{neighborhoods}{A list of neighborhoods for each vertex. 
  The default value is  \code{NULL}.  If \code{graph.them} is NULL, and
  \code{neighborhoods} is not \code{NULL}, computes the \sQuote{them} variant
  of the scan statistics using the list of vertex neighborhoods in
  \code{neighborhoods} }
  \item{\dots}{Arguments passed to \code{FUN}, the function that computes the local statistics}
}
\details{
  If both \code{graph.them} and  \code{neighborhoods} are \code{NULL}, 
  \code{local.scan} computes the \sQuote{us} variant of the scan statistics. 
  Otherwise,  \code{graph.them} should be an igraph object and the  \sQuote{them} variant is computed using \code{graph.them} to extract   the neighborhood information.
  If \code{graph.them} is   \code{NULL}, and   \code{neighborhoods} is not 
  \code{NULL}, \code{neighborhoods} should be a list of vertex neighborhoods, and
  the \sQuote{them}  variant is computed using \code{neighborhoods}. 
  \code{neighborhoods} should a list of length \code{vcount(graph.us)}. 
}
\value{The computed scan statistics. TODO}
\references{ Priebe, C. E., Conroy, J. M., Marchette, D. J., Park, Y. (2005).
Scan Statistics on Enron Graphs. \emph{Computational and Mathematical Organization Theory}.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{TODO}
\examples{
##TODO}
\keyword{graphs}