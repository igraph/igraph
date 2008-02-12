\name{components}
\alias{subcomponent}
\concept{Subcomponent}
\title{In- or out- component of a vertex}
\description{Finds all vertices reachable from a given vertex, or the
  opposite: all vertices from which a given vertex is reachable via a
  directed path.}
\usage{
subcomponent(graph, v, mode = c("all", "out", "in"))
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{v}{The vertex to start the search from.}
  \item{mode}{Character string, either \dQuote{in}, \dQuote{out} or
    \dQuote{all}. If \dQuote{in} all vertices from which \code{v} is
    reachable are listed. If \dQuote{out} all vertices reachable from
    \code{v} are returned. If \dQuote{all} returns the union of
    these. It is ignored for undirected graphs.} 
}
\details{A breadh-first search is conducted starting from vertex
  \code{v}.}
\value{Numeric vector, the ids of the vertices in the same component as
  \code{v}.} 
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{clusters}}}
\examples{
g <- erdos.renyi.game(100, 1/200)
subcomponent(g, 1, "in")
subcomponent(g, 1, "out")
subcomponent(g, 1, "all")
}
\keyword{graphs}
