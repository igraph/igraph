\name{all.st.mincuts}
\alias{all.st.mincuts}
\concept{Edge cuts}
\concept{(s,t)-cuts}
\concept{Minimum cuts}
\concept{Minimum (s,t)-cuts}
\title{List all minimum (s,t)-cuts of a graph}
\description{
}
\usage{
all.st.mincuts(graph, source, target, capacity = NULL)
}
\arguments{
  \item{graph}{The input graph. It must be directed.}
  \item{source}{The id of the source vertex.}
  \item{target}{The id of the target vertex.}
  \item{capacity}{Numeric vector giving the edge capacities. If this is
    \code{NULL} and the graph has a \sQuote{\code{weight}} edge
    attribute, then this attribute defines the edge capacities. For
    forcing unit edge capacities, even for graphs that have a
    \sQuote{\code{weight}} edge attribute, supply \code{NA} here.}
}
\details{
  Given a \eqn{G} directed graph and two, different and non-ajacent
  vertices, \eqn{s} and \eqn{t}, an \eqn{(s,t)}-cut is a set of edges,
  such that after removing this edges from \eqn{G} there is no directed
  path from \eqn{s} to \eqn{t}.

  The size of an \eqn{(s,t)}-cut is defined as the sum of their
  capacities (or weights). For unweighed (=equally weighted) graphs,
  this is simply the number of edges.
  
  An \eqn{(s,t)}-cut is minimum if it is of the smallest possible
  size.
}
\value{
  A list with entries:
  \item{value}{Numeric scalar, the size of the minimum cut(s).}
  \item{cuts}{A list of numeric vectors containing edge ids. Each vector
    is a minimum \eqn{(s,t)}-cut.}
  \item{partition1s}{A list of numeric vectors containing vertex
    ids, they correspond to the edge cuts. Each vertex set is a 
    generator of the corresponding cut, i.e. the graph \eqn{G=(V,E)} the
    vertex set \eqn{X} and its complementer \eqn{V-X}, generates the cut
    that contains exactly the edges that go from \eqn{X} to \eqn{V-X}.}
}
\references{
  JS Provan and DR Shier: A Paradigm for listing
  (s,t)-cuts in graphs, \emph{Algorithmica} 15, 351--372, 1996.  
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{
  \code{\link{all.st.cuts}}, \code{\link{minimum.size.separators}}
}
\examples{
# A difficult graph, from the Provan-Shier paper
g <- graph.formula(s --+ a:b, a:b --+ t,
                   a --+ 1:2:3:4:5, 1:2:3:4:5 --+ b)
all.st.mincuts(g, source="s", target="t")
}
\keyword{graphs}
