\name{all.st.cuts}
\alias{all.st.cuts}
\concept{Edge cuts}
\concept{(s,t)-cuts}
\title{List all (s,t)-cuts of a graph}
\description{
  List all (s,t)-cuts in a directed graph.
}
\usage{
all.st.cuts(graph, source, target)
}
\arguments{
  \item{graph}{The input graph. It must be directed.}
  \item{source}{The source vertex.}
  \item{target}{The target vertex.}
}
\details{
  Given a \eqn{G} directed graph and two, different and non-ajacent
  vertices, \eqn{s} and \eqn{t}, an \eqn{(s,t)}-cut is a set of edges,
  such that after removing these edges from \eqn{G} there is no directed
  path from \eqn{s} to \eqn{t}.  
}
\value{
  A list with entries:
  \item{cuts}{A list of numeric vectors containing edge ids. Each vector
    is an \eqn{(s,t)}-cut.}
  \item{partition1s}{A list of numeric vectors containing vertex
    ids, they correspond to the edge cuts. Each vertex set is a 
    generator of the corresponding cut, i.e. in the graph \eqn{G=(V,E)}, the
    vertex set \eqn{X} and its complementer \eqn{V-X}, generates the cut
    that contains exactly the edges that go from \eqn{X} to \eqn{V-X}.}
}
\references{
  JS Provan and DR Shier: A Paradigm for listing
  (s,t)-cuts in graphs, \emph{Algorithmica} 15, 351--372, 1996.    
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{
  \code{\link{all.st.mincuts}} to list all minimum cuts.
}
\examples{
# A very simple graph
g <- graph.formula(a -+ b -+ c -+ d -+ e)
all.st.cuts(g, source="a", target="e")

# A somewhat more difficult graph
g2 <- graph.formula(s --+ a:b, a:b --+ t,
                   a --+ 1:2:3, 1:2:3 --+ b)
all.st.cuts(g2, source="s", target="t")
}
\keyword{graphs}
