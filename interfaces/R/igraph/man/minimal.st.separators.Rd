\name{minimal.st.separators}
\alias{minimal.st.separators}
\concept{Minimal (s,t) separators}
\concept{Vertex separator}
\title{Minimum size vertex separators}
\description{List all vertex sets that are minimal (s,t) separators for
  some s and t, in an undirected graph.}
\usage{
minimal.st.separators(graph)
}
\arguments{
  \item{graph}{The input graph. It may be directed, but edge directions
    are ignored.}
}
\details{
  A \eqn{(s,t)} vertex separator is a set of vertices, such that after
  their removal from the graph, there is no path between \eqn{s} and
  \eqn{t} in the graph.

  A \eqn{(s,t)} vertex separator is minimal if none of its subsets is an
  \eqn{(s,t)} vertex separator.
}
\value{
  A list of numeric vectors. Each vector contains a vertex set (defined
  by vertex ids), each vector is an (s,t) separator of the input graph,
  for some \eqn{s} and \eqn{t}.
}
\references{
  Anne Berry, Jean-Paul Bordat and Olivier Cogis: Generating All the
  Minimal Separators of a Graph, In: Peter Widmayer, Gabriele Neyer
  and Stephan Eidenbenz (editors): \emph{Graph-theoretic concepts in
    computer science}, 1665, 167--172, 1999. Springer. 
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
%\seealso{}
\examples{
ring <- graph.ring(4)
minimal.st.separators(ring)

chvatal <- graph.famous("chvatal")
minimal.st.separators(chvatal)
}
\keyword{graphs}
