\name{is_dag}
\alias{is.dag}
\alias{is_dag}
\concept{DAG}
\title{Directed acyclic graphs}
\description{This function tests whether the given graph is a DAG, a
  directed acyclic graph.}
\usage{
is_dag(graph)
}
\arguments{
  \item{graph}{The input graph. It may be undirected, in which case
    \code{FALSE} is reported.}
}
\details{
  \code{is_dag} checks whether there is a directed cycle in the
  graph. If not, the graph is a DAG.
}
\value{A logical vector of length one.}
\author{Tamas Nepusz \email{ntamas@gmail.com} for the C code,
  Gabor Csardi \email{csardi.gabor@gmail.com} for the R interface.
}
\examples{
g <- make_tree(10)
is_dag(g)
g2 <- g + edge(5,1)
is_dag(g2)
}
\keyword{graphs}
