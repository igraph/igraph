\name{bipartite.random.game}
\alias{bipartite.random.game}
\concept{Random graph model}
\concept{Bipartite graph}
\title{Bipartite random graphs}
\description{Generate bipartite graphs using the Erdos-Renyi model}
\usage{
bipartite.random.game(n1, n2, type = c("gnp", "gnm"), p, m,
    directed = FALSE, mode = c("out", "in", "all")) 
}
\arguments{
  \item{n1}{Integer scalar, the number of bottom vertices.}
  \item{n2}{Integer scalar, the number of top vertices.}
  \item{type}{Character scalar, the type of the graph, \sQuote{gnp}
    creates a $G(n,p)$ graph, \sQuote{gnm} creates a $G(n,m)$ graph. See
    details below.} 
  \item{p}{Real scalar, connection probability for $G(n,p)$
    graphs. Should not be given for $G(n,m)$ graphs.}
  \item{m}{Integer scalar, the number of edges for $G(n,p)$
    graphs. Should not be given for $G(n,p)$ graphs.}
  \item{directed}{Logical scalar, whether to create a directed
    graph. See also the \code{mode} argument.}
  \item{mode}{Character scalar, specifies how to direct the edges in
    directed graphs. If it is \sQuote{out}, then directed edges point
    from bottom vertices to top vertices. If it is \sQuote{in}, edges
    point from top vertices to bottom vertices. \sQuote{out} and
    \sQuote{in} do not generate mutual edges. If this argument is 
    \sQuote{all}, then each edge direction is considered
    independently and mutual edges might be generated. This
    argument is ignored for undirected graphs.
  }
}
\details{
  Similarly to unipartite (one-mode) networks, we can define the 
  $G(n,p)$, and $G(n,m)$ graph classes for bipartite graphs, via their 
  generating process. In $G(n,p)$ every possible edge between top and
  bottom vertices is realized with probablity $p$, independently of the 
  rest of the edges. In $G(n,m)$, we uniformly choose $m$ edges to
  realize.
}
\value{A bipartite igraph graph.}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{erdos.renyi.game}} for the unipartite version.}
\examples{
## empty graph
bipartite.random.game(10, 5, p=0)

## full graph
bipartite.random.game(10, 5, p=1)

## random bipartite graph
bipartite.random.game(10, 5, p=.1)

## directed bipartite graph, G(n,m)
bipartite.random.game(10, 5, type="Gnm", m=20, directed=TRUE, mode="all")
}
\keyword{graphs}