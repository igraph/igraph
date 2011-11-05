\name{static.fitness.game}
\alias{static.fitness.game}
\concept{Random graph model}
\title{Random graphs from vertex fitness scores}
\description{
  This function generates a non-growing random graph with edge probabilities
  proportional to node fitness scores.
}
\usage{
static.fitness.game (no.of.edges, fitness.out, fitness.in,
                     loops = FALSE, multiple = FALSE) 

}
\arguments{
  \item{no.of.edges}{The number of edges in the generated graph.}
  \item{fitness.out}{A numeric vector containing the fitness of each
    vertex. For directed graphs, this specifies the out-fitness
    of each vertex.}
  \item{fitness.in}{If \code{NULL} (the default), the generated graph
    will be undirected. If not \code{NULL}, then it should be a numeric
    vector and it specifies the in-fitness of each vertex.

    If this argument is not \code{NULL}, then a directed graph is
    generated, otherwise an undirected one.
  }
  \item{loops}{Logical scalar, whether to allow loop edges in the
    graph.}
  \item{multiple}{Logical scalar, whether to allow multiple edges in the
    graph.}
}
\details{
  This game generates a directed or undirected random graph where the
  probability of an edge between vertices \eqn{i} and \eqn{j} depends on
  the fitness scores of the two vertices involved. For undirected graphs,
  each vertex has a single fitness score. For directed graphs, each
  vertex has an out- and an in-fitness, and the probability of an edge
  from \eqn{i} to \eqn{j} depends on the out-fitness of vertex \eqn{i}
  and the in-fitness of vertex \eqn{j}.
  
  The generation process goes as follows. We start from \eqn{N}
  disconnected nodes (where \eqn{N} is given by the length of the fitness
  vector). Then we randomly select two vertices \eqn{i} and \eqn{j}, with
  probabilities proportional to their fitnesses. (When the generated
  graph is directed, \eqn{i} is selected according to the out-fitnesses
  and \eqn{j} is selected according to the in-fitnesses). If the
  vertices are not connected yet (or if multiple edges are allowed), we 
  connect them; otherwise we select a new pair. This is repeated until
  the desired number of links are created.
  
  It can be shown that the \emph{expected} degree of each vertex will be
  proportional to its fitness, although the actual, observed degree will not
  be. If you need to generate a graph with an exact degree sequence, consider
  \code{\link{degree.sequence.game}} instead.
  
  This model is commonly used to generate static scale-free networks. To
  achieve this, you have to draw the fitness scores from the desired power-law
  distribution. Alternatively, you may use
  \code{\link{static.power.law.game}} which generates the fitnesses for
  you with a given exponent.
}
\value{An igraph graph, directed or undirected.}
\references{
  Goh K-I, Kahng B, Kim D: Universal behaviour of load distribution
  in scale-free networks. \emph{Phys Rev Lett} 87(27):278701, 2001.
}
\author{Tamas Nepusz \email{ntamas@gmail.com}}
\examples{
N <- 10000
g <- static.fitness.game(5*N, sample((1:50)^-2, N, replace=TRUE))
degree.distribution(g)
\dontrun{plot(degree.distribution(g, cumulative=TRUE), log="xy")}
}
\keyword{graphs}
