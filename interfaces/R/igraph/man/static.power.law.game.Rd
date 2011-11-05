\name{static.power.law.game}
\alias{static.power.law.game}
\concept{Random graph model}
\title{Scale-free random graphs, from vertex fitness scores}
\description{
  This function generates a non-growing random graph with expected
  power-law degree distributions.
}
\usage{
static.power.law.game (no.of.nodes, no.of.edges, exponent.out,
                       exponent.in = -1, loops = FALSE,
                       multiple = FALSE, finite.size.correction = TRUE) 
}
\arguments{
  \item{no.of.nodes}{The number of vertices in the generated graph.}
  \item{no.of.edges}{The number of edges in the generated graph.}
  \item{exponent.out}{Numeric scalar, the power law exponent of the
    degree distribution. For directed graphs, this specifies the
    exponent of the out-degree distribution. It must be greater than or
    equal to 2. If you pass \code{Inf} here, you will get back an
    Erdos-Renyi random network.}
  \item{exponent.in}{Numeric scalar. If negative, the generated graph
    will be undirected. If greater than or equal to 2, this argument
    specifies the exponent of the in-degree distribution. If
    non-negative but less than 2, an error will be generated.}
  \item{loops}{Logical scalar, whether to allow loop edges in the
    generated graph.}
  \item{multiple}{Logical scalar, whether to allow multiple edges in the
    generated graph.}
  \item{finite.size.correction}{Logical scalar, whether to use the
    proposed finite size correction of Cho et al., see references
    below.}
}
\details{
 This game generates a directed or undirected random graph where the
 degrees of vertices follow power-law distributions with prescribed
 exponents. For directed graphs, the exponents of the in- and out-degree
 distributions may be specified separately.

 The game simply uses \code{\link{static.fitness.game}} with appropriately
 constructed fitness vectors. In particular, the fitness of vertex
 \eqn{i} is \eqn{i^{-alpha}}{i^(-alpha)}, where \eqn{alpha = 1/(gamma-1)}
 and gamma is the exponent given in the arguments.

 To remove correlations between in- and out-degrees in case of directed
 graphs, the in-fitness vector will be shuffled after it has been set up
 and before \code{\link{static.fitness.game}} is called.

 Note that significant finite size effects may be observed for exponents
 smaller than 3 in the original formulation of the game. This function
 provides an argument that lets you remove the finite size effects by
 assuming that the fitness of vertex \eqn{i} is 
 \eqn{(i+i_0-1)^{-alpha}}{(i+i0-1)^(-alpha)}
 where \eqn{i_0}{i0} is a constant chosen appropriately to ensure that
 the maximum degree is less than the square root of the number of edges
 times the average degree; see the paper of Chung and Lu, and Cho et al
 for more details.
}
\value{An igraph graph, directed or undirected.}
\references{
  Goh K-I, Kahng B, Kim D: Universal behaviour of load distribution
  in scale-free networks. \emph{Phys Rev Lett} 87(27):278701, 2001.

 Chung F and Lu L: Connected components in a random graph with given
 degree sequences. \emph{Annals of Combinatorics} 6, 125-145, 2002.

 Cho YS, Kim JS, Park J, Kahng B, Kim D: Percolation transitions in
 scale-free networks under the Achlioptas process. \emph{Phys Rev Lett}
 103:135702, 2009.
}
\author{Tamas Nepusz \email{ntamas@gmail.com}}
\examples{
g <- static.power.law.game(10000, 30000, 2.2, 2.3)
\dontrun{plot(degree.distribution(g, cumulative=TRUE, mode="out"), log="xy")}
}
\keyword{graphs}
 
  
