\name{transitivity}
\alias{transitivity}
\concept{Transitivity}
\concept{Clustering coefficient}
\title{Transitivity of a graph}
\description{Transitivity measures the probability that the adjacent
  vertices of a vertex are connected. This is sometimes also called the
  clustering coefficient. }
\usage{
transitivity(graph, type=c("undirected", "global", "globalundirected",
       "localundirected", "local", "average", "localaverage",
       "localaverageundirected", "barrat", "weighted"), vids=NULL,
       weights=NULL, isolates=c("NaN", "zero"))
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{type}{The type of the transitivity to calculate. Possible
    values:
    \describe{
      \item{\code{global}}{The global transitivity of an
	undirected graph (directed graphs are considered as undirected ones
	as well). This is simply the ratio of the triangles and the
	connected triples in the graph. For directed graph the direction
	of the edges is ignored. }
      \item{\code{local}}{The local transitivity of an undirected graph,
	this is calculated for each vertex given in the \code{vids}
	argument. The local transitivity of a vertex is the ratio of the
	triangles connected to the vertex and the triples centered on
	the vertex. For directed graph the direction of the edges is
	ignored. }
      \item{\code{undirected}}{This is the same as \code{global}.}
      \item{\code{globalundirected}}{This is the same as \code{global}.}
      \item{\code{localundirected}}{This is the same as \code{local}.}
      \item{\code{barrat}}{The weighted transitivity as defined
	A. Barrat. See details below.}
      \item{\code{weighted}}{The same as \code{barrat}.}
    }
  }
  \item{vids}{The vertex ids for the local transitivity will be
    calculated. This will be ignored for global transitivity types.
    The default value is \code{NULL}, in this case all vertices are
    considered. It is slightly faster to supply \code{NULL} here than
    \code{V(graph)}.    
  }
  \item{weights}{Optional weights for weighted transitivity. It is
    ignored for other transitivity measures. If it is \code{NULL} (the
    default) and the graph has a \code{weight} edge attribute, then it
    is used automatically.
  }
  \item{isolates}{Character scalar, defines how to treat vertices with
    degree zero and one. If it is \sQuote{\code{NaN}} then they local
    transitivity is reported as \code{NaN} and they are not included in
    the averaging, for the transitivity types that calculate an
    average. If there are no vertices with degree
    two or higher, then the averaging will still result \code{NaN}. If
    it is \sQuote{\code{zero}}, then we report 0 transitivity for them,
    and they are included in the averaging, if an average is
    calculated.
  }
}
\details{
   Note that there are essentially two classes of transitivity measures,
   one is a vertex-level, the other a graph level property.

   There are several generalizations of transitivity to weighted graphs,
   here we use the definition by A. Barrat, this is a local vertex-level
   quantity, its formula is

   \deqn{C_i^w=\frac{1}{s_i(k_i-1)}\sum_{j,h}\frac{w_{ij}+w_{ih}}{2}a_{ij}a_{ih}a_{jh}}{
     weighted C_i = 1/s_i 1/(k_i-1) sum( (w_ij+w_ih)/2 a_ij a_ih a_jh, j, h)}

   \eqn{s_i}{s_i} is the strength of vertex \eqn{i}{i}, see
   \code{\link{graph.strength}}, \eqn{a_{ij}}{a_ij} are elements of the
   adjacency matrix, \eqn{k_i}{k_i} is the vertex degree,
   \eqn{w_{ij}}{w_ij} are the weights.
   
   This formula gives back the normal not-weighted local transitivity if
   all the edge weights are the same.

   The \code{barrat} type of transitivity does not work for graphs with
   multiple and/or loop edges. If you want to calculate it for a
   directed graph, call \code{\link{as.undirected}} with the
   \code{collapse} mode first.
}
\value{
  For \sQuote{\code{global}} a single number, or \code{NaN} if there are no
  connected triples in the graph.

  For \sQuote{\code{local}} a vector of transitivity scores, one for
  each vertex in \sQuote{\code{vids}}.
}
\references{ Wasserman, S., and Faust, K. (1994). \emph{Social Network
    Analysis: Methods and Applications.} Cambridge: Cambridge University
  Press.

  Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro
  Vespignani: The architecture of complex weighted networks,
  Proc. Natl. Acad. Sci. USA 101, 3747 (2004)
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
transitivity(g)
g2 <- erdos.renyi.game(1000, 10/1000)
transitivity(g2)   # this is about 10/1000

# Weighted version, the figure from the Barrat paper
gw <- graph.formula(A-B:C:D:E, B-C:D, C-D)
E(gw)$weight <- 1
E(gw)[ V(gw)[name == "A"] \%--\% V(gw)[name == "E" ] ]$weight <- 5
transitivity(gw, vids=V(gw)[name=="A"], type="weighted")

# Weighted reduces to "local" if weights are the same
gw2 <- erdos.renyi.game(1000, 10/1000)
E(gw2)$weight <- 1
all(transitivity(gw2, type="local") == transitivity(gw2, type="weighted"))
all( transitivity(gw2, vids=2:vcount(g)-1, type="local") ==
     transitivity(gw2, vids=2:vcount(g)-1, type="weighted") )
}
\keyword{graphs}

