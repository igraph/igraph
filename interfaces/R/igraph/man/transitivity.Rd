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
       "localaverageundirected"), vids=NULL)
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{type}{The type of the transitivity to calculate. Possible
    values:
    \itemize{
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
    }
  }
  \item{vids}{The vertex ids for the local transitivity will be
    calculated. This will be ignored for global transitivity types.
    The default value is \code{NULL}, in this case all vertices are
    considered. It is slightly faster to supply \code{NULL} here than
    \code{V(graph)}.
  }
}
%\details{}
\value{
  For \sQuote{\code{global}} a single number, or \code{NaN} if there are no
  connected triples in the graph.

  For \sQuote{\code{local}} a vector of transitivity scores, one for
  each vertex in \sQuote{\code{vids}}.
}
\references{ Wasserman, S., and Faust, K. (1994). \emph{Social Network
  Analysis: Methods and Applications.} Cambridge: Cambridge University
  Press. }
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
transitivity(g)
g2 <- erdos.renyi.game(1000, 10/1000)
transitivity(g2)   # this is about 10/1000
}
\keyword{graphs}
