\name{constraint}
\alias{constraint}
\concept{Burt's constraint}
\title{Burt's constraint}
\description{Given a graph, \code{constraint} calculates Burt's constraint for
  each vertex.
}
\usage{
constraint(graph, nodes=V(graph), weights=NULL)
}
\arguments{
  \item{graph}{A graph object, the input graph.}
  \item{nodes}{The vertices for which the constraint will be
    calculated. Defaults to all vertices.}
  \item{weights}{The weights of the edges. If this is \code{NULL} and
    there is a \code{weight} edge attribute this is used. If there is no
    such edge attribute all edges will have the same weight.}
}
\details{Burt's constraint is higher if ego has less, or mutually
  stronger related (i.e. more redundant) contacts. Burt's measure of
  constraint, \eqn{C_i}{C[i]}, of vertex \eqn{i}'s ego network
  \eqn{V_i}{V[i]}, is defined for directed and valued graphs,
  \deqn{C_i=\sum_{j \in V_i \setminus \{i\}} (p_{ij}+\sum_{q \in V_i
      \setminus \{i,j\}} p_{iq} p_{qj})^2}{%
    C[i] = sum( sum( p[i,j] + p[i,q] p[q,j])^2, q in V[i], q != i,j ), j in
    V[], j != i).
  }
  for a graph of order (ie. number of vertices) \eqn{N}, where
  proportional tie strengths are defined as 
  \deqn{p_{ij} = \frac{a_{ij}+a_{ji}}{\sum_{k \in V_i \setminus \{i\}}(a_{ik}+a_{ki})},}{%
    p[i,j]=(a[i,j]+a[j,i]) / sum(a[i,k]+a[k,i], k in V[i], k != i),
  }
  \eqn{a_{ij}}{a[i,j]} are elements of \eqn{A} and the latter being the
  graph adjacency matrix. For isolated vertices, constraint is undefined.
}
\value{A numeric vector of constraint scores}
\author{Jeroen Bruggeman \email{j.p.bruggeman@uva.nl} and Gabor Csardi
  \email{csardi@rmki.kfki.hu}
}
\references{Burt, R.S. (2004). Structural holes and good ideas. \emph{American
    Journal of Sociology} 110, 349-399.
}
\examples{
g <- erdos.renyi.game(20, 5/20)
constraint(g)
}
\keyword{graphs}

