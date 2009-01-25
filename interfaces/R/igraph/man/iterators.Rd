\name{iterators}
\alias{iterators}
\alias{V}
\alias{E}
\alias{V<-}
\alias{E<-}
\alias{\%--\%}
\alias{\%->\%}
\alias{\%<-\%}
\alias{[<-.igraph.es}
\alias{[.igraph.es}
\alias{$<-.igraph.es}
\alias{$.igraph.es}
\alias{[<-.igraph.vs}
\alias{[.igraph.vs}
\alias{$<-.igraph.vs}
\alias{$.igraph.vs}
\alias{print.igraph.es}
\alias{print.igraph.vs}
\title{Vertex and edge sequences and iterators}
\description{Vertex and edge sequences are central concepts of
  igraph.}
\usage{
V(graph)
E(graph, P=NULL, path=NULL, directed=TRUE)
}
\arguments{
  \item{graph}{A graph object.}
  \item{P}{Numeric vector for selecting edges by giving their end
    points. See details below.}
  \item{path}{Numeric vector, this is for selecting all edges along a
    path. See also details below.}
  \item{directed}{Logcal constant, can be supplied only if either
    \code{P} or \code{path} is also present and gives whether the pairs
    or the path are directed or not.} 
}
\details{
  It is often needed to perform an operation on a subset of vertices of
  edges in a graph.

  A vertex sequence is simply a vector containing vertex ids, but it has
  a special class attribute which makes it possible to perform graph
  specific operations on it, like selecting a subset of the vertices
  based on some vertex attributes.

  A vertex sequence is created by \code{V(g)} this selects are vertices
  in increasing vertex id order. A vertex sequence can be indexed by a
  numeric vector, and a subset of all vertices can be
  selected.

  Vertex sequences provide powerful operations for dealing with vertex
  attributes. A vertex sequence can be indexed with the
  \sQuote{\code{$}} operator to select (or modify) the attributes of a
  subset of vertices. A vertex sequence can be indexed by a logical
  expression, and this expression may contain the names of the vertex
  attributes and ordinary variables as well. The return value of such a
  construct (ie. a vertex sequence indexed by a logical expression) is
  another vertex sequence containing only vertices from the original
  sequence for which the expression evaluates to TRUE.

  Let us see an example to make everything clear. We assign random
  numbers between 1 and 100 to the vertices, and select those vertices
  for which the number is less than 50. We set the color of these
  vertices to red.
  \preformatted{
    g <- graph.ring(10)
    V(g)$number <- sample(1:100, vcount(g), replace=TRUE)
    V(g)$color <- "grey"
    V(g)[ number < 50 ]$color <- "red"
    plot(g, layout=layout.circle, vertex.color=V(g)$color,
         vertex.label=V(g)$number)
  }

  There is a similar notation for edges. \code{E(g)} selects all edges
  from the \sQuote{\code{g}} graph. Edge sequences can be also indexed
  with logical expressions containing edge attributes:
  \preformatted{
    g <- graph.ring(10)
    E(g)$weight <- runif(ecount(g))
    E(g)$width <- 1
    E(g)[ weight >= 0.5 ]$width <- 3
    plot(g, layout=layout.circle, edge.width=E(g)$width, edge.color="black")
  }

  It is important to note that, whenever we use iterators to assign new
  attribute values, the new values are recycled. So in the following
  example half of the vertices will be black, the other half red, in an
  alternated way.
  \preformatted{
    g <- graph.ring(10)
    V(g)$color <- c("black", "red")
    plot(g, layout=layout.circle)
  }

  For the recycling, the standard R rules apply and a warning is given
  if the number of items to replace is not a multiple of the replacement
  length. E.g. the following code gives a warning, because we set the
  attribute for three vertices, but supply only two values:
  \preformatted{
    g <- graph.tree(10)
    V(g)$color <- "grey"
    V(g)[0:2]$color <- c("green", "blue")
  }
  
  If a new vertex/edge attribute is created with an assignment, but only
  a subset of of vertices are specified, then the rest is set to
  \code{NA} if the new values are in a vector and to \code{NULL} if they
  are a list. Try the following:
  \preformatted{
    V(g)[5]$foo <- "foo"
    V(g)$foo
    V(g)[5]$bar <- list(bar="bar")
    V(g)$bar
  }
  
  There are some special functions which are only defined in the
  indexing expressions of vertex and edge sequences. For vertex
  sequences these are: \code{nei}, \code{adj}, \code{from} and
  \code{to}, \code{innei} and \code{outnei}.

  \code{nei} has a mandatory and an optional argument, the first is
  another vertex sequence, the second is a mode argument similar to that
  of the \code{\link{neighbors}} function. \code{nei} returns a logical
  vector of the same length as the indexed vertex sequence and evaluates
  to \code{TRUE} for those vertices only which have a neighbor vertex in
  the vertex sequence supplied as a parameter. Thus for selecting all
  neighbors of vertices 0 and 1 one can write:
  \preformatted{
    V(g) [ nei( 0:1 ) ]
  }
  The mode argument (just like for \code{\link{neighbors}}) gives the
  type of the neighbors to be included, it is interpreted only in
  directed graphs, and defaults to all types of neighbors. See the
  example below. \code{innei(v)} is a shorthand for the \sQuote{incoming}
  neighbors (\code{nei(v, mode="in")}), and \code{outnei(v)} is a
  shorthand for the \sQuote{outgoing} neighbors
  (\code{nei(v,mode="out")}).

  \code{adj} takes an edge sequence as an argument and returns
  \code{TRUE} for vertices which have at least one adjacent edge in it.

  \code{from} and \code{to} are similar to \code{adj} but only edges
  originated at (\code{from}) or pointing to (\code{to}) are taken into
  account.

  For edge sequences the special functions are: \code{adj}, \code{from},
  \code{to}, \code{\%--\%}, \code{\%->\%} and \code{\%<-\%}.

  \code{adj} takes a vertex sequence as an argument and returns
  \code{NULL} for edges which have an adjacent vertex in it.

  \code{from} and \code{to} are similar to \code{adj}, but only vertices
  at the source (\code{from}) or target (\code{to}) of the edge.

  The \code{\%--\%} operator selects edges connecting two vertex
  sequences, the direction of the edges is ignored. The \code{\%->\%} is
  different only for directed graphs and only edges pointing from the
  left hand side argument to the right hand side argument are selected.
  \code{\%<-\%} is exactly the opposite, it selects edges pointing from
  the right hand side to the left hand side.  

  \code{E} has two optional arguments: \code{P} and \code{path}. If
  given \code{P} can be used to select edges based on their end points,
  eg. \code{E(g, P=c(0,1))} selects edge \code{0->1}.

  \code{path} can be used to select all edges along a path. The path
  should be given with the visited vertex ids in the appropriate order.
  
  See also the examples below.  
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
# mean degree of vertices in the largest cluster in a random graph
g <- erdos.renyi.game(100, 2/100)
c <- clusters(g)
vsl <- which(which.max(c$csize)-1==c$membership)-1
mean(degree(g, vsl))

# set the color of these vertices to red, others greens
V(g)$color <- "green"
V(g)[vsl]$color <- "red"
\dontrun{plot(g, vertex.size=3, labels=NA, vertex.color="a:color",
              layout=layout.fruchterman.reingold)}

# the longest geodesic within the largest cluster
long <- numeric()
for (v in vsl) {
  paths <- get.shortest.paths(g, from=v, to=vsl)
  fl <- paths[[ which.max(sapply(paths, length)) ]]
  if (length(fl) > length(long)) {
    long <- fl
  }
}

# the mode argument of the nei() function
g <- graph( c(0,1, 1,2, 1,3, 3,1) )
V(g)[ nei( c(1,3) ) ]
V(g)[ nei( c(1,3), "in") ]
V(g)[ nei( c(1,3), "out") ]

# operators for edge sequences
g <- barabasi.game(100, power=0.3)
E(g) [ 0:2 \%--\% 1:5 ]
E(g) [ 0:2 \%->\% 1:5 ]
E(g) [ 0:2 \%<-\% 1:5 ]

# the edges along the diameter
g <- barabasi.game(100, directed=FALSE)
d <- get.diameter(g)
E(g, path=d)
}
\keyword{graphs}
