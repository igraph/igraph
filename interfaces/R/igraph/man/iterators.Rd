\name{iterators}
\alias{V}
\alias{E}
\title{Vertex and edge sequences and iterators}
\description{Vertex and edge sequences are central concepts of
  igraph.}
\usage{
V(graph)
E(graph)
}
\arguments{
  \item{graph}{A graph object.}
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
  numberic vector, this way only a subset of all vertices can be
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
    plot(g, layout=layout.circle, vertex.color="a:color")
  }

  There is a similar notation for edges. \code{E(g)} selects all edges
  from the \sQuote{\code{g}} graph. Edge sequences can be also indexed
  with logical expressions containing edge attributes:
  \preformatted{
    g <- graph.ring(10)
    E(g)$weight <- runif(ecount(g))
    E(g)$width <- 1
    E(g)[ weight >= 0.5 ]$width <- 2
    plot(g, layout=layout.circle, edge.width="a:width", edge.color="black")
  }

  Edge sequence indexing is a bit different than vertex sequence
  indexing because it specifies two special variable and three special
  operators to be used in the logical indexing expression.
  The two variables are called \sQuote{\code{from}} and
  \sQuote{\code{to}} and they contain the vertex ids of the two end
  points of the edge. For undirected edges \sQuote{\code{from}} is the
  source vertex and \sQuote{\code{to}} is the target vertex; they are
  the same for undirected graphs. Ie. to select all edges adjacent to a
  given vertex you can write
  \preformatted{
    g <- barabasi.game(100)
    E(g)[ from==2 | to==2 ]    
  }

  The special operators preform similar function but they allow a more
  convenient notation. \sQuote{\code{\%->\%}} is a binary operator and
  evaluates to   TRUE for all edges pointing from vertices in the first
  argument to those in the second.

  The arguments of this operator can be numberic vectors or vertex
  sequences. 
  
  There is also a
  \sQuote{\code{\%<-\%}} operator for completeness. This is the same as
  the previous one with the arguments exchanged.

  The undirected version of the \dQuote{connects} operator is
  \sQuote{\code{\%--\%}}. In directed graphs this selects edges pointing
  from the first argument to the second or backwards.

  For undirected graphs \sQuote{\code{\%->\%}}, \sQuote{\code{\%<-\%}}
  and \sQuote{\code{\%--\%}} are the same.

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
}
\keyword{graphs}
