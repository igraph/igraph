\docType{package}
\name{igraph-package}
\alias{igraph-package}
\alias{igraph}
\title{The igraph package}
\description{igraph is a library for network analysis.}
\section{Introduction}{
  The main goals of the igraph library is 
  to provide a set of data types and functions for 1) pain-free
  implementation of graph algorithms, 2) fast handling of large graphs,
  with millions of vertices and edges, 3) allowing rapid prototyping via
  high level languages like R.
}
\section{Creating graphs}{
  There are many functions in igraph for creating graphs, both
  deterministic and stochastic; stochastic graph constructors are called
  \sQuote{games} in igraph.

  To create small graphs with a given structure probably the
  \code{\link{graph.formula}} function is easiest. It uses R's formula
  interface, its manual page contains many examples. Another option is
  \code{\link{graph}}, which takes numeric vertex ids directly.
  \code{\link{graph.atlas}} creates graph from the Graph Atlas,
  \code{\link{graph.famous}} can create some special graphs.

  To create graphs from field data, \code{\link{graph.edgelist}},
  \code{\link{graph.data.frame}} and \code{\link{graph.adjacency}} are
  probably the best choices.

  The igraph include some classic random graphs like the Erdos-Renyi GNP
  and GNM graphs (\code{\link{erdos.renyi.game}}) and some recent
  popular models, like preferential attachment
  (\code{\link{barabasi.game}}) and the small-world model
  (\code{\link{watts.strogatz.game}}). 
}
\section{Vertex and edge IDs}{
  Vertices and edges have numerical vertex ids in igraph. Vertex ids are
  always consecutive and they start with zero. I.e. for a graph with
  \sQuote{n} vertices the vertex ids are between \sQuote{0} and
  \sQuote{n-1}. If some operation changes the number of vertices in the
  graphs, e.g. a subgraph is created via \code{\link{subgraph}}, then
  the vertices are renumbered to satisfty this criteria.

  The same is true for the edges as well, edge ids are always between
  zero and \sQuote{m-1} where \sQuote{m} is the total number of edges in
  the graph.

  It is often desirable to follow vertices along a number of graph
  operations, and vertex ids don't allow this because of the
  renumbering. The solution is to assign attributes to the
  vertices. These are kept by all operations, if possible. See more
  about attributes in the next section.
}
\section{Attributes}{
  In igraph it is possible to assign attributes to the vertices or edges
  of a graph, or to the graph itself. igraph provides flexible
  constructs for selecting a set of vertices or edges based on their
  attribute values, see \code{\link{get.vertex.attribute}} and
  \code{\link{iterators}} for details.

  Some vertex/edge/graph attributes are treated specially. One of them
  is the \sQuote{name} attribute. This is used for printing the graph
  instead of the numerical ids, if it exists. Vertex names can also be
  used to specify a vector or set of vertices, in all igraph
  functions. E.g. \code{\link{degree}} has a \sQuote{\code{v}} argument
  that gives the vertices for which the degree is calculated. This
  argument can be given as character vector of vertex names as well.

  Other attributes define visualization parameters, see
  \code{\link{igraph.plotting}} for details.

  Attribute values can be set to any R object, but note that storing the
  graph in some file formats might result the loss of complex attribute
  values. All attribute values are preserved if you use
  \code{\link[base]{save}} and \code{\link[base]{load}} to store/retrieve your
  graphs.
}
\section{Visualization}{
  igraph provides three different ways for visualization. The first is
  the \sQuote{\code{plot.igraph}} function. (Actually you don't need to
  write \sQuote{\code{plot.igraph}}, \sQuote{\code{plot}} is
  enough. This function uses base R graphic and can be used with any R
  device.

  The second function is \code{\link{tkplot}}, which uses a Tk GUI for
  basic interactive graph manipulation. (Tk is quite resource hungry, so
  don't try this for very large graphs.)

  The third way requires the \code{rgl} package and uses OpenGL. See the
  \code{\link{rglplot}} function for the details.

  Make sure you read \code{\link{igraph.plotting}} before you start
  plotting your graphs.
}
\section{File formats}{
  igraph can handle various graph file formats, usually both for reading
  and writing. We suggest that you use the GraphML file format for your
  graphs, except if the graphs are too big. For big graphs a simpler
  format is recommended. See \code{\link{read.graph}} and
  \code{\link{write.graph}} for details.
}
\section{Further information}{
  The igraph homepage is at \url{http://cneurocvs.rmki.kfki.hu/igraph}.
  See especially the documentation section. Join the igraph-help mailing
  list if you have questions or comments.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
