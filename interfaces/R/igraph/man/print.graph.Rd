\name{print.igraph}
\alias{print.igraph}
\alias{str.igraph}
\alias{summary.igraph}
\title{Print graphs to the terminal}
\description{These functions attempt to print a graph to the terminal in
  a human readable form.}
\usage{
\method{print}{igraph}(x, full=getIgraphOpt("print.full"),
  graph.attributes=getIgraphOpt("print.graph.attributes"),
  vertex.attributes=getIgraphOpt("print.vertex.attributes"),
  edge.attributes=getIgraphOpt("print.edge.attributes"),
  names=TRUE, \dots)
\method{summary}{igraph}(object, \dots)
\method{str}{igraph}(object, \dots)
}
\arguments{
  \item{x}{The graph to print.}
  \item{full}{Logical scalar, whether to print the graph structure
    itself as well.}
  \item{graph.attributes}{Logical constant, whether to print graph attributes.}
  \item{vertex.attributes}{Logical constant, whether to print vertex attributes.}
  \item{edge.attributes}{Logical constant, whether to print edge
    attributes.}
  \item{names}{Logical constant, whether to print symbolic vertex names
    (ie. the \code{name} vertex attribute) or vertex ids.}
  \item{object}{The graph of which the summary will be printed.}
  \item{\dots}{Additional agruments.}
}
\details{
  \code{summary.igraph} prints the number of vertices, edges and whether the
  graph is directed.

  \code{str.igraph} prints the same information, and also lists
  the edges, and optionally graph, vertex and/or edge attributes.

  \code{print.igraph} behaves either as \code{summary.igraph} or
  \code{str.igraph} depending on the \code{full} argument. See also the
  \sQuote{print.full} igraph option and \code{\link{getIgraphOpt}}.

  The graph summary printed by \code{summary.igraph} (and
  \code{print.igraph} and \code{str.igraph}) consists one or more
  lines. The first line contains the basic properties of the graph, and
  the rest contains its attributes. Here is an example, a small star
  graph with weighed directed edges and named vertices: \preformatted{    IGRAPH DNW- 10 9 -- In-star
    + attr: name (g/c), mode (g/c), center (g/n), name (v/c),
      weight (e/n) }
  The first line always starts with \code{IGRAPH}, showing you that the
  object is an igraph graph. Then a four letter long code string is
  printed. The first letter distinguishes between directed
  (\sQuote{\code{D}}) and undirected (\sQuote{\code{U}}) graphs. The
  second letter is \sQuote{\code{N}} for named graphs, i.e. graphs with
  the \code{name} vertex attribute set. The third letter is
  \sQuote{\code{W}} for weighted graphs, i.e. graphs with the
  \code{weight} edge attribute set. The fourth letter is
  \sQuote{\code{B}} for bipartite graphs, i.e. for graphs with the
  \code{type} vertex attribute set.

  Then, after two dashes, the name of the graph is printed, if it has
  one, i.e. if the \code{name} graph attribute is set.

  From the second line, the attributes of the graph are listed,
  separated by a comma. After the attribute names, the kind of the
  attribute -- graph (\sQuote{\code{g}}), vertex (\sQuote{\code{v}}) or
  edge (\sQuote{\code{e}}) -- is denoted, and the type of the 
  attribute as well, character (\sQuote{\code{c}}), numeric
  (\sQuote{\code{n}}), logical (\sQuote{\code{l}}), or
  other (\sQuote{\code{x}}).  
  
  As of igraph 0.4 \code{str.igraph} and \code{print.igraph} use the
  \code{max.print} option, see \code{\link[base]{options}} for details.
}
\value{All these functions return the graph invisibly.}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
%\seealso{}
\examples{
g <- graph.ring(10)
g
summary(g)
}
\keyword{graphs}
