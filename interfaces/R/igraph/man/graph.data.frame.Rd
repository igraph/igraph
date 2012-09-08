\name{graph.data.frame}
\alias{graph.data.frame}
\alias{get.data.frame}
\concept{Data frame}
\title{Creating igraph graphs from data frames or vice-versa}
\description{
  This function creates an igraph graph from one or two data frames
  containing the (symbolic) edge list and edge/vertex attributes.
}
\usage{
graph.data.frame(d, directed=TRUE, vertices=NULL)
get.data.frame(x, what=c("edges", "vertices", "both"))
}
\arguments{
  \item{d}{A data frame containing a symbolic edge list in the first two
    columns. Additional columns are considered as edge attributes.}
  \item{directed}{Logical scalar, whether or not to create a directed graph.}
  \item{vertices}{A data frame with vertex metadata, or \code{NULL}. See
    details below.}
  \item{x}{An igraph object.}
  \item{what}{Character constant, whether to return info about vertices,
    edges, or both. The default is \sQuote{edges}.}
}
\details{
  \code{graph.data.frame} creates igraph graphs from one or two data
  frames. It has two modes of operatation, depending whether the
  \code{vertices} argument is \code{NULL} or not.

  If \code{vertices} is \code{NULL}, then the first two columns of
  \code{d} are used as a symbolic edge list and additional columns as
  edge attributes. The names of the attributes are taken from the names
  of the columns.

  If \code{vertices} is not \code{NULL}, then it must be a data frame
  giving vertex metadata. The first column of \code{vertices} is assumed
  to contain symbolic vertex names, this will be added to the graphs as
  the \sQuote{\code{name}} vertex attribute. Other columns will be added
  as additional vertex attributes. If \code{vertices} is not \code{NULL}
  then the symbolic edge list given in \code{d} is checked to contain
  only vertex names listed in \code{vertices}.

  Typically, the data frames are exported from some speadsheat software
  like Excel and are imported into R via \code{\link{read.table}},
  \code{\link{read.delim}} or \code{\link{read.csv}}.

  \code{get.data.frame} converts the igraph graph into one or more data
  frames, depending on the \code{what} argument.

  If the \code{what} argument is \code{edges} (the default), then the
  edges of the graph and also the edge attributes are returned. The
  edges will be in the first two columns, named \code{from} and
  \code{to}. (This also denotes edge direction for directed graphs.)
  For named graphs, the vertex names will be included in these columns,
  for other graphs, the numeric vertex ids. The edge attributes will be
  in the other columns. It is not a good idea to have an edge attribute
  named \code{from} or \code{to}, because then the column named in the
  data frame will not be unique. The edges are listed in the order of
  their numeric ids.

  If the \code{what} argument is \code{vertices}, then vertex attributes
  are returned. Vertices are listed in the order of their numeric vertex
  ids.

  If the \code{what} argument is \code{both}, then both vertex and edge
  data is returned, in a list with named entries \code{vertices} and
  \code{edges}.
}
\note{
  For \code{graph.data.frame} \code{NA} elements in the first two
  columns \sQuote{d} are replaced by the string \dQuote{NA} before
  creating the graph. This means that all \code{NA}s will correspond to
  a single vertex. 

  \code{NA} elements in the first column of \sQuote{vertices} are also
  replaced by the string \dQuote{NA}, but the rest of \sQuote{vertices}
  is not touched. In other words, vertex names (=the first column)
  cannot be \code{NA}, but other vertex attributes can.
}
\value{
  An igraph graph object for \code{graph.data.frame}, and either a data
  frame or a list of two data frames named \code{edges} and
  \code{vertices} for \code{as.data.frame}.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{ \code{\link{graph.constructors}} and
  \code{\link{graph.formula}} for other ways to create graphs,
  \code{\link{read.table}} to read in tables from files.}
\examples{
## A simple example with a couple of actors
## The typical case is that these tables are read in from files....
actors <- data.frame(name=c("Alice", "Bob", "Cecil", "David",
                            "Esmeralda"),
                     age=c(48,33,45,34,21),
                     gender=c("F","M","F","M","F"))
relations <- data.frame(from=c("Bob", "Cecil", "Cecil", "David",
                               "David", "Esmeralda"),
                        to=c("Alice", "Bob", "Alice", "Alice", "Bob", "Alice"),
                        same.dept=c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),
                        friendship=c(4,5,5,2,1,1), advice=c(4,5,5,4,2,3))
g <- graph.data.frame(relations, directed=TRUE, vertices=actors)
print(g, e=TRUE, v=TRUE)

## The opposite operation
get.data.frame(g, what="vertices")
get.data.frame(g, what="edges")
}
\keyword{graphs}
