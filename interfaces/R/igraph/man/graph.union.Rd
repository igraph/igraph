\name{graph.union}
\alias{graph.union}
\alias{\%u\%}
\concept{Graph operators}
\title{Union of graphs}
\description{The union of two or more graphs are created. The graphs may
  have identical or overlapping vertex sets.}
\usage{
graph.union(\dots, byname = "auto")
}
\arguments{
  \item{\dots}{Graph objects or lists of graph objects.}
  \item{byname}{A logical scalar, or the character scalar
    \code{auto}. Whether to perform the operation based on symbolic
    vertex names. If it is \code{auto}, that means \code{TRUE} if all
    graphs are named and \code{FALSE} otherwise. A warning is generated
    if \code{auto} and some (but not all) graphs are named.}
}
\details{
  \code{graph.union} creates the union of two or more graphs.
  Edges which are included in at least one graph will be part of the new
  graph. This function can be also used via the \%u\% operator.

  If the \code{byname} argument is \code{TRUE} (or \code{auto} and all
  graphs are named), then the operation is performed on symbolic vertex
  names instead of the internal numeric vertex ids.

  \code{graph.union} keeps the attributes of all graphs. All
  graph, vertex and edge attributes are copied to the result. If an
  attribute is present in multiple graphs and would result a name clash,
  then this attribute is renamed by adding suffixes: _1, _2, etc.

  The \code{name} vertex attribute is treated specially if the operation
  is performed based on symbolic vertex names. In this case \code{name}
  must be present in all graphs, and it is not renamed in the result
  graph.

  An error is generated if some input graphs are directed and others are
  undirected.  
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
## Union of two social networks with overlapping sets of actors
net1 <- graph.formula(D-A:B:F:G, A-C-F-A, B-E-G-B, A-B, F-G,
                      H-F:G, H-I-J)
net2 <- graph.formula(D-A:F:Y, B-A-X-F-H-Z, F-Y)
str(net1 \%u\% net2)
}
\keyword{graphs}
