\name{graph.intersection}
\alias{graph.intersection}
\alias{\%s\%}
\concept{Graph operators}
\title{Intersection of graphs}
\description{The intersection of two or more graphs are created.
  The graphs may have identical or overlapping vertex sets.}
\usage{
graph.intersection(\dots, byname = "auto", keep.all.vertices = TRUE)
}
\arguments{
  \item{\dots}{Graph objects or lists of graph objects.}
  \item{byname}{A logical scalar, or the character scalar
    \code{auto}. Whether to perform the operation based on symbolic
    vertex names. If it is \code{auto}, that means \code{TRUE} if all
    graphs are named and \code{FALSE} otherwise. A warning is generated
    if \code{auto} and some (but not all) graphs are named.}
  \item{keep.all.vertices}{Logical scalar, whether to keep vertices that
    only appear in a subset of the input graphs.}
}
\details{
  \code{graph.intersection} creates the intersection of two or more
  graphs: only edges present in all graphs will be included.
  The corresponding operator is \%s\%.

  If the \code{byname} argument is \code{TRUE} (or \code{auto} and all
  graphs are named), then the operation is performed on symbolic vertex
  names instead of the internal numeric vertex ids.

  \code{graph.intersection} keeps the attributes of all graphs. All
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
## Common part of two social networks
net1 <- graph.formula(D-A:B:F:G, A-C-F-A, B-E-G-B, A-B, F-G,
                      H-F:G, H-I-J)
net2 <- graph.formula(D-A:F:Y, B-A-X-F-H-Z, F-Y)
str(net1 \%s\% net2)
}
\keyword{graphs}
