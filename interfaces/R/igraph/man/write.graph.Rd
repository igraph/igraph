\name{write.graph}
\alias{write.graph}
\alias{write.graph.edgelist}
\alias{write.graph.lgl}
\alias{write.graph.ncol}
\title{Writing the graph to a file in some format}
\description{\code{write.graph} is a general function for exporting
  graphs to foreign file formats, however not many formats are
  implemented right now.}
\usage{
write.graph(graph, file, format="edgelist", \dots)
}
\arguments{
  \item{graph}{The graph to export.}
  \item{file}{A connection or a string giving the file name to write the
    graph to.}
  \item{format}{Character string giving the file format. Right now 
    \code{edgelist}, \code{lgl} and \code{ncol} are implemented.}
  \item{\dots}{Other, format specific arguments, see below.}
}
\details{
  The \code{edgelist} format is a simple text file, with one edge in a
  line, the two vertex ids separated by a space character. The file is
  sorted by the first and the second column.

  The \code{lgl} format is also a simple text file, this is the
  format expected by the 'Large Graph Layout' layout generator software.
  See \link{read.graph.lgl} for details.
  Additional arguments:
  \itemize{
    \item{names}{If you want to write symbolic vertex names instead of
      vertex ids, supply the name of the vertex attribute containing the
      symbolic names here. By default the \sQuote{name} attribute is
      used if there is one. Supply \code{NULL} if you want to numeric use
      vertex ids.} 
    \item{weights}{If you want to write edge weights to the file, supply
      the name of the edge attribute here. Give \code{NULL} here if you
      want to omit the weights.}
    \item{isolates}{Logical, if \code{TRUE} the isolate vertices are
      also written to the file, they are omitted by default.}
  }

  The \code{ncol} format is also used by LGL, it is a text file, see
  \link{read.graph.ncol} for details.
  Additional arguments:
  \itemize{
    \item{names}{If you want to write symbolic vertex names instead of
      vertex ids, supply the name of the vertex attribute containing the
      symbolic names here. By default the \sQuote{name} attribute is
      used if there is one. Supply \code{NULL} if you want to numeric use
      vertex ids.} 
    \item{weights}{If you want to write edge weights to the file, supply
      the name of the edge attribute here. Give \code{NULL} here if you
      want to omit the weights.}
  }
} 
\value{A NULL, invisibly.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{read.graph}} }
\examples{
g <- graph.ring(10)
\dontrun{write.graph(g, "/tmp/g.txt", "edgelist")}
}
\keyword{graphs}
