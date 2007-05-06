\name{write.graph}
\alias{write.graph}
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
    \code{pajek}, \code{graphml}, \code{edgelist}, \code{lgl} and
    \code{ncol} are implemented.} 
  \item{\dots}{Other, format specific arguments, see below.}
}
\details{
  The \code{edgelist} format is a simple text file, with one edge in a
  line, the two vertex ids separated by a space character. The file is
  sorted by the first and the second column.

  The Pajek format is a text file, see \code{\link{read.graph}} for
  details. Appropriate vertex and edge attributes are also written to
  the file. This format has no additional arguments.

  The GraphML format is a flexible XML based format. See
  \code{\link{read.graph}} for GraphML details. Vertex and edge
  attributes are also written to the file. This format has no additional
  arguments.
  
  The \code{lgl} format is also a simple text file, this is the
  format expected by the 'Large Graph Layout' layout generator software.
  See \link{read.graph} for details.
  Additional arguments:
  \itemize{
    \item{names}{If you want to write symbolic vertex names instead of
      vertex ids, supply the name of the vertex attribute containing the
      symbolic names here. By default the \sQuote{name} attribute is
      used if there is one. Supply \code{NULL} if you want to use numeric 
      vertex ids even if there is a \sQuote{name} vertex attribute.} 
    \item{weights}{If you want to write edge weights to the file, supply
      the name of the edge attribute here. By defaults the vertex
      attribute \sQuote{weights} are used if they are installed.
      Supply \code{NULL} here if you want to omit the weights.}
    \item{isolates}{Logical, if \code{TRUE} the isolate vertices are
      also written to the file, they are omitted by default.}
  }

  The \code{ncol} format is also used by LGL, it is a text file, see
  \link{read.graph} for details.
  Additional arguments:
  \itemize{
    \item{names}{If you want to write symbolic vertex names instead of
      vertex ids, supply the name of the vertex attribute containing the
      symbolic names here. By default the \sQuote{name} attribute is
      used if there is one. Supply \code{NULL} if you want to use numeric 
      vertex ids even if there is a \sQuote{name} vertex attribute.} 
    \item{weights}{If you want to write edge weights to the file, supply
      the name of the edge attribute here. By defaults the vertex
      attribute \sQuote{weights} are used if they are installed.
      Supply \code{NULL} here if you want to omit the weights.}
  }

  The \code{dimacs} file format, more specifically the 
  version for network flow problems, see the files at
  \url{ftp://dimacs.rutgers.edu/pub/netflow/general-info/}
  
  This is a line-oriented text file (ASCII) format. The first
  character of each line defines the type of the line. If the first
  character is \code{c} the line is a comment line and it is
  ignored. There is one problem line (\code{p} in the file, it
  must appear before any node and arc descriptor lines. The problem
  line has three fields separated by spaces: the problem type
  (\code{min}, \code{max} or \code{asn}), the
  number of vertices and number of edges in the graph.
  Exactly two node identification lines are expected
  (\code{n}), one for the source, one for the target vertex.
  These have two fields: the id of the vertex and the type of the
  vertex, either \code{s} (=source) or \code{t}
  (=target). Arc lines start with \code{a} and have three
  fields: the source vertex, the target vertex and the edge capacity.
  
  Vertex ids are numbered from 1.
  
  Additional arguments:
  \itemize{
    \item{source}{The id of the source vertex, if \code{NULL} (the
      default) then it is taken from the \code{source} graph attribute.}
    \item{target}{The id of the target vertex, if \code{NULL} (the
      default) then it is taken from the \code{target} graph attribute.}
    \item{capacity}{A numeric vector giving the edge capacities. If
      \code{NULL} (the default) then it is taken from the
      \code{capacity} edge attribute.}
  }
}
\value{A NULL, invisibly.}
\references{Adai AT, Date SV, Wieland S, Marcotte EM. LGL: creating a
  map of protein function with an algorithm for visualizing very large
  biological networks. \emph{J Mol Biol.} 2004 Jun 25;340(1):179-90.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{read.graph}} }
\examples{
g <- graph.ring(10)
\dontrun{write.graph(g, "/tmp/g.txt", "edgelist")}
}
\keyword{graphs}
