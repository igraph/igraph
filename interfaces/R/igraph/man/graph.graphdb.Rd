\name{graph.graphdb}
\alias{graph.graphdb}
\concept{Graph database}
\title{Load a graph from the graph database for testing graph
  isomorphism.}
\description{This function downloads a graph from a database created for
  the evaluation of graph isomorphism testing algothitms.}
\usage{
graph.graphdb (url = NULL, prefix = "iso", type = "r001", nodes = NULL, 
    pair = "A", which = 0, base = "http://cneurocvs.rmki.kfki.hu/graphdb/gzip", 
    compressed = TRUE, directed = TRUE) 
}
\arguments{
  \item{url}{If not \code{NULL} it is a complete URL with the file to import.}
  \item{prefix}{Gives the prefix. See details below. Possible values:
    \code{iso}, \code{i2}, \code{si4}, \code{si6}, \code{mcs10},
    \code{mcs30}, \code{mcs50}, \code{mcs70}, \code{mcs90}.
  }
  \item{type}{Gives the graph type identifier. See details
    below. Possible values: \code{r001}, \code{r005}, \code{r01},
    \code{r02}, \code{m2D}, \code{m2Dr2}, \code{m2Dr4}, \code{m2Dr6}
    \code{m3D}, \code{m3Dr2}, \code{m3Dr4}, \code{m3Dr6}, \code{m4D},
    \code{m4Dr2}, \code{m4Dr4}, \code{m4Dr6}, \code{b03}, 
    \code{b03m}, \code{b06}, \code{b06m}, \code{b09}, \code{b09m}.
  }
  \item{nodes}{The number of vertices in the graph.}
  \item{pair}{Specifies which graph of the pair to read. Possible
    values: \code{A} and \code{B}.}
  \item{which}{Gives the number of the graph to read. For every graph
    type there are a number of actual graphs in the database. This
    argument specifies which one to read.}
  \item{base}{The base address of the database. See details below.}
  \item{compressed}{Logical constant, if TRUE than the file is expected
    to be compressed by gzip. If \code{url} is \code{NULL} then a
    \sQuote{\code{.gz}} suffix is added to the filename.}
  \item{directed}{Logical constant, whether to create a directed graph.}
}
\details{
  \code{graph.graphdb} reads a graph from the graph database from an FTP
  or HTTP server or from a local copy. It has two modes of operation:

  If the \code{url} argument is specified then it should the complete
  path to a local or remote graph database file. In this case
  we simply call \code{\link{read.graph}} with the proper arguments to
  read the file.

  If \code{url} is \code{NULL}, and this is the default, then the
  filename is assembled from the \code{base}, \code{prefix}, \code{type},
  \code{nodes}, \code{pair} and \code{which} arguments.
  
  See the documentation for the graph database at
  \url{http://amalfi.dis.unina.it/graph/db/doc/graphdbat.html} for
  the actual format of a graph database file and other information.
}
\value{
  A new graph object.
}
\references{M. De Santo, P. Foggia, C. Sansone, M. Vento:
  A large database of graphs and its use for benchmarking graph
  isomorphism algorithms, \emph{Pattern Recognition Letters},
  Volume 24, Issue 8 (May 2003)
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{read.graph}}, \code{\link{graph.isomorphic.vf2}}}
\examples{
\dontrun{
g <- graph.graphdb(prefix="iso", type="r001", nodes=20, pair="A",
  which=10, compressed=TRUE)
g2 <- graph.graphdb(prefix="iso", type="r001", nodes=20, pair="B",
  which=10, compressed=TRUE)
graph.isomorphic.vf2(g, g2)	% should be TRUE
g3 <- graph.graphdb(url="http://cneurocvs.rmki.kfki.hu/graphdb/gzip/iso/bvg/b06m/iso_b06m_m200.A09.gz")
}}
\keyword{graphs}
