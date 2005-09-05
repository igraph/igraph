\name{read.graph}
\alias{read.graph}
\alias{read.graph.pajek}
\alias{read.graph.edgelist}
\title{Reading foreign file formats}
\description{The \code{read.graph} function is able to read graphs in
  various representations from a file, or from a http
  connection. Currently the Pajek and edgelist forms are supported, the
  default is being the first.}
\usage{
read.graph(file, format = "pajek", \dots)
}
\arguments{
  \item{file}{The connection to read from. This can be a local file, or
    a \code{http} or \code{ftp} URI.}
  \item{format}{Character constant giving the file format. Right now
    \dQuote{pajek} and \dQuote{edgelist} are supported. }
  \item{\dots}{Additional arguments, these are used as graph attributes.}
}
\details{
  The \code{read.graph} function may have additional arguments depending
  on the file format (the \code{format} argument).
  \itemize{
    \item{pajek}{
      The pajek format is used by the popular network analysis program
      with the same name, see references. Right now igraph supports most
      of the Pajek format, and even if it does not support some feature
      like vertex shapes, it silently ignores that piece of information,
      and does not fail. 
      
      Additional arguments:
      \itemize{
	\item{attributes}{ A logical constant, giving whether the vertex
	  attributes are also to be read from the file. Defaults to TRUE.}
      }
    }
    \item{edgelist}{
      
      This format is a simple text file with numeric vertex ids defining
      the edges. There is no need to have newline characters between the
      edges, a simple space will also do.
      
      There are no additional arguments for this format.
    }
  }
  
  Other additional arguments will be used as graph attributes.
}
\value{A graph object.}
\references{W. de Nooy, A. Mrvar, V. Batagelj:
  Exploratory Social Network Analysis with Pajek, Structural Analysis in
  the Social Sciences 27, Cambridge University Press, 2005. ISBN:0521602629.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{write.graph}} }
\examples{
\dontrun{g <- read.graph(format="pajek",
         "http://vlado.fmf.uni-lj.si/pub/networks/data/Erdos/Erdos02.net")
diameter(g)
}
}
\keyword{graphs}
