\name{read.graph}
\alias{read.graph}
\alias{read.graph.edgelist}
\title{Reading foreign file formats}
\description{The \code{read.graph} function is able to read graphs in
  various representations from a file, or from a http
  connection. Currently the Pajek and edgelist forms are supported, the
  default is being the first.}
\usage{
read.graph(file, format = "edgelist", \dots)
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
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{write.graph}} }
% \examples{}
\keyword{graphs}
