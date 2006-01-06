\name{read.graph}
\alias{read.graph}
\alias{read.graph.edgelist}
\alias{read.graph.ncol}
\alias{read.graph.lgl}
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
    \dQuote{edgelist}, \dQuote{ncol} and \dQuote{lgl} are supported. }
  \item{\dots}{Additional arguments, see below.}
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
    \item{ncol}{
      This format is used by the Large Graph Layout program
      (\url{http://bioinformatics.icmb.utexas.edu/lgl}), and it is
      simply a symbolic weighted edge list. It is a simple text file
      with one edge per line. An edge is defined by two symbolic vertex
      names separated by whitespace. (The symbolic vertex names
      themselves cannot contain whitespace. They might follow by an
      optional number, this will be the weight of the edge; the number
      can be negative and can be in scientific notation. If there is no
      weight specified to an edge it is assumed to be zero. 

      The resulting graph is always undirected. LGL cannot deal with
      files which contain multiple or loop edges, this is however not
      checked here, as \a igraph is happy with these.

      Additional arguments:
      \itemize{
	\item{names}{Logical constant, whether to add the symbolic names
	  as vertex attributes to the graph. If TRUE the name of the
	  vertex attribute will be \sQuote{name}.}
	\item{weights}{Logical constant, whether to add the weights of
	  the edges as edge attribute \sQuote{weight}.
	}
      }
    }
    \item{lgl}{The \code{lgl} format is used by the Large Graph Layout
      visualization software
      (\url{http://bioinformatics.icmb.utexas.edu/lgl}), it can describe
      undirected optionally weighted graphs. From the LGL manual:
      \dQuote{The second format is the LGL file format (.lgl file
      suffix). This is yet another graph file format that tries to be
      as stingy as possible with space, yet keeping the edge file in a
      human readable (not binary) format. The format itself is like
      the following:
      \preformatted{
	# vertex1name 
	vertex2name [optionalWeight]
	vertex3name [optionalWeight]
      }
      Here, the first vertex of an edge is preceded with a pound sign
      '\#'.  Then each vertex that shares an edge with that vertex is
      listed one per line on subsequent lines.}
      
      LGL cannot handle loop and multiple edges or directed graphs, but
      in \a igraph it is not an error to have multiple and loop edges.
      
      Additional arguments:
      \itemize{
	\item{names}{Logical constant, whether to add the symbolic names
	  as vertex attributes to the graph. If TRUE the name of the
	  vertex attribute will be \sQuote{name}.}
	\item{weights}{Logical constant, whether to add the weights of
	  the edges as edge attribute \sQuote{weight}.
	}
      }      
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
