\name{read.graph}
\alias{read.graph}
\alias{LGL}
\alias{Pajek}
\alias{GraphML}
\alias{GML}
\alias{DL}
\alias{UCINET}
\title{Reading foreign file formats}
\description{The \code{read.graph} function is able to read graphs in
  various representations from a file, or from a http
  connection. Currently some simple formats are supported.}
\usage{
read.graph(file, format = c("edgelist", "pajek", "ncol", "lgl",
        "graphml", "dimacs", "graphdb", "gml", "dl"), \dots)
}
\arguments{
  \item{file}{The connection to read from. This can be a local file, or
    a \code{http} or \code{ftp} connection. It can also be a character
    string  with the file name or URI.}
  \item{format}{Character constant giving the file format. Right now
    \code{edgelist}, \code{pajek}, \code{graphml}, \code{gml},
    \code{ncol}, \code{lgl}, \code{dimacs} and \code{graphdb}
    are supported, the default is \code{edgelist}. As of igraph 0.4 this
    argument is case insensitive.
  }
  \item{\dots}{Additional arguments, see below.}
}
\details{
  The \code{read.graph} function may have additional arguments depending
  on the file format (the \code{format} argument). See the details
  separately for each file format, below.
}
\section{Edge list format}{
      This format is a simple text file with numeric vertex ids defining
      the edges. There is no need to have newline characters between the
      edges, a simple space will also do.
      
      Additional arguments:
      \describe{
	\item{n}{The number of vertices in the graph. If it is smaller
	  than or equal to the largest integer in the file, then it is
	  ignored; so it is safe to set it to zero (the default).}
	\item{directed}{Logical scalar, whether to create a directed
	  graph. The default value is \code{TRUE}.}
      }
}
\section{Pajek format}{
      Pajek it a popular network analysis program for Windows. (See the
      Pajek homepage at
      \url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}.)
      It has a quite flexible but not very well
      documented file format, see the Pajek manual on the Pajek homepage
      for some information about the file format.

      \code{igraph} implements only a subset of the Pajek format:
      \itemize{
	\item Only .net files are supported, Pajek project files (which
	  can contain many graph and also other type of data) are
	  not. Poject files might be supported in a forthcoming igraph
	  release if they turned out to be needed.
	\item Time events networks are not supported.
	\item Hypergraphs (graphs with non-binary edges) are not
	  supported as igraph cannot handle them.
	\item Graphs containing both directed and undirected edges are
	  not supported as igraph cannot represent them.
	\item Graph with multiple edge sets are not supported.
      }

      From version 0.6.1 igraph supports reading bipartite (two-mode)
      graphs from Pajek files and adds the \code{type} vertex
      attribute. A warning is given if invalid edges (edges connecting
      vertices of the same type) are present in the file.
      
      Vertex and edge attributes defined in the Pajek file will be also
      read and assigned to the graph object to be created. These are
      mainly parameters for graph visualization, but not exclusively,
      eg. the file might contain edge weights as well.

      The following vertex attributes might be added:
      \tabular{ll}{
	igraph name \tab description, Pajek attribute \cr
	id 	    \tab Vertex id \cr
	x, y, z     \tab The \sQuote{x}, \sQuote{y} and \sQuote{z} coordinate of the vertex \cr
	vertexsize  \tab The size of the vertex when plotted (\code{size} in Pajek). \cr
	shape 	    \tab The shape of the vertex when plotted. \cr
	color 	    \tab Vertex color (\code{ic} in Pajek) if given with symbolic name \cr
	framecolor  \tab Border color (\code{bc} in Pajek) if given with symbolic name \cr
	labelcolor  \tab Label color (\code{lc} in Pajek) if given with symbolic name \cr
	xfact, yfact \tab The \code{x_fact} and \code{y_fact} Pajek attributes. \cr
	labeldist \tab The distance of the label from the vertex. (\code{lr} in Pajek.) \cr
	labeldegree, \tab \cr
	labeldegree2 \tab The \code{la} and \code{lphi} Pajek attributes \cr
	framewidth \tab The width of the border (\code{bw} in Pajek). \cr
	fontsize  \tab Size of the label font (\code{fos} in Pajek.) \cr
	rotation \tab The rotation of the vertex (\code{phi} in Pajek). \cr
	radius \tab Radius, for some vertex shapes (\code{r} in Pajek). \cr
	diamondratio \tab For the diamond shape (\code{q} in Pajek). \cr
	type \tab vertex types in bipartite (two-mode) graphs. \cr
      }

      These igraph attributes are only created if there is at least one vertex
      in the Pajek file which has the corresponding associated
      information. Eg. if there are vertex coordinates for at least one vertex
      then the \sQuote{x}, \sQuote{y} and possibly also \sQuote{z}
      vertex attributes will be created. For those vertices for which the
      attribute is not defined, \code{NaN} is assigned. 
      
      The following edge attributes might be added:
      \tabular{ll}{
	igraph name \tab description, Pajek attribute \cr
	weight 	\tab Edge weights. \cr
	label 	\tab \code{l} in Pajek. \cr
	color 	\tab Edge color, if the color is given with a symbolic name, \code{c} in Pajek. \cr
	color-red, \tab \cr
	color-green, \tab \cr
	color-blue \tab Edge color if it was given in RGB notation, \code{c} in Pajek. \cr
	edgewidth \tab \code{w} in Pajek. \cr
	arrowsize \tab \code{s} in Pajek. \cr
	hook1, hook2 \tab \code{h1} and \code{h2} in Pajek. \cr
	angle1, angle2 \tab \code{a1} and \code{a2} in Pajek, Bezier curve parameters. \cr
	velocity1, \tab \cr
	velocity2 \tab \code{k1} and \code{k2} in Pajek, Bezier curve parameter. \cr
	arrowpos \tab \code{ap} in Pajek. \cr
	labelpos \tab \code{lp} in Pajek. \cr
	labelangle, \tab \cr
	labelangle2 \tab \code{lr} and \code{lphi} in Pajek. \cr
	labeldegree \tab \code{la} in Pajek. \cr
	fontsize \tab \code{fos} in Pajek. \cr
	arrowtype \tab \code{a} in Pajek. \cr
	linepattern \tab \code{p} in Pajek. \cr
	labelcolor \tab \code{lc} in Pajek. \cr
      }
      
      There are no additional arguments for this format.
}
\section{GraphML file format}{
      GraphML is an XML-based file format (an XML application in the XML
      terminology) to describe graphs. It is a modern format, and can
      store graphs with an extensible set of vertex and edge attributes,
      and generalized graphs which igraph cannot handle. Thus igraph
      supports only a subset of the GraphML language: 
      \itemize{
	\item Hypergraphs are not supported.
	\item Nested graphs are not supported.
	\item Mixed graphs, ie. graphs with both directed and
	  undirected edges are not supported. read.graph() sets the
	  graph directed if this is the default in the GraphML file,
	  even if all the edges are in fact undirected.
      }

      See the GraphML homepage at
      \url{http://graphml.graphdrawing.org} for more information about
      the GraphML format.       
      
      Additional arguments:
      \describe{
	\item{index}{If the GraphML file contains more than one graphs,
	  this argument can be used to select the graph to read. By
	  default the first graph is read (index 0).}
      }
}
\section{GML file format}{
      GML is a simple textual format, see
      \url{http://www.infosun.fim.uni-passau.de/Graphlet/GML/} for details.
 
      Although all syntactically correct GML can be parsed, 
      we implement only a subset of this format, some attributes might be
      ignored. Here is a list of all the differences:
      \itemize{
	\item Only \code{node} and \code{edge} attributes are 
	  used, and only if they have a simple type: integer, real or
	  string. So if an attribute is an array or a record, then it is 
	  ignored. This is also true if only some values of the
	  attribute are complex.
	\item Top level attributes except for \code{Version} and the
	  first \code{graph} attribute are completely ignored.
	\item Graph attributes except for \code{node} and
	  \code{edge} are completely ignored.
	\item There is no maximum line length.
	\item There is no maximum keyword length.
	\item Character entities in strings are not interpreted.
	\item We allow \code{inf} (infinity) and \code{nan}
	  (not a number) as a real number. This is case insensitive, so
	  \code{nan}, \code{NaN} and \code{NAN} are equal.
      }
      
      Please contact us if you cannot live with these
      limitations of the GML parser.

      There are not additional argument for this format.
}
\section{DL file format}{
      The DL format is a simple textual file format used by the UCINET
      software. See
      \url{http://www.analytictech.com/networks/dataentry.htm} for
      examples. All formats mentioned here is supported by igraph.

      Note the specification does not mention whether the
      format is case sensitive or not. For igraph DL files are case
      sensitive, i.e. \sQuote{Larry} and \sQuote{larry} are not the
      same.
      
      Additional arguments:
      \describe{
	\item{directed}{Logical scalar, whether to create a directed
	  graph. The default is to make the graph directed.}
      }
}
\section{NCOL format}{
      This format is used by the Large Graph Layout program
      (\url{http://bioinformatics.icmb.utexas.edu/lgl}), and it is
      simply a symbolic weighted edge list. It is a simple text file
      with one edge per line. An edge is defined by two symbolic vertex
      names separated by whitespace. (The symbolic vertex names
      themselves cannot contain whitespace.) They might followed by an
      optional number, this will be the weight of the edge; the number
      can be negative and can be in scientific notation. If there is no
      weight specified to an edge it is assumed to be zero. 

      The resulting graph is always undirected. LGL cannot deal with
      files which contain multiple or loop edges, this is however not
      checked here, as igraph is happy with these.

      Additional arguments:
      \describe{
	\item{names}{Logical constant, whether to add the symbolic names
	  as vertex attributes to the graph. If TRUE the name of the
	  vertex attribute will be \sQuote{name}.}
	\item{weights}{Character scalar, specifies whether edge weights
	  should be added to the graph. Possible values are and their
	  meaning are: \sQuote{no}, edge weights will not be added;
	  \sQuote{yes}, edge weights will be added, if they are not
	  present in the file, then all edges get zero weight;
	  \sQuote{auto}, edge weights will added if they are present in
	  the file, otherwise not. The default is \sQuote{auto}.
	}
	\item{directed}{Logical constant, whether to create a directed
	  graph. The default is undirected.} 
      }
}
\section{LGL file format}{
      The \code{lgl} format is used by the Large Graph Layout
      visualization software
      (\url{http://bioinformatics.icmb.utexas.edu/lgl}), it can describe
      undirected optionally weighted graphs. From the LGL manual:
      \dQuote{The second format is the LGL file format (.lgl file
      suffix). This is yet another graph file format that tries to be
      as stingy as possible with space, yet keeping the edge file in a
      human readable (not binary) format. The format itself is like
      the following: \preformatted{	 # vertex1name 
	vertex2name [optionalWeight]
	vertex3name [optionalWeight]
      } Here, the first vertex of an edge is preceded with a pound sign
      '\#'.  Then each vertex that shares an edge with that vertex is
      listed one per line on subsequent lines.}
      
      LGL cannot handle loop and multiple edges or directed graphs, but
      in igraph it is not an error to have multiple and loop edges.
      
      Additional arguments:
      \describe{
	\item{names}{Logical constant, whether to add the symbolic names
	  as vertex attributes to the graph. If TRUE the name of the
	  vertex attribute will be \sQuote{name}.}
	\item{weights}{Character scalar, specifies whether edge weights
	  should be added to the graph. Possible values are and their
	  meaning are: \sQuote{no}, edge weights will not be added;
	  \sQuote{yes}, edge weights will be added, if they are not
	  present in the file, then all edges get zero weight;
	  \sQuote{auto}, edge weights will added if they are present in
	  the file, otherwise not. The default is \sQuote{auto}.
	}
      }      
}
\section{DIMACS file format}{
      The DIMACS file format, more specifically the 
      version for network flow problems, see the files at
      \url{ftp://dimacs.rutgers.edu/pub/netflow/general-info/}

      This is a line-oriented text file (ASCII) format. The first
      character of each line defines the type of the line. If the first
      character is \code{c} the line is a comment line and it is
      ignored. There is one problem line (\code{p}) in the file, it
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

      The source vertex is assigned to the \code{source}, the target
      vertex to the \code{target} graph attribute. The edge capacities
      are assigned to the \code{capacity} edge attribute.

      Additional arguments:
      \describe{
	\item{directed}{Logical scalar, whether to create a directed
	  graph. By default a directed graph is created.
	}
      }
}
\section{GraphDB format}{
      This is a binary format, used in the graph database
      for isomorphism testing (\url{http://amalfi.dis.unina.it/graph/})
      From the graph database homepage
      (\url{http://amalfi.dis.unina.it/graph/db/doc/graphdbat-2.html}):

      \emph{The graphs are stored in a compact binary format, one graph per
	file. The file is composed of 16 bit words, which are represented
	using the so-called little-endian convention, i.e. the least
	significant byte of the word is stored first.}
      
      \emph{Then, for each node, the file contains the list of edges coming
	out of the node itself. The list is represented by a word encoding
	its length, followed by a word for each edge, representing the
	destination node of the edge. Node numeration is 0-based, so the
	first node of the graph has index 0.}

      See also \code{\link{graph.graphdb}}.
      
      Only unlabelled graphs are implemented.

      Additional attributes:
      \describe{
	\item{directed}{Logical scalar. Whether to create a directed
	  graph.}
      }
}
\value{A graph object.}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{ \code{\link{write.graph}} }
% \examples{}
\keyword{graphs}
