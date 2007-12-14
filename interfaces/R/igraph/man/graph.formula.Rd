\name{graph.formula}
\alias{graph.formula}
\title{Creating (small) graphs via a simple interface}
\description{
  This function is useful if you want to create a small (named) graph
  quickly, it works for both directed and undirected graphs.
}
\usage{
graph.formula(...)
}
\arguments{
  \item{...}{The formulae giving the structure of the graph, see details
    below.}
}
\details{
  \code{graph.formula} is very handy for creating small graphs quickly.
  You need to supply one or more R expressions giving the structure of
  the graph. The expressions consist of vertex names and edge
  operators. An edge operator is a sequence of \sQuote{\code{-}} and
  \sQuote{\code{+}} characters, the former is for the edges and the
  latter is used for arrow heads. The edges can be arbitrarily long,
  ie. you may use as many \sQuote{\code{-}} characters to \dQuote{draw}
  them as you like.

  If all edge operators consist of only \sQuote{\code{-}} characters
  then the graph will be undirected, whereas a single \sQuote{\code{+}}
  character implies a directed graph.

  Let us see some simple examples. Without arguments the function
  creates an empty graph:
  \preformatted{
    graph.formula()
  }

  A simple undirected graph with two vertices called \sQuote{A} and
  \sQuote{B} and one edge only:
  \preformatted{
    graph.formula(A-B)
  }

  Remember that the length of the edges does not matter, so we could
  have written the following, this creates the same graph:
  \preformatted{
    graph.formula( A-----B )
  }

  If you have many disconnected components in the graph, separate them
  with commas. You can also give isolate vertices.
  \preformatted{
    graph.formula( A--B, C--D, E--F, G--H, I, J, K )
  }

  The \sQuote{\code{:}} operator can be used to define vertex sets. If
  an edge operator connects two vertex sets then every edge from the
  first set will be connected to every edge in the second set. The
  following form creates a full graph, including loop edges:
  \preformatted{
    graph.formula( A:B:C:D -- A:B:C:D )
  }

  In directed graphs, edges will be created only if the edge operator
  includes a arrow head (\sQuote{+}) \emph{at the end} of the edge:
  \preformatted{
    graph.formula( A -+ B -+ C )
    graph.formula( A +- B -+ C )
    graph.formula( A +- B -- C )
  }
  Thus in the third example no edge is created between vertices \code{B}
  and \code{C}.

  Mutual edges can be also created with a simple edge operator:
  \preformatted{
    graph.formula( A +-+ B +---+ C ++ D + E)
  }
  Note again that the length of the edge operators is arbitrary,
  \sQuote{\code{+}}, \sQuote{\code{++}} and \sQuote{\code{+-----+}} have
  exactly the same meaning.

  If the vertex names include spaces or other special characters then
  you need to quote them:
  \preformatted{
    graph.formula( "this is" +- "a silly" -+ "graph here" )
  }
  You can include any character in the vertex names this way, even
  \sQuote{+} and \sQuote{-} characters.

  See more examples below.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph}} for more general graph creation methods.}
\examples{
# A simple undirected graph
g <- graph.formula( Alice-Bob-Cecil-Alice, Daniel-Cecil-Eugene, Cecil-Gordon )
g

# Another undirected graph, ":" notation
g2 <- graph.formula( Alice-Bob:Cecil:Daniel, Cecil:Daniel-Eugene:Gordon )
g2

# A directed graph
g3 <- graph.formula( Alice +-+ Bob --+ Cecil +-- Daniel,
                     Eugene --+ Gordon:Helen )
g3

# A graph with isolate vertices
g4 <- graph.formula( Alice -- Bob -- Daniel, Cecil:Gordon, Helen )
g4
V(g4)$name

# "Arrows" can be arbitrarily long
g5 <- graph.formula( Alice +---------+ Bob )
g5

# Special vertex names
g6 <- graph.formula( "+" -- "-", "*" -- "/", "\%\%" -- "\%/\%" )
g6
}
\keyword{graphs}
