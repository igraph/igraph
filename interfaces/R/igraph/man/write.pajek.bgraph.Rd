\name{write.pajek.bgraph}
\alias{write.pajek.bgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Write graphs and their cohesive block hierarchy as Pajek files}
\description{
  Creates a series of Pajek-compatible files from a \code{bgraph} object
  as output by \code{\link{cohesive.blocks}}.
}
\usage{
write.pajek.bgraph(graph, filename, hierarchy = FALSE)
}
\arguments{
  \item{graph}{A \code{bgraph} object as output by
    \code{\link{cohesive.blocks}}.} 
  \item{filename}{The filename \emph{without file extension} to use.}
  \item{hierarchy}{Logical. Should a separate \code{.clu} file be
    created for each block in the graph? (see details)} 
}
\details{
This function writes at least two files to disk. First, the file
\code{filename.net} contains the basic graph structure---vertices and
edges. The second is the file \code{filename.clu}, which clusters the
vertices according to vertex connectivity of their maxmally cohesive
subgraph. This uses the same logic as the default vertex coloring in
\code{\link{plot.bgraph}}.

If \code{hierarchy} is \code{TRUE}, more detailed and complete
information is written to a further series of files.  
For each block a seperate \code{.clu} file is written, each named like
\dQuote{\code{filename_block2(1).clu}}, where the first number refers to the
block number and the number in parenthesis refers to that block's
cohesion. Finally, if \code{hierarchy} is \code{TRUE}, a tree
representing the block hierarchy is written to
\code{filename_blocktree.net}, allowing the subset-structure of the
individual block structures to be determined. 
}
\references{Pajek website:
  \url{http://vlado.fmf.uni-lj.si/pub/networks/pajek/}}
\author{Peter McMahan \email{peter.mcmahan@gmail.com}}
\note{This function has not been thoroughly tested, and may contain bugs}
\seealso{\code{\link{cohesive.blocks}}, \code{\link{plot.bgraph}},
  \code{\link{write.graph}}}
\keyword{graphs}
