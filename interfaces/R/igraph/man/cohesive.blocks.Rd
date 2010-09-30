\name{cohesive.blocks}
\alias{cohesive.blocks}
\alias{cohesiveBlocks}
\alias{blocks}
\alias{blockGraphs}
\alias{cohesion}
\alias{hierarchy}
\alias{parent}
\alias{plotHierarchy}
\alias{exportPajek}
\alias{maxcohesion}
\concept{Structurally cohesive blocks}
\title{Calculate Cohesive Blocks}
\description{Calculates cohesive blocks for objects of class \code{igraph}.}
\usage{
cohesive.blocks(graph)

blocks(blocks)
blockGraphs(blocks, graph)
cohesion(blocks)
hierarchy(blocks)
parent(blocks)
plotHierarchy(blocks,
   layout=layout.reingold.tilford(hierarchy(blocks), root=0), \dots)
exportPajek(blocks, graph, file, project.file = TRUE)
maxcohesion(blocks)

\method{print}{cohesiveBlocks}(x, \dots)
\method{summary}{cohesiveBlocks}(object, \dots)
\method{length}{cohesiveBlocks}(x)
\method{plot}{cohesiveBlocks}(x, y, colbar = rainbow(max(cohesion(x))+1),
     col = colbar[maxcohesion(x)+1], mark.groups = blocks(x)[-1], \dots)
}
\arguments{
    \item{graph}{For \code{cohesive.blocks} a graph object of class
      \code{igraph}. It must be undirected and simple. (See
      \code{\link{is.simple}}.)

      For \code{blockGraphs} and \code{exportPajek} the same graph must
      be supplied whose cohesive block structure is given in the
      \code{blocks} argument.
    }
    \item{blocks,x,object}{A \code{cohesiveBlocks} object, created with
      the \code{cohesive.blocks} function.}
    \item{file}{Defines the file (or connection) the Pajek file is
      written to.

      If the \code{project.file} argument is \code{TRUE},
      then it can be a filename (with extension), a file object, or in
      general any king of connection object. The file/connection will be
      opened if it wasn't already.

      If the \code{project.file} argument is \code{FALSE}, then several
      files are created and \code{file} must be a character scalar
      containing the base name of the files, without extension. (But it
      can contain the path to the files.)

      See also details below.
    }
    \item{project.file}{Logical scalar, whether to create a single Pajek
      project file containing all the data, or to create separated files
      for each item. See details below.}
    \item{y}{The graph whose cohesive blocks are supplied in the
      \code{x} argument.}
    \item{colbar}{Color bar for the vertex colors. Its length should be
      at least \eqn{m+1}, where \eqn{m} is the maximum cohesion in the
      graph. Alternatively, the vertex colors can also be directly
    specified via the \code{col} argument.}
    \item{col}{A vector of vertex colors, in any of the usual
      formats. (Symbolic color names (e.g. \sQuote{red}, \sQuote{blue},
      etc.) , RGB colors (e.g. \sQuote{#FF9900FF}), integer numbers
      referring to the current palette. By default the given
      \code{colbar} is used and vertices with the same maximal cohesion
      will have the same color.}
    \item{mark.groups}{A list of vertex sets to mark on the plot by
      circling them. By default all cohesive blocks are marked, except
      the one corresponding to the all vertices.}
    \item{layout}{The layout of a plot, it is simply passed on to
      \code{plot.igraph}, see the possible formats there. By default the
      Reingold-Tilford layout generator is used.}
    \item{\dots}{Additional arguments. \code{plotHierarchy} and
      \code{plot} pass them to \code{plot.igraph}.
      \code{print} and \code{summary} ignore them.
    }
}
\details{
  Cohesive blocking is a method of determining hierarchical subsets of
  graph vertices based on their structural cohesion (or vertex
  connectivity). For a given graph \eqn{G}, a subset of its vertices
  \eqn{S\subset V(G)}{S} is said to be maximally \eqn{k}-cohesive if there is
  no superset of \eqn{S} with vertex connectivity greater than or equal to
  \eqn{k}. Cohesive blocking is a process through which, given a
  \eqn{k}-cohesive set of vertices, maximally \eqn{l}-cohesive subsets are
  recursively identified with \eqn{l>k}. Thus a hiearchy of vertex subsets
  is found, whith the entire graph \eqn{G} at its root. 
  
  The function \code{cohesive.blocks} implements cohesive blocking.
  It returns a \code{cohesiveBlocks} object. \code{cohesiveBlocks}
  should be handled as an opaque class, i.e. its internal structure
  should not be accessed directly, but through the functions listed
  here.

  The function \code{length} can be used on \code{cohesiveBlocks}
  objects and it gives the number of blocks.
  
  The function \code{blocks} returns the actual blocks stored in the
  \code{cohesiveBlocks} object. They are returned in a list of numeric
  vectors, each containing vertex ids.

  The function \code{blockGraphs} is similar, but returns the blocks as
  (induced) subgraphs of the input graph. The various (graph, vertex and
  edge) attributes are kept in the subgraph.

  The function \code{cohesion} returns a numeric vector, the cohesion of
  the different blocks. The order of the blocks is the same as for
  the \code{blocks} and \code{blockGraphs} functions.

  The block hierarchy can be queried using the \code{hierarchy}
  function. It returns an igraph graph, its vertex ids are ordered
  according the order of the blocks in the \code{blocks} and
  \code{blockGraphs}, \code{cohesion}, etc. functions.

  \code{parent} gives the parent vertex of each block, in the block
  hierarchy, for the root vertex it gives 0.

  \code{plotHierarchy} plots the hierarchy tree of the cohesive blocks
  on the active graphics device, by calling \code{igraph.plot}.

  The \code{exportPajek} function can be used to export the graph and
  its cohesive blocks in Pajek format. It can either export a single
  Pajek project file with all the information, or a set of files,
  depending on its \code{project.file} argument. If \code{project.file}
  is \code{TRUE}, then the following information is written to the file
  (or connection) given in the \code{file} argument: (1) the input
  graph, together with its attributes, see \code{\link{write.graph}} for
  details; (2) the hierarchy graph; and (3) one binary partition for
  each cohesive block. If \code{project.file} is \code{FALSE}, then the
  \code{file} argument must be a character scalar and it is used as the
  base name for the generated files. If \code{file} is
  \sQuote{basename}, then the following files are created: (1)
  \sQuote{basename.net} for the original graph; (2)
  \sQuote{basename_hierarchy.net} for the hierarchy graph; (3)
  \sQuote{basename_block_x.net} for each cohesive block, where
  \sQuote{x} is the number of the block, starting with one.

  \code{maxcohesion} returns the maximal cohesion of each vertex,
  i.e. the cohesion of the most cohesive block of the vertex.

  The generic function \code{summary} works on \code{cohesiveBlocks}
  objects and it prints a one line summary to the terminal.

  The generic function \code{print} is also defined on
  \code{cohesiveBlocks} objects and it is invoked automatically if the
  name of the \code{cohesiveBlocks} object is typed in. It produces
  an output like this: \preformatted{
    Cohesive block structure:
    B-1         c 1, n 23
    '- B-2      c 2, n 14   oooooooo.. .o......oo ooo 
       '- B-4   c 5, n  7   ooooooo... .......... ... 
    '- B-3      c 2, n 10   ......o.oo o.oooooo.. ... 
       '- B-5   c 3, n  4   ......o.oo o......... ...     
  }
  The left part shows the block structure, in this case for five
  blocks. The first block always corresponds to the whole graph, even if
  its cohesion is zero. Then cohesion of the block and the number of
  vertices in the block are shown. The last part is only printed if the
  display is wide enough and shows the vertices in the blocks,
  ordered by vertex ids. \sQuote{o} means that the vertex is included,
  a dot means that it is not, and the vertices are shown in groups of
  ten.
  
  The generic function \code{plot} plots the graph, showing one or more
  cohesive blocks in it.  
}
\value{
  \code{cohesive.blocks} returns a \code{cohesiveBlocks}
  object.

  \code{blocks} returns a list of numeric vectors, containing vertex
  ids.

  \code{blockGraphs} returns a list of igraph graphs, corresponding to
  the cohesive blocks.

  \code{cohesion} returns a numeric vector, the cohesion of each block.

  \code{hierarchy} returns an igraph graph, the representation of the
  cohesive block hierarchy.

  \code{parent} returns a numeric vector giving the parent block of each
  cohesive block, in the block hierarchy. The block at the root of the
  hierarchy has no parent and \code{0} is returned for it.

  \code{plotHierarchy}, \code{plot} and \code{exportPajek} return
  \code{NULL}, invisibly.

  \code{maxcohesion} returns a numeric vector with one entry for each
  vertex, giving the cohesion of its most cohesive block.

  \code{print} and \code{summary} return the \code{cohesiveBlocks}
  object itself, invisibly.

  \code{length} returns a numeric scalar, the number of blocks.  
}
\references{
  J. Moody and D. R. White. Structural cohesion and embeddedness: A
  hierarchical concept of social groups. \emph{American Sociological
    Review}, 68(1):103--127, Feb 2003. 
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com} for the current
  implementation, Peter McMahan \email{peter.mcmahan@gmail.com} wrote
  the first version in R.}
\seealso{\code{\link{graph.cohesion}}}
\examples{
## The graph from the Moody-White paper
mw <- graph.formula(1-2:3:4:5:6, 2-3:4:5:7, 3-4:6:7, 4-5:6:7,
                    5-6:7:21, 6-7, 7-8:11:14:19, 8-9:11:14, 9-10,
                    10-12:13, 11-12:14, 12-16, 13-16, 14-15, 15-16,
                    17-18:19:20, 18-20:21, 19-20:22:23, 20-21,
                    21-22:23, 22-23)

mwBlocks <- cohesive.blocks(mw)

# Inspect block membership and cohesion
mwBlocks
blocks(mwBlocks)
cohesion(mwBlocks)

# Save results in a Pajek file
\dontrun{
exportPajek(mwBlocks, mw, file="/tmp/mwBlocks.paj")
}

# Plot the results
if (interactive()) {
  plot(mwBlocks, mw)
}

## The science camp network
camp <- graph.formula(Harry:Steve:Don:Bert - Harry:Steve:Don:Bert,
                      Pam:Brazey:Carol:Pat - Pam:Brazey:Carol:Pat,
                      Holly   - Carol:Pat:Pam:Jennie:Bill,
                      Bill    - Pauline:Michael:Lee:Holly,
                      Pauline - Bill:Jennie:Ann,
                      Jennie  - Holly:Michael:Lee:Ann:Pauline,
                      Michael - Bill:Jennie:Ann:Lee:John,
                      Ann     - Michael:Jennie:Pauline,
                      Lee     - Michael:Bill:Jennie,
                      Gery    - Pat:Steve:Russ:John,
                      Russ    - Steve:Bert:Gery:John,
                      John    - Gery:Russ:Michael)
campBlocks <- cohesive.blocks(camp)
campBlocks

if (interactive()) {
  plot(campBlocks, camp, vertex.label=V(camp)$name, margin=-0.2,
       vertex.shape="rectangle", vertex.size=24, vertex.size2=8,
       mark.border=1, colbar=c(NA, NA,"cyan","orange") )
}
}
\keyword{graphs}
