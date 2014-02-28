\name{layout.reingold.tilford}
\alias{layout.raingold.tilford}
\concept{Graph layout}
\title{The Reingold-Tilford graph layout algorithm}
\description{A tree-like layout, it is perfect for trees, acceptable for
  graphs with not too many cycles.
}
\usage{
layout.reingold.tilford (graph, root=numeric(), circlular=FALSE,
     rootlevel=numeric(), mode="out", flip.y=TRUE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{root}{The index of the root vertex or root vertices.
    If this is a non-empty vector then the supplied vertex ids are used 
    as the roots of the trees (or a single tree if the graph is connected).
    If it is an empty vector, then the root vertices are automatically
    calculated based on topological sorting, performed with the opposite
    mode than the \code{mode} argument. After the vertices have been
    sorted, one is selected from each component.}
  \item{circular}{Logical scalar, whether to plot the tree in a circular
    fashion. Defaults to \code{FALSE}, so the tree branches are going
    bottom-up (or top-down, see the \code{flip.y} argument.}
  \item{rootlevel}{This argument can be useful when drawing forests
    which are not trees (i.e. they are unconnected and have tree
    components). It specifies the level of the root vertices for every
    tree in the forest. It is only considered if the \code{roots}
    argument is not an empty vector.}
  \item{mode}{Specifies which edges to consider when building the tree.
    If it is \sQuote{out}, then only the outgoing, if it is \sQuote{in},
    then only the incoming edges of a parent are considered. If it is 
    \sQuote{all} then all edges are used (this was the behavior in 
    igraph 0.5 and before). This parameter also influences how the root 
    vertices are calculated, if they are not given. See the \code{roots}
    parameter.}
  \item{flip.y}{Logical scalar, whether to flip the \sQuote{y}
    coordinates. The default is flipping because that puts the root
    vertex on the top.}
}
\details{
 Arranges the nodes in a tree where the given node is used as the root.
 The tree is directed downwards and the parents are centered above its
 children. For the exact algorithm, the refernce below.

 If the given graph is not a tree, a breadth-first search is executed
 first to obtain a possible spanning tree.
}
\value{
  A numeric matrix with two columns, and one row for each vertex.
}
\references{
  Reingold, E and Tilford, J (1981). Tidier drawing of trees.
  \emph{IEEE Trans. on Softw. Eng.}, SE-7(2):223--228.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
tree <- graph.tree(20, 3)
plot(tree, layout=layout.reingold.tilford)
plot(tree, layout=layout.reingold.tilford(tree, flip.y=FALSE))
plot(tree, layout=layout.reingold.tilford(tree, circular=TRUE))

tree2 <- graph.tree(10, 3) + graph.tree(10, 2)
plot(tree2, layout=layout.reingold.tilford)
plot(tree2, layout=layout.reingold.tilford(tree2, root=c(1,11),
                                           rootlevel=c(2,1))
}
\keyword{graphs}
