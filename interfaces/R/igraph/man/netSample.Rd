\name{netSample}

\alias{alanine}
\alias{immuno}
\alias{cellCycle}
\alias{synSmall}
\alias{synLarge}

\docType{data}
\title{Real-World Networks}

\description{
A sample of networks taken from real-world data (some published in peer-reviewed journals).
}

\usage{
data(alanine)
data(immuno)
data(cellCycle)
data(synSmall)
data(synLarge)
}

\format{
Networks of the class \code{igraph} (see the doc of \pkg{igraph}).
}

\details{
\code{alanine} is the directed and strongly connected network of configurations of the di-alanine peptide. It consists of 1832 vertices representing configurations visited by the peptide while folded by Molecular Dynamics. The configurations are characterized by the dihedral angles (vertex attributes \code{phi} and \code{psi}). Edges represent transitions between configurations and the number of transitions is given by the edge attribute \code{weight}.

\code{immuno} is the undirected and connected network of interactions in the immunoglobulin protein. It is made up of 1316 vertices representing amino-acids and an edge is drawn between two amino-acids if the shortest distance between their \eqn{C_\alpha} atoms is smaller than the threshold value \eqn{\theta = 8 Angstrom}.

\code{cellCycle} is the cell cycle network of the yeast fission \emph{Schizosaccharomices pombe}. This complete and undirected network is composed of 692 vertices representing periodically expressed genes. The gene name, phase, and position in the cell cycle are available through the vertex attributes \code{name}, \code{phase} and \code{pos}. Edge weights represent the phase difference between two genes, say \eqn{i} and \eqn{j}, and are set to \deqn{\textrm{exp}[\beta \cos(\phi_i -\phi_j)],} where \eqn{\beta = 10}.

\code{synSmall} and \code{synLarge} are two connected components of a linguistic network based on the relation of synonymy in French. The nodes are the words in a given sense (vertex attribute \code{name}). Two nodes are connected if they are considered synonyms. The small and large components have respective size of 185 and 9812 vertices.
}

\source{
\code{alanine}

D. Gfeller, P. De Los Rios, A. Caﬂisch, and F. Rao, Complex network analysis of free-energy landscapes, \emph{Proc. Natl. Acad. Sci. U.S.A.}, \bold{104}(1817), 2007.
\url{http://www.pnas.org/content/104/6/1817.abstract?etoc}

D. Gfeller, and P. De Los Rios, Spectral Coarse Graining of Complex Networks, \emph{Physical Review Letters}, \bold{99}(3), 2007.
\url{http://arxiv.org/abs/0706.0812}


\code{immuno}

D. Gfeller, Simplifying complex networks: from a clustering to a coarse graining strategy, \emph{PhD Thesis EPFL}, \bold{no 3888}, 2007.
\url{http://library.epfl.ch/theses/?nr=3888}


\code{cellCycle}

C. Caretta Cartozo, P. De Los Rios, F. Piazza, and P. Li\`o, Bottleneck genes and community structure in the cell cycle network of \emph{S. pombe}, \emph{PLoS computational Biology}, \bold{3}(6), 2007.
\url{http://www.ploscompbiol.org}

A. Oliva, A. Rosebrock, F. Ferrezuelo, S. Pyne, H. Chen, S. Skiena, B. Futcher, and J. Leatherwood, The cell cycle–regulated genes of \emph{schizosaccharomyces pombe}, \emph{PLoS Biology}, \bold{3}(7), 2005.
\url{http://biology.plosjournals.org}


\code{synSmall} and \code{synLarge}

D. Gfeller, J.-C. Chappelier, and P. De Los Rios, Finding instabilities in the community structure of complex networks, \emph{Physical Review E}, \bold{72}(056135), 2005.
\url{http://arxiv.org/abs/cond-mat/0503593}
}

\examples{

# alanine is transformed into an undirected network
# while preserving the edge weights
data(alanine)
g <- as.undirected(alanine, mode = "each")
g <- simplify(g, remove.loops = FALSE)
# computing the sum of weights of all edges adjacent to a vertex  
l <- get.adjedgelist(g)
we <- E(g)$weight
n <- vcount(g)
weight.sum <- numeric(n)
for(i in 1:n)
	weight.sum[i] <- sum(we[ l[[i]]+1 ])
weight.sum <- 20*weight.sum / max(weight.sum) + 1

g <- simplify(g) # remove self-loops for better visualization
layout <- layout.lgl(g) #can take some time
plot(g, layout = layout, vertex.size = weight.sum, vertex.label = NA)
#TODO: improve the layout


## working with cellCycle
## remove edges with phase difference smaller than 0.2
data(cellCycle)
g <- delete.edges(cellCycle, which(E(cellCycle)$weight < 18045.76)-1)
layout <- layout.kamada.kawai(g)
col <- c("red", "yellow", "blue", "green")
plot(g, layout = layout, vertex.color = col[V(g)$pos],
	vertex.size = 5, vertex.label = NA)
## should end up with a circular network
## colors display the gene position along the cell cycle
}
\keyword{datasets}

