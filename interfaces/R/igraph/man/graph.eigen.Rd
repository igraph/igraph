\name{spectrum}
\alias{graph.eigen}
\alias{spectrum}
\alias{igraph.eigen.default}
\alias{eigen_defaults}
\concept{Eigenvalues}
\concept{Eigenvectors}
\title{Eigenvalues and eigenvectors of the adjacency matrix of a graph}
\description{
  Calculate selected eigenvalues and eigenvectors of a (supposedly
  sparse) graph.
}
\usage{
spectrum (graph, algorithm = c("arpack", "auto", "lapack", "comp_auto", 
    "comp_lapack", "comp_arpack"), which = list(),
    options = arpack_defaults)
eigen_defaults
}
\arguments{
  \item{graph}{The input graph, can be directed or undirected.}
  \item{algorithm}{The algorithm to use. Currently only \code{arpack} is
    implemented, which uses the ARPACK solver. See also
    \code{\link{arpack}}.}
  \item{which}{A list to specify which eigenvalues and eigenvectors to
    calculate. By default the leading (i.e. largest magnitude)
    eigenvalue and the corresponding eigenvector is calculated.}
  \item{options}{Options for the ARPACK solver. See
    \code{\link{arpack_defaults}}.}
}
\details{
  The \code{which} argument is a list and it specifies which eigenvalues
  and corresponding eigenvectors to calculate: There are eight options:
  \enumerate{
    \item Eigenvalues with the largest magnitude. Set \code{pos} to
    \code{LM}, and \code{howmany} to the number of eigenvalues you
    want.
    \item Eigenvalues with the smallest magnitude. Set \code{pos} to
    \code{SM} and \code{howmany} to the number of eigenvalues you want.
    \item Largest eigenvalues. Set \code{pos} to \code{LA} and
    \code{howmany} to the number of eigenvalues you want.
    \item Smallest eigenvalues. Set \code{pos} to \code{SA} and
    \code{howmany} to the number of eigenvalues you want.
    \item Eigenvalues from both ends of the spectrum. Set \code{pos} to
    \code{BE} and \code{howmany} to the number of eigenvalues you
    want. If \code{howmany} is odd, then one more eigenvalue is returned
    from the larger end.
    \item Selected eigenvalues. This is not (yet) implemented currently.
    \item Eigenvalues in an interval. This is not (yet) implemented.
    \item All eigenvalues. This is not implemented yet. The standard
    \code{eigen} function does a better job at this, anyway.
  }
  
  Note that ARPACK might be unstable for graphs with multiple
  components, e.g. graphs with isolate vertices.
}
\value{
  Depends on the algorithm used.

  For \code{arpack} a list with three entries is returned:
  \item{options}{See the return value for \code{arpack} for a complete
    description.}
  \item{values}{Numeric vector, the eigenvalues.}
  \item{vectors}{Numeric matrix, with the eigenvectors as columns.}
}
% \references{}
\author{ Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{\code{\link{as_adj}} to create a (sparse) adjacency
  matrix.}
\examples{
## Small example graph, leading eigenvector by default
kite <- make_graph("Krackhardt_kite")
spectrum(kite)[c("values", "vectors")]

## Double check
eigen(as_adj(kite, sparse=FALSE))$vectors[,1]

## Should be the same as 'eigen_centrality' (but rescaled)
cor(eigen_centrality(kite)$vector, spectrum(kite)$vectors)

## Smallest eigenvalues
spectrum(kite, which=list(pos="SM", howmany=2))$values
}
\keyword{graphs}
