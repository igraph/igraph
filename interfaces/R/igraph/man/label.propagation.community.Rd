\name{label.propagation.community}
\alias{label.propagation.community}
\concept{Community structure}
\title{Finding communities based on propagating labels}
\description{
  This is a fast, nearly linear time algorithm for detecting community
  structure in networks. In works by labeling the vertices with unique
  labels and then updating the labels by majority voting in the
  neighborhood of the vertex.
}
\usage{
label.propagation.community (graph, weights = NULL,
            initial = NULL, fixed = NULL)
}
\arguments{
  \item{graph}{The input graph, should be undirected to make sense.}
  \item{weights}{An optional weight vector. It should contain a positive
    weight for all the edges. The \sQuote{weight} edge attribute is used
    if present. Supply \sQuote{\code{NA}} here if you want to ignore the
    \sQuote{weight} edge attribute.}
  \item{initial}{The initial state. If \code{NULL}, every vertex will have
    a different label at the beginning. Otherwise it must be a vector
    with an entry for each vertex. Non-negative values denote different
    labels, negative entries denote vertices without labels.}
  \item{fixed}{Logical vector denoting which labels are fixed. Of course
    this makes sense only if you provided an initial state, otherwise
    this element will be ignored. Also note that vertices without labels
    cannot be fixed.}
}
\details{
  This function implements the community detection method described in:
  Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm
  to detect community structures in large-scale networks. Phys Rev E
  76, 036106. (2007). This version extends the original method by
  the ability to take edge weights into consideration and also
  by allowing some labels to be fixed.

  From the abstract of the paper:
  \dQuote{In our algorithm every node is initialized with a unique label
    and at every step each node adopts the label that most of its
    neighbors currently have. In this iterative process densely
    connected groups of  nodes form a consensus on a unique label to form
    communities.}
}
\value{
  A numeric vector giving the final label (=community) of each vertex.
}
\references{
  Raghavan, U.N. and Albert, R. and Kumara, S.: Near linear time algorithm
  to detect community structures in large-scale networks. \emph{Phys Rev E}
  76, 036106. (2007)
}
\author{
  Tamas Nepusz \email{ntamas@rmki.kfki.hu} for the C implementation,
  Gabor Csardi \email{Gabor.Csardi@unil.ch} for this manual page.
}
\seealso{\code{\link{fastgreedy.community}},
  \code{\link{walktrap.community}} and \code{\link{spinglass.community}}
  for other community detection methods.}
\examples{
  g <- erdos.renyi.game(10, 5/10) \%du\% erdos.renyi.game(9, 5/9)
  g <- add.edges(g, c(0, 11))
  label.propagation.community(g)
}
\keyword{graphs}
