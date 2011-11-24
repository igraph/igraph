\name{Spectral coarse graining}
\alias{SCG}
\concept{Spectral coarse graining}
\title{Spectral Coarse Graining}
\description{
  Functions to perform the Spectral Coarse Graining (SCG) of matrices
  and graphs.
}
\section{Introduction}{
  The SCG functions provide a framework, called Spectral Coarse Graining
  (SCG), for reducing large graphs while preserving their
  \emph{spectral-related features}, that is features closely related
  with the eigenvalues and eigenvectors of a graph matrix (which for now
  can be the adjacency, the stochastic, or the Laplacian matrix).

  Common examples of such features comprise the first-passage-time of
  random walkers on Markovian graphs, thermodynamic properties of
  lattice models in statistical physics (e.g. Ising model), and the
  epidemic threshold of epidemic network models (SIR and SIS models).

  SCG differs from traditional clustering schemes by producing a
  \emph{coarse-grained graph} (not just a partition of the vertices),
  representative of the original one. As shown in [1], Principal
  Component Analysis can be viewed as a particular SCG, called
  \emph{exact SCG}, where the matrix to be coarse-grained is the
  covariance matrix of some data set. 

  SCG should be of interest to practitioners of various
  fields dealing with problems where matrix eigenpairs play an important
  role, as for instance is the case of dynamical processes on networks.
}
\section{SCG in brief}{
  The main idea of SCG is to operate on a matrix a shrinkage operation
  specifically designed to preserve some of the matrix eigenpairs while
  not altering other important matrix features (such as its structure). 
  Mathematically, this idea was expressed as follows. Consider a
  (complex) \eqn{n \times n}{n x n} matrix \eqn{M} and form the product
  \deqn{\widetilde{M} = LMR^*,}{M'=LMR*,}
  where \eqn{\tilde n < n}{n' < n} and
  \eqn{L,R \in \mathbf{C}^{\tilde n \times n}}{L,R in C[n'xn]} are such
  that \eqn{LR^*=I_{\tilde n}}{LR*=I[n']} (\eqn{R^*}{R*} denotes the
  conjugate transpose of \eqn{R}). Under these assumptions, it can be
  shown that \eqn{P=R^*L}{P=R*L} is an \eqn{\tilde n}{n'}-rank projector
  and that, if \eqn{(\lambda, v)}{(lambda, v)} is a (right) 
  eigenpair of \eqn{M} (i.e. \eqn{Mv=\lambda v}{Mv=lambda v}) and
  \eqn{P} is orthogonal, there exists an eigenvalue \eqn{\tilde
    \lambda}{lambda'} of \eqn{\widetilde M}{M'} such that
  \deqn{|\lambda-\tilde{\lambda}| \le \textrm{const} \Vert e_P(v)\Vert
    [1 + O(\Vert e_P(v)\Vert^2)],}{%
  |lambda-lambda'| <= const ||e[P](v)|| [1+O(||e[P](v)||^2)],}
  where \eqn{\Vert e_P(v)\Vert = \Vert
    v-Pv\Vert}{||e[P](v)||=||v-Pv||}. 
  Hence, if \eqn{P} (or equivalently \eqn{L}, \eqn{R}) is chosen so as
  to make \eqn{\Vert e_P(v)\Vert}{||e[P](v)||} as small as possible, one
  can preserve to any desired level the original eigenvalue
  \eqn{\lambda}{lambda} in the coarse-grained matrix \eqn{\widetilde M}{M'};
  under extra assumptions on \eqn{M}, this result can be generalized to
  eigenvectors [1]. This leads to the following generic definition of a
  SCG problem.

  Given \eqn{M\in\mathbf{C}^{n\times n}}{M in C[nxn]} and
  \eqn{(\lambda,v)}{(lambda,v)} a (right) eigenpair of \eqn{M} to be
  preserved by the coarse graining, the problem is to find a projector
  \eqn{\widehat{P}}{P'} solving
  \deqn{\min_{P\in \Omega} \Vert e_{P}(v)\Vert,}{%
  min(||e[P](v)||, p in Omega),}
  where \eqn{\Omega}{Omega} is a set of projectors in
  \eqn{\mathbf{C}^{n\times n}}{C[nxn]} described by some 
  ad hoc constraints \eqn{c_{1},\dots c_{r}}{c[1], ..., c[r]}
  (e.g. \eqn{c_{1}:P\in\mathbf{R}^{n\times n}, c_{2}:P=P^{t},
    c_{3}:P_{ij}\ge 0}{c[1]: P in R[nxn], c[2]: P=t(P), c[3]: P[i,j] >=
    0}, etc).

  Choosing pertinent constraints to solve the SCG problem is of great
  importance in applications. For instance, in the absence of
  constraints the SCG problem is solved trivially by
  \eqn{\widehat{P}=vv^*}{P'=vv*} (\eqn{v} is assumed normalized). We have
  designed a particular constraint, called \emph{homogeneous mixing},
  which ensures that vertices belonging to the same group are merged
  consistently from a physical point of view (see [1] for
  details). Under this constraint the SCG problem reduces to finding the
  partition of \eqn{\{1,\dots,n\}}{1, ..., n} (labeling the original
  vertices) minimizing
  \deqn{\Vert e_P(v)\Vert^{2} =
    \sum_{\alpha=1}^{\tilde{n}}\sum_{i\in\alpha}[v(i)-(Pv)(i)]^{2},}{%
    ||e[P](v)||^2 = sum([v(i)-(Pv)(i)]^2; alpha=1,...,n', i in alpha),} 
  where \eqn{\alpha}{alpha} denotes a group (i.e. a block) in a partition of
  \eqn{ \{1,\dots,n\} }{{1, ..., n}}, and \eqn{|\alpha|}{|alpha|} is the
  number of elements in \eqn{\alpha}{alpha}. 

  If \eqn{M} is symmetric or stochastic, for instance, then it may be
  desirable (or mandatory) to choose \eqn{L}, \eqn{R} so that
  \eqn{\widetilde M}{M'} is symmetric or stochastic as well. This
  \emph{structural constraint} has led to the construction of particular
  semi-projectors for symmetric [1], stochastic [3] and Laplacian [2]
  matrices, that are made available.

  In short, the coarse graining of matrices and graphs involves: 
  \enumerate{
    \item Retrieving a matrix or a graph matrix \eqn{M} from the
    problem. 
    \item Computing the eigenpairs of \eqn{M} to be preserved in the
    coarse-grained graph or matrix. 
    \item Setting some problem-specific constraints (e.g. dimension of
    the coarse-grained object). 
    \item Solving the constrained SCG problem, that is finding
    \eqn{\widehat{P}}{P'}. 
    \item Computing from \eqn{\widehat{P}}{P'} two semi-projectors
    \eqn{\widehat{L}}{L'} and \eqn{\widehat{R}}{R'} (e.g. following the method
    proposed in [1]). 
    \item Working out the product \eqn{\widetilde{M}=\widehat{L} M
      \widehat{R}^*}{M'=L'MR'*} and, if needed, defining from
    \eqn{\widetilde{M}}{M'} a coarse-grained graph. 
  }
}
\section{Functions for performing SCG}{
  The main function is the \dQuote{all-in-one}
  \code{\link{scg}}. This function handles all the steps involved in the
  Spectral Coarse Graining (SCG) of some particular matrices and graphs
  as described above and in reference [1]. In more details,
  \code{\link{scg}} computes some prescribed eigenpairs of a matrix or a
  graph matrix (for now adjacency, Laplacian and stochastic matrices are
  available), works out an optimal partition to preserve the eigenpairs,
  and finally outputs a coarse-grained matrix or graph along with other
  useful information. 
 
  These steps can also be carried out independently: (1) Use
  \code{\link{get.adjacency}}, \code{\link{graph.laplacian}} or
  \code{\link{stochasticMatrix}} to compute a matrix \eqn{M}. (2) Work
  out some prescribed eigenpairs of \eqn{M} e.g. by 
  means of \code{eigen} or \code{arpack}. (3) Invoke one the four
  algorithms of the function \code{\link{scgGrouping}} to get a
  partition that will preserve the eigenpairs in the coarse-grained
  matrix. (4) Compute the semi-projectors \eqn{L} and \eqn{R} using
  \code{\link{scgSemiProjectors}} and from there the coarse-grained
  matrix \eqn{\widetilde{M}=LMR^*}{M'=LMR*}. If necessary, construct a
  coarse-grained graph from \eqn{\widetilde{M}}{M'} (e.g. as in [1]). 
}
\references{
  D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios,
  Shrinking Matrices while Preserving their Eigenpairs with Application
  to the Spectral Coarse Graining of Graphs. Submitted to \emph{SIAM
    Journal on Matrix Analysis and Applications}, 2008. 
  \url{http://people.epfl.ch/david.morton}
  
  D. Gfeller, and P. De Los Rios, Spectral Coarse Graining and
  Synchronization in Oscillator Networks. \emph{Physical Review
    Letters}, \bold{100}(17), 2008. 
  \url{http://arxiv.org/abs/0708.2055}

  D. Gfeller, and P. De Los Rios, Spectral Coarse Graining of Complex
  Networks, \emph{Physical Review Letters}, \bold{99}(3), 2007. 
  \url{http://arxiv.org/abs/0706.0812}
}
\author{David Morton de Lachapelle \email{david.morton@epfl.ch},
  \email{david.mortondelachapelle@swissquote.ch}}
\keyword{graphs}
