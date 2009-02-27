\name{arpack}
\alias{arpack}
\alias{arpack-options}
\alias{igraph.arpack.default}
\alias{arpack.unpack.complex}
\concept{Eigenvalues}
\concept{Eigenvectors}
\concept{ARPACK}
\title{ARPACK eigenvector calculation}
\description{Interface to the ARPACK library for calculating
  eigenvectors of sparse matrices}
\usage{
arpack(func, extra = NULL, sym = FALSE, options = igraph.arpack.default, 
    env = parent.frame(), complex=!sym)
arpack.unpack.complex(vectors, values, nev)
}
\arguments{
  \item{func}{The function to perform the matrix-vector
    multiplication. ARPACK requires to perform these by the user. The
    function gets the vector \eqn{x} as the first argument, and it should
    return \eqn{Ax}, where \eqn{A} is the \dQuote{input matrix}. (The
    input matrix is never given explicitly.) The second argument is
    \code{extra}.}
  \item{extra}{Extra argument to supply to \code{func}.}
  \item{sym}{Logical scalar, whether the input matrix is
    symmetric. Always supply \code{TRUE} here if it is, since it can
    speed up the computation. }
  \item{options}{Options to ARPACK, a named list to overwrite some of
    the default option values. See details below.}
  \item{env}{The environment in which \code{func} will be evaluated.}
  \item{complex}{Whether to convert the eigenvectors returned by ARPACK
    into R complex vectors. By default this is not done for symmetric
    problems (these only have real eigenvectors/values), but only
    non-symmetric ones. If you have a non-symmetric problem, but you're
    sure that the results will be real, then supply \code{FALSE} here.
    The conversion is done by calling \code{arpack.unpack.complex}.
  }
  \item{vectors}{Eigenvectors, as returned by ARPACK.}
  \item{values}{Eigenvalues, as returned by ARPACK}
  \item{nev}{The number of eigenvectors/values to extract. This can be
    less than or equal to the number of eigenvalues requested in the
    original ARPACK call.}
}
\details{
  ARPACK is a library for solving large scale eigenvalue problems.
  The package is designed to compute a few eigenvalues and corresponding
  eigenvectors of a general \eqn{n} by \eqn{n} matrix \eqn{A}. It is
  most appropriate for large sparse or structured matrices \eqn{A} where
  structured means that a matrix-vector product \code{w <- Av} requires
  order \eqn{n} rather than the usual order \eqn{n^2} floating point
  operations. Please see 
  \url{http://www.caam.rice.edu/software/ARPACK/} for details.

  This function is an interface to ARPACK. igraph does not contain all
  ARPACK routines, only the ones dealing with symmetric and
  non-symmetric eigenvalue problems using double precision real
  numbers.

  The eigenvalue calculation in ARPACK (in the simplest
  case) involves the calculation of the \eqn{Av} product where \eqn{A}
  is the matrix we work with and \eqn{v} is an arbitrary vector. The
  function supplied in the \code{fun} argument is expected to perform
  this product. If the product can be done efficiently, e.g. if the
  matrix is sparse, then \code{arpack} is usually able to calculate the
  eigenvalues very quickly.

  The \code{options} argument specifies what kind of calculation to
  perform. It is a list with the following members, they correspond
  directly to ARPACK parameters. On input it has the following fields:
  \itemize{
    \item{bmat}{Character constant, possible values: \sQuote{\code{I}},
      stadard eigenvalue problem, \eqn{Ax=\lambda x}{A*x=lambda*x}; and
      \sQuote{\code{G}}, generalized eigenvalue problem,
      \eqn{Ax=\lambda B x}{A*x=lambda B*x}. Currently only
      \sQuote{\code{I}} is supported.}
    \item{n}{Numeric scalar. The dimension of the eigenproblem. You only
      need to set this if you call \code{\link{arpack}}
      directly. (I.e. not needed for \code{\link{evcent}},
      \code{\link{page.rank}}, etc.)}
    \item{which}{Specify which eigenvalues/vectors to compute, character
      constant with exactly two characters.

      Possible values for symmetric input matrices:
      \itemize{
	\item{\sQuote{\code{LA}}}{Compute \code{nev} largest (algebraic)
	  eigenvalues.}
	\item{\sQuote{\code{SA}}}{Compute \code{nev} smallest (algebraic)
	  eigenvalues.}
	\item{\sQuote{\code{LM}}}{Compute \code{nev} largest (in magnitude)
	  eigenvalues.}
	\item{\sQuote{\code{SM}}}{Compute \code{nev} smallest (in magnitude)
	  eigenvalues.}
	\item{\sQuote{\code{BE}}}{Compute \code{nev} eigenvalues, half
	  from each end of the spectrum. When \code{nev} is odd, compute
	  one more from the high end than from the low end.}
      }

      Possible values for non-symmetric input matrices:
      \itemize{
	\item{\sQuote{\code{LM}}}{Compute \code{nev} eigenvalues of
	  largest magnitude.}
	\item{\sQuote{\code{SM}}}{Compute \code{nev} eigenvalues of
	  smallest magnitude.}
	\item{\sQuote{\code{LR}}}{Compute \code{nev} eigenvalues of
	  largest real part.}
	\item{\sQuote{\code{SR}}}{Compute \code{nev} eigenvalues of
	  smallest real part.}
	\item{\sQuote{\code{LI}}}{Compute \code{nev} eigenvalues of
	  largest imaginary part.}
	\item{\sQuote{\code{SI}}}{Compute \code{nev} eigenvalues of
	  smallest imaginary part.}
      }	

      This parameter is sometimes overwritten by the various functions,
      e.g. \code{\link{page.rank}} always sets \sQuote{\code{LM}}.
    }
    \item{nev}{Numeric scalar. The number of eigenvalues to be computed.}
    \item{tol}{Numeric scalar. Stopping criterion: the relative accuracy
      of the Ritz value is considered acceptable if its error is less
      than \code{tol} times its estimated value. If this is set to zero
      then machine precision is used.}
    \item{ncv}{Number of Lanczos vectors to be generated.}
    \item{ldv}{Numberic scalar. It should be set to
      zero in the current implementation.}
    \item{ishift}{Either zero or one. If zero then the shifts are
      provided by the user via reverse communication. If one then exact
      shifts with respect to the reduced tridiagonal matrix \eqn{T}.
      Please always set this to one.}
    \item{maxiter}{Maximum number of Arnoldi update iterations allowed. }
    \item{nb}{Blocksize to be used in the recurrence. Please always
      leave this on the default value, one.}
    \item{mode}{The type of the eigenproblem to be solved.
      Possible values if the input matrix is symmetric:
      \itemize{
	\item{1}{\eqn{Ax=\lambda x}{A*x=lambda*x}, \eqn{A} is symmetric.}
	\item{2}{\eqn{Ax=\lambda Mx}{A*x=lambda*M*x}, \eqn{A} is
	  symmetric, \eqn{M} is symmetric positive definite.}
	\item{3}{\eqn{Kx=\lambda Mx}{K*x=lambda*M*x}, \eqn{K} is
	  symmetric, \eqn{M} is symmetric positive semi-definite.}
	\item{4}{\eqn{Kx=\lambda KGx}{K*x=lambda*KG*x}, \eqn{K} is
	  symmetric positive semi-definite, \eqn{KG} is symmetric
	  indefinite.}
	\item{5}{\eqn{Ax=\lambda Mx}{A*x=lambda*M*x}, \eqn{A} is
	  symmetric, \eqn{M} is symmetric positive
	  semi-definite. (Cayley transformed mode.)}
      }
      Please note that only \code{mode==1} was tested and other values
      might not work properly.

      Possible values if the input matrix is not symmetric:
      \itemize{
	\item{1}{\eqn{Ax=\lambda x}{A*x=lambda*x}.}
	\item{2}{\eqn{Ax=\lambda Mx}{A*x=lambda*M*x}, \eqn{M} is
	  symmetric positive definite.}
	\item{3}{\eqn{Ax=\lambda Mx}{A*x=lambda*M*x}, \eqn{M} is
	  symmetric semi-definite.}
	\item{4}{\eqn{Ax=\lambda Mx}{A*x=lambda*M*x}, \eqn{M} is
	      symmetric semi-definite.}
      }
      Please note that only \code{mode==1} was tested and other values
      might not work properly.
    }
    \item{start}{Not used currently. Later it be used to set a starting
      vector.}
    \item{sigma}{Not used currently.}
    \item{sigmai}{Not use currently.}

    On output the following additional fields are added:
    \itemize{
      \item{info}{Error flag of ARPACK. Possible values:
	\itemize{
	  \item{0}{Normal exit.}
	  \item{1}{Maximum number of iterations taken.}
	  \item{3}{No shifts could be applied during a cycle of the
	    Implicitly restarted Arnoldi iteration. One possibility
	    is to increase the size of \code{ncv} relative to \code{nev}.}
	}

	ARPACK can return more error conditions than these, but they are
	converted to regular igraph errors.
      }
      \item{iter}{Number of Arnoldi iterations taken.}
      \item{nconv}{Number of \dQuote{converged} Ritz values. This
	represents the number of Ritz values that satisfy the
	convergence critetion. }
      \item{numop}{Total number of matrix-vector multiplications.}
      \item{numopb}{Not used currently.}
      \item{numreo}{Total number of steps of re-orthogonalization.}
    }
  }
  Please see the ARPACK documentation for additional details.

  \code{arpack.unpack.complex} is a (semi-)internal function that
  converts the output of the non-symmetric ARPACK solver to a more
  readable format. It is called internally by \code{arpack}.
  
}
\value{
  A named list with the following members:
  \item{values}{Numeric vector, the desired eigenvalues.}
  \item{vectors}{Numeric matrix, the desired eigenvectors as
    columns. If \code{complex=TRUE} (the default for non-symmetric
    problems), then the matrix is complex.}
  \item{options}{A named list with the supplied \code{options} and
    some information about the performed calculation, including an
    ARPACK exit code. See the details above.
  }
}
\references{
  D.C. Sorensen, Implicit Application of Polynomial Filters in
  a k-Step Arnoldi Method. \emph{SIAM J. Matr. Anal. Apps.}, 13 (1992),
  pp 357-385.
  
  R.B. Lehoucq, Analysis and Implementation of an Implicitly
  Restarted Arnoldi Iteration. \emph{Rice University Technical Report}
  TR95-13, Department of Computational and Applied Mathematics.
  
  B.N. Parlett & Y. Saad, Complex Shift and Invert Strategies for
  Real Matrices. \emph{Linear Algebra and its Applications}, vol 88/89,
  pp 575-595, (1987).
}
\author{Rich Lehoucq, Kristi Maschhoff, Danny Sorensen, Chao Yang for
  ARPACK, Gabor Csardi \email{csardi@rmki.kfki.hu} for the R interface.}
\seealso{\code{\link{evcent}}, \code{\link{page.rank}},
  \code{\link{hub.score}}, \code{\link{leading.eigenvector.community}}
  are some of the functions in igraph which use ARPACK. The ARPACK
  homepage is at \url{http://www.caam.rice.edu/software/ARPACK/}.
}
\examples{
# Identity matrix
f <- function(x, extra=NULL) x
arpack(f, options=list(n=10, nev=2, ncv=4), sym=TRUE)

# Graph laplacian of a star graph (undirected), n>=2
# Note that this is a linear operation
f <- function(x, extra=NULL) {
  y <- x
  y[1] <- (length(x)-1)*x[1] - sum(x[-1])
  for (i in 2:length(x)) {
    y[i] <- x[i] - x[1]
  }
  y
}

arpack(f, options=list(n=10, nev=1, ncv=3), sym=TRUE)

# double check
eigen(graph.laplacian(graph.star(10, mode="undirected")))

}
\keyword{graphs}
