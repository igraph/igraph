/*  -- translated by f2c (version 20240504).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* > \brief <b> DGESV computes the solution to system of linear equations A * X = B for GE matrices</b>   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DGESV + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesv.f
">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesv.f
">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesv.f
">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )   

         INTEGER            INFO, LDA, LDB, N, NRHS   
         INTEGER            IPIV( * )   
         DOUBLE PRECISION   A( LDA, * ), B( LDB, * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DGESV computes the solution to a real system of linear equations   
   >    A * X = B,   
   > where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   
   >   
   > The LU decomposition with partial pivoting and row interchanges is   
   > used to factor A as   
   >    A = P * L * U,   
   > where P is a permutation matrix, L is unit lower triangular, and U is   
   > upper triangular.  The factored form of A is then used to solve the   
   > system of equations A * X = B.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The number of linear equations, i.e., the order of the   
   >          matrix A.  N >= 0.   
   > \endverbatim   
   >   
   > \param[in] NRHS   
   > \verbatim   
   >          NRHS is INTEGER   
   >          The number of right hand sides, i.e., the number of columns   
   >          of the matrix B.  NRHS >= 0.   
   > \endverbatim   
   >   
   > \param[in,out] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          On entry, the N-by-N coefficient matrix A.   
   >          On exit, the factors L and U from the factorization   
   >          A = P*L*U; the unit diagonal elements of L are not stored.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.  LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[out] IPIV   
   > \verbatim   
   >          IPIV is INTEGER array, dimension (N)   
   >          The pivot indices that define the permutation matrix P;   
   >          row i of the matrix was interchanged with row IPIV(i).   
   > \endverbatim   
   >   
   > \param[in,out] B   
   > \verbatim   
   >          B is DOUBLE PRECISION array, dimension (LDB,NRHS)   
   >          On entry, the N-by-NRHS matrix of right hand side matrix B.   
   >          On exit, if INFO = 0, the N-by-NRHS solution matrix X.   
   > \endverbatim   
   >   
   > \param[in] LDB   
   > \verbatim   
   >          LDB is INTEGER   
   >          The leading dimension of the array B.  LDB >= max(1,N).   
   > \endverbatim   
   >   
   > \param[out] INFO   
   > \verbatim   
   >          INFO is INTEGER   
   >          = 0:  successful exit   
   >          < 0:  if INFO = -i, the i-th argument had an illegal value   
   >          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization   
   >                has been completed, but the factor U is exactly   
   >                singular, so the solution could not be computed.   
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup doubleGEsolve   

    =====================================================================   
   Subroutine */ int igraphdgesv_(integer *n, integer *nrhs, doublereal *a, integer 
	*lda, integer *ipiv, doublereal *b, integer *ldb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int igraphdgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), igraphxerbla_(char *, integer *, 
	    ftnlen), igraphdgetrs_(char *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *);


/*  -- LAPACK driver routine (version 3.4.0) --   
    -- LAPACK is a software package provided by Univ. of Tennessee,    --   
    -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--   
       November 2011   


    =====================================================================   


       Test the input parameters.   

       Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;
    b_dim1 = *ldb;
    b_offset = 1 + b_dim1;
    b -= b_offset;

    /* Function Body */
    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*lda < max(1,*n)) {
	*info = -4;
    } else if (*ldb < max(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DGESV ", &i__1, (ftnlen)6);
	return 0;
    }

/*     Compute the LU factorization of A. */

    igraphdgetrf_(n, n, &a[a_offset], lda, &ipiv[1], info);
    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */

	igraphdgetrs_("No transpose", n, nrhs, &a[a_offset], lda, &ipiv[1], &b[
		b_offset], ldb, info);
    }
    return 0;

/*     End of DGESV */

} /* igraphdgesv_ */

