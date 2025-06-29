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

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = 1.;
static integer c_n1 = -1;

/* > \brief \b DGETRS   

    =========== DOCUMENTATION ===========   

   Online html documentation available at   
              http://www.netlib.org/lapack/explore-html/   

   > \htmlonly   
   > Download DGETRS + dependencies   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgetrs.
f">   
   > [TGZ]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgetrs.
f">   
   > [ZIP]</a>   
   > <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgetrs.
f">   
   > [TXT]</a>   
   > \endhtmlonly   

    Definition:   
    ===========   

         SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )   

         CHARACTER          TRANS   
         INTEGER            INFO, LDA, LDB, N, NRHS   
         INTEGER            IPIV( * )   
         DOUBLE PRECISION   A( LDA, * ), B( LDB, * )   


   > \par Purpose:   
    =============   
   >   
   > \verbatim   
   >   
   > DGETRS solves a system of linear equations   
   >    A * X = B  or  A**T * X = B   
   > with a general N-by-N matrix A using the LU factorization computed   
   > by DGETRF.   
   > \endverbatim   

    Arguments:   
    ==========   

   > \param[in] TRANS   
   > \verbatim   
   >          TRANS is CHARACTER*1   
   >          Specifies the form of the system of equations:   
   >          = 'N':  A * X = B  (No transpose)   
   >          = 'T':  A**T* X = B  (Transpose)   
   >          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)   
   > \endverbatim   
   >   
   > \param[in] N   
   > \verbatim   
   >          N is INTEGER   
   >          The order of the matrix A.  N >= 0.   
   > \endverbatim   
   >   
   > \param[in] NRHS   
   > \verbatim   
   >          NRHS is INTEGER   
   >          The number of right hand sides, i.e., the number of columns   
   >          of the matrix B.  NRHS >= 0.   
   > \endverbatim   
   >   
   > \param[in] A   
   > \verbatim   
   >          A is DOUBLE PRECISION array, dimension (LDA,N)   
   >          The factors L and U from the factorization A = P*L*U   
   >          as computed by DGETRF.   
   > \endverbatim   
   >   
   > \param[in] LDA   
   > \verbatim   
   >          LDA is INTEGER   
   >          The leading dimension of the array A.  LDA >= max(1,N).   
   > \endverbatim   
   >   
   > \param[in] IPIV   
   > \verbatim   
   >          IPIV is INTEGER array, dimension (N)   
   >          The pivot indices from DGETRF; for 1<=i<=N, row i of the   
   >          matrix was interchanged with row IPIV(i).   
   > \endverbatim   
   >   
   > \param[in,out] B   
   > \verbatim   
   >          B is DOUBLE PRECISION array, dimension (LDB,NRHS)   
   >          On entry, the right hand side matrix B.   
   >          On exit, the solution matrix X.   
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
   > \endverbatim   

    Authors:   
    ========   

   > \author Univ. of Tennessee   
   > \author Univ. of California Berkeley   
   > \author Univ. of Colorado Denver   
   > \author NAG Ltd.   

   > \date November 2011   

   > \ingroup doubleGEcomputational   

    =====================================================================   
   Subroutine */ int igraphdgetrs_(char *trans, integer *n, integer *nrhs, 
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    extern logical igraphlsame_(char *, char *);
    extern /* Subroutine */ int igraphdtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *), igraphxerbla_(
	    char *, integer *, ftnlen), igraphdlaswp_(integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, integer *);
    logical notran;


/*  -- LAPACK computational routine (version 3.4.0) --   
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
    notran = igraphlsame_(trans, "N");
    if (! notran && ! igraphlsame_(trans, "T") && ! igraphlsame_(
	    trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < max(1,*n)) {
	*info = -5;
    } else if (*ldb < max(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	igraphxerbla_("DGETRS", &i__1, (ftnlen)6);
	return 0;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return 0;
    }

    if (notran) {

/*        Solve A * X = B.   

          Apply row interchanges to the right hand sides. */

	igraphdlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);

/*        Solve L*X = B, overwriting B with X. */

	igraphdtrsm_("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

/*        Solve U*X = B, overwriting B with X. */

	igraphdtrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, &
		a[a_offset], lda, &b[b_offset], ldb);
    } else {

/*        Solve A**T * X = B.   

          Solve U**T *X = B, overwriting B with X. */

	igraphdtrsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

/*        Solve L**T *X = B, overwriting B with X. */

	igraphdtrsm_("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &a[
		a_offset], lda, &b[b_offset], ldb);

/*        Apply row interchanges to the solution vectors. */

	igraphdlaswp_(nrhs, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
    }

    return 0;

/*     End of DGETRS */

} /* igraphdgetrs_ */

