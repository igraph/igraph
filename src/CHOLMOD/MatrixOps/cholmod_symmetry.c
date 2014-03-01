/* ========================================================================== */
/* === MatrixOps/cholmod_symmetry =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Determines if a sparse matrix is rectangular, unsymmetric, symmetric,
 * skew-symmetric, or Hermitian.  It does so by looking at its numerical values
 * of both upper and lower triangular parts of a CHOLMOD "unsymmetric"
 * matrix, where A->stype == 0.  The transpose of A is NOT constructed.
 *
 * If not unsymmetric, it also determines if the matrix has a diagonal whose
 * entries are all real and positive (and thus a candidate for sparse Cholesky
 * if A->stype is changed to a nonzero value).
 *
 * Note that a Matrix Market "general" matrix is either rectangular or
 * unsymmetric.
 *
 * The row indices in the column of each matrix MUST be sorted for this function
 * to work properly (A->sorted must be TRUE).  This routine returns EMPTY if
 * A->stype is not zero, or if A->sorted is FALSE.  The exception to this rule
 * is if A is rectangular.
 *
 * If option == 0, then this routine returns immediately when it finds a
 * non-positive diagonal entry (or one with nonzero imaginary part).   If the
 * matrix is not a candidate for sparse Cholesky, it returns the value
 * CHOLMOD_MM_UNSYMMETRIC, even if the matrix might in fact be symmetric or
 * Hermitian.
 *
 * This routine is useful inside the MATLAB backslash, which must look at an
 * arbitrary matrix (A->stype == 0) and determine if it is a candidate for
 * sparse Cholesky.  In that case, option should be 0.
 *
 * This routine is also useful when writing a MATLAB matrix to a file in
 * Rutherford/Boeing or Matrix Market format.  Those formats require a
 * determination as to the symmetry of the matrix, and thus this routine should
 * not return upon encountering the first non-positive diagonal.  In this case,
 * option should be 1.
 *
 * If option is 2, this function can be used to compute the numerical and
 * pattern symmetry, where 0 is a completely unsymmetric matrix, and 1 is a
 * perfectly symmetric matrix.  This option is used when computing the following
 * statistics for the matrices in the UF Sparse Matrix Collection.
 *
 *	numerical symmetry: number of matched offdiagonal nonzeros over
 *	the total number of offdiagonal entries.  A real entry A(i,j), i ~= j,
 *	is matched if A (j,i) == A (i,j), but this is only counted if both
 *	A(j,i) and A(i,j) are nonzero.  This does not depend on Z.
 *	(If A is complex, then the above test is modified; A (i,j) is matched
 *	if conj (A (j,i)) == A (i,j)).
 *
 *	Then numeric symmetry = xmatched / nzoffdiag, or 1 if nzoffdiag = 0.
 *  
 *	pattern symmetry: number of matched offdiagonal entries over the
 *	total number of offdiagonal entries.  An entry A(i,j), i ~= j, is
 *	matched if A (j,i) is also an entry.
 *
 *	Then pattern symmetry = pmatched / nzoffdiag, or 1 if nzoffdiag = 0.
 *  
 * The symmetry of a matrix with no offdiagonal entries is equal to 1.
 *
 * A workspace of size ncol integers is allocated; EMPTY is returned if this
 * allocation fails.
 *
 * Summary of return values:
 *
 *  EMPTY (-1)			    out of memory, stype not zero, A not sorted
 *  CHOLMOD_MM_RECTANGULAR 1	    A is rectangular
 *  CHOLMOD_MM_UNSYMMETRIC 2	    A is unsymmetric
 *  CHOLMOD_MM_SYMMETRIC 3	    A is symmetric, but with non-pos. diagonal
 *  CHOLMOD_MM_HERMITIAN 4	    A is Hermitian, but with non-pos. diagonal
 *  CHOLMOD_MM_SKEW_SYMMETRIC 5	    A is skew symmetric
 *  CHOLMOD_MM_SYMMETRIC_POSDIAG 6  A is symmetric with positive diagonal
 *  CHOLMOD_MM_HERMITIAN_POSDIAG 7  A is Hermitian with positive diagonal
 *
 * See also the spsym mexFunction, which is a MATLAB interface for this code.
 *
 * If the matrix is a candidate for sparse Cholesky, it will return a result
 * CHOLMOD_MM_SYMMETRIC_POSDIAG if real, or CHOLMOD_MM_HERMITIAN_POSDIAG if
 * complex.  Otherwise, it will return a value less than this.  This is true
 * regardless of the value of the option parameter.
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"


/* ========================================================================== */
/* === get_value ============================================================ */
/* ========================================================================== */

/* Get the pth value in the matrix. */

static void get_value
(
    double *Ax,	    /* real values, or real/imag. for CHOLMOD_COMPLEX type */
    double *Az,	    /* imaginary values for CHOLMOD_ZOMPLEX type */
    Int p,	    /* get the pth entry */
    Int xtype,	    /* A->xtype: pattern, real, complex, or zomplex */
    double *x,	    /* the real part */
    double *z	    /* the imaginary part */
)
{
    switch (xtype)
    {
	case CHOLMOD_PATTERN:
	    *x = 1 ;
	    *z = 0 ;
	    break ;

	case CHOLMOD_REAL:
	    *x = Ax [p] ;
	    *z = 0 ;
	    break ;

	case CHOLMOD_COMPLEX:
	    *x = Ax [2*p] ;
	    *z = Ax [2*p+1] ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    *x = Ax [p] ;
	    *z = Az [p] ;
	    break ;
    }
}


/* ========================================================================== */
/* === cholmod_symmetry ===================================================== */
/* ========================================================================== */

/* Determine the symmetry of a matrix, and check its diagonal.
 *
 * option 0:  Do not count # of matched pairs.  Quick return if the
 *	      the matrix has a zero, negative, or imaginary diagonal entry.
 *
 * option 1:  Do not count # of matched pairs.  Do not return quickly if
 *	      the matrix has a zero, negative, or imaginary diagonal entry.
 *	The result 1 to 7 is accurately computed:
 *
 *	EMPTY (-1)		out of memory, stype not zero, A not sorted
 *	CHOLMOD_MM_RECTANGULAR 1	A is rectangular
 *	CHOLMOD_MM_UNSYMMETRIC 2	A is unsymmetric
 *	CHOLMOD_MM_SYMMETRIC 3		A is symmetric, with non-pos. diagonal
 *	CHOLMOD_MM_HERMITIAN 4		A is Hermitian, with non-pos. diagonal
 *	CHOLMOD_MM_SKEW_SYMMETRIC 5	A is skew symmetric
 *	CHOLMOD_MM_SYMMETRIC_POSDIAG 6  is symmetric with positive diagonal
 *	CHOLMOD_MM_HERMITIAN_POSDIAG 7  A is Hermitian with positive diagonal
 *
 *	The routine returns as soon as the above is determined (that is, it
 *	can return as soon as it determines the matrix is unsymmetric).
 *
 * option 2:  All of the above, but also compute the number of matched off-
 *	diagonal entries (of two types).  xmatched is the number of 
 *	nonzero entries for which A(i,j) = conj(A(j,i)).  pmatched is
 *	the number of entries (i,j) for which A(i,j) and A(j,i) are both in
 *	the pattern of A (the value doesn't matter).  nzoffdiag is the total
 *	number of off-diagonal entries in the pattern.  nzdiag is the number of
 *	diagonal entries in the pattern.
 *
 * With option 0 or 1, or if the matrix is rectangular, xmatched, pmatched,
 * nzoffdiag, and nzdiag are not computed.
 *
 * Note that a matched pair, A(i,j) and A(j,i) for i != j, is counted twice
 * (once per entry).
 */

int CHOLMOD(symmetry)
(
    /* ---- input ---- */
    cholmod_sparse *A,
    int option,			/* option 0, 1, or 2 (see above) */
    /* ---- output --- */	/* outputs ignored if any are NULL */
    Int *p_xmatched,		/* # of matched numerical entries */
    Int *p_pmatched,		/* # of matched entries in pattern */
    Int *p_nzoffdiag,		/* # of off diagonal entries */
    Int *p_nzdiag,		/* # of diagonal entries */
    /* --------------- */
    cholmod_common *Common
)
{
    double aij_real = 0, aij_imag = 0, aji_real = 0, aji_imag = 0 ;
    double *Ax, *Az ;
    Int *Ap, *Ai, *Anz, *munch ;
    Int packed, nrow, ncol, xtype, is_symmetric, is_skew, is_hermitian, posdiag,
	j, p, pend, i, piend, result, xmatched, pmatched, nzdiag, i2, found ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;
    ASSERT (CHOLMOD(dump_sparse) (A, "cholmod_symmetry", Common) >= 0) ;

    if (p_xmatched == NULL || p_pmatched == NULL
	|| p_nzoffdiag == NULL || p_nzdiag == NULL)
    {
	/* option 2 is not performed if any output parameter is NULL */
	option = MAX (option, 1) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    packed = A->packed ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    xtype = A->xtype ;

    /* ---------------------------------------------------------------------- */
    /* check if rectangular, unsorted, or stype is not zero */
    /* ---------------------------------------------------------------------- */

    if (nrow != ncol)
    {
	/* matrix is rectangular */
	return (CHOLMOD_MM_RECTANGULAR) ;
    }

    if (!(A->sorted) || A->stype != 0)
    {
	/* this function cannot determine the type or symmetry */
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* this function requires uninitialized Int workspace of size ncol */
    CHOLMOD(allocate_work) (0, ncol, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (EMPTY) ;
    }

    munch = Common->Iwork ;	    /* the munch array is size ncol */

    /* ---------------------------------------------------------------------- */
    /* determine symmetry of a square matrix */
    /* ---------------------------------------------------------------------- */

    /* a complex or zomplex matrix is Hermitian until proven otherwise */
    is_hermitian = (xtype >= CHOLMOD_COMPLEX) ;

    /* any matrix is symmetric until proven otherwise */
    is_symmetric = TRUE ;

    /* a non-pattern matrix is skew-symmetric until proven otherwise */
    is_skew = (xtype != CHOLMOD_PATTERN) ;

    /* a matrix has positive diagonal entries until proven otherwise */
    posdiag = TRUE ;

    /* munch pointers start at the top of each column */
    for (j = 0 ; j < ncol ; j++)
    {
	munch [j] = Ap [j] ;
    }

    xmatched = 0 ;
    pmatched = 0 ;
    nzdiag = 0 ;

    for (j = 0 ; j < ncol ; j++)	/* examine each column of A */
    {

	/* ------------------------------------------------------------------ */
	/* look at the entire munch column j */
	/* ------------------------------------------------------------------ */

	/* start at the munch point of column j, and go to end of the column */
	p = munch [j] ;
	pend = (packed) ? (Ap [j+1]) : (Ap [j] + Anz [j]) ;

	for ( ; p < pend ; p++)
	{
	    /* get the row index of A(i,j) */
	    i = Ai [p] ;

	    if (i < j)
	    {

		/* ---------------------------------------------------------- */
		/* A(i,j) in triu(A), but matching A(j,i) not in tril(A) */
		/* ---------------------------------------------------------- */

		/* entry A(i,j) is unmatched; it appears in the upper triangular
		 * part, but not the lower triangular part.  The matrix is
		 * unsymmetric. */
		is_hermitian = FALSE ;
		is_symmetric = FALSE ;
		is_skew = FALSE ;

	    }
	    else if (i == j)
	    {

		/* ---------------------------------------------------------- */
		/* the diagonal A(j,j) is present; check its value */
		/* ---------------------------------------------------------- */

		get_value (Ax, Az, p, xtype, &aij_real, &aij_imag) ;
		if (aij_real != 0. || aij_imag != 0.)
		{
		    /* diagonal is nonzero; matrix is not skew-symmetric */
		    nzdiag++ ;
		    is_skew = FALSE ;
		}
		if (aij_real <= 0. || aij_imag != 0.)
		{
		    /* diagonal negative or imaginary; not chol candidate */
		    posdiag = FALSE ;
		}
		if (aij_imag != 0.)
		{
		    /* imaginary part is present; not Hermitian */
		    is_hermitian = FALSE ;
		}

	    }
	    else /* i > j */
	    {

		/* ---------------------------------------------------------- */
		/* consider column i, up to and including row j */
		/* ---------------------------------------------------------- */

		/* munch the entry at top of column i up to and incl row j */
		piend = (packed) ? (Ap [i+1]) : (Ap [i] + Anz [i]) ;

		found = FALSE ;

		for ( ; munch [i] < piend ; munch [i]++)
		{

		    i2 = Ai [munch [i]] ;

		    if (i2 < j)
		    {

			/* -------------------------------------------------- */
			/* A(i2,i) in triu(A) but A(i,i2) not in tril(A) */
			/* -------------------------------------------------- */

			/* The matrix is unsymmetric. */
			is_hermitian = FALSE ;
			is_symmetric = FALSE ;
			is_skew = FALSE ;

		    }
		    else if (i2 == j)
		    {

			/* -------------------------------------------------- */
			/* both A(i,j) and A(j,i) exist in the matrix */
			/* -------------------------------------------------- */

			/* this is one more matching entry in the pattern */
			pmatched += 2 ;
			found = TRUE ;

			/* get the value of A(i,j) */
			get_value (Ax, Az, p, xtype, &aij_real, &aij_imag) ;

			/* get the value of A(j,i) */
			get_value (Ax, Az, munch [i],
			    xtype, &aji_real, &aji_imag) ;

			/* compare A(i,j) with A(j,i) */
			if (aij_real != aji_real || aij_imag != aji_imag)
			{
			    /* the matrix cannot be symmetric */
			    is_symmetric = FALSE ;
			}
			if (aij_real != -aji_real || aij_imag != aji_imag)
			{
			    /* the matrix cannot be skew-symmetric */
			    is_skew = FALSE ;
			}
			if (aij_real != aji_real || aij_imag != -aji_imag)
			{
			    /* the matrix cannot be Hermitian */
			    is_hermitian = FALSE ;
			}
			else
			{
			    /* A(i,j) and A(j,i) are numerically matched */
			    xmatched += 2 ;
			}

		    }
		    else /* i2 > j */
		    {

			/* -------------------------------------------------- */
			/* entry A(i2,i) is not munched; consider it later */
			/* -------------------------------------------------- */

			break ;
		    }
		}

		if (!found)
		{
		    /* A(i,j) in tril(A) but A(j,i) not in triu(A).
		     * The matrix is unsymmetric. */
		    is_hermitian = FALSE ;
		    is_symmetric = FALSE ;
		    is_skew = FALSE ;
		}
	    }

	    if (option < 2 && !(is_symmetric || is_skew || is_hermitian))
	    {
		/* matrix is unsymmetric; terminate the test */
		return (CHOLMOD_MM_UNSYMMETRIC) ;
	    }
	}

	/* ------------------------------------------------------------------ */
	/* quick return if not Cholesky candidate */
	/* ------------------------------------------------------------------ */

	if (option < 1 && (!posdiag || nzdiag < ncol))
	{
	    /* Diagonal entry not present, or present but negative or with
	     * nonzero imaginary part.  Quick return for option 0. */
	    return (CHOLMOD_MM_UNSYMMETRIC) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the results */
    /* ---------------------------------------------------------------------- */

    if (nzdiag < ncol)
    {
        /* not all diagonal entries are present */
        posdiag = FALSE ;
    }

    if (option >= 2)
    {
	*p_xmatched = xmatched ;
	*p_pmatched = pmatched ;
	*p_nzoffdiag = CHOLMOD(nnz) (A, Common) - nzdiag ;
	*p_nzdiag = nzdiag ;
    }

    result = CHOLMOD_MM_UNSYMMETRIC ;
    if (is_hermitian)
    {
	/* complex Hermitian matrix, with either pos. or non-pos. diagonal */
	result = posdiag ? CHOLMOD_MM_HERMITIAN_POSDIAG : CHOLMOD_MM_HERMITIAN ;
    }
    else if (is_symmetric)
    {
	/* real or complex symmetric matrix, with pos. or non-pos. diagonal */
	result = posdiag ? CHOLMOD_MM_SYMMETRIC_POSDIAG : CHOLMOD_MM_SYMMETRIC ;
    }
    else if (is_skew)
    {
	/* real or complex skew-symmetric matrix */
	result = CHOLMOD_MM_SKEW_SYMMETRIC ;
    }
    return (result) ;
}
#endif
