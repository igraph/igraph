/* ========================================================================== */
/* === MatrixOps/cholmod_drop =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Drop small entries from A, and entries in the ignored part of A if A
 * is symmetric.  None of the matrix operations drop small numerical entries
 * from a matrix, except for this one.  NaN's and Inf's are kept.
 *
 * workspace: none
 *
 * Supports pattern and real matrices, complex and zomplex not supported.
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"


/* ========================================================================== */
/* === cholmod_drop ========================================================= */
/* ========================================================================== */

int CHOLMOD(drop)
(
    /* ---- input ---- */
    double tol,		/* keep entries with absolute value > tol */
    /* ---- in/out --- */
    cholmod_sparse *A,	/* matrix to drop entries from */
    /* --------------- */
    cholmod_common *Common
)
{
    double aij ;
    double *Ax ;
    Int *Ap, *Ai, *Anz ;
    Int packed, i, j, nrow, ncol, p, pend, nz, values ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_REAL, FALSE) ;
    Common->status = CHOLMOD_OK ;
    ASSERT (CHOLMOD(dump_sparse) (A, "A predrop", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Anz = A->nz ;
    packed = A->packed ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    values = (A->xtype != CHOLMOD_PATTERN) ;
    nz = 0 ;

    if (values)
    {

	/* ------------------------------------------------------------------ */
	/* drop small numerical entries from A, and entries in ignored part */
	/* ------------------------------------------------------------------ */

	if (A->stype > 0)
	{

	    /* -------------------------------------------------------------- */
	    /* A is symmetric, with just upper triangular part stored */
	    /* -------------------------------------------------------------- */

	    for (j = 0 ; j < ncol ; j++)
	    {
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Ap [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    aij = Ax [p] ;
		    if (i <= j && (fabs (aij) > tol || IS_NAN (aij)))
		    {
			Ai [nz] = i ;
			Ax [nz] = aij ;
			nz++ ;
		    }
		}
	    }

	}
	else if (A->stype < 0)
	{

	    /* -------------------------------------------------------------- */
	    /* A is symmetric, with just lower triangular part stored */
	    /* -------------------------------------------------------------- */

	    for (j = 0 ; j < ncol ; j++)
	    {
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Ap [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    aij = Ax [p] ;
		    if (i >= j && (fabs (aij) > tol || IS_NAN (aij)))
		    {
			Ai [nz] = i ;
			Ax [nz] = aij ;
			nz++ ;
		    }
		}
	    }
	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* both parts of A present, just drop small entries */
	    /* -------------------------------------------------------------- */

	    for (j = 0 ; j < ncol ; j++)
	    {
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Ap [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    aij = Ax [p] ;
		    if (fabs (aij) > tol || IS_NAN (aij))
		    {
			Ai [nz] = i ;
			Ax [nz] = aij ;
			nz++ ;
		    }
		}
	    }
	}
	Ap [ncol] = nz ;

	/* reduce A->i and A->x in size */
	ASSERT (MAX (1,nz) <= A->nzmax) ;
	CHOLMOD(reallocate_sparse) (nz, A, Common) ;
	ASSERT (Common->status >= CHOLMOD_OK) ;

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* consider only the pattern of A */
	/* ------------------------------------------------------------------ */

	/* Note that cholmod_band_inplace calls cholmod_reallocate_sparse */
	if (A->stype > 0)
	{
	    CHOLMOD(band_inplace) (0, ncol, 0, A, Common) ;
	}
	else if (A->stype < 0)
	{
	    CHOLMOD(band_inplace) (-nrow, 0, 0, A, Common) ;
	}
    }

    ASSERT (CHOLMOD(dump_sparse) (A, "A dropped", Common) >= 0) ;
    return (TRUE) ;
}
#endif
