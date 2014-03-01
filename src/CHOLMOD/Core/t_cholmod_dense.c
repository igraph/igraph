/* ========================================================================== */
/* === Core/t_cholmod_dense ================================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Template routine for cholmod_dense.  All xtypes supported, except that there
 * are no dense matrices with an xtype of pattern. */

#include "cholmod_template.h"

/* ========================================================================== */
/* === t_cholmod_sparse_to_dense ============================================ */
/* ========================================================================== */

static cholmod_dense *TEMPLATE (cholmod_sparse_to_dense)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to copy */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Xx, *Az, *Xz ;
    Int *Ap, *Ai, *Anz ;
    cholmod_dense *X ;
    Int i, j, p, pend, nrow, ncol, packed ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;
    packed = A->packed ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;

    /* ---------------------------------------------------------------------- */
    /* allocate result */
    /* ---------------------------------------------------------------------- */

    X = CHOLMOD(zeros) (nrow, ncol, XTYPE2, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    Xx = X->x ;
    Xz = X->z ;

    /* ---------------------------------------------------------------------- */
    /* copy into dense matrix */
    /* ---------------------------------------------------------------------- */

    if (A->stype < 0)
    {
	/* A is symmetric with lower stored, but both parts of X are present */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i >= j)
		{
		    ASSIGN2 (Xx, Xz, i+j*nrow, Ax, Az, p) ;
		    ASSIGN2_CONJ (Xx, Xz, j+i*nrow, Ax, Az, p) ;
		}
	    }
	}
    }
    else if (A->stype > 0)
    {
	/* A is symmetric with upper stored, but both parts of X are present */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i <= j)
		{
		    ASSIGN2 (Xx, Xz, i+j*nrow, Ax, Az, p) ;
		    ASSIGN2_CONJ (Xx, Xz, j+i*nrow, Ax, Az, p) ;
		}
	    }
	}
    }
    else
    {
	/* both parts of A and X are present */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		ASSIGN2 (Xx, Xz, i+j*nrow, Ax, Az, p) ;
	    }
	}
    }

    return (X) ;
}


#ifndef PATTERN

/* There are no dense matrices of xtype CHOLMOD_PATTERN */

/* ========================================================================== */
/* === t_cholmod_dense_to_sparse ============================================ */
/* ========================================================================== */

static cholmod_sparse *TEMPLATE (cholmod_dense_to_sparse)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* matrix to copy */
    int values,		/* TRUE if values to be copied, FALSE otherwise */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Xx, *Cx, *Xz, *Cz ;
    Int *Ci, *Cp ;
    cholmod_sparse *C ;
    Int i, j, p, d, nrow, ncol, nz ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = X->nrow ;
    ncol = X->ncol ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;

    /* ---------------------------------------------------------------------- */
    /* count the number of nonzeros in the result */
    /* ---------------------------------------------------------------------- */

    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    if (ENTRY_IS_NONZERO (Xx, Xz, i+j*d))
	    {
		nz++ ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* allocate the result C */
    /* ---------------------------------------------------------------------- */

    C = CHOLMOD(allocate_sparse) (nrow, ncol, nz, TRUE, TRUE, 0,
	    values ? XTYPE : CHOLMOD_PATTERN, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;
    Cz = C->z ;

    /* ---------------------------------------------------------------------- */
    /* copy the dense matrix X into the sparse matrix C */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = p ;
	for (i = 0 ; i < nrow ; i++)
	{
	    if (ENTRY_IS_NONZERO (Xx, Xz, i+j*d))
	    {
		Ci [p] = i ;
		if (values)
		{
		    ASSIGN (Cx, Cz, p, Xx, Xz, i+j*d) ;
		}
		p++ ;
	    }
	}
    }
    ASSERT (p == nz) ;
    Cp [ncol] = nz ;

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (C, "C", Common) >= 0) ;
    return (C) ;
}


/* ========================================================================== */
/* === t_cholmod_copy_dense2 ================================================ */
/* ========================================================================== */

/* Y = X, where X and Y are both already allocated.  */

static int TEMPLATE (cholmod_copy_dense2)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* matrix to copy */
    /* ---- output --- */
    cholmod_dense *Y	/* copy of matrix X */
)
{
    double *Xx, *Xz, *Yx, *Yz ;
    Int i, j, nrow, ncol, dy, dx ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Xx = X->x ;
    Xz = X->z ;
    Yx = Y->x ;
    Yz = Y->z ;
    dx = X->d ;
    dy = Y->d ;
    nrow = X->nrow ;
    ncol = X->ncol ;

    /* ---------------------------------------------------------------------- */
    /* copy */
    /* ---------------------------------------------------------------------- */

    CLEAR (Yx, Yz, 0) ;
    for (j = 0 ; j < ncol ; j++)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    ASSIGN (Yx, Yz, i+j*dy, Xx, Xz, i+j*dx) ;
	}
    }
    return (TRUE) ;
}

#endif

#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
