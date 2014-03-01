/* ========================================================================== */
/* === MatrixOps/cholmod_norm =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* r = norm (A), compute the infinity-norm, 1-norm, or 2-norm of a sparse or
 * dense matrix.  Can compute the 2-norm only for a dense column vector.
 * Returns -1 if an error occurs.
 *
 * Pattern, real, complex, and zomplex sparse matrices are supported.
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"


/* ========================================================================== */
/* === abs_value ============================================================ */
/* ========================================================================== */

/* Compute the absolute value of a real, complex, or zomplex value */

static double abs_value
(
    int xtype,
    double *Ax,
    double *Az,
    Int p,
    cholmod_common *Common
)
{
    double s = 0 ;
    switch (xtype)
    {
	case CHOLMOD_PATTERN:
	    s = 1 ;
	    break ;

	case CHOLMOD_REAL:
	    s = fabs (Ax [p]) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    s = Common->hypotenuse (Ax [2*p], Ax [2*p+1]) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    s = Common->hypotenuse (Ax [p], Az [p]) ;
	    break ;
    }
    return (s) ;
}


/* ========================================================================== */
/* === cholmod_norm_dense =================================================== */
/* ========================================================================== */

double CHOLMOD(norm_dense)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* matrix to compute the norm of */
    int norm,		/* type of norm: 0: inf. norm, 1: 1-norm, 2: 2-norm */
    /* --------------- */
    cholmod_common *Common
)
{
    double xnorm, s, x, z ;
    double *Xx, *Xz, *W ;
    Int nrow, ncol, d, i, j, use_workspace, xtype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (X, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;
    ncol = X->ncol ;
    if (norm < 0 || norm > 2 || (norm == 2 && ncol > 1))
    {
	ERROR (CHOLMOD_INVALID, "invalid norm") ;
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = X->nrow ;
    d = X->d ;
    Xx = X->x ;
    Xz = X->z ;
    xtype = X->xtype ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace, if needed */
    /* ---------------------------------------------------------------------- */

    W = NULL ;
    use_workspace = (norm == 0 && ncol > 4) ;
    if (use_workspace)
    {
	CHOLMOD(allocate_work) (0, 0, nrow, Common) ;
	W = Common->Xwork ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* oops, no workspace */
	    use_workspace = FALSE ;
	}
    }


    /* ---------------------------------------------------------------------- */
    /* compute the norm */
    /* ---------------------------------------------------------------------- */

    xnorm = 0 ;

    if (use_workspace)
    {

	/* ------------------------------------------------------------------ */
	/* infinity-norm = max row sum, using stride-1 access of X */
	/* ------------------------------------------------------------------ */

	DEBUG (for (i = 0 ; i < nrow ; i++) ASSERT (W [i] == 0)) ;

	/* this is faster than stride-d, but requires O(nrow) workspace */
	for (j = 0 ; j < ncol ; j++)
	{
	    for (i = 0 ; i < nrow ; i++)
	    {
		W [i] += abs_value (xtype, Xx, Xz, i+j*d, Common) ;
	    }
	}
	for (i = 0 ; i < nrow ; i++)
	{
	    s = W [i] ;
	    if ((IS_NAN (s) || s > xnorm) && !IS_NAN (xnorm))
	    {
		xnorm = s ;
	    }
	    W [i] = 0 ;
	}

    }
    else if (norm == 0)
    {

	/* ------------------------------------------------------------------ */
	/* infinity-norm = max row sum, using stride-d access of X */
	/* ------------------------------------------------------------------ */

	for (i = 0 ; i < nrow ; i++)
	{
	    s = 0 ;
	    for (j = 0 ; j < ncol ; j++)
	    {
		s += abs_value (xtype, Xx, Xz, i+j*d, Common) ;
	    }
	    if ((IS_NAN (s) || s > xnorm) && !IS_NAN (xnorm))
	    {
		xnorm = s ;
	    }
	}

    }
    else if (norm == 1)
    {

	/* ------------------------------------------------------------------ */
	/* 1-norm = max column sum */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < ncol ; j++)
	{
	    s = 0 ;
	    for (i = 0 ; i < nrow ; i++)
	    {
		s += abs_value (xtype, Xx, Xz, i+j*d, Common) ;
	    }
	    if ((IS_NAN (s) || s > xnorm) && !IS_NAN (xnorm))
	    {
		xnorm = s ;
	    }
	}
    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* 2-norm = sqrt (sum (X.^2)) */
	/* ------------------------------------------------------------------ */

	switch (xtype)
	{

	    case CHOLMOD_REAL:
		for (i = 0 ; i < nrow ; i++)
		{
		    x = Xx [i] ;
		    xnorm += x*x ;
		}
		break ; 

	    case CHOLMOD_COMPLEX:
		for (i = 0 ; i < nrow ; i++)
		{
		    x = Xx [2*i  ] ;
		    z = Xx [2*i+1] ;
		    xnorm += x*x + z*z ;
		}
		break ; 

	    case CHOLMOD_ZOMPLEX:
		for (i = 0 ; i < nrow ; i++)
		{
		    x = Xx [i] ;
		    z = Xz [i] ;
		    xnorm += x*x + z*z ;
		}
		break ; 
	}

	xnorm = sqrt (xnorm) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    return (xnorm) ;
}


/* ========================================================================== */
/* === cholmod_norm_sparse ================================================== */
/* ========================================================================== */

double CHOLMOD(norm_sparse)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to compute the norm of */
    int norm,		/* type of norm: 0: inf. norm, 1: 1-norm */
    /* --------------- */
    cholmod_common *Common
)
{
    double anorm, s ;
    double *Ax, *Az, *W ;
    Int *Ap, *Ai, *Anz ;
    Int i, j, p, pend, nrow, ncol, packed, xtype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;
    ncol = A->ncol ;
    nrow = A->nrow ;
    if (norm < 0 || norm > 1)
    {
	ERROR (CHOLMOD_INVALID, "invalid norm") ;
	return (EMPTY) ;
    }
    if (A->stype && nrow != ncol)
    {
	ERROR (CHOLMOD_INVALID, "matrix invalid") ;
	return (EMPTY) ;
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
    xtype = A->xtype ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace, if needed */
    /* ---------------------------------------------------------------------- */

    W = NULL ;
    if (A->stype || norm == 0)
    {
	CHOLMOD(allocate_work) (0, 0, nrow, Common) ;
	W = Common->Xwork ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    return (EMPTY) ;
	}
	DEBUG (for (i = 0 ; i < nrow ; i++) ASSERT (W [i] == 0)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute the norm */
    /* ---------------------------------------------------------------------- */

    anorm = 0 ;

    if (A->stype > 0)
    {

	/* ------------------------------------------------------------------ */
	/* A is symmetric with upper triangular part stored */
	/* ------------------------------------------------------------------ */

	/* infinity-norm = 1-norm = max row/col sum */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		s = abs_value (xtype, Ax, Az, p, Common) ;
		if (i == j)
		{
		    W [i] += s ;
		}
		else if (i < j)
		{
		    W [i] += s ;
		    W [j] += s ;
		}
	    }
	}

    }
    else if (A->stype < 0)
    {

	/* ------------------------------------------------------------------ */
	/* A is symmetric with lower triangular part stored */
	/* ------------------------------------------------------------------ */

	/* infinity-norm = 1-norm = max row/col sum */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		s = abs_value (xtype, Ax, Az, p, Common) ;
		if (i == j)
		{
		    W [i] += s ;
		}
		else if (i > j)
		{
		    W [i] += s ;
		    W [j] += s ;
		}
	    }
	}

    }
    else if (norm == 0)
    {

	/* ------------------------------------------------------------------ */
	/* A is unsymmetric, compute the infinity-norm */
	/* ------------------------------------------------------------------ */

	/* infinity-norm = max row sum */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		W [Ai [p]] += abs_value (xtype, Ax, Az, p, Common) ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* A is unsymmetric, compute the 1-norm */
	/* ------------------------------------------------------------------ */

	/* 1-norm = max column sum */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    if (xtype == CHOLMOD_PATTERN)
	    {
		s = pend - p ;
	    }
	    else
	    {
		s = 0 ;
		for ( ; p < pend ; p++)
		{
		    s += abs_value (xtype, Ax, Az, p, Common) ;
		}
	    }
	    if ((IS_NAN (s) || s > anorm) && !IS_NAN (anorm))
	    {
		anorm = s ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* compute the max row sum */
    /* ---------------------------------------------------------------------- */

    if (A->stype || norm == 0)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    s = W [i] ;
	    if ((IS_NAN (s) || s > anorm) && !IS_NAN (anorm))
	    {
		anorm = s ;
	    }
	    W [i] = 0 ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    return (anorm) ;
}
#endif
