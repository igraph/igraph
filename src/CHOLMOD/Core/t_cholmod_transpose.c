/* ========================================================================== */
/* === Core/t_cholmod_transpose ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Template routine for cholmod_transpose.  All xtypes are supported.  For
 * complex matrices, either the array tranpose or complex conjugate transpose
 * can be computed. */

#include "cholmod_template.h"

/* ========================================================================== */
/* === t_cholmod_transpose_unsym ============================================ */
/* ========================================================================== */

/* Compute F = A', A (:,f)', or A (p,f)', where A is unsymmetric and F is
 * already allocated.  The complex case performs either the array transpose
 * or complex conjugate transpose.
 *
 * workspace:
 * Iwork (MAX (nrow,ncol)) if fset is present
 * Iwork (nrow) if fset is NULL
 */

static int TEMPLATE (cholmod_transpose_unsym)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to transpose */
    Int *Perm,		/* size nrow, if present (can be NULL) */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    Int nf,		/* size of fset */
    /* ---- output --- */
    cholmod_sparse *F,	/* F = A', A(:,f)', or A(p,f)' */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Az, *Fx, *Fz ;
    Int *Ap, *Anz, *Ai, *Fp, *Fnz, *Fj, *Wi, *Iwork ;
    Int j, p, pend, nrow, ncol, Apacked, use_fset, fp, Fpacked, jj, permute ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    /* ensure the xtype of A and F match (ignored if this is pattern version) */
    if (!XTYPE_OK (A->xtype))
    {
	ERROR (CHOLMOD_INVALID, "real/complex mismatch") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    use_fset = (fset != NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;

    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Ax = A->x ;		/* size nz, real values of A */
    Az = A->z ;		/* size nz, imag values of A */
    Anz = A->nz ;
    Apacked = A->packed ;
    ASSERT (IMPLIES (!Apacked, Anz != NULL)) ;

    permute = (Perm != NULL) ;

    Fp = F->p ;		/* size A->nrow+1, row pointers of F */
    Fj = F->i ;		/* size nz, column indices of F */
    Fx = F->x ;		/* size nz, real values of F */
    Fz = F->z ;		/* size nz, imag values of F */
    Fnz = F->nz ;
    Fpacked = F->packed ;
    ASSERT (IMPLIES (!Fpacked, Fnz != NULL)) ;

    nf = (use_fset) ? nf : ncol ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wi = Iwork ;		/* size nrow (i/l/l) */

    /* ---------------------------------------------------------------------- */
    /* construct the transpose */
    /* ---------------------------------------------------------------------- */

    for (jj = 0 ; jj < nf ; jj++)
    {
	j = (use_fset) ? (fset [jj]) : jj ;
	p = Ap [j] ;
	pend = (Apacked) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    fp = Wi [Ai [p]]++ ;
	    Fj [fp] = j ;
#ifdef NCONJUGATE
	    ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
#else
	    ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
#endif
	}
    }

    return (TRUE) ;
}


/* ========================================================================== */
/* === t_cholmod_transpose_sym ============================================== */
/* ========================================================================== */

/* Compute F = A' or A (p,p)', where A is symmetric and F is already allocated.
 * The complex case performs either the array transpose or complex conjugate
 * transpose.
 *
 * workspace:  Iwork (nrow) if Perm NULL, Iwork (2*nrow) if Perm non-NULL.
 */

static int TEMPLATE (cholmod_transpose_sym)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to transpose */
    Int *Perm,		/* size n, if present (can be NULL) */
    /* ---- output --- */
    cholmod_sparse *F,	/* F = A' or A(p,p)' */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Az, *Fx, *Fz ;
    Int *Ap, *Anz, *Ai, *Fp, *Fj, *Wi, *Pinv, *Iwork ;
    Int p, pend, packed, fp, upper, permute, jold, n, i, j, iold ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    /* ensure the xtype of A and F match (ignored if this is pattern version) */
    if (!XTYPE_OK (A->xtype))
    {
	ERROR (CHOLMOD_INVALID, "real/complex mismatch") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    permute = (Perm != NULL) ;
    n = A->nrow ;
    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Ax = A->x ;		/* size nz, real values of A */
    Az = A->z ;		/* size nz, imag values of A */
    Anz = A->nz ;
    packed = A->packed ;
    ASSERT (IMPLIES (!packed, Anz != NULL)) ;
    upper = (A->stype > 0) ;

    Fp = F->p ;		/* size A->nrow+1, row pointers of F */
    Fj = F->i ;		/* size nz, column indices of F */
    Fx = F->x ;		/* size nz, real values of F */
    Fz = F->z ;		/* size nz, imag values of F */

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wi = Iwork ;	/* size n (i/l/l) */
    Pinv = Iwork + n ;	/* size n (i/i/l) , unused if Perm NULL */

    /* ---------------------------------------------------------------------- */
    /* construct the transpose */
    /* ---------------------------------------------------------------------- */

    if (permute)
    {
	if (upper)
	{
	    /* permuted, upper */
	    for (j = 0 ; j < n ; j++)
	    {
		jold = Perm [j] ;
		p = Ap [jold] ;
		pend = (packed) ? Ap [jold+1] : p + Anz [jold] ;
		for ( ; p < pend ; p++)
		{
		    iold = Ai [p] ;
		    if (iold <= jold)
		    {
			i = Pinv [iold] ;
			if (i < j)
			{
			    fp = Wi [i]++ ;
			    Fj [fp] = j ;
#ifdef NCONJUGATE
			    ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
#else
			    ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
#endif
			}
			else
			{
			    fp = Wi [j]++ ;
			    Fj [fp] = i ;
			    ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
			}
		    }
		}
	    }
	}
	else
	{
	    /* permuted, lower */
	    for (j = 0 ; j < n ; j++)
	    {
		jold = Perm [j] ;
		p = Ap [jold] ;
		pend = (packed) ? Ap [jold+1] : p + Anz [jold] ;
		for ( ; p < pend ; p++)
		{
		    iold = Ai [p] ;
		    if (iold >= jold)
		    {
			i = Pinv [iold] ;
			if (i > j)
			{
			    fp = Wi [i]++ ;
			    Fj [fp] = j ;
#ifdef NCONJUGATE
			    ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
#else
			    ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
#endif
			}
			else
			{
			    fp = Wi [j]++ ;
			    Fj [fp] = i ;
			    ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
			}
		    }
		}
	    }
	}
    }
    else
    {
	if (upper)
	{
	    /* unpermuted, upper */
	    for (j = 0 ; j < n ; j++)
	    {
		p = Ap [j] ;
		pend = (packed) ? Ap [j+1] : p + Anz [j] ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i <= j)
		    {
			fp = Wi [i]++ ;
			Fj [fp] = j ;
#ifdef NCONJUGATE
			ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
#else
			ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
#endif
		    }
		}
	    }
	}
	else
	{
	    /* unpermuted, lower */
	    for (j = 0 ; j < n ; j++)
	    {
		p = Ap [j] ;
		pend = (packed) ? Ap [j+1] : p + Anz [j] ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= j)
		    {
			fp = Wi [i]++ ;
			Fj [fp] = j ;
#ifdef NCONJUGATE
			ASSIGN (Fx, Fz, fp, Ax, Az, p) ;
#else
			ASSIGN_CONJ (Fx, Fz, fp, Ax, Az, p) ;
#endif
		    }
		}
	    }
	}
    }

    return (TRUE) ;
}

#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
#undef NCONJUGATE
