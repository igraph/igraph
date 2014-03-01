/* ========================================================================== */
/* === Cholesky/cholmod_rcond =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Return a rough estimate of the reciprocal of the condition number:
 * the minimum entry on the diagonal of L (or absolute entry of D for an LDL'
 * factorization) divided by the maximum entry (squared for LL').  L can be
 * real, complex, or zomplex.  Returns -1 on error, 0 if the matrix is singular
 * or has a zero entry on the diagonal of L, 1 if the matrix is 0-by-0, or
 * min(diag(L))/max(diag(L)) otherwise.  Never returns NaN; if L has a NaN on
 * the diagonal it returns zero instead.
 *
 * For an LL' factorization,  (min(diag(L))/max(diag(L)))^2 is returned.
 * For an LDL' factorization, (min(diag(D))/max(diag(D))) is returned.
 */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"

/* ========================================================================== */
/* === LMINMAX ============================================================== */
/* ========================================================================== */

/* Update lmin and lmax for one entry L(j,j) */

#define FIRST_LMINMAX(Ljj,lmin,lmax) \
{ \
    double ljj = Ljj ; \
    if (IS_NAN (ljj)) \
    { \
	return (0) ; \
    } \
    lmin = ljj ; \
    lmax = ljj ; \
}

#define LMINMAX(Ljj,lmin,lmax) \
{ \
    double ljj = Ljj ; \
    if (IS_NAN (ljj)) \
    { \
	return (0) ; \
    } \
    if (ljj < lmin) \
    { \
	lmin = ljj ; \
    } \
    else if (ljj > lmax) \
    { \
	lmax = ljj ; \
    } \
}

/* ========================================================================== */
/* === cholmod_rcond ======================================================== */
/* ========================================================================== */

double CHOLMOD(rcond)	    /* return min(diag(L)) / max(diag(L)) */
(
    /* ---- input ---- */
    cholmod_factor *L,
    /* --------------- */
    cholmod_common *Common
)
{
    double lmin, lmax, rcond ;
    double *Lx ;
    Int *Lpi, *Lpx, *Super, *Lp ;
    Int n, e, nsuper, s, k1, k2, psi, psend, psx, nsrow, nscol, jj, j ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (L, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    n = L->n ;
    if (n == 0)
    {
	return (1) ;
    }
    if (L->minor < L->n)
    {
	return (0) ;
    }

    e = (L->xtype == CHOLMOD_COMPLEX) ? 2 : 1 ;

    if (L->is_super)
    {
	/* L is supernodal */
	nsuper = L->nsuper ;	/* number of supernodes in L */
	Lpi = L->pi ;		/* column pointers for integer pattern */
	Lpx = L->px ;		/* column pointers for numeric values */
	Super = L->super ;	/* supernode sizes */
	Lx = L->x ;		/* numeric values */
	FIRST_LMINMAX (Lx [0], lmin, lmax) ;	/* first diagonal entry of L */
	for (s = 0 ; s < nsuper ; s++)
	{
	    k1 = Super [s] ;		/* first column in supernode s */
	    k2 = Super [s+1] ;		/* last column in supernode is k2-1 */
	    psi = Lpi [s] ;		/* first row index is L->s [psi] */
	    psend = Lpi [s+1] ;		/* last row index is L->s [psend-1] */
	    psx = Lpx [s] ;		/* first numeric entry is Lx [psx] */
	    nsrow = psend - psi ;	/* supernode is nsrow-by-nscol */
	    nscol = k2 - k1 ;
	    for (jj = 0 ; jj < nscol ; jj++)
	    {
		LMINMAX (Lx [e * (psx + jj + jj*nsrow)], lmin, lmax) ;
	    }
	}
    }
    else
    {
	/* L is simplicial */
	Lp = L->p ;
	Lx = L->x ;
	if (L->is_ll)
	{
	    /* LL' factorization */
	    FIRST_LMINMAX (Lx [Lp [0]], lmin, lmax) ;
	    for (j = 1 ; j < n ; j++)
	    {
		LMINMAX (Lx [e * Lp [j]], lmin, lmax) ;
	    }
	}
	else
	{
	    /* LDL' factorization, the diagonal might be negative */
	    FIRST_LMINMAX (fabs (Lx [Lp [0]]), lmin, lmax) ;
	    for (j = 1 ; j < n ; j++)
	    {
		LMINMAX (fabs (Lx [e * Lp [j]]), lmin, lmax) ;
	    }
	}
    }
    rcond = lmin / lmax ;
    if (L->is_ll)
    {
	rcond = rcond*rcond ;
    }
    return (rcond) ;
}
#endif
