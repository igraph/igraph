/* ========================================================================== */
/* === Core/t_cholmod_change_factor ========================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Template routine for cholmod_change_factor.  All xtypes supported. */

#include "cholmod_template.h"

/* ========================================================================== */
/* === t_change_simplicial_numeric ========================================== */
/* ========================================================================== */

static void TEMPLATE (change_simplicial_numeric)
(
    cholmod_factor *L,
    Int to_ll,
    Int to_packed,
    Int *newLi,
    double *newLx,
    double *newLz,
    Int lnz,
    Int grow,
    double grow1,
    Int grow2,
    Int make_ll,
    Int make_monotonic,
    Int make_ldl,
    cholmod_common *Common
)
{
    double xlen, dj [1], ljj [1], lj2 [1] ;
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz ;
    Int n, j, len, pnew, pold, k, p, pend ;

    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lnz = L->nz ;

    if (make_ll)
    {
	L->minor = n ;
    }

    if (make_monotonic)
    {

	/* ------------------------------------------------------------------ */
	/* reorder the columns to make them monotonic */
	/* ------------------------------------------------------------------ */

	pnew = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    /* copy and pack column j */
	    len = Lnz [j] ;
	    PRINT2 (("j: "ID" Lnz[j] "ID" len "ID" p "ID"\n",
			j, Lnz [j], len, pnew)) ;
	    pold = Lp [j] ;
	    ASSERT (Li [pold] == j) ;

	    if (make_ll)
	    {

		/* ---------------------------------------------------------- */
		/* copy and convert LDL' to LL' */
		/* ---------------------------------------------------------- */

		/* dj = Lx [pold] ; */
		ASSIGN_REAL (dj,0, Lx,pold) ;

		if (IS_LE_ZERO (dj [0]))
		{
		    /* Conversion has failed; matrix is not positive definite.
		     * Do not modify the column so that the LDL' factorization
		     * can be restored if desired, by converting back to LDL'.
		     * Continue the conversion, but flag the error. */
		    if (L->minor == (size_t) n)
		    {
			ERROR (CHOLMOD_NOT_POSDEF, "L not positive definite") ;
			L->minor = j ;
		    }
		    for (k = 0 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (newLx, newLz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    ljj [0] = sqrt (dj [0]) ;
		    newLi [pnew] = j ;
		    /* newLx [pnew] = ljj ; */
		    ASSIGN_REAL (newLx, pnew, ljj, 0) ;
		    CLEAR_IMAG (newLx, newLz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] * ljj ; */
			MULT_REAL (newLx, newLz, pnew+k, Lx, Lz, pold+k, ljj,0);
		    }
		}

	    }
	    else if (make_ldl)
	    {

		/* ---------------------------------------------------------- */
		/* copy and convert LL' to LDL' */
		/* ---------------------------------------------------------- */

		/* ljj = Lx [pold] ; */
		ASSIGN_REAL (ljj, 0, Lx, pold) ;

		if (ljj [0] <= 0)
		{
		    /* matrix is not positive-definite; copy column as-is */
		    for (k = 0 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (newLx, newLz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    newLi [pnew] = j ;
		    /* newLx [pnew] = ljj*ljj ; */
		    lj2 [0] = ljj [0] * ljj [0] ;
		    ASSIGN_REAL (newLx, pnew, lj2, 0) ;
		    CLEAR_IMAG (newLx, newLz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			newLi [pnew + k] = Li [pold + k] ;
			/* newLx [pnew + k] = Lx [pold + k] / ljj ; */
			DIV_REAL (newLx, newLz, pnew+k, Lx, Lz, pold+k, ljj,0) ;
		    }
		}

	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* copy and leave LL' or LDL' as-is */
		/* ---------------------------------------------------------- */

		for (k = 0 ; k < len ; k++)
		{
		    newLi [pnew + k] = Li [pold + k] ;
		    /* newLx [pnew + k] = Lx [pold + k] ; */
		    ASSIGN (newLx, newLz, pnew+k, Lx, Lz, pold+k) ;
		}
	    }

	    Lp [j] = pnew ;

	    /* compute len in double to avoid integer overflow */
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (Int) xlen ;
	    }
	    ASSERT (len >= Lnz [j] && len <= n-j) ;
	    pnew += len ;
	    ASSERT (pnew > 0) ;	    /* integer overflow case already covered */
	}
	Lp [n] = pnew ;
	PRINT1 (("final pnew = "ID", lnz "ID" lnzmax %g\n",
		    pnew, lnz, (double) L->nzmax)) ;
	ASSERT (pnew <= lnz) ;

	/* free the old L->i and L->x and replace with the new ones */
	CHOLMOD(free) (L->nzmax, sizeof (Int), L->i, Common) ;

#ifdef REAL
	CHOLMOD(free) (L->nzmax, sizeof (double), L->x, Common) ;
#elif defined (COMPLEX)
	CHOLMOD(free) (L->nzmax, 2*sizeof (double), L->x, Common) ;
#else
	CHOLMOD(free) (L->nzmax, sizeof (double), L->x, Common) ;
	CHOLMOD(free) (L->nzmax, sizeof (double), L->z, Common) ;
#endif

	L->i = newLi ;
	L->x = newLx ;
	L->z = newLz ;
	L->nzmax = lnz ;

	/* reconstruct the link list */
	natural_list (L) ;

    }
    else if (to_packed)
    {

	/* ------------------------------------------------------------------ */
	/* already monotonic, just pack the columns of L */
	/* ------------------------------------------------------------------ */

	pnew = 0 ;

	if (make_ll)
	{

	    /* -------------------------------------------------------------- */
	    /* pack and convert LDL' to LL' */
	    /* -------------------------------------------------------------- */

	    for (j = 0 ; j < n ; j++)
	    {
		/* pack column j */
		pold = Lp [j] ;
		len = Lnz [j] ;
		ASSERT (len > 0) ;
		ASSERT (Li [pold] == j) ;
		PRINT2 (("col "ID" pnew "ID" pold "ID"\n", j, pnew, pold)) ;

		/* dj = Lx [pold] ; */
		ASSIGN_REAL (dj,0, Lx,pold) ;

		if (IS_LE_ZERO (dj [0]))
		{
		    /* Conversion has failed; matrix is not positive definite.
		     * Do not modify the column so that the LDL' factorization
		     * can be restored if desired, by converting back to LDL'.
		     * Continue the conversion, but flag the error. */
		    if (L->minor == (size_t) n)
		    {
			ERROR (CHOLMOD_NOT_POSDEF, "L not positive definite") ;
			L->minor = j ;
		    }
		    for (k = 0 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (Lx, Lz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    ljj [0] = sqrt (dj [0]) ;
		    Li [pnew] = j ;

		    /* Lx [pnew] = ljj ; */
		    ASSIGN_REAL (Lx, pnew, ljj, 0) ;
		    CLEAR_IMAG (Lx, Lz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] * ljj ; */
			MULT_REAL (Lx, Lz, pnew+k, Lx, Lz, pold+k, ljj,0) ;
		    }
		}
		Lp [j] = pnew ;
		pnew += len ;
	    }

	}
	else if (make_ldl)
	{

	    /* -------------------------------------------------------------- */
	    /* pack and convert LL' to LDL' */
	    /* -------------------------------------------------------------- */

	    for (j = 0 ; j < n ; j++)
	    {
		/* pack column j */
		pold = Lp [j] ;
		len = Lnz [j] ;

		/* ljj = Lx [pold] ; */
		ASSIGN_REAL (ljj, 0, Lx, pold) ;

		ASSERT (len > 0) ;
		PRINT2 (("col "ID" pnew "ID" pold "ID"\n", j, pnew, pold)) ;
		if (ljj [0] <= 0)
		{
		    /* matrix is not positive-definite; pack column as-is */
		    for (k = 0 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (Lx, Lz, pnew+k, Lx, Lz, pold+k) ;
		    }
		}
		else
		{
		    Li [pnew] = Li [pold] ;

		    /* Lx [pnew] = ljj*ljj ; */
		    lj2 [0] = ljj [0] * ljj [0] ;
		    ASSIGN_REAL (Lx, pnew, lj2, 0) ;
		    CLEAR_IMAG (Lx, Lz, pnew) ;

		    for (k = 1 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] / ljj ; */
			DIV_REAL (Lx, Lz, pnew+k, Lx, Lz, pold+k, ljj,0) ;
		    }
		}
		Lp [j] = pnew ;
		pnew += len ;
	    }

	}
	else
	{

	    /* ---------------------------------------------------------- */
	    /* pack and leave LL' or LDL' as-is */
	    /* ---------------------------------------------------------- */

	    for (j = 0 ; j < n ; j++)
	    {
		/* pack column j */
		pold = Lp [j] ;
		len = Lnz [j] ;
		ASSERT (len > 0) ;
		PRINT2 (("col "ID" pnew "ID" pold "ID"\n", j, pnew, pold)) ;
		if (pnew < pold)
		{
		    PRINT2 (("    pack this column\n")) ;
		    for (k = 0 ; k < len ; k++)
		    {
			Li [pnew + k] = Li [pold + k] ;
			/* Lx [pnew + k] = Lx [pold + k] ; */
			ASSIGN (Lx, Lz, pnew+k, Lx, Lz, pold+k) ;
		    }
		    Lp [j] = pnew ;
		}
		pnew += len ;
	    }
	}

	Lp [n] = pnew ;
	PRINT2 (("Lp [n] = "ID"\n", pnew)) ;

    }
    else if (make_ll)
    {

	/* ------------------------------------------------------------------ */
	/* convert LDL' to LL', but do so in-place */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < n ; j++)
	{
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;

	    /* dj = Lx [p] ; */
	    ASSIGN_REAL (dj,0, Lx,p) ;

	    if (IS_LE_ZERO (dj [0]))
	    {
		/* Conversion has failed; matrix is not positive definite.
		 * Do not modify the column so that the LDL' factorization
		 * can be restored if desired, by converting back to LDL'.
		 * Continue the conversion, but flag the error. */
		if (L->minor == (size_t) n)
		{
		    ERROR (CHOLMOD_NOT_POSDEF, "L not positive definite") ;
		    L->minor = j ;
		}
	    }
	    else
	    {
		ljj [0] = sqrt (dj [0]) ;
		/* Lx [p] = ljj ; */
		ASSIGN_REAL (Lx,p, ljj,0) ;
		CLEAR_IMAG (Lx, Lz, p) ;

		for (p++ ; p < pend ; p++)
		{
		    /* Lx [p] *= ljj ; */
		    MULT_REAL (Lx,Lz,p, Lx,Lz,p, ljj,0) ;
		}
	    }
	}

    }
    else if (make_ldl)
    {

	/* ------------------------------------------------------------------ */
	/* convert LL' to LDL', but do so in-place */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < n ; j++)
	{
	    p = Lp [j] ;
	    pend = p + Lnz [j] ;

	    /* ljj = Lx [p] ; */
	    ASSIGN_REAL (ljj, 0, Lx, p) ;

	    if (ljj [0] > 0)
	    {
		/* Lx [p] = ljj*ljj ; */
		lj2 [0] = ljj [0] * ljj [0] ;
		ASSIGN_REAL (Lx, p, lj2, 0) ;
		CLEAR_IMAG (Lx, Lz, p) ;

		for (p++ ; p < pend ; p++)
		{
		    /* Lx [p] /= ljj ; */
		    DIV_REAL (Lx,Lz,p, Lx,Lz,p, ljj,0) ;
		}
	    }
	}
    }

    L->is_ll = to_ll ;

    DEBUG (CHOLMOD(dump_factor) (L, "done change simplicial numeric", Common)) ;
}


/* ========================================================================== */
/* === t_ll_super_to_simplicial_numeric ===================================== */
/* ========================================================================== */

/* A supernodal L can only be real or complex, not zomplex */

#ifndef ZOMPLEX

static void TEMPLATE (ll_super_to_simplicial_numeric)
(
    cholmod_factor *L,
    Int to_packed,
    Int to_ll,
    cholmod_common *Common
)
{
    double ljj [1], lj2 [1] ;
    double *Lx ;
    Int *Ls, *Lpi, *Lpx, *Super, *Lp, *Li, *Lnz ;
    Int n, lnz, s, nsuper, p, psi, psx, psend, nsrow, nscol, ii, jj, j, k1, k2,
	q ;

    L->is_ll = to_ll ;

    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lnz = L->nz ;
    lnz = L->nzmax ;

    n = L->n ;
    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;

    p = 0 ;

    for (s = 0 ; s < nsuper ; s++)
    {
	k1 = Super [s] ;
	k2 = Super [s+1] ;
	psi = Lpi [s] ;
	psend = Lpi [s+1] ;
	psx = Lpx [s] ;
	nsrow = psend - psi ;
	nscol = k2 - k1 ;

	for (jj = 0 ; jj < nscol ; jj++)
	{
	    /* column j of L starts here */
	    j = jj + k1 ;

	    if (to_ll)
	    {
		if (to_packed)
		{

		    /* ------------------------------------------------------ */
		    /* convert to LL' packed */
		    /* ------------------------------------------------------ */

		    Lp [j] = p ;
		    PRINT2 (("Col j "ID" p "ID"\n", j, p)) ;
		    for (ii = jj ; ii < nsrow ; ii++)
		    {
			/* get L(i,j) from supernode and store in column j */
			ASSERT (p < (Int) (L->xsize) && p <= psx+ii+jj*nsrow) ;
			Li [p] = Ls [psi + ii] ;
			/* Lx [p] = Lx [psx + ii + jj*nsrow] ; */
			q = psx + ii + jj*nsrow ;
			ASSIGN (Lx,-,p, Lx,-,q) ;
			PRINT2 (("  i "ID" ", Li [p])) ;
			XPRINT2 (Lx,-,q) ;
			PRINT2 (("\n")) ;
			p++ ;
		    }
		    Lnz [j] = p - Lp [j] ;

		}
		else
		{

		    /* ------------------------------------------------------ */
		    /* convert to LL' unpacked */
		    /* ------------------------------------------------------ */

		    p = psx + jj + jj*nsrow ;
		    Lp [j] = p ;
		    Li [p] = j ;
		    Lnz [j] = nsrow - jj ;
		    p++ ;
		    for (ii = jj + 1 ; ii < nsrow ; ii++)
		    {
			/* get L(i,j) from supernode and store in column j */
			Li [psx + ii + jj*nsrow] = Ls [psi + ii] ;
		    }

		}
	    }
	    else
	    {
		if (to_packed)
		{

		    /* ------------------------------------------------------ */
		    /* convert to LDL' packed */
		    /* ------------------------------------------------------ */

		    Lp [j] = p ;
		    PRINT2 (("Col j "ID" p "ID"\n", Lp [j], p)) ;
		    /* ljj = Lx [psx + jj + jj*nsrow] ; */
		    ASSIGN_REAL (ljj, 0, Lx, psx + jj + jj*nsrow) ;

		    if (ljj [0] <= 0)
		    {
			/* the matrix is not positive definite; do not divide */
			/* Lx [p] = ljj ; */
			ASSIGN_REAL (Lx, p, ljj, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
			ljj [0] = 1 ;
		    }
		    else
		    {
			lj2 [0] = ljj [0] * ljj [0] ;
			/* Lx [p] = ljj*ljj ; */
			ASSIGN_REAL (Lx, p, lj2, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
		    }
		    Li [p] = j ;
		    p++ ;
		    for (ii = jj + 1 ; ii < nsrow ; ii++)
		    {
			/* get L(i,j) from supernode and store in column j */
			ASSERT (p < (Int) (L->xsize) && p <= psx+ii+jj*nsrow) ;
			Li [p] = Ls [psi + ii] ;

			/* Lx [p] = Lx [psx + ii + jj*nsrow] / ljj ; */
			q = psx + ii + jj*nsrow ;
			DIV_REAL (Lx, Lz, p, Lx, Lz, q, ljj,0) ;

			PRINT2 (("  i "ID" %g\n", Li [p], Lx [p])) ;
			p++ ;
		    }
		    Lnz [j] = p - Lp [j] ;

		}
		else
		{

		    /* ------------------------------------------------------ */
		    /* convert to LDL' unpacked */
		    /* ------------------------------------------------------ */

		    p = psx + jj + jj*nsrow ;
		    Lp [j] = p ;

		    /* ljj = Lx [p] ; */
		    ASSIGN_REAL (ljj,0, Lx,p) ;

		    if (ljj [0] <= 0)
		    {
			/* the matrix is not positive definite; do not divide */
			/* Lx [p] = ljj ; */
			ASSIGN_REAL (Lx, p, ljj, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
			ljj [0] = 1 ;
		    }
		    else
		    {
			lj2 [0] = ljj [0] * ljj [0] ;
			/* Lx [p] = ljj*ljj ; */
			ASSIGN_REAL (Lx, p, lj2, 0) ;
			CLEAR_IMAG (Lx, Lz, p) ;
		    }
		    Li [p] = j ;
		    Lnz [j] = nsrow - jj ;
		    p++ ;
		    for (ii = jj + 1 ; ii < nsrow ; ii++)
		    {
			/* get L(i,j) from supernode and store in column j */
			Li [psx + ii + jj*nsrow] = Ls [psi + ii] ;

			/* Lx [psx + ii + jj*nsrow] /= ljj ; */
			q = psx + ii + jj*nsrow ;
			DIV_REAL (Lx, Lz, q, Lx, Lz, q, ljj,0) ;
		    }
		}
	    }
	}
    }

    if (to_packed)
    {
	Lp [n] = p ;
	PRINT1 (("Final Lp "ID" n "ID" lnz "ID"\n", p, n, lnz)) ;
	ASSERT (Lp [n] == lnz) ;
	ASSERT (lnz <= (Int) (L->xsize)) ;
	/* reduce size of L->x to match L->i.  This cannot fail. */
	L->x = CHOLMOD(realloc) (lnz, 
#ifdef COMPLEX
		2 *
#endif
		sizeof (double), L->x, &(L->xsize), Common) ;
	ASSERT (lnz == (Int) (L->xsize)) ;
	Common->status = CHOLMOD_OK ;
    }
    else
    {
	Lp [n] = Lpx [nsuper] ;
	ASSERT (MAX (1,Lp [n]) == (Int) (L->xsize)) ;
	ASSERT (MAX (1,Lp [n]) == (Int) (L->nzmax)) ;
    }
}

#endif

#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
