/* ========================================================================== */
/* === Cholesky/t_cholmod_lsolve ============================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2013, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Template routine to solve Lx=b with unit or non-unit diagonal, or
 * solve LDx=b.
 *
 * The numeric xtype of L and Y must match.  Y contains b on input and x on
 * output, stored in row-form.  Y is nrow-by-n, where nrow must equal 1 for the 
 * complex or zomplex cases, and nrow <= 4 for the real case.
 *
 * This file is not compiled separately.  It is included in t_cholmod_solve.c
 * instead.  It contains no user-callable routines.
 *
 * workspace: none
 *
 * Supports real, complex, and zomplex factors.
 */

/* undefine all prior definitions */
#undef FORM_NAME
#undef LSOLVE

/* -------------------------------------------------------------------------- */
/* define the method */
/* -------------------------------------------------------------------------- */

#ifdef LL
/* LL': solve Lx=b with non-unit diagonal */
#define FORM_NAME(prefix,rank) prefix ## ll_lsolve_ ## rank

#elif defined (LD)
/* LDL': solve LDx=b */
#define FORM_NAME(prefix,rank) prefix ## ldl_ldsolve_ ## rank

#else
/* LDL': solve Lx=b with unit diagonal */
#define FORM_NAME(prefix,rank) prefix ## ldl_lsolve_ ## rank

#endif

/* LSOLVE(k) defines the name of a routine for an n-by-k right-hand-side. */

#define LSOLVE(prefix,rank) FORM_NAME(prefix,rank)

#ifdef REAL

/* ========================================================================== */
/* === LSOLVE (1) =========================================================== */
/* ========================================================================== */

/* Solve Lx=b, where b has 1 column  */

static void LSOLVE (PREFIX,1)
(
    cholmod_factor *L,
    double X [ ]                        /* n-by-1 in row form */
)
{
    double *Lx = L->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int j, n = L->n ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y = X [j] ;
#ifdef LL
	    y /= Lx [p] ;
	    X [j] = y ;
#elif defined (LD)
	    X [j] = y / Lx [p] ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		X [Li [p]] -= Lx [p] * y ;
	    }
	    j++ ;	/* advance to next column of L */

	}
	else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2] ;
	    Int q = Lp [j+1] ;
#ifdef LL
	    y [0] = X [j] / Lx [p] ;
	    y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
	    X [j  ] = y [0] ;
	    X [j+1] = y [1] ;
#elif defined (LD)
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    X [j  ] = y [0] / Lx [p] ;
	    X [j+1] = y [1] / Lx [q] ;
#else
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    X [j+1] = y [1] ;
#endif
	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] ;
	    }
	    j += 2 ;	    /* advance to next column of L */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3] ;
	    Int q = Lp [j+1] ;
	    Int r = Lp [j+2] ;
#ifdef LL
	    y [0] = X [j] / Lx [p] ;
	    y [1] = (X [j+1] - Lx [p+1] * y [0]) / Lx [q] ;
	    y [2] = (X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1]) / Lx [r] ;
	    X [j  ] = y [0] ;
	    X [j+1] = y [1] ;
	    X [j+2] = y [2] ;
#elif defined (LD)
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
	    X [j  ] = y [0] / Lx [p] ;
	    X [j+1] = y [1] / Lx [q] ;
	    X [j+2] = y [2] / Lx [r] ;
#else
	    y [0] = X [j] ;
	    y [1] = X [j+1] - Lx [p+1] * y [0] ;
	    y [2] = X [j+2] - Lx [p+2] * y [0] - Lx [q+1] * y [1] ;
	    X [j+1] = y [1] ;
	    X [j+2] = y [2] ;
#endif
	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		X [Li [p]] -= Lx [p] * y [0] + Lx [q] * y [1] + Lx [r] * y [2] ;
	    }
	    j += 3 ;	    /* advance to next column of L */
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (2) =========================================================== */
/* ========================================================================== */

/* Solve Lx=b, where b has 2 columns */

static void LSOLVE (PREFIX,2)
(
    cholmod_factor *L,
    double X [ ][2]		/* n-by-2 in row form */
)
{
    double *Lx = L->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int j, n = L->n ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [2] ;
	    y [0] = X [j][0] ;
	    y [1] = X [j][1] ;
#ifdef LL
	    y [0] /= Lx [p] ;
	    y [1] /= Lx [p] ;
	    X [j][0] = y [0] ;
	    X [j][1] = y [1] ;
#elif defined (LD)
	    X [j][0] = y [0] / Lx [p] ;
	    X [j][1] = y [1] / Lx [p] ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		Int i = Li [p] ;
		X [i][0] -= Lx [p] * y [0] ;
		X [i][1] -= Lx [p] * y [1] ;
	    }
	    j++ ;	/* advance to next column of L */

	}
	else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][2] ;
	    Int q = Lp [j+1] ;
	    y [0][0] = X [j][0] ;
	    y [0][1] = X [j][1] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [1][0] = (X [j+1][0] - Lx [p+1] * y [0][0]) / Lx [q] ;
	    y [1][1] = (X [j+1][1] - Lx [p+1] * y [0][1]) / Lx [q] ;
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
#elif defined (LD)
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    X [j  ][0] = y [0][0] / Lx [p] ;
	    X [j  ][1] = y [0][1] / Lx [p] ;
	    X [j+1][0] = y [1][0] / Lx [q] ;
	    X [j+1][1] = y [1][1] / Lx [q] ;
#else
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
#endif
	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		Int i = Li [p] ;
		X [i][0] -= Lx [p] * y [0][0] + Lx [q] * y [1][0] ;
		X [i][1] -= Lx [p] * y [0][1] + Lx [q] * y [1][1] ;
	    }
	    j += 2 ;	    /* advance to next column of L */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][2] ;
	    Int q = Lp [j+1] ;
	    Int r = Lp [j+2] ;
	    y [0][0] = X [j][0] ;
	    y [0][1] = X [j][1] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [1][0] = (X [j+1][0] - Lx[p+1] * y[0][0]) / Lx [q] ;
	    y [1][1] = (X [j+1][1] - Lx[p+1] * y[0][1]) / Lx [q] ;
	    y [2][0] = (X [j+2][0] - Lx[p+2] * y[0][0] - Lx[q+1]*y[1][0])/Lx[r];
	    y [2][1] = (X [j+2][1] - Lx[p+2] * y[0][1] - Lx[q+1]*y[1][1])/Lx[r];
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+2][0] = y [2][0] ;
	    X [j+2][1] = y [2][1] ;
#elif defined (LD)
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [2][0] = X [j+2][0] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [j+2][1] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    X [j  ][0] = y [0][0] / Lx [p] ;
	    X [j  ][1] = y [0][1] / Lx [p] ;
	    X [j+1][0] = y [1][0] / Lx [q] ;
	    X [j+1][1] = y [1][1] / Lx [q] ;
	    X [j+2][0] = y [2][0] / Lx [r] ;
	    X [j+2][1] = y [2][1] / Lx [r] ;
#else
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [2][0] = X [j+2][0] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [j+2][1] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+2][0] = y [2][0] ;
	    X [j+2][1] = y [2][1] ;
#endif
	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		Int i = Li [p] ;
		X[i][0] -= Lx[p] * y[0][0] + Lx[q] * y[1][0] + Lx[r] * y[2][0] ;
		X[i][1] -= Lx[p] * y[0][1] + Lx[q] * y[1][1] + Lx[r] * y[2][1] ;
	    }
	    j += 3 ;	    /* advance to next column of L */
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (3) =========================================================== */
/* ========================================================================== */

/* Solve Lx=b, where b has 3 columns */

static void LSOLVE (PREFIX,3)
(
    cholmod_factor *L,
    double X [ ][3]			/* n-by-3 in row form */
)
{
    double *Lx = L->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int j, n = L->n ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [3] ;
	    y [0] = X [j][0] ;
	    y [1] = X [j][1] ;
	    y [2] = X [j][2] ;
#ifdef LL
	    y [0] /= Lx [p] ;
	    y [1] /= Lx [p] ;
	    y [2] /= Lx [p] ;
	    X [j][0] = y [0] ;
	    X [j][1] = y [1] ;
	    X [j][2] = y [2] ;
#elif defined (LD)
	    X [j][0] = y [0] / Lx [p] ;
	    X [j][1] = y [1] / Lx [p] ;
	    X [j][2] = y [2] / Lx [p] ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		Int i = Li [p] ;
		double lx = Lx [p] ;
		X [i][0] -= lx * y [0] ;
		X [i][1] -= lx * y [1] ;
		X [i][2] -= lx * y [2] ;
	    }
	    j++ ;	/* advance to next column of L */

	}
	else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][3] ;
	    Int q = Lp [j+1] ;
	    y [0][0] = X [j][0] ;
	    y [0][1] = X [j][1] ;
	    y [0][2] = X [j][2] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [1][0] = (X [j+1][0] - Lx [p+1] * y [0][0]) / Lx [q] ;
	    y [1][1] = (X [j+1][1] - Lx [p+1] * y [0][1]) / Lx [q] ;
	    y [1][2] = (X [j+1][2] - Lx [p+1] * y [0][2]) / Lx [q] ;
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j  ][2] = y [0][2] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
#elif defined (LD)
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    X [j  ][0] = y [0][0] / Lx [p] ;
	    X [j  ][1] = y [0][1] / Lx [p] ;
	    X [j  ][2] = y [0][2] / Lx [p] ;
	    X [j+1][0] = y [1][0] / Lx [q] ;
	    X [j+1][1] = y [1][1] / Lx [q] ;
	    X [j+1][2] = y [1][2] / Lx [q] ;
#else
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
#endif
	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		Int i = Li [p] ;
		double lx [2] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		X [i][0] -= lx [0] * y [0][0] + lx [1] * y [1][0] ;
		X [i][1] -= lx [0] * y [0][1] + lx [1] * y [1][1] ;
		X [i][2] -= lx [0] * y [0][2] + lx [1] * y [1][2] ;
	    }
	    j += 2 ;	    /* advance to next column of L */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][3] ;
	    Int q = Lp [j+1] ;
	    Int r = Lp [j+2] ;
	    y [0][0] = X [j][0] ;
	    y [0][1] = X [j][1] ;
	    y [0][2] = X [j][2] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [1][0] = (X [j+1][0] - Lx[p+1] * y[0][0]) / Lx [q] ;
	    y [1][1] = (X [j+1][1] - Lx[p+1] * y[0][1]) / Lx [q] ;
	    y [1][2] = (X [j+1][2] - Lx[p+1] * y[0][2]) / Lx [q] ;
	    y [2][0] = (X [j+2][0] - Lx[p+2] * y[0][0] - Lx[q+1]*y[1][0])/Lx[r];
	    y [2][1] = (X [j+2][1] - Lx[p+2] * y[0][1] - Lx[q+1]*y[1][1])/Lx[r];
	    y [2][2] = (X [j+2][2] - Lx[p+2] * y[0][2] - Lx[q+1]*y[1][2])/Lx[r];
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j  ][2] = y [0][2] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
	    X [j+2][0] = y [2][0] ;
	    X [j+2][1] = y [2][1] ;
	    X [j+2][2] = y [2][2] ;
#elif defined (LD)
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    y [2][0] = X [j+2][0] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [j+2][1] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [j+2][2] - Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    X [j  ][0] = y [0][0] / Lx [p] ;
	    X [j  ][1] = y [0][1] / Lx [p] ;
	    X [j  ][2] = y [0][2] / Lx [p] ;
	    X [j+1][0] = y [1][0] / Lx [q] ;
	    X [j+1][1] = y [1][1] / Lx [q] ;
	    X [j+1][2] = y [1][2] / Lx [q] ;
	    X [j+2][0] = y [2][0] / Lx [r] ;
	    X [j+2][1] = y [2][1] / Lx [r] ;
	    X [j+2][2] = y [2][2] / Lx [r] ;
#else
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    y [2][0] = X [j+2][0] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [j+2][1] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [j+2][2] - Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
	    X [j+2][0] = y [2][0] ;
	    X [j+2][1] = y [2][1] ;
	    X [j+2][2] = y [2][2] ;
#endif
	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		Int i = Li [p] ;
		double lx [3] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		lx [2] = Lx [r] ;
		X [i][0] -= lx[0] * y[0][0] + lx[1] * y[1][0] + lx[2] * y[2][0];
		X [i][1] -= lx[0] * y[0][1] + lx[1] * y[1][1] + lx[2] * y[2][1];
		X [i][2] -= lx[0] * y[0][2] + lx[1] * y[1][2] + lx[2] * y[2][2];
	    }
	    j += 3 ;	    /* advance to next column of L */
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (4) =========================================================== */
/* ========================================================================== */

/* Solve Lx=b, where b has 4 columns */

static void LSOLVE (PREFIX,4)
(
    cholmod_factor *L,
    double X [ ][4]			    /* n-by-4 in row form */
)
{
    double *Lx = L->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int j, n = L->n ;

    for (j = 0 ; j < n ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j+1, and j+2) */
	if (lnz < 4 || lnz != Lnz [j+1] + 1 || Li [p+1] != j+1)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [4] ;
	    y [0] = X [j][0] ;
	    y [1] = X [j][1] ;
	    y [2] = X [j][2] ;
	    y [3] = X [j][3] ;
#ifdef LL
	    y [0] /= Lx [p] ;
	    y [1] /= Lx [p] ;
	    y [2] /= Lx [p] ;
	    y [3] /= Lx [p] ;
	    X [j][0] = y [0] ;
	    X [j][1] = y [1] ;
	    X [j][2] = y [2] ;
	    X [j][3] = y [3] ;
#elif defined (LD)
	    X [j][0] = y [0] / Lx [p] ;
	    X [j][1] = y [1] / Lx [p] ;
	    X [j][2] = y [2] / Lx [p] ;
	    X [j][3] = y [3] / Lx [p] ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		Int i = Li [p] ;
		double lx = Lx [p] ;
		X [i][0] -= lx * y [0] ;
		X [i][1] -= lx * y [1] ;
		X [i][2] -= lx * y [2] ;
		X [i][3] -= lx * y [3] ;
	    }
	    j++ ;	/* advance to next column of L */

	}
	else if (lnz != Lnz [j+2] + 2 || Li [p+2] != j+2)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][4] ;
	    Int q = Lp [j+1] ;
	    y [0][0] = X [j][0] ;
	    y [0][1] = X [j][1] ;
	    y [0][2] = X [j][2] ;
	    y [0][3] = X [j][3] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [0][3] /= Lx [p] ;
	    y [1][0] = (X [j+1][0] - Lx [p+1] * y [0][0]) / Lx [q] ;
	    y [1][1] = (X [j+1][1] - Lx [p+1] * y [0][1]) / Lx [q] ;
	    y [1][2] = (X [j+1][2] - Lx [p+1] * y [0][2]) / Lx [q] ;
	    y [1][3] = (X [j+1][3] - Lx [p+1] * y [0][3]) / Lx [q] ;
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j  ][2] = y [0][2] ;
	    X [j  ][3] = y [0][3] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
	    X [j+1][3] = y [1][3] ;
#elif defined (LD)
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [j+1][3] - Lx [p+1] * y [0][3] ;
	    X [j  ][0] = y [0][0] / Lx [p] ;
	    X [j  ][1] = y [0][1] / Lx [p] ;
	    X [j  ][2] = y [0][2] / Lx [p] ;
	    X [j  ][3] = y [0][3] / Lx [p] ;
	    X [j+1][0] = y [1][0] / Lx [q] ;
	    X [j+1][1] = y [1][1] / Lx [q] ;
	    X [j+1][2] = y [1][2] / Lx [q] ;
	    X [j+1][3] = y [1][3] / Lx [q] ;
#else
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [j+1][3] - Lx [p+1] * y [0][3] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
	    X [j+1][3] = y [1][3] ;
#endif
	    for (p += 2, q++ ; p < pend ; p++, q++)
	    {
		Int i = Li [p] ;
		double lx [2] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		X [i][0] -= lx [0] * y [0][0] + lx [1] * y [1][0] ;
		X [i][1] -= lx [0] * y [0][1] + lx [1] * y [1][1] ;
		X [i][2] -= lx [0] * y [0][2] + lx [1] * y [1][2] ;
		X [i][3] -= lx [0] * y [0][3] + lx [1] * y [1][3] ;
	    }
	    j += 2 ;	    /* advance to next column of L */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][4] ;
	    Int q = Lp [j+1] ;
	    Int r = Lp [j+2] ;
	    y [0][0] = X [j][0] ;
	    y [0][1] = X [j][1] ;
	    y [0][2] = X [j][2] ;
	    y [0][3] = X [j][3] ;
#ifdef LL
	    y [0][0] /= Lx [p] ;
	    y [0][1] /= Lx [p] ;
	    y [0][2] /= Lx [p] ;
	    y [0][3] /= Lx [p] ;
	    y [1][0] = (X [j+1][0] - Lx[p+1] * y[0][0]) / Lx [q] ;
	    y [1][1] = (X [j+1][1] - Lx[p+1] * y[0][1]) / Lx [q] ;
	    y [1][2] = (X [j+1][2] - Lx[p+1] * y[0][2]) / Lx [q] ;
	    y [1][3] = (X [j+1][3] - Lx[p+1] * y[0][3]) / Lx [q] ;
	    y [2][0] = (X [j+2][0] - Lx[p+2] * y[0][0] - Lx[q+1]*y[1][0])/Lx[r];
	    y [2][1] = (X [j+2][1] - Lx[p+2] * y[0][1] - Lx[q+1]*y[1][1])/Lx[r];
	    y [2][2] = (X [j+2][2] - Lx[p+2] * y[0][2] - Lx[q+1]*y[1][2])/Lx[r];
	    y [2][3] = (X [j+2][3] - Lx[p+2] * y[0][3] - Lx[q+1]*y[1][3])/Lx[r];
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j  ][2] = y [0][2] ;
	    X [j  ][3] = y [0][3] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
	    X [j+1][3] = y [1][3] ;
	    X [j+2][0] = y [2][0] ;
	    X [j+2][1] = y [2][1] ;
	    X [j+2][2] = y [2][2] ;
	    X [j+2][3] = y [2][3] ;
#elif defined (LD)
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [j+1][3] - Lx [p+1] * y [0][3] ;
	    y [2][0] = X [j+2][0] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [j+2][1] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [j+2][2] - Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    y [2][3] = X [j+2][3] - Lx [p+2] * y [0][3] - Lx [q+1] * y [1][3] ;
	    X [j  ][0] = y [0][0] / Lx [p] ;
	    X [j  ][1] = y [0][1] / Lx [p] ;
	    X [j  ][2] = y [0][2] / Lx [p] ;
	    X [j  ][3] = y [0][3] / Lx [p] ;
	    X [j+1][0] = y [1][0] / Lx [q] ;
	    X [j+1][1] = y [1][1] / Lx [q] ;
	    X [j+1][2] = y [1][2] / Lx [q] ;
	    X [j+1][3] = y [1][3] / Lx [q] ;
	    X [j+2][0] = y [2][0] / Lx [r] ;
	    X [j+2][1] = y [2][1] / Lx [r] ;
	    X [j+2][2] = y [2][2] / Lx [r] ;
	    X [j+2][3] = y [2][3] / Lx [r] ;
#else
	    y [1][0] = X [j+1][0] - Lx [p+1] * y [0][0] ;
	    y [1][1] = X [j+1][1] - Lx [p+1] * y [0][1] ;
	    y [1][2] = X [j+1][2] - Lx [p+1] * y [0][2] ;
	    y [1][3] = X [j+1][3] - Lx [p+1] * y [0][3] ;
	    y [2][0] = X [j+2][0] - Lx [p+2] * y [0][0] - Lx [q+1] * y [1][0] ;
	    y [2][1] = X [j+2][1] - Lx [p+2] * y [0][1] - Lx [q+1] * y [1][1] ;
	    y [2][2] = X [j+2][2] - Lx [p+2] * y [0][2] - Lx [q+1] * y [1][2] ;
	    y [2][3] = X [j+2][3] - Lx [p+2] * y [0][3] - Lx [q+1] * y [1][3] ;
	    X [j+1][0] = y [1][0] ;
	    X [j+1][1] = y [1][1] ;
	    X [j+1][2] = y [1][2] ;
	    X [j+1][3] = y [1][3] ;
	    X [j+2][0] = y [2][0] ;
	    X [j+2][1] = y [2][1] ;
	    X [j+2][2] = y [2][2] ;
	    X [j+2][3] = y [2][3] ;
#endif
	    for (p += 3, q += 2, r++ ; p < pend ; p++, q++, r++)
	    {
		Int i = Li [p] ;
		double lx [3] ;
		lx [0] = Lx [p] ;
		lx [1] = Lx [q] ;
		lx [2] = Lx [r] ;
		X [i][0] -= lx[0] * y[0][0] + lx[1] * y[1][0] + lx[2] * y[2][0];
		X [i][1] -= lx[0] * y[0][1] + lx[1] * y[1][1] + lx[2] * y[2][1];
		X [i][2] -= lx[0] * y[0][2] + lx[1] * y[1][2] + lx[2] * y[2][2];
		X [i][3] -= lx[0] * y[0][3] + lx[1] * y[1][3] + lx[2] * y[2][3];
	    }
	    j += 3 ;	    /* advance to next column of L */
	}
    }
}

#endif


/* ========================================================================== */
/* === LSOLVE (k) =========================================================== */
/* ========================================================================== */

static void LSOLVE (PREFIX,k)
(
    cholmod_factor *L,
    cholmod_dense *Y,		    /* nr-by-n where nr is 1 to 4 */
    Int *Yseti, Int ysetlen
)
{

    double yx [2] ;
#ifdef ZOMPLEX
    double yz [1] ;
    double *Lz = L->z ;
    double *Xz = Y->z ;
#endif
    double *Lx = L->x ;
    double *Xx = Y->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int n = L->n, jj, jjiters ;

    ASSERT (L->xtype == Y->xtype) ; /* L and Y must have the same xtype */
    ASSERT (L->n == Y->ncol) ;	    /* dimensions must match */
    ASSERT (Y->nrow == Y->d) ;	    /* leading dimension of Y = # rows of Y */
    ASSERT (L->xtype != CHOLMOD_PATTERN) ;  /* L is not symbolic */
    ASSERT (!(L->is_super)) ;	    /* L is simplicial LL' or LDL' */

#ifdef REAL

    if (Yseti == NULL)
    {

        /* ------------------------------------------------------------------ */
        /* real case, no Yseti, with 1 to 4 RHS's and dynamic supernodes */
        /* ------------------------------------------------------------------ */

        ASSERT (Y->nrow <= 4) ;

        switch (Y->nrow)
        {
            case 1: LSOLVE (PREFIX,1) (L, Y->x) ; break ;
            case 2: LSOLVE (PREFIX,2) (L, Y->x) ; break ;
            case 3: LSOLVE (PREFIX,3) (L, Y->x) ; break ;
            case 4: LSOLVE (PREFIX,4) (L, Y->x) ; break ;
        }

    }
    else
#endif
    {

        /* ------------------------------------------------------------------ */
        /* solve a complex linear system or solve with Yseti */
        /* ------------------------------------------------------------------ */

        ASSERT (Y->nrow == 1) ;

        jjiters = Yseti ? ysetlen : n ;

        for (jj = 0 ; jj < jjiters ; jj++)
        {
            Int j = Yseti ? Yseti [jj] : jj ;

            /* get the start, end, and length of column j */
            Int p = Lp [j] ;
            Int lnz = Lnz [j] ;
            Int pend = p + lnz ;

            /* y = X [j] ; */
            ASSIGN (yx,yz,0, Xx,Xz,j) ;

#ifdef LL
            /* y /= Lx [p] ; */
            /* X [j] = y ; */
            DIV_REAL (yx,yz,0, yx,yz,0, Lx,p) ;
            ASSIGN (Xx,Xz,j, yx,yz,0) ;
#elif defined (LD)
            /* X [j] = y / Lx [p] ; */
            DIV_REAL (Xx,Xz,j, yx,yz,0, Lx,p) ;
#endif

            for (p++ ; p < pend ; p++)
            {
                /* X [Li [p]] -= Lx [p] * y ; */
                Int i = Li [p] ;
                MULTSUB (Xx,Xz,i, Lx,Lz,p, yx,yz,0) ;
            }
        }
    }
}

/* prepare for the next inclusion of this file in cholmod_solve.c */
#undef LL
#undef LD
