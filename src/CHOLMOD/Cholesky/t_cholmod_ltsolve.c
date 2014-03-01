/* ========================================================================== */
/* === Cholesky/t_cholmod_ltsolve =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2013, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Template routine to solve L'x=b with unit or non-unit diagonal, or
 * solve DL'x=b.
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
#undef DIAG

/* -------------------------------------------------------------------------- */
/* define the method */
/* -------------------------------------------------------------------------- */

#ifdef LL
/* LL': solve Lx=b with non-unit diagonal */
#define FORM_NAME(prefix,rank) prefix ## ll_ltsolve_ ## rank
#define DIAG

#elif defined (LD)
/* LDL': solve LDx=b */
#define FORM_NAME(prefix,rank) prefix ## ldl_dltsolve_ ## rank
#define DIAG

#else
/* LDL': solve Lx=b with unit diagonal */
#define FORM_NAME(prefix,rank) prefix ## ldl_ltsolve_ ## rank

#endif

/* LSOLVE(k) defines the name of a routine for an n-by-k right-hand-side. */
#define LSOLVE(prefix,rank) FORM_NAME(prefix,rank)

#ifdef REAL

/* ========================================================================== */
/* === LSOLVE (1) =========================================================== */
/* ========================================================================== */

/* Solve L'x=b, where b has 1 column  */

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

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j-1, and j-2) */
	if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y = X [j] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y /= d ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		y -= Lx [p] * X [Li [p]] ;
	    }
#ifdef LL
	    X [j] = y / d ;
#else
	    X [j] = y ;
#endif
	    j-- ;	/* advance to the next column of L */

	}
	else if (lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2], t ;
	    Int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0] = X [j  ] / d [0] ;
	    y [1] = X [j-1] / d [1] ;
#else
	    y [0] = X [j  ] ;
	    y [1] = X [j-1] ;
#endif
	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		Int i = Li [p] ;
		y [0] -= Lx [p] * X [i] ;
		y [1] -= Lx [q] * X [i] ;
	    }
#ifdef LL
	    y [0] /= d [0] ;
	    y [1] = (y [1] - t * y [0]) / d [1] ;
#else
	    y [1] -= t * y [0] ;
#endif
	    X [j  ] = y [0] ;
	    X [j-1] = y [1] ;
	    j -= 2 ;	    /* advance to the next column of L */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3], t [3] ;
	    Int q = Lp [j-1] ;
	    Int r = Lp [j-2] ;
#ifdef DIAG
	    double d [3] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
	    d [2] = Lx [r] ;
#endif
	    t [0] = Lx [q+1] ;
	    t [1] = Lx [r+1] ;
	    t [2] = Lx [r+2] ;
#ifdef LD
	    y [0] = X [j]   / d [0] ;
	    y [1] = X [j-1] / d [1] ;
	    y [2] = X [j-2] / d [2] ;
#else
	    y [0] = X [j] ;
	    y [1] = X [j-1] ;
	    y [2] = X [j-2] ;
#endif
	    for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
	    {
		Int i = Li [p] ;
		y [0] -= Lx [p] * X [i] ;
		y [1] -= Lx [q] * X [i] ;
		y [2] -= Lx [r] * X [i] ;
	    }
#ifdef LL
	    y [0] /= d [0] ;
	    y [1] = (y [1] - t [0] * y [0]) / d [1] ;
	    y [2] = (y [2] - t [2] * y [0] - t [1] * y [1]) / d [2] ;
#else
	    y [1] -= t [0] * y [0] ;
	    y [2] -= t [2] * y [0] + t [1] * y [1] ;
#endif
	    X [j-2] = y [2] ;
	    X [j-1] = y [1] ;
	    X [j  ] = y [0] ;
	    j -= 3 ;	    /* advance to the next column of L */
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (2) =========================================================== */
/* ========================================================================== */

/* Solve L'x=b, where b has 2 columns */

static void LSOLVE (PREFIX,2)
(
    cholmod_factor *L,
    double X [ ][2]		    /* n-by-2 in row form */
)
{
    double *Lx = L->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int j, n = L->n ;

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j-1, and j-2) */
	if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [2] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y [0] = X [j][0] / d ;
	    y [1] = X [j][1] / d ;
#else
	    y [0] = X [j][0] ;
	    y [1] = X [j][1] ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		Int i = Li [p] ;
		y [0] -= Lx [p] * X [i][0] ;
		y [1] -= Lx [p] * X [i][1] ;
	    }
#ifdef LL
	    X [j][0] = y [0] / d ;
	    X [j][1] = y [1] / d ;
#else
	    X [j][0] = y [0] ;
	    X [j][1] = y [1] ;
#endif
	    j-- ;	/* advance to the next column of L */

	}
	else if (lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][2], t ;
	    Int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0][0] = X [j  ][0] / d [0] ;
	    y [0][1] = X [j  ][1] / d [0] ;
	    y [1][0] = X [j-1][0] / d [1] ;
	    y [1][1] = X [j-1][1] / d [1] ;
#else
	    y [0][0] = X [j  ][0] ;
	    y [0][1] = X [j  ][1] ;
	    y [1][0] = X [j-1][0] ;
	    y [1][1] = X [j-1][1] ;
#endif
	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		Int i = Li [p] ;
		y [0][0] -= Lx [p] * X [i][0] ;
		y [0][1] -= Lx [p] * X [i][1] ;
		y [1][0] -= Lx [q] * X [i][0] ;
		y [1][1] -= Lx [q] * X [i][1] ;
	    }
#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [1][0] = (y [1][0] - t * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t * y [0][1]) / d [1] ;
#else
	    y [1][0] -= t * y [0][0] ;
	    y [1][1] -= t * y [0][1] ;
#endif
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j-1][0] = y [1][0] ;
	    X [j-1][1] = y [1][1] ;
	    j -= 2 ;	    /* advance to the next column of L */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][2], t [3] ;
	    Int q = Lp [j-1] ;
	    Int r = Lp [j-2] ;
#ifdef DIAG
	    double d [3] ; 
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
	    d [2] = Lx [r] ;
#endif
	    t [0] = Lx [q+1] ;
	    t [1] = Lx [r+1] ;
	    t [2] = Lx [r+2] ;
#ifdef LD
	    y [0][0] = X [j  ][0] / d [0] ;
	    y [0][1] = X [j  ][1] / d [0] ;
	    y [1][0] = X [j-1][0] / d [1] ;
	    y [1][1] = X [j-1][1] / d [1] ;
	    y [2][0] = X [j-2][0] / d [2] ;
	    y [2][1] = X [j-2][1] / d [2] ;
#else
	    y [0][0] = X [j  ][0] ;
	    y [0][1] = X [j  ][1] ;
	    y [1][0] = X [j-1][0] ;
	    y [1][1] = X [j-1][1] ;
	    y [2][0] = X [j-2][0] ;
	    y [2][1] = X [j-2][1] ;
#endif
	    for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
	    {
		Int i = Li [p] ;
		y [0][0] -= Lx [p] * X [i][0] ;
		y [0][1] -= Lx [p] * X [i][1] ;
		y [1][0] -= Lx [q] * X [i][0] ;
		y [1][1] -= Lx [q] * X [i][1] ;
		y [2][0] -= Lx [r] * X [i][0] ;
		y [2][1] -= Lx [r] * X [i][1] ;
	    }
#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [1][0] = (y [1][0] - t [0] * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t [0] * y [0][1]) / d [1] ;
	    y [2][0] = (y [2][0] - t [2] * y [0][0] - t [1] * y [1][0]) / d [2];
	    y [2][1] = (y [2][1] - t [2] * y [0][1] - t [1] * y [1][1]) / d [2];
#else
	    y [1][0] -= t [0] * y [0][0] ;
	    y [1][1] -= t [0] * y [0][1] ;
	    y [2][0] -= t [2] * y [0][0] + t [1] * y [1][0] ;
	    y [2][1] -= t [2] * y [0][1] + t [1] * y [1][1] ;
#endif
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j-1][0] = y [1][0] ;
	    X [j-1][1] = y [1][1] ;
	    X [j-2][0] = y [2][0] ;
	    X [j-2][1] = y [2][1] ;
	    j -= 3 ;	    /* advance to the next column of L */
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (3) =========================================================== */
/* ========================================================================== */

/* Solve L'x=b, where b has 3 columns */

static void LSOLVE (PREFIX,3)
(
    cholmod_factor *L,
    double X [ ][3]		    /* n-by-3 in row form */
)
{
    double *Lx = L->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int j, n = L->n ;

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j-1, and j-2) */
	if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [3] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y [0] = X [j][0] / d ;
	    y [1] = X [j][1] / d ;
	    y [2] = X [j][2] / d ;
#else
	    y [0] = X [j][0] ;
	    y [1] = X [j][1] ;
	    y [2] = X [j][2] ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		Int i = Li [p] ;
		y [0] -= Lx [p] * X [i][0] ;
		y [1] -= Lx [p] * X [i][1] ;
		y [2] -= Lx [p] * X [i][2] ;
	    }
#ifdef LL
	    X [j][0] = y [0] / d ;
	    X [j][1] = y [1] / d ;
	    X [j][2] = y [2] / d ;
#else
	    X [j][0] = y [0] ;
	    X [j][1] = y [1] ;
	    X [j][2] = y [2] ;
#endif
	    j-- ;	/* advance to the next column of L */

	}
	else if (lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][3], t ;
	    Int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0][0] = X [j  ][0] / d [0] ;
	    y [0][1] = X [j  ][1] / d [0] ;
	    y [0][2] = X [j  ][2] / d [0] ;
	    y [1][0] = X [j-1][0] / d [1] ;
	    y [1][1] = X [j-1][1] / d [1] ;
	    y [1][2] = X [j-1][2] / d [1] ;
#else
	    y [0][0] = X [j  ][0] ;
	    y [0][1] = X [j  ][1] ;
	    y [0][2] = X [j  ][2] ;
	    y [1][0] = X [j-1][0] ;
	    y [1][1] = X [j-1][1] ;
	    y [1][2] = X [j-1][2] ;
#endif
	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		Int i = Li [p] ;
		y [0][0] -= Lx [p] * X [i][0] ;
		y [0][1] -= Lx [p] * X [i][1] ;
		y [0][2] -= Lx [p] * X [i][2] ;
		y [1][0] -= Lx [q] * X [i][0] ;
		y [1][1] -= Lx [q] * X [i][1] ;
		y [1][2] -= Lx [q] * X [i][2] ;
	    }
#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [0][2] /= d [0] ;
	    y [1][0] = (y [1][0] - t * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t * y [0][1]) / d [1] ;
	    y [1][2] = (y [1][2] - t * y [0][2]) / d [1] ;
#else
	    y [1][0] -= t * y [0][0] ;
	    y [1][1] -= t * y [0][1] ;
	    y [1][2] -= t * y [0][2] ;
#endif
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j  ][2] = y [0][2] ;
	    X [j-1][0] = y [1][0] ;
	    X [j-1][1] = y [1][1] ;
	    X [j-1][2] = y [1][2] ;
	    j -= 2 ;	    /* advance to the next column of L */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of three columns of L */
	    /* -------------------------------------------------------------- */

	    double y [3][3], t [3] ;
	    Int q = Lp [j-1] ;
	    Int r = Lp [j-2] ;
#ifdef DIAG
	    double d [3] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
	    d [2] = Lx [r] ;
#endif
	    t [0] = Lx [q+1] ;
	    t [1] = Lx [r+1] ;
	    t [2] = Lx [r+2] ;
#ifdef LD
	    y [0][0] = X [j  ][0] / d [0] ;
	    y [0][1] = X [j  ][1] / d [0] ;
	    y [0][2] = X [j  ][2] / d [0] ;
	    y [1][0] = X [j-1][0] / d [1] ;
	    y [1][1] = X [j-1][1] / d [1] ;
	    y [1][2] = X [j-1][2] / d [1] ;
	    y [2][0] = X [j-2][0] / d [2] ;
	    y [2][1] = X [j-2][1] / d [2] ;
	    y [2][2] = X [j-2][2] / d [2] ;
#else
	    y [0][0] = X [j  ][0] ;
	    y [0][1] = X [j  ][1] ;
	    y [0][2] = X [j  ][2] ;
	    y [1][0] = X [j-1][0] ;
	    y [1][1] = X [j-1][1] ;
	    y [1][2] = X [j-1][2] ;
	    y [2][0] = X [j-2][0] ;
	    y [2][1] = X [j-2][1] ;
	    y [2][2] = X [j-2][2] ;
#endif
	    for (p++, q += 2, r += 3 ; p < pend ; p++, q++, r++)
	    {
		Int i = Li [p] ;
		y [0][0] -= Lx [p] * X [i][0] ;
		y [0][1] -= Lx [p] * X [i][1] ;
		y [0][2] -= Lx [p] * X [i][2] ;
		y [1][0] -= Lx [q] * X [i][0] ;
		y [1][1] -= Lx [q] * X [i][1] ;
		y [1][2] -= Lx [q] * X [i][2] ;
		y [2][0] -= Lx [r] * X [i][0] ;
		y [2][1] -= Lx [r] * X [i][1] ;
		y [2][2] -= Lx [r] * X [i][2] ;
	    }
#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [0][2] /= d [0] ;
	    y [1][0] = (y [1][0] - t [0] * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t [0] * y [0][1]) / d [1] ;
	    y [1][2] = (y [1][2] - t [0] * y [0][2]) / d [1] ;
	    y [2][0] = (y [2][0] - t [2] * y [0][0] - t [1] * y [1][0]) / d [2];
	    y [2][1] = (y [2][1] - t [2] * y [0][1] - t [1] * y [1][1]) / d [2];
	    y [2][2] = (y [2][2] - t [2] * y [0][2] - t [1] * y [1][2]) / d [2];
#else
	    y [1][0] -= t [0] * y [0][0] ;
	    y [1][1] -= t [0] * y [0][1] ;
	    y [1][2] -= t [0] * y [0][2] ;
	    y [2][0] -= t [2] * y [0][0] + t [1] * y [1][0] ;
	    y [2][1] -= t [2] * y [0][1] + t [1] * y [1][1] ;
	    y [2][2] -= t [2] * y [0][2] + t [1] * y [1][2] ;
#endif
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j  ][2] = y [0][2] ;
	    X [j-1][0] = y [1][0] ;
	    X [j-1][1] = y [1][1] ;
	    X [j-1][2] = y [1][2] ;
	    X [j-2][0] = y [2][0] ;
	    X [j-2][1] = y [2][1] ;
	    X [j-2][2] = y [2][2] ;
	    j -= 3 ;	    /* advance to the next column of L */
	}
    }
}


/* ========================================================================== */
/* === LSOLVE (4) =========================================================== */
/* ========================================================================== */

/* Solve L'x=b, where b has 4 columns */

static void LSOLVE (PREFIX,4)
(
    cholmod_factor *L,
    double X [ ][4]		    /* n-by-4 in row form */
)
{
    double *Lx = L->x ;
    Int *Li = L->i ;
    Int *Lp = L->p ;
    Int *Lnz = L->nz ;
    Int j, n = L->n ;

    for (j = n-1 ; j >= 0 ; )
    {
	/* get the start, end, and length of column j */
	Int p = Lp [j] ;
	Int lnz = Lnz [j] ;
	Int pend = p + lnz ;

	/* find a chain of supernodes (up to j, j-1, and j-2) */
	if (j < 4 || lnz != Lnz [j-1] - 1 || Li [Lp [j-1]+1] != j)
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a single column of L */
	    /* -------------------------------------------------------------- */

	    double y [4] ;
#ifdef DIAG
	    double d = Lx [p] ;
#endif
#ifdef LD
	    y [0] = X [j][0] / d ;
	    y [1] = X [j][1] / d ;
	    y [2] = X [j][2] / d ;
	    y [3] = X [j][3] / d ;
#else
	    y [0] = X [j][0] ;
	    y [1] = X [j][1] ;
	    y [2] = X [j][2] ;
	    y [3] = X [j][3] ;
#endif
	    for (p++ ; p < pend ; p++)
	    {
		Int i = Li [p] ;
		y [0] -= Lx [p] * X [i][0] ;
		y [1] -= Lx [p] * X [i][1] ;
		y [2] -= Lx [p] * X [i][2] ;
		y [3] -= Lx [p] * X [i][3] ;
	    }
#ifdef LL
	    X [j][0] = y [0] / d ;
	    X [j][1] = y [1] / d ;
	    X [j][2] = y [2] / d ;
	    X [j][3] = y [3] / d ;
#else
	    X [j][0] = y [0] ;
	    X [j][1] = y [1] ;
	    X [j][2] = y [2] ;
	    X [j][3] = y [3] ;
#endif
	    j-- ;	/* advance to the next column of L */

	}
	else /* if (j == 1 || lnz != Lnz [j-2]-2 || Li [Lp [j-2]+2] != j) */
	{

	    /* -------------------------------------------------------------- */
	    /* solve with a supernode of two columns of L */
	    /* -------------------------------------------------------------- */

	    double y [2][4], t ;
	    Int q = Lp [j-1] ;
#ifdef DIAG
	    double d [2] ;
	    d [0] = Lx [p] ;
	    d [1] = Lx [q] ;
#endif
	    t = Lx [q+1] ;
#ifdef LD
	    y [0][0] = X [j  ][0] / d [0] ;
	    y [0][1] = X [j  ][1] / d [0] ;
	    y [0][2] = X [j  ][2] / d [0] ;
	    y [0][3] = X [j  ][3] / d [0] ;
	    y [1][0] = X [j-1][0] / d [1] ;
	    y [1][1] = X [j-1][1] / d [1] ;
	    y [1][2] = X [j-1][2] / d [1] ;
	    y [1][3] = X [j-1][3] / d [1] ;
#else
	    y [0][0] = X [j  ][0] ;
	    y [0][1] = X [j  ][1] ;
	    y [0][2] = X [j  ][2] ;
	    y [0][3] = X [j  ][3] ;
	    y [1][0] = X [j-1][0] ;
	    y [1][1] = X [j-1][1] ;
	    y [1][2] = X [j-1][2] ;
	    y [1][3] = X [j-1][3] ;
#endif
	    for (p++, q += 2 ; p < pend ; p++, q++)
	    {
		Int i = Li [p] ;
		y [0][0] -= Lx [p] * X [i][0] ;
		y [0][1] -= Lx [p] * X [i][1] ;
		y [0][2] -= Lx [p] * X [i][2] ;
		y [0][3] -= Lx [p] * X [i][3] ;
		y [1][0] -= Lx [q] * X [i][0] ;
		y [1][1] -= Lx [q] * X [i][1] ;
		y [1][2] -= Lx [q] * X [i][2] ;
		y [1][3] -= Lx [q] * X [i][3] ;
	    }
#ifdef LL
	    y [0][0] /= d [0] ;
	    y [0][1] /= d [0] ;
	    y [0][2] /= d [0] ;
	    y [0][3] /= d [0] ;
	    y [1][0] = (y [1][0] - t * y [0][0]) / d [1] ;
	    y [1][1] = (y [1][1] - t * y [0][1]) / d [1] ;
	    y [1][2] = (y [1][2] - t * y [0][2]) / d [1] ;
	    y [1][3] = (y [1][3] - t * y [0][3]) / d [1] ;
#else
	    y [1][0] -= t * y [0][0] ;
	    y [1][1] -= t * y [0][1] ;
	    y [1][2] -= t * y [0][2] ;
	    y [1][3] -= t * y [0][3] ;
#endif
	    X [j  ][0] = y [0][0] ;
	    X [j  ][1] = y [0][1] ;
	    X [j  ][2] = y [0][2] ;
	    X [j  ][3] = y [0][3] ;
	    X [j-1][0] = y [1][0] ;
	    X [j-1][1] = y [1][1] ;
	    X [j-1][2] = y [1][2] ;
	    X [j-1][3] = y [1][3] ;
	    j -= 2 ;	    /* advance to the next column of L */
	}

	/* NOTE: with 4 right-hand-sides, it suffices to exploit dynamic
	 * supernodes of just size 1 and 2.  3-column supernodes are not
	 * needed. */
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

#ifdef DIAG
    double d [1] ;
#endif
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

        for (jj = jjiters-1 ; jj >= 0 ; jj--)
        {

            Int j = Yseti ? Yseti [jj] : jj ;

            /* get the start, end, and length of column j */
            Int p = Lp [j] ;
            Int lnz = Lnz [j] ;
            Int pend = p + lnz ;

            /* y = X [j] ; */
            ASSIGN (yx,yz,0, Xx,Xz,j) ;

#ifdef DIAG
            /* d = Lx [p] ; */
            ASSIGN_REAL (d,0, Lx,p) ;
#endif
#ifdef LD
            /* y /= d ; */
            DIV_REAL (yx,yz,0, yx,yz,0, d,0) ;
#endif

            for (p++ ; p < pend ; p++)
            {
                /* y -= conj (Lx [p]) * X [Li [p]] ; */
                Int i = Li [p] ;
                MULTSUBCONJ (yx,yz,0, Lx,Lz,p, Xx,Xz,i) ;
            }

#ifdef LL
            /* X [j] = y / d ; */
            DIV_REAL (Xx,Xz,j, yx,yz,0, d,0) ;
#else
            /* X [j] = y ; */
            ASSIGN (Xx,Xz,j, yx,yz,0) ;
#endif

        }
    }
}

/* prepare for the next inclusion of this file in cholmod_solve.c */
#undef LL
#undef LD
