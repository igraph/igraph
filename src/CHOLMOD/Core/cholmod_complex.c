/* ========================================================================== */
/* === Core/cholmod_complex ================================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* If you convert a matrix that contains uninitialized data, valgrind will
 * complain.  This can occur in a factor L which has gaps (a partial
 * factorization, or after updates that change the nonzero pattern), an
 * unpacked sparse matrix, a dense matrix with leading dimension d > # of rows,
 * or any matrix (dense, sparse, triplet, or factor) with more space allocated
 * than is used.  You can safely ignore any of these complaints by valgrind. */

#include "cholmod_internal.h"
#include "cholmod_core.h"

/* ========================================================================== */
/* === cholmod_hypot ======================================================== */
/* ========================================================================== */

/* There is an equivalent routine called hypot in <math.h>, which conforms
 * to ANSI C99.  However, CHOLMOD does not assume that ANSI C99 is available.
 * You can use the ANSI C99 hypot routine with:
 *
 *	#include <math.h>
 *	Common->hypotenuse = hypot ;
 *
 * Default value of the Common->hypotenuse pointer is cholmod_hypot.
 *
 * s = hypot (x,y) computes s = sqrt (x*x + y*y) but does so more accurately.
 * The NaN cases for the double relops x >= y and x+y == x are safely ignored.
 * 
 * Source: Algorithm 312, "Absolute value and square root of a complex number,"
 * P. Friedland, Comm. ACM, vol 10, no 10, October 1967, page 665.
 */

double CHOLMOD(hypot) (double x, double y)
{
    double s, r ;
    x = fabs (x) ;
    y = fabs (y) ;
    if (x >= y)
    {
	if (x + y == x)
	{
	    s = x ;
	}
	else
	{
	    r = y / x ;
	    s = x * sqrt (1.0 + r*r) ;
	}
    }
    else
    {
	if (y + x == y)
	{
	    s = y ;
	}
	else
	{
	    r = x / y ;
	    s = y * sqrt (1.0 + r*r) ;
	}
    } 
    return (s) ;
}


/* ========================================================================== */
/* === cholmod_divcomplex =================================================== */
/* ========================================================================== */

/* c = a/b where c, a, and b are complex.  The real and imaginary parts are
 * passed as separate arguments to this routine.  The NaN case is ignored
 * for the double relop br >= bi.  Returns 1 if the denominator is zero,
 * 0 otherwise.  Note that this return value is the single exception to the
 * rule that all CHOLMOD routines that return int return TRUE if successful
 * or FALSE otherise.
 *
 * This uses ACM Algo 116, by R. L. Smith, 1962, which tries to avoid
 * underflow and overflow.
 *
 * c can be the same variable as a or b.
 *
 * Default value of the Common->complex_divide pointer is cholmod_divcomplex.
 */

int CHOLMOD(divcomplex)
(
    double ar, double ai,	/* real and imaginary parts of a */
    double br, double bi,	/* real and imaginary parts of b */
    double *cr, double *ci	/* real and imaginary parts of c */
)
{
    double tr, ti, r, den ;
    if (fabs (br) >= fabs (bi))
    {
	r = bi / br ;
	den = br + r * bi ;
	tr = (ar + ai * r) / den ;
	ti = (ai - ar * r) / den ;
    }
    else
    {
	r = br / bi ;
	den = r * br + bi ;
	tr = (ar * r + ai) / den ;
	ti = (ai * r - ar) / den ;
    }
    *cr = tr ;
    *ci = ti ;
    return (IS_ZERO (den)) ;
}


/* ========================================================================== */
/* === change_complexity ==================================================== */
/* ========================================================================== */

/* X and Z represent an array of size nz, with numeric xtype given by xtype_in.
 *
 * If xtype_in is:
 * CHOLMOD_PATTERN: X and Z must be NULL.
 * CHOLMOD_REAL:    X is of size nz, Z must be NULL.
 * CHOLMOD_COMPLEX: X is of size 2*nz, Z must be NULL.
 * CHOLMOD_ZOMPLEX: X is of size nz, Z is of size nz.
 *
 * The array is changed into the numeric xtype given by xtype_out, with the
 * same definitions of X and Z above.  Note that the input conditions, above,
 * are not checked.  These are checked in the caller routine.
 *
 * Returns TRUE if successful, FALSE otherwise.  X and Z are not modified if
 * not successful.
 */

static int change_complexity
(
    /* ---- input ---- */
    Int nz,		/* size of X and/or Z */
    int xtype_in,	/* xtype of X and Z on input */
    int xtype_out,	/* requested xtype of X and Z on output */
    int xtype1,		/* xtype_out must be in the range [xtype1 .. xtype2] */
    int xtype2,
    /* ---- in/out --- */
    void **XX,		/* old X on input, new X on output */
    void **ZZ,		/* old Z on input, new Z on output */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Xold, *Zold, *Xnew, *Znew ;
    Int k ;
    size_t nz2 ;

    if (xtype_out < xtype1 || xtype_out > xtype2)
    {
	ERROR (CHOLMOD_INVALID, "invalid xtype") ;
	return (FALSE) ;
    }

    Common->status = CHOLMOD_OK ;
    Xold = *XX ;
    Zold = *ZZ ;

    switch (xtype_in)
    {

	/* ------------------------------------------------------------------ */
	/* converting from pattern */
	/* ------------------------------------------------------------------ */

	case CHOLMOD_PATTERN:

	    switch (xtype_out)
	    {

		/* ---------------------------------------------------------- */
		/* pattern -> real */
		/* ---------------------------------------------------------- */

		case CHOLMOD_REAL:
		    /* allocate X and set to all ones */
		    Xnew = CHOLMOD(malloc) (nz, sizeof (double), Common) ;
		    if (Common->status < CHOLMOD_OK)
		    {
			return (FALSE) ;
		    }
		    for (k = 0 ; k < nz ; k++)
		    {
			Xnew [k] = 1 ;
		    }
		    *XX = Xnew ;
		    break ;

		/* ---------------------------------------------------------- */
		/* pattern -> complex */
		/* ---------------------------------------------------------- */

		case CHOLMOD_COMPLEX:
		    /* allocate X and set to all ones */
		    Xnew = CHOLMOD(malloc) (nz, 2*sizeof (double), Common) ;
		    if (Common->status < CHOLMOD_OK)
		    {
			return (FALSE) ;
		    }
		    for (k = 0 ; k < nz ; k++)
		    {
			Xnew [2*k  ] = 1 ;
			Xnew [2*k+1] = 0 ;
		    }
		    *XX = Xnew ;
		    break ;

		/* ---------------------------------------------------------- */
		/* pattern -> zomplex */
		/* ---------------------------------------------------------- */

		case CHOLMOD_ZOMPLEX:
		    /* allocate X and Z and set to all ones */
		    Xnew = CHOLMOD(malloc) (nz, sizeof (double), Common) ;
		    Znew = CHOLMOD(malloc) (nz, sizeof (double), Common) ;
		    if (Common->status < CHOLMOD_OK)
		    {
			CHOLMOD(free) (nz, sizeof (double), Xnew, Common) ;
			CHOLMOD(free) (nz, sizeof (double), Znew, Common) ;
			return (FALSE) ;
		    }
		    for (k = 0 ; k < nz ; k++)
		    {
			Xnew [k] = 1 ;
			Znew [k] = 0 ;
		    }
		    *XX = Xnew ;
		    *ZZ = Znew ;
		    break ;
	    }
	    break ;

	/* ------------------------------------------------------------------ */
	/* converting from real */
	/* ------------------------------------------------------------------ */

	case CHOLMOD_REAL:

	    switch (xtype_out)
	    {

		/* ---------------------------------------------------------- */
		/* real -> pattern */
		/* ---------------------------------------------------------- */

		case CHOLMOD_PATTERN:
		    /* free X */
		    *XX = CHOLMOD(free) (nz, sizeof (double), *XX, Common) ;
		    break ;

		/* ---------------------------------------------------------- */
		/* real -> complex */
		/* ---------------------------------------------------------- */

		case CHOLMOD_COMPLEX:
		    /* allocate a new X and copy the old X */
		    Xnew = CHOLMOD(malloc) (nz, 2*sizeof (double), Common) ;
		    if (Common->status < CHOLMOD_OK)
		    {
			return (FALSE) ;
		    }
		    for (k = 0 ; k < nz ; k++)
		    {
			Xnew [2*k  ] = Xold [k] ;
			Xnew [2*k+1] = 0 ;
		    }
		    CHOLMOD(free) (nz, sizeof (double), *XX, Common) ;
		    *XX = Xnew ;
		    break ;

		/* ---------------------------------------------------------- */
		/* real -> zomplex */
		/* ---------------------------------------------------------- */

		case CHOLMOD_ZOMPLEX:
		    /* allocate a new Z and set it to zero */
		    Znew = CHOLMOD(malloc) (nz, sizeof (double), Common) ;
		    if (Common->status < CHOLMOD_OK)
		    {
			return (FALSE) ;
		    }
		    for (k = 0 ; k < nz ; k++)
		    {
			Znew [k] = 0 ;
		    }
		    *ZZ = Znew ;
		    break ;
	    }
	    break ;

	/* ------------------------------------------------------------------ */
	/* converting from complex */
	/* ------------------------------------------------------------------ */

	case CHOLMOD_COMPLEX:

	    switch (xtype_out)
	    {

		/* ---------------------------------------------------------- */
		/* complex -> pattern */
		/* ---------------------------------------------------------- */

		case CHOLMOD_PATTERN:
		    /* free X */
		    *XX = CHOLMOD(free) (nz, 2*sizeof (double), *XX, Common) ;
		    break ;

		/* ---------------------------------------------------------- */
		/* complex -> real */
		/* ---------------------------------------------------------- */

		case CHOLMOD_REAL:
		    /* pack the real part of X, discarding the imaginary part */
		    for (k = 0 ; k < nz ; k++)
		    {
			Xold [k] = Xold [2*k] ;
		    }
		    /* shrink X in half (this cannot fail) */
		    nz2 = 2*nz ;
		    *XX = CHOLMOD(realloc) (nz, sizeof (double), *XX, &nz2,
			    Common) ;
		    break ;

		/* ---------------------------------------------------------- */
		/* complex -> zomplex */
		/* ---------------------------------------------------------- */

		case CHOLMOD_ZOMPLEX:
		    /* allocate X and Z and copy the old X into them */
		    Xnew = CHOLMOD(malloc) (nz, sizeof (double), Common) ;
		    Znew = CHOLMOD(malloc) (nz, sizeof (double), Common) ;
		    if (Common->status < CHOLMOD_OK)
		    {
			CHOLMOD(free) (nz, sizeof (double), Xnew, Common) ;
			CHOLMOD(free) (nz, sizeof (double), Znew, Common) ;
			return (FALSE) ;
		    }
		    for (k = 0 ; k < nz ; k++)
		    {
			Xnew [k] = Xold [2*k  ] ;
			Znew [k] = Xold [2*k+1] ;
		    }
		    CHOLMOD(free) (nz, 2*sizeof (double), *XX, Common) ;
		    *XX = Xnew ;
		    *ZZ = Znew ;
		    break ;
	    }
	    break ;

	/* ------------------------------------------------------------------ */
	/* converting from zomplex */
	/* ------------------------------------------------------------------ */

	case CHOLMOD_ZOMPLEX:

	    switch (xtype_out)
	    {

		/* ---------------------------------------------------------- */
		/* zomplex -> pattern */
		/* ---------------------------------------------------------- */

		case CHOLMOD_PATTERN:
		    /* free X and Z */
		    *XX = CHOLMOD(free) (nz, sizeof (double), *XX, Common) ;
		    *ZZ = CHOLMOD(free) (nz, sizeof (double), *ZZ, Common) ;
		    break ;

		/* ---------------------------------------------------------- */
		/* zomplex -> real */
		/* ---------------------------------------------------------- */

		case CHOLMOD_REAL:
		    /* free the imaginary part */
		    *ZZ = CHOLMOD(free) (nz, sizeof (double), *ZZ, Common) ;
		    break ;

		/* ---------------------------------------------------------- */
		/* zomplex -> complex */
		/* ---------------------------------------------------------- */

		case CHOLMOD_COMPLEX:
		    Xnew = CHOLMOD(malloc) (nz, 2*sizeof (double), Common) ;
		    if (Common->status < CHOLMOD_OK)
		    {
			return (FALSE) ;
		    }
		    for (k = 0 ; k < nz ; k++)
		    {
			Xnew [2*k  ] = Xold [k] ;
			Xnew [2*k+1] = Zold [k] ;
		    }
		    CHOLMOD(free) (nz, sizeof (double), *XX, Common) ;
		    CHOLMOD(free) (nz, sizeof (double), *ZZ, Common) ;
		    *XX = Xnew ;
		    *ZZ = NULL ;
		    break ;

	    }
	    break ;
    }

    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_sparse_xtype ================================================= */
/* ========================================================================== */

/* Change the numeric xtype of a sparse matrix.  Supports any type on input
 * and output (pattern, real, complex, or zomplex). */

int CHOLMOD(sparse_xtype)
(
    /* ---- input ---- */
    int to_xtype,	/* requested xtype */
    /* ---- in/out --- */
    cholmod_sparse *A,	/* sparse matrix to change */
    /* --------------- */
    cholmod_common *Common
)
{
    Int ok ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;

    ok = change_complexity (A->nzmax, A->xtype, to_xtype,
	    CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, &(A->x), &(A->z), Common) ;
    if (ok)
    {
	A->xtype = to_xtype ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_triplet_xtype ================================================ */
/* ========================================================================== */

/* Change the numeric xtype of a triplet matrix.  Supports any type on input
 * and output (pattern, real, complex, or zomplex). */

int CHOLMOD(triplet_xtype)
(
    /* ---- input ---- */
    int to_xtype,	/* requested xtype */
    /* ---- in/out --- */
    cholmod_triplet *T,	/* triplet matrix to change */
    /* --------------- */
    cholmod_common *Common
)
{
    Int ok ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (T, FALSE) ;
    RETURN_IF_XTYPE_INVALID (T, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    ok = change_complexity (T->nzmax, T->xtype, to_xtype,
	    CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, &(T->x), &(T->z), Common) ;
    if (ok)
    {
	T->xtype = to_xtype ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_dense_xtype ================================================= */
/* ========================================================================== */

/* Change the numeric xtype of a dense matrix.  Supports real, complex or
 * zomplex on input and output */

int CHOLMOD(dense_xtype)
(
    /* ---- input ---- */
    int to_xtype,	/* requested xtype */
    /* ---- in/out --- */
    cholmod_dense *X,	/* dense matrix to change */
    /* --------------- */
    cholmod_common *Common
)
{
    Int ok ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (X, FALSE) ;
    RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    ok = change_complexity (X->nzmax, X->xtype, to_xtype,
	    CHOLMOD_REAL, CHOLMOD_ZOMPLEX, &(X->x), &(X->z), Common) ;
    if (ok)
    {
	X->xtype = to_xtype ;
    }
    return (ok) ;
}


/* ========================================================================== */
/* === cholmod_factor_xtype ================================================= */
/* ========================================================================== */

/* Change the numeric xtype of a factor.  Supports real, complex or zomplex on
 * input and output.   Supernodal zomplex factors are not supported. */

int CHOLMOD(factor_xtype)
(
    /* ---- input ---- */
    int to_xtype,	/* requested xtype */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to change */
    /* --------------- */
    cholmod_common *Common
)
{
    Int ok ;
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    if (L->is_super &&
	    (L->xtype == CHOLMOD_ZOMPLEX || to_xtype == CHOLMOD_ZOMPLEX))
    {
	ERROR (CHOLMOD_INVALID, "invalid xtype for supernodal L") ;
	return (FALSE) ;
    }
    ok = change_complexity ((L->is_super ? L->xsize : L->nzmax), L->xtype,
	    to_xtype, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, &(L->x), &(L->z), Common) ;
    if (ok)
    {
	L->xtype = to_xtype ;
    }
    return (ok) ;
}
