/* ========================================================================== */
/* === Core/cholmod_change_factor =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Change the numeric/symbolic, LL/LDL, simplicial/super, packed/unpacked,
 * monotonic/non-monotonic status of a cholmod_factor object.
 *
 * There are four basic classes of factor types:
 *
 * (1) simplicial symbolic:  Consists of two size-n arrays: the fill-reducing
 *	permutation (L->Perm) and the nonzero count for each column of L
 *	(L->ColCount).  All other factor types also include this information.
 *	L->ColCount may be exact (obtained from the analysis routines), or
 *	it may be a guess.  During factorization, and certainly after update/
 *	downdate, the columns of L can have a different number of nonzeros.
 *	L->ColCount is used to allocate space.  L->ColCount is exact for the
 *	supernodal factorizations.  The nonzero pattern of L is not kept.
 *
 * (2) simplicial numeric:  These represent L in a compressed column form.  The
 *	variants of this type are:
 *
 *	LDL':	L is unit diagonal.  Row indices in column j are located in
 *	    L->i [L->p [j] ... L->p [j] + L->nz [j]], and corresponding numeric
 *	    values are in the same locations in L->x.  The total number of
 *	    entries is the sum of L->nz [j].  The unit diagonal is not stored;
 *	    D is stored on the diagonal of L instead.  L->p may or may not be
 *	    monotonic.  The order of storage of the columns in L->i and L->x is
 *	    given by a doubly-linked list (L->prev and L->next).  L->p is of
 *	    size n+1, but only the first n entries are used (it is used if L
 *	    is converted to a sparse matrix via cholmod_factor_to_sparse).
 *
 *	    For the complex case, L->x is stored interleaved with real/imag
 *	    parts, and is of size 2*lnz*sizeof(double).  For the zomplex case,
 *	    L->x is of size lnz*sizeof(double) and holds the real part; L->z
 *	    is the same size and holds the imaginary part.
 *
 *	LL':  This is identical to the LDL' form, except that the non-unit
 *	    diagonal of L is stored as the first entry in each column of L.
 *
 * (3) supernodal symbolic:  A representation of the nonzero pattern of the
 *	supernodes for a supernodal factorization.  There are L->nsuper
 *	supernodes.  Columns L->super [k] to L->super [k+1]-1 are in the kth
 *	supernode.  The row indices for the kth supernode are in
 *	L->s [L->pi [k] ... L->pi [k+1]-1].  The numerical values are not
 *	allocated (L->x), but when they are they will be located in
 *	L->x [L->px [k] ... L->px [k+1]-1], and the L->px array is defined
 *	in this factor type.
 *
 *	For the complex case, L->x is stored interleaved with real/imag parts,
 *	and is of size 2*L->xsize*sizeof(double).  The zomplex supernodal case
 *	is not supported, since it is not compatible with LAPACK and the BLAS.
 *
 * (4) supernodal numeric:  Always an LL' factorization.  L is non-unit
 *      diagonal.  L->x contains the numerical values of the supernodes, as
 *      described above for the supernodal symbolic factor.
 *	For the complex case, L->x is stored interleaved, and is of size
 *	2*L->xsize*sizeof(double).  The zomplex supernodal case is not
 *	supported, since it is not compatible with LAPACK and the BLAS.
 *
 *      FUTURE WORK: support a supernodal LDL' factor.
 *
 *
 * In all cases, the row indices in each column (L->i for simplicial L and
 * L->s for supernodal L) are kept sorted from low indices to high indices.
 * This means the diagonal of L (or D for LDL' factors) is always kept as the
 * first entry in each column.
 *
 * The cholmod_change_factor routine can do almost all possible conversions.
 * It cannot do the following conversions:
 *
 *	(1) Simplicial numeric types cannot be converted to a supernodal
 *	    symbolic type.  This would simultaneously deallocate the
 *	    simplicial pattern and numeric values and reallocate uninitialized
 *	    space for the supernodal pattern.  This isn't useful for the user,
 *	    and not needed by CHOLMOD's own routines either.
 *
 *	(2) Only a symbolic factor (simplicial to supernodal) can be converted
 *	    to a supernodal numeric factor.
 *
 * Some conversions are meant only to be used internally by other CHOLMOD
 * routines, and should not be performed by the end user.  They allocate space
 * whose contents are undefined:
 *
 *	(1) converting from simplicial symbolic to supernodal symbolic.
 *	(2) converting any factor to supernodal numeric.
 *
 * workspace: no conversion routine uses workspace in Common.  No temporary
 *	workspace is allocated.
 *
 * Supports all xtypes, except that there is no supernodal zomplex L.
 *
 * The to_xtype parameter is used only when converting from symbolic to numeric
 * or numeric to symbolic.  It cannot be used to convert a numeric xtype (real,
 * complex, or zomplex) to a different numeric xtype.  For that conversion,
 * use cholmod_factor_xtype instead.
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"

static void natural_list (cholmod_factor *L) ;

/* ========================================================================== */
/* === TEMPLATE ============================================================= */
/* ========================================================================== */

#define REAL
#include "t_cholmod_change_factor.c"
#define COMPLEX
#include "t_cholmod_change_factor.c"
#define ZOMPLEX
#include "t_cholmod_change_factor.c"


/* ========================================================================== */
/* === L_is_packed ========================================================== */
/* ========================================================================== */

/* Return TRUE if the columns of L are packed, FALSE otherwise.  For debugging
 * only. */

#ifndef NDEBUG
static int L_is_packed (cholmod_factor *L, cholmod_common *Common)
{
    Int j ;
    Int *Lnz = L->nz ;
    Int *Lp = L->p ;
    Int n = L->n ;

    if (L->xtype == CHOLMOD_PATTERN || L->is_super)
    {
	return (TRUE) ;
    }

    if (Lnz == NULL || Lp == NULL)
    {
	return (TRUE) ;
    }

    for (j = 0 ; j < n ; j++)
    {
	PRINT3 (("j: "ID" Lnz "ID" Lp[j+1] "ID" Lp[j] "ID"\n", j, Lnz [j],
		Lp [j+1], Lp [j])) ;
	if (Lnz [j] != (Lp [j+1] - Lp [j]))
	{
	    PRINT2 (("L is not packed\n")) ;
	    return (FALSE) ;
	}
    }
    return (TRUE) ;
}
#endif


/* ========================================================================== */
/* === natural_list ========================================================= */
/* ========================================================================== */

/* Create a naturally-ordered doubly-linked list of columns. */

static void natural_list (cholmod_factor *L)
{
    Int head, tail, n, j ;
    Int *Lnext, *Lprev ;
    Lnext = L->next ;
    Lprev = L->prev ;
    ASSERT (Lprev != NULL && Lnext != NULL) ;
    n = L->n ;
    head = n+1 ;
    tail = n ;
    Lnext [head] = 0 ;
    Lprev [head] = EMPTY ;
    Lnext [tail] = EMPTY ;
    Lprev [tail] = n-1 ;
    for (j = 0 ; j < n ; j++)
    {
	Lnext [j] = j+1 ;
	Lprev [j] = j-1 ;
    }
    Lprev [0] = head ;
    L->is_monotonic = TRUE ;
}


/* ========================================================================== */
/* === allocate_simplicial_numeric ========================================== */
/* ========================================================================== */

/* Allocate O(n) arrays for simplicial numeric factorization.  Initializes
 * the link lists only.  Does not allocate the L->i, L->x, or L->z arrays. */

static int allocate_simplicial_numeric
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    Int n ;
    Int *Lp, *Lnz, *Lprev, *Lnext ;
    size_t n1, n2 ;

    PRINT1 (("Allocate simplicial\n")) ;

    ASSERT (L->xtype == CHOLMOD_PATTERN || L->is_super) ;
    ASSERT (L->p == NULL) ;
    ASSERT (L->nz == NULL) ;
    ASSERT (L->prev == NULL) ;
    ASSERT (L->next == NULL) ;

    n = L->n ;

    /* this cannot cause size_t overflow */
    n1 = ((size_t) n) + 1 ;
    n2 = ((size_t) n) + 2 ;

    Lp = CHOLMOD(malloc) (n1, sizeof (Int), Common) ;
    Lnz = CHOLMOD(malloc) (n, sizeof (Int), Common) ;
    Lprev = CHOLMOD(malloc) (n2, sizeof (Int), Common) ;
    Lnext = CHOLMOD(malloc) (n2, sizeof (Int), Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free) (n1, sizeof (Int), Lp,    Common) ;
	CHOLMOD(free) (n,   sizeof (Int), Lnz,   Common) ;
	CHOLMOD(free) (n2, sizeof (Int), Lprev, Common) ;
	CHOLMOD(free) (n2, sizeof (Int), Lnext, Common) ;
	PRINT1 (("Allocate simplicial failed\n")) ;
	return (FALSE) ;	/* out of memory */
    }

    /* ============================================== commit the changes to L */

    L->p = Lp ;
    L->nz = Lnz ;
    L->prev = Lprev ;
    L->next = Lnext ;
    /* initialize a doubly linked list for columns in natural order */
    natural_list (L) ;
    PRINT1 (("Allocate simplicial done\n")) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === simplicial_symbolic_to_super_symbolic ================================ */
/* ========================================================================== */

/* Convert a simplicial symbolic factor supernodal symbolic factor.  Does not
 * initialize the new space. */

static int simplicial_symbolic_to_super_symbolic
(
    cholmod_factor *L,
    cholmod_common *Common
)
{
    Int nsuper, xsize, ssize ;
    Int *Lsuper, *Lpi, *Lpx, *Ls ;
    size_t nsuper1 ;

    ASSERT (L->xtype == CHOLMOD_PATTERN && !(L->is_super)) ;

    xsize  = L->xsize ;
    ssize  = L->ssize ;
    nsuper = L->nsuper ;
    nsuper1 = ((size_t) nsuper) + 1 ;

    PRINT1 (("simple sym to super sym: ssize "ID" xsize "ID" nsuper "ID""
	" status %d\n", ssize, xsize, nsuper, Common->status)) ;

    /* O(nsuper) arrays, where nsuper <= n */
    Lsuper = CHOLMOD(malloc) (nsuper1, sizeof (Int), Common) ;
    Lpi    = CHOLMOD(malloc) (nsuper1, sizeof (Int), Common) ;
    Lpx    = CHOLMOD(malloc) (nsuper1, sizeof (Int), Common) ;

    /* O(ssize) array, where ssize <= nnz(L), and usually much smaller */
    Ls = CHOLMOD(malloc) (ssize, sizeof (Int), Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free) (nsuper1, sizeof (Int), Lsuper, Common) ;
	CHOLMOD(free) (nsuper1, sizeof (Int), Lpi,    Common) ;
	CHOLMOD(free) (nsuper1, sizeof (Int), Lpx,    Common) ;
	CHOLMOD(free) (ssize,    sizeof (Int), Ls,     Common) ;
	return (FALSE) ;	/* out of memory */
    }

    /* ============================================== commit the changes to L */

    ASSERT (Lsuper != NULL && Lpi != NULL && Lpx != NULL && Ls != NULL) ;

    L->maxcsize = 0 ;
    L->maxesize = 0 ;

    L->super = Lsuper ;
    L->pi = Lpi ;
    L->px = Lpx ;
    L->s  = Ls ;
    Ls [0] = EMPTY ;	    /* supernodal pattern undefined */

    L->is_super = TRUE ;
    L->is_ll = TRUE ;	    /* supernodal LDL' not supported */
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = DTYPE ;
    L->minor = L->n ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === any_to_simplicial_symbolic =========================================== */
/* ========================================================================== */

/* Convert any factor L to a simplicial symbolic factor, leaving only L->Perm
 * and L->ColCount.  Cannot fail.  Any of the components of L (except Perm and
 * ColCount) may already be free'd.  */

static void any_to_simplicial_symbolic
(
    cholmod_factor *L,
    int to_ll,
    cholmod_common *Common
)
{
    Int n, lnz, xs, ss, s, e ;
    size_t n1, n2 ;

    /* ============================================== commit the changes to L */

    n = L->n ;
    lnz = L->nzmax ;
    s = L->nsuper + 1 ;
    xs = (L->is_super) ? ((Int) (L->xsize)) : (lnz) ;
    e = (L->xtype == CHOLMOD_COMPLEX ? 2 : 1) ;
    ss = L->ssize ;

    /* this cannot cause size_t overflow */
    n1 = ((size_t) n) + 1 ;
    n2 = ((size_t) n) + 2 ;

    /* free all but the symbolic analysis (Perm and ColCount) */
    L->p     = CHOLMOD(free) (n1,  sizeof (Int),      L->p,     Common) ;
    L->i     = CHOLMOD(free) (lnz, sizeof (Int),      L->i,     Common) ;
    L->x     = CHOLMOD(free) (xs,  e*sizeof (double), L->x,     Common) ;
    L->z     = CHOLMOD(free) (lnz, sizeof (double),   L->z,     Common) ;
    L->nz    = CHOLMOD(free) (n,   sizeof (Int),      L->nz,    Common) ;
    L->next  = CHOLMOD(free) (n2,  sizeof (Int),      L->next,  Common) ;
    L->prev  = CHOLMOD(free) (n2,  sizeof (Int),      L->prev,  Common) ;
    L->super = CHOLMOD(free) (s,   sizeof (Int),      L->super, Common) ;
    L->pi    = CHOLMOD(free) (s,   sizeof (Int),      L->pi,    Common) ;
    L->px    = CHOLMOD(free) (s,   sizeof (Int),      L->px,    Common) ;
    L->s     = CHOLMOD(free) (ss,  sizeof (Int),      L->s,     Common) ;
    L->nzmax = 0 ;
    L->is_super = FALSE ;
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = DTYPE ;
    L->minor = n ;
    L->is_ll = to_ll ;
}


/* ========================================================================== */
/* === ll_super_to_super_symbolic =========================================== */
/* ========================================================================== */

/* Convert a numerical supernodal L to symbolic supernodal.  Cannot fail. */

static void ll_super_to_super_symbolic
(
    cholmod_factor *L,
    cholmod_common *Common
)
{

    /* ============================================== commit the changes to L */

    /* free all but the supernodal numerical factor */
    ASSERT (L->xtype != CHOLMOD_PATTERN && L->is_super && L->is_ll) ;
    DEBUG (CHOLMOD(dump_factor) (L, "start to super symbolic", Common)) ;
    L->x = CHOLMOD(free) (L->xsize,
	    (L->xtype == CHOLMOD_COMPLEX ? 2 : 1) * sizeof (double), L->x,
	    Common) ;
    L->xtype = CHOLMOD_PATTERN ;
    L->dtype = DTYPE ;
    L->minor = L->n ;
    L->is_ll = TRUE ;	    /* supernodal LDL' not supported */
    DEBUG (CHOLMOD(dump_factor) (L, "done  to super symbolic", Common)) ;
}


/* ========================================================================== */
/* === simplicial_symbolic_to_simplicial_numeric ============================ */
/* ========================================================================== */

/* Convert a simplicial symbolic L to a simplicial numeric L; allocate space
 * for L using L->ColCount from symbolic analysis, and set L to identity.
 *
 * If packed < 0, then this routine is creating a copy of another factor
 * (via cholmod_copy_factor).  In this case, the space is not initialized. */

static void simplicial_symbolic_to_simplicial_numeric
(
    cholmod_factor *L,
    int to_ll,
    int packed,
    int to_xtype,
    cholmod_common *Common
)
{
    double grow0, grow1, xlen, xlnz ;
    double *Lx, *Lz ;
    Int *Li, *Lp, *Lnz, *ColCount ;
    Int n, grow, grow2, p, j, lnz, len, ok, e ;

    ASSERT (L->xtype == CHOLMOD_PATTERN && !(L->is_super)) ;
    if (!allocate_simplicial_numeric (L, Common))
    {
	PRINT1 (("out of memory, allocate simplicial numeric\n")) ;
	return ;	/* out of memory */
    }
    ASSERT (L->ColCount != NULL && L->nz != NULL && L->p != NULL) ;
    ASSERT (L->x == NULL && L->z == NULL && L->i == NULL) ;

    ColCount = L->ColCount ;
    Lnz = L->nz ;
    Lp = L->p ;
    ok = TRUE ;
    n = L->n ;

    if (packed < 0)
    {

	/* ------------------------------------------------------------------ */
	/* used by cholmod_copy_factor to allocate a copy of a factor object */
	/* ------------------------------------------------------------------ */

	lnz = L->nzmax ;
	L->nzmax = 0 ;

    }
    else if (packed)
    {

	/* ------------------------------------------------------------------ */
	/* LDL' or LL' packed */
	/* ------------------------------------------------------------------ */

	PRINT1 (("convert to packed LL' or LDL'\n")) ;
	lnz = 0 ;
	for (j = 0 ; ok && j < n ; j++)
	{
	    /* ensure len is in the range 1 to n-j */
	    len = ColCount [j] ;
	    len = MAX (1, len) ;
	    len = MIN (len, n-j) ;
	    lnz += len ;
	    ok = (lnz >= 0) ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Lp [j] = j ;
	}
	for (j = 0 ; j < n ; j++)
	{
	    Lnz [j] = 1 ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* LDL' unpacked */
	/* ------------------------------------------------------------------ */

	PRINT1 (("convert to unpacked\n")) ;
	/* compute new lnzmax */
	/* if any parameter is NaN, grow is false */
	grow0 = Common->grow0 ;
	grow1 = Common->grow1 ;
	grow2 = Common->grow2 ;
	grow0 = IS_NAN (grow0) ? 1 : grow0 ;
	grow1 = IS_NAN (grow1) ? 1 : grow1 ;
	/* fl.pt. compare, but no NaN's: */
	grow = (grow0 >= 1.0) && (grow1 >= 1.0) && (grow2 > 0) ;
	PRINT1 (("init, grow1 %g grow2 "ID"\n", grow1, grow2)) ;
	/* initialize Lp and Lnz for each column */
	lnz = 0 ;
	for (j = 0 ; ok && j < n ; j++)
	{
	    Lp [j] = lnz ;
	    Lnz [j] = 1 ;

	    /* ensure len is in the range 1 to n-j */
	    len = ColCount [j] ;
	    len = MAX (1, len) ;
	    len = MIN (len, n-j) ;

	    /* compute len in double to avoid integer overflow */
	    PRINT1 (("ColCount ["ID"] = "ID"\n", j, len)) ;
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (Int) xlen ;
	    }
	    ASSERT (len >= 1 && len <= n-j) ;
	    lnz += len ;
	    ok = (lnz >= 0) ;
	}
	if (ok)
	{
	    Lp [n] = lnz ;
	    if (grow)
	    {
		/* add extra space */
		xlnz = (double) lnz ;
		xlnz *= grow0 ;
		xlnz = MIN (xlnz, Size_max) ;
		xlnz = MIN (xlnz, ((double) n * (double) n + (double) n) / 2) ;
		lnz = (Int) xlnz ;
	    }
	}
    }

    lnz = MAX (1, lnz) ;

    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
    }

    /* allocate L->i, L->x, and L->z */
    PRINT1 (("resizing from zero size to lnz "ID"\n", lnz)) ;
    ASSERT (L->nzmax == 0) ;
    e = (to_xtype == CHOLMOD_COMPLEX ? 2 : 1) ;
    if (!ok || !CHOLMOD(realloc_multiple) (lnz, 1, to_xtype, &(L->i), NULL,
		&(L->x), &(L->z), &(L->nzmax), Common))
    {
	L->p    = CHOLMOD(free) (n+1, sizeof (Int),      L->p, Common) ;
	L->nz   = CHOLMOD(free) (n,   sizeof (Int),      L->nz, Common) ;
	L->prev = CHOLMOD(free) (n+2, sizeof (Int),      L->prev, Common) ;
	L->next = CHOLMOD(free) (n+2, sizeof (Int),      L->next, Common) ;
	L->i    = CHOLMOD(free) (lnz, sizeof (Int),      L->i, Common) ;
	L->x    = CHOLMOD(free) (lnz, e*sizeof (double), L->x, Common) ;
	L->z    = CHOLMOD(free) (lnz, sizeof (double),   L->z, Common) ;
	PRINT1 (("cannot realloc simplicial numeric\n")) ;
	return ;	/* out of memory */
    }

    /* ============================================== commit the changes to L */

    /* initialize L to be the identity matrix */
    L->xtype = to_xtype ;
    L->dtype = DTYPE ;
    L->minor = n ;

    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;

#if 0
    if (lnz == 1)
    {
	/* the user won't expect to access this entry, but some CHOLMOD
	 * routines may.  Set it to zero so that valgrind doesn't complain. */
	switch (to_xtype)
	{
	    case CHOLMOD_REAL:
		Lx [0] = 0 ;
		break ;

	    case CHOLMOD_COMPLEX:
		Lx [0] = 0 ;
		Lx [1] = 0 ;
		break ;

	    case CHOLMOD_ZOMPLEX:
		Lx [0] = 0 ;
		Lz [0] = 0 ;
		break ;
	}
    }
#endif

    if (packed >= 0)
    {
	/* create the unit diagonal for either the LL' or LDL' case */

	switch (L->xtype)
	{
	    case CHOLMOD_REAL:
		for (j = 0 ; j < n ; j++)
		{
		    ASSERT (Lp [j] < Lp [j+1]) ;
		    p = Lp [j] ;
		    Li [p] = j ;
		    Lx [p] = 1 ;
		}
		break ;

	    case CHOLMOD_COMPLEX:
		for (j = 0 ; j < n ; j++)
		{
		    ASSERT (Lp [j] < Lp [j+1]) ;
		    p = Lp [j] ;
		    Li [p] = j ;
		    Lx [2*p  ] = 1 ;
		    Lx [2*p+1] = 0 ;
		}
		break ;

	    case CHOLMOD_ZOMPLEX:
		for (j = 0 ; j < n ; j++)
		{
		    ASSERT (Lp [j] < Lp [j+1]) ;
		    p = Lp [j] ;
		    Li [p] = j ;
		    Lx [p] = 1 ;
		    Lz [p] = 0 ;
		}
		break ;
	}
    }

    L->is_ll = to_ll ;

    PRINT1 (("done convert simplicial symbolic to numeric\n")) ;
}


/* ========================================================================== */
/* === change_simplicial_numeric ============================================ */
/* ========================================================================== */

/* Change LL' to LDL', LDL' to LL', or leave as-is.
 *
 * If to_packed is TRUE, then the columns of L are packed and made monotonic
 * (to_monotonic is ignored; it is implicitly TRUE).
 *
 * If to_monotonic is TRUE but to_packed is FALSE, the columns of L are made
 * monotonic but not packed.
 *
 * If both to_packed and to_monotonic are FALSE, then the columns of L are
 * left as-is, and the conversion is done in place.
 *
 * If L is already monotonic, or if it is to be left non-monotonic, then this
 * conversion always succeeds.
 *
 * When converting an LDL' to LL' factorization, any column with a negative
 * or zero diagonal entry is not modified so that conversion back to LDL' will
 * succeed.  This can result in a matrix L with a negative entry on the diagonal
 * If the kth entry on the diagonal of D is negative, it and the kth column of
 * L are left unchanged.  A subsequent conversion back to an LDL' form will also
 * leave the column unchanged, so the correct LDL' factorization will be
 * restored.  L->minor is set to the smallest k for which D (k,k) is negative.
 */

static void change_simplicial_numeric
(
    cholmod_factor *L,
    int to_ll,
    int to_packed,
    int to_monotonic,
    cholmod_common *Common
)
{
    double grow0, grow1, xlen, xlnz ;
    void *newLi, *newLx, *newLz ;
    double *Lx, *Lz ;
    Int *Lp, *Li, *Lnz ;
    Int make_monotonic, grow2, n, j, lnz, len, grow, ok, make_ll, make_ldl ;
    size_t nzmax0 ;

    PRINT1 (("\n===Change simplicial numeric: %d %d %d\n",
	    to_ll, to_packed, to_monotonic)) ;
    DEBUG (CHOLMOD(dump_factor) (L, "change simplicial numeric", Common)) ;
    ASSERT (L->xtype != CHOLMOD_PATTERN && !(L->is_super)) ;

    make_monotonic = ((to_packed || to_monotonic) && !(L->is_monotonic)) ;
    make_ll  = (to_ll && !(L->is_ll)) ;
    make_ldl = (!to_ll && L->is_ll) ;

    n = L->n ;
    Lp = L->p ;
    Li = L->i ;
    Lx = L->x ;
    Lz = L->z ;
    Lnz = L->nz ;

    grow = FALSE ;
    grow0 = Common->grow0 ;
    grow1 = Common->grow1 ;
    grow2 = Common->grow2 ;
    grow0 = IS_NAN (grow0) ? 1 : grow0 ;
    grow1 = IS_NAN (grow1) ? 1 : grow1 ;
    ok = TRUE ;
    newLi = NULL ;
    newLx = NULL ; 
    newLz = NULL ; 
    lnz = 0 ;

    if (make_monotonic)
    {

	/* ------------------------------------------------------------------ */
	/* Columns out of order, but will be reordered and optionally packed. */
	/* ------------------------------------------------------------------ */

	PRINT1 (("L is non-monotonic\n")) ;

	/* compute new L->nzmax */
	if (!to_packed)
	{
	    /* if any parameter is NaN, grow is false */
	    /* fl.pt. comparisons below are false if any parameter is NaN */
	    grow = (grow0 >= 1.0) && (grow1 >= 1.0) && (grow2 > 0) ;
	}
	for (j = 0 ; ok && j < n ; j++)
	{
	    len = Lnz [j] ;
	    ASSERT (len >= 1 && len <= n-j) ;

	    /* compute len in double to avoid integer overflow */
	    if (grow)
	    {
		xlen = (double) len ;
		xlen = grow1 * xlen + grow2 ;
		xlen = MIN (xlen, n-j) ;
		len = (Int) xlen ;
	    }
	    ASSERT (len >= Lnz [j] && len <= n-j) ;

	    PRINT2 (("j: "ID" Lnz[j] "ID" len "ID" p "ID"\n",
			j, Lnz [j], len, lnz)) ;

	    lnz += len ;
	    ok = (lnz >= 0) ;
	}

	if (!ok)
	{
	    ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	    return ;
	}

	if (grow)
	{
	    xlnz = (double) lnz ;
	    xlnz *= grow0 ;
	    xlnz = MIN (xlnz, Size_max) ;
	    xlnz = MIN (xlnz, ((double) n * (double) n + (double) n) / 2) ;
	    lnz = (Int) xlnz ;
	}

	lnz = MAX (1, lnz) ;
	PRINT1 (("final lnz "ID"\n", lnz)) ;
	nzmax0 = 0 ;

	CHOLMOD(realloc_multiple) (lnz, 1, L->xtype, &newLi, NULL,
		&newLx, &newLz, &nzmax0, Common) ;

	if (Common->status < CHOLMOD_OK)
	{
	    return ;	    /* out of memory */
	}
    }

    /* ============================================== commit the changes to L */

    /* ---------------------------------------------------------------------- */
    /* convert the simplicial L, using template routine */
    /* ---------------------------------------------------------------------- */

    switch (L->xtype)
    {

	case CHOLMOD_REAL:
	    r_change_simplicial_numeric (L, to_ll, to_packed,
		    newLi, newLx, newLz, lnz, grow, grow1, grow2,
		    make_ll, make_monotonic, make_ldl, Common) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    c_change_simplicial_numeric (L, to_ll, to_packed,
		    newLi, newLx, newLz, lnz, grow, grow1, grow2,
		    make_ll, make_monotonic, make_ldl, Common) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    z_change_simplicial_numeric (L, to_ll, to_packed,
		    newLi, newLx, newLz, lnz, grow, grow1, grow2,
		    make_ll, make_monotonic, make_ldl, Common) ;
	    break ;
    }

    DEBUG (CHOLMOD(dump_factor) (L, "L simplicial changed", Common)) ;
}


/* ========================================================================== */
/* === ll_super_to_simplicial_numeric ======================================= */
/* ========================================================================== */

/* Convert a supernodal numeric factorization to any simplicial numeric one.
 * Leaves L->xtype unchanged (real or complex, not zomplex since there is
 * no supernodal zomplex L). */

static void ll_super_to_simplicial_numeric
(
    cholmod_factor *L,
    int to_packed,
    int to_ll,
    cholmod_common *Common
)
{
    Int *Ls, *Lpi, *Lpx, *Super, *Li ;
    Int n, lnz, s, nsuper, psi, psend, nsrow, nscol, k1, k2, erows ;

    DEBUG (CHOLMOD(dump_factor) (L, "start LL super to simplicial", Common)) ;
    PRINT1 (("super -> simplicial (%d %d)\n", to_packed, to_ll)) ;
    ASSERT (L->xtype != CHOLMOD_PATTERN && L->is_ll && L->is_super) ;
    ASSERT (L->x != NULL && L->i == NULL) ;

    n = L->n ;
    nsuper = L->nsuper ;
    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Super = L->super ;

    /* Int overflow cannot occur since supernodal L already exists */

    if (to_packed)
    {
	/* count the number of nonzeros in L.  Each supernode is of the form
	 *
	 *    l	. . .	    For this example, nscol = 4 (# columns). nsrow = 9.
	 *    l l . .	    The "." entries are allocated in the supernodal
	 *    l l l .	    factor, but not used.  They are not copied to the
	 *    l l l l	    simplicial factor.  Some "l" and "e" entries may be
	 *    e e e e	    numerically zero and even symbolically zero if a
	 *    e e e e	    tight simplicial factorization or resymbol were
	 *    e e e e	    done, because of numerical cancellation and relaxed
	 *    e e e e	    supernode amalgamation, respectively.
	 *    e e e e
	 */
	lnz = 0 ;
	for (s = 0 ; s < nsuper ; s++)
	{
	    k1 = Super [s] ;
	    k2 = Super [s+1] ;
	    psi = Lpi [s] ;
	    psend = Lpi [s+1] ;
	    nsrow = psend - psi ;
	    nscol = k2 - k1 ;
	    ASSERT (nsrow >= nscol) ;
	    erows = nsrow - nscol ;

	    /* lower triangular part, including the diagonal,
	     * counting the "l" terms in the figure above. */
	    lnz += nscol * (nscol+1) / 2 ;

	    /* rectangular part, below the diagonal block (the "e" terms) */
	    lnz += nscol * erows ;
	}
	ASSERT (lnz <= (Int) (L->xsize)) ;
    }
    else
    {
	/* Li will be the same size as Lx */
	lnz = L->xsize ;
    }
    ASSERT (lnz >= 0) ;
    PRINT1 (("simplicial lnz = "ID"  to_packed: %d  to_ll: %d L->xsize %g\n",
		lnz, to_ll, to_packed, (double) L->xsize)) ;

    Li = CHOLMOD(malloc) (lnz, sizeof (Int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return ;	/* out of memory */
    }

    if (!allocate_simplicial_numeric (L, Common))
    {
	CHOLMOD(free) (lnz, sizeof (Int), Li, Common) ;
	return ;	/* out of memory */
    }

    /* ============================================== commit the changes to L */

    L->i = Li ;
    L->nzmax = lnz ;

    /* ---------------------------------------------------------------------- */
    /* convert the supernodal L, using template routine */
    /* ---------------------------------------------------------------------- */

    switch (L->xtype)
    {

	case CHOLMOD_REAL:
	    r_ll_super_to_simplicial_numeric (L, to_packed, to_ll, Common) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    c_ll_super_to_simplicial_numeric (L, to_packed, to_ll, Common) ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* free unused parts of L */
    /* ---------------------------------------------------------------------- */

    L->super = CHOLMOD(free) (nsuper+1, sizeof (Int), L->super, Common) ;
    L->pi    = CHOLMOD(free) (nsuper+1, sizeof (Int), L->pi, Common) ;
    L->px    = CHOLMOD(free) (nsuper+1, sizeof (Int), L->px, Common) ;
    L->s     = CHOLMOD(free) (L->ssize, sizeof (Int), L->s, Common) ;

    L->ssize = 0 ;
    L->xsize = 0 ;
    L->nsuper = 0 ;
    L->maxesize = 0 ;
    L->maxcsize = 0 ;

    L->is_super = FALSE ;

    DEBUG (CHOLMOD(dump_factor) (L, "done  LL super to simplicial", Common)) ;
}


/* ========================================================================== */
/* === super_symbolic_to_ll_super =========================================== */
/* ========================================================================== */

/* Convert a supernodal symbolic factorization to a supernodal numeric
 * factorization by allocating L->x.  Contents of L->x are undefined.
 */

static int super_symbolic_to_ll_super
(
    int to_xtype,
    cholmod_factor *L,
    cholmod_common *Common
)
{
    double *Lx ;
    Int wentry = (to_xtype == CHOLMOD_REAL) ? 1 : 2 ;
    PRINT1 (("convert super sym to num\n")) ;
    ASSERT (L->xtype == CHOLMOD_PATTERN && L->is_super) ;
    Lx = CHOLMOD(malloc) (L->xsize, wentry * sizeof (double), Common) ;
    PRINT1 (("xsize %g\n", (double) L->xsize)) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    /* ============================================== commit the changes to L */

    if (L->xsize == 1)
    {
	/* the caller won't expect to access this entry, but some CHOLMOD
	 * routines may.  Set it to zero so that valgrind doesn't complain. */
	switch (to_xtype)
	{
	    case CHOLMOD_REAL:
		Lx [0] = 0 ;
		break ;

	    case CHOLMOD_COMPLEX:
		Lx [0] = 0 ;
		Lx [1] = 0 ;
		break ;
	}
    }

    L->x = Lx ;
    L->xtype = to_xtype ;
    L->dtype = DTYPE ;
    L->minor = L->n ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_change_factor ================================================ */
/* ========================================================================== */

/* Convert a factor L.  Some conversions simply allocate uninitialized space
 * that meant to be filled later.
 *
 * If the conversion fails, the factor is left in its original form, with one
 * exception.  Converting a supernodal symbolic factor to a simplicial numeric
 * one (with L=D=I) may leave the factor in simplicial symbolic form.
 *
 * Memory allocated for each conversion is listed below.
 */

int CHOLMOD(change_factor)
(
    /* ---- input ---- */
    int to_xtype,	/* convert to CHOLMOD_PATTERN, _REAL, _COMPLEX, or
			 * _ZOMPLEX */
    int to_ll,		/* TRUE: convert to LL', FALSE: LDL' */
    int to_super,	/* TRUE: convert to supernodal, FALSE: simplicial */
    int to_packed,	/* TRUE: pack simplicial columns, FALSE: do not pack */
    int to_monotonic,	/* TRUE: put simplicial columns in order, FALSE: not */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    if (to_xtype < CHOLMOD_PATTERN || to_xtype > CHOLMOD_ZOMPLEX)
    {
	ERROR (CHOLMOD_INVALID, "xtype invalid") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    PRINT1 (("-----convert from (%d,%d,%d,%d,%d) to (%d,%d,%d,%d,%d)\n",
    L->xtype, L->is_ll, L->is_super, L_is_packed (L, Common), L->is_monotonic,
    to_xtype, to_ll,    to_super,    to_packed,               to_monotonic)) ;

    /* ensure all parameters are TRUE/FALSE */
    to_ll = BOOLEAN (to_ll) ;
    to_super = BOOLEAN (to_super) ;

    ASSERT (BOOLEAN (L->is_ll) == L->is_ll) ;
    ASSERT (BOOLEAN (L->is_super) == L->is_super) ;

    if (to_super && to_xtype == CHOLMOD_ZOMPLEX)
    {
	ERROR (CHOLMOD_INVALID, "supernodal zomplex L not supported") ;
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* convert */
    /* ---------------------------------------------------------------------- */

    if (to_xtype == CHOLMOD_PATTERN)
    {

	/* ------------------------------------------------------------------ */
	/* convert to symbolic */
	/* ------------------------------------------------------------------ */

	if (!to_super)
	{

	    /* -------------------------------------------------------------- */
	    /* convert any factor into a simplicial symbolic factor */
	    /* -------------------------------------------------------------- */

	    any_to_simplicial_symbolic (L, to_ll, Common) ;    /* cannot fail */

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* convert to a supernodal symbolic factor */
	    /* -------------------------------------------------------------- */

	    if (L->xtype != CHOLMOD_PATTERN && L->is_super)
	    {
		/* convert from supernodal numeric to supernodal symbolic.
		 * this preserves symbolic pattern of L, discards numeric
		 * values */
		ll_super_to_super_symbolic (L, Common) ;       /* cannot fail */
	    }
	    else if (L->xtype == CHOLMOD_PATTERN && !(L->is_super))
	    {
		/* convert from simplicial symbolic to supernodal symbolic.
		 * contents of supernodal pattern are uninitialized.  Not meant
		 * for the end user. */
		simplicial_symbolic_to_super_symbolic (L, Common) ;
	    }
	    else
	    {
		/* cannot convert from simplicial numeric to supernodal
		 * symbolic */
		ERROR (CHOLMOD_INVALID,
			"cannot convert L to supernodal symbolic") ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* convert to numeric */
	/* ------------------------------------------------------------------ */
	    
	if (to_super)
	{

	    /* -------------------------------------------------------------- */
	    /* convert to supernodal numeric factor */
	    /* -------------------------------------------------------------- */

	    if (L->xtype == CHOLMOD_PATTERN)
	    {
		if (L->is_super)
		{
		    /* Convert supernodal symbolic to supernodal numeric.
		     * Contents of supernodal numeric values are uninitialized.
		     * This is used by cholmod_super_numeric.  Not meant for
		     * the end user. */
		    super_symbolic_to_ll_super (to_xtype, L, Common) ;
		}
		else
		{
		    /* Convert simplicial symbolic to supernodal numeric.
		     * Contents not defined.  This is used by
		     * Core/cholmod_copy_factor only.  Not meant for the end
		     * user. */
		    if (!simplicial_symbolic_to_super_symbolic (L, Common))
		    {
			/* failure, convert back to simplicial symbolic */
			any_to_simplicial_symbolic (L, to_ll, Common) ;
		    }
		    else
		    {
			/* conversion to super symbolic OK, allocate numeric
			 * part */
			super_symbolic_to_ll_super (to_xtype, L, Common) ;
		    }
		}
	    }
	    else
	    {
		/* nothing to do if L is already in supernodal numeric form */
		if (!(L->is_super))
		{
		    ERROR (CHOLMOD_INVALID,
			"cannot convert simplicial L to supernodal") ;
		}
		/* FUTURE WORK: convert to/from supernodal LL' and LDL' */
	    }

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* convert any factor to simplicial numeric */
	    /* -------------------------------------------------------------- */

	    if (L->xtype == CHOLMOD_PATTERN && !(L->is_super))
	    {

		/* ---------------------------------------------------------- */
		/* convert simplicial symbolic to simplicial numeric (L=I,D=I)*/
		/* ---------------------------------------------------------- */

		simplicial_symbolic_to_simplicial_numeric (L, to_ll, to_packed,
			to_xtype, Common) ;

	    }
	    else if (L->xtype != CHOLMOD_PATTERN && L->is_super)
	    {

		/* ---------------------------------------------------------- */
		/* convert a supernodal LL' to simplicial numeric */
		/* ---------------------------------------------------------- */

		ll_super_to_simplicial_numeric (L, to_packed, to_ll, Common) ;

	    }
	    else if (L->xtype == CHOLMOD_PATTERN && L->is_super)
	    {

		/* ---------------------------------------------------------- */
		/* convert a supernodal symbolic to simplicial numeric (L=D=I)*/
		/* ---------------------------------------------------------- */

		any_to_simplicial_symbolic (L, to_ll, Common) ;
		/* if the following fails, it leaves the factor in simplicial
		 * symbolic form */
		simplicial_symbolic_to_simplicial_numeric (L, to_ll, to_packed,
			to_xtype, Common) ;

	    }
	    else
	    {

		/* ---------------------------------------------------------- */
		/* change a simplicial numeric factor */
		/* ---------------------------------------------------------- */

		/* change LL' to LDL', LDL' to LL', or leave as-is.  pack the
		 * columns of L, or leave as-is.  Ensure the columns are
		 * monotonic, or leave as-is. */

		change_simplicial_numeric (L, to_ll, to_packed, to_monotonic,
			Common) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    return (Common->status >= CHOLMOD_OK) ;
}
