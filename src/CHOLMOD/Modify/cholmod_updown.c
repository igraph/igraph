/* ========================================================================== */
/* === Modify/cholmod_updown ================================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Modify Module.
 * Copyright (C) 2005-2006, Timothy A. Davis and William W. Hager.
 * The CHOLMOD/Modify Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Updates/downdates the LDL' factorization (symbolic, then numeric), by
 * computing a new factorization of
 *
 *  	Lnew * Dnew * Lnew' = Lold * Dold * Lold' +/- C*C'
 *
 * C must be sorted.  It can be either packed or unpacked.  As in all CHOLMOD
 * routines, the columns of L are sorted on input, and also on output.
 *
 * If the factor is not an unpacked LDL' or dynamic LDL', it is converted
 * to an LDL' dynamic factor.  An unpacked LDL' factor may be updated, but if
 * any one column runs out of space, the factor is converted to an LDL'
 * dynamic one.  If the initial conversion fails, the factor is returned
 * unchanged.
 *
 * If memory runs out during the update, the factor is returned as a simplicial
 * symbolic factor.  That is, everything is freed except for the fill-reducing
 * ordering and its corresponding column counts (typically computed by
 * cholmod_analyze).
 *
 * Note that the fill-reducing permutation L->Perm is NOT used.  The row
 * indices of C refer to the rows of L, not A.  If your original system is
 * LDL' = PAP' (where P = L->Perm), and you want to compute the LDL'
 * factorization of A+CC', then you must permute C first.  That is:
 *
 *	PAP' = LDL'
 *	P(A+CC')P' = PAP'+PCC'P' = LDL' + (PC)(PC)' = LDL' + Cnew*Cnew'
 *	where Cnew = P*C.
 *
 * You can use the cholmod_submatrix routine in the MatrixOps module
 * to permute C, with:
 *
 * Cnew = cholmod_submatrix (C, L->Perm, L->n, NULL, -1, TRUE, TRUE, Common) ;
 *
 * Note that the sorted input parameter to cholmod_submatrix must be TRUE,
 * because cholmod_updown requires C with sorted columns.
 *
 * The system Lx=b can also be updated/downdated.  The old system was Lold*x=b.
 * The new system is Lnew*xnew = b + deltab.  The old solution x is overwritten
 * with xnew.  Note that as in the update/downdate of L itself, the fill-
 * reducing permutation L->Perm is not used.  x and b are in the permuted
 * ordering, not your original ordering.  x and b are n-by-1; this routine
 * does not handle multiple right-hand-sides.
 *
 * workspace: Flag (nrow), Head (nrow+1), W (maxrank*nrow), Iwork (nrow),
 * where maxrank is 2, 4, or 8.
 *
 * Only real matrices are supported.  A symbolic L is converted into a
 * numeric identity matrix.
 */

#ifndef NMODIFY

#include "cholmod_internal.h"
#include "cholmod_modify.h"


/* ========================================================================== */
/* === cholmod_updown ======================================================= */
/* ========================================================================== */

/* Compute the new LDL' factorization of LDL'+CC' (an update) or LDL'-CC'
 * (a downdate).  The factor object L need not be an LDL' factorization; it
 * is converted to one if it isn't. */

int CHOLMOD(updown)
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(updown_mask) (update, C, NULL, NULL, L, NULL, NULL,
	Common)) ;
}


/* ========================================================================== */
/* === cholmod_updown_solve ================================================= */
/* ========================================================================== */

/* Does the same as cholmod_updown, except that it also updates/downdates the
 * solution to Lx=b+DeltaB.  x and b must be n-by-1 dense matrices.  b is not
 * need as input to this routine, but a sparse change to b is (DeltaB).  Only
 * entries in DeltaB corresponding to columns modified in L are accessed; the
 * rest are ignored.
 */

int CHOLMOD(updown_solve)
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(updown_mask) (update, C, NULL, NULL, L, X, DeltaB,
	Common)) ;
}


/* ========================================================================== */
/* === Power2 =============================================================== */
/* ========================================================================== */

/* Power2 [i] is smallest power of 2 that is >= i (for i in range 0 to 8) */

static Int Power2 [ ] =
{
/*  0  1  2  3  4  5  6  7  8 */
    0, 1, 2, 4, 4, 8, 8, 8, 8
} ;

/* ========================================================================== */
/* === debug routines ======================================================= */
/* ========================================================================== */

#ifndef NDEBUG

static void dump_set (Int s, Int **Set_ps1, Int **Set_ps2, Int j, Int n,
	cholmod_common *Common)
{
    Int *p, len, i, ilast ;

    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return ;
    }

    len = Set_ps2 [s] - Set_ps1 [s] ;
    PRINT2 (("Set s: "ID" len: "ID":", s, len)) ;
    ASSERT (len > 0) ;
    ilast = j ;
    for (p = Set_ps1 [s] ; p < Set_ps2 [s] ; p++)
    {
	i = *p ;
	PRINT3 ((" "ID"", i)) ;
	ASSERT (i > ilast && i < n) ;
	ilast = i ;
    }
    PRINT3 (("\n")) ;
}

static void dump_col
(
    char *w, Int j, Int p1, Int p2, Int *Li, double *Lx, Int n,
    cholmod_common *Common
)
{
    Int p, row, lastrow ;

    if (CHOLMOD(dump) < -1)
    {
	/* no checks if debug level is -2 or less */
	return ;
    }

    PRINT3 (("\n\nDUMP COL==== j = "ID"  %s: p1="ID" p2="ID" \n", j, w, p1,p2));
    lastrow = -1 ;
    for (p = p1 ; p < p2 ; p++)
    {
	PRINT3 (("   "ID": ", p)) ;
	row = Li [p] ;
	PRINT3 ((""ID"  ", Li [p])) ;
	PRINT3 (("%g ", Lx [p])) ;
	PRINT3 (("\n")) ;
	ASSERT (row > lastrow && row < n) ;
	lastrow = row ;
    }
    ASSERT (p1 < p2) ;
    ASSERT (Li [p1] == j) ;
    PRINT3 (("\n")) ;
}
#endif


/* ========================================================================== */
/* === a path =============================================================== */
/* ========================================================================== */

/* A path is a set of nodes of the etree which are all affected by the same
 * columns of C. */

typedef struct Path_struct
{
    Int start ;		/* column at which to start, or EMPTY if initial */
    Int end ;		/* column at which to end, or EMPTY if initial */
    Int ccol ;		/* column of C to which path refers */
    Int parent ;	/* parent path */
    Int c ;		/* child of j along this path */
    Int next ;		/* next path in link list */
    Int rank ;		/* number of rank-1 paths merged onto this path */
    Int order ;		/* dfs order of this path */
    Int wfirst ;	/* first column of W to affect this path */
    Int pending ;	/* column at which the path is pending */
    Int botrow ;	/* for partial update/downdate of solution to Lx=b */ 

} Path_type ;


/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* Compute the DFS order of the set of paths.  This can be recursive because
 * there are at most 23 paths to sort: one for each column of C (8 at most),
 * and one for each node in a balanced binary tree with 8 leaves (15).
 * Stack overflow is thus not a problem.  */

static void dfs
(
    Path_type *Path,	/* the set of Paths */
    Int k,		/* the rank of the update/downdate */
    Int path,		/* which path to work on */
    Int *path_order,	/* the current path order */
    Int *w_order,	/* the current order of the columns of W */
    Int depth,
    Int npaths		/* total number of paths */
)
{
    Int c ;		/* child path */

    ASSERT (path >= 0 && path < npaths) ;
    if (path < k)
    {
	/* this is a leaf node, corresponding to column W (:,path) */
	/* and column C (:, Path [path].ccol) */
	ASSERT (Path [path].ccol >= 0) ;
	Path [path].wfirst = *w_order ;
	Path [path].order = *w_order ;
	(*w_order)++ ;
    }
    else
    {
	/* this is a non-leaf path, within the tree */
	ASSERT (Path [path].c != EMPTY) ;
	ASSERT (Path [path].ccol == EMPTY) ;
	/* order each child path */
	for (c = Path [path].c ; c != EMPTY ; c = Path [c].next)
	{
	    dfs (Path, k, c, path_order, w_order, depth+1, npaths) ;
	    if (Path [path].wfirst == EMPTY)
	    {
		Path [path].wfirst = Path [c].wfirst ;
	    }
	}
	/* order this path next */
	Path [path].order = (*path_order)++ ;
    }
}


/* ========================================================================== */
/* === numeric update/downdate routines ===================================== */
/* ========================================================================== */

#define WDIM 1
#include "t_cholmod_updown.c"
#define WDIM 2
#include "t_cholmod_updown.c"
#define WDIM 4
#include "t_cholmod_updown.c"
#define WDIM 8
#include "t_cholmod_updown.c"


/* ========================================================================== */
/* === cholmod_updown_mark ================================================== */
/* ========================================================================== */

/* Update/downdate LDL' +/- C*C', and update/downdate selected portions of the
 * solution to Lx=b.
 *
 * The original system is L*x = b.  The new system is Lnew*xnew = b + deltab.
 * deltab(i) can be nonzero only if column i of L is modified by the update/
 * downdate.  If column i is not modified, the deltab(i) is not accessed.
 *
 * The solution to Lx=b is not modified if either X or DeltaB are NULL.
 *
 * Rowmark and colmark:
 * --------------------
 *
 * rowmark and colmark affect which portions of L take part in the update/
 * downdate of the solution to Lx=b.  They do not affect how L itself is
 * updated/downdated.  They are both ignored if X or DeltaB are NULL.
 *
 * If not NULL, rowmark is an integer array of size n where L is n-by-n.
 * rowmark [j] defines the part of column j of L that takes part in the update/
 * downdate of the forward solve, Lx=b.  Specifically, if i = rowmark [j],
 * then L(j:i-1,j) is used, and L(i:end,j) is ignored.
 *
 * If not NULL, colmark is an integer array of size C->ncol.  colmark [ccol]
 * for a column C(:,ccol) redefines those parts of L that take part in the
 * update/downdate of Lx=b.  Each column of C affects a set of columns of L.
 * If column ccol of C affects column j of L, then the new rowmark [j] of
 * column j of L is defined as colmark [ccol].  In a multiple-rank update/
 * downdate, if two or more columns of C affect column j, its new rowmark [j]
 * is the colmark of the least-numbered column of C.  colmark is ignored if
 * it is NULL, in which case rowmark is not modified.  If colmark [ccol] is
 * EMPTY (-1), then rowmark is not modified for that particular column of C.
 * colmark is ignored if it is NULL, or rowmark, X, or DeltaB are NULL.
 *
 * The algorithm for modifying the solution to Lx=b when rowmark and colmark
 * are NULL is as follows:
 *
 *	for each column j of L that is modified:
 *	    deltab (j:end) += L (j:end,j) * x(j)
 *	modify L
 *	for each column j of L that is modified:
 *	    x (j) = deltab (j)
 *	    deltab (j) = 0
 *	    deltab (j+1:end) -= L (j+1:end,j) * x(j)
 *
 * If rowmark is non-NULL but colmark is NULL:
 *
 *	for each column j of L that is modified:
 *	    deltab (j:rowmark(j)-1) += L (j:rowmark(j)-1,j) * x(j)
 *	modify L
 *	for each column j of L that is modified:
 *	    x (j) = deltab (j)
 *	    deltab (j) = 0
 *	    deltab (j+1:rowmark(j)-1) -= L (j+1:rowmark(j)-1,j) * x(j)
 *
 * If both rowmark and colmark are non-NULL:
 *
 *	for each column j of L that is modified:
 *	    deltab (j:rowmark(j)-1) += L (j:rowmark(j)-1,j) * x(j)
 *	modify L
 *	for each column j of L that is modified:
 *	    modify rowmark (j) according to colmark
 *	for each column j of L that is modified:
 *	    x (j) = deltab (j)
 *	    deltab (j) = 0
 *	    deltab (j+1:rowmark(j)-1) -= L (j+1:rowmark(j)-1,j) * x(j)
 *
 * Note that if the rank of C exceeds k = Common->maxrank (which is 2, 4, or 8),
 * then the update/downdate is done as a series of rank-k updates.  In this
 * case, the above algorithm is repeated for each block of k columns of C.
 *
 * Unless it leads to no changes in rowmark, colmark should be used only if
 * C->ncol <= Common->maxrank, because the update/downdate is done with maxrank
 * columns at a time.  Otherwise, the results are undefined.
 *
 * This routine is an "expert" routine.  It is meant for use in LPDASA only.
 */

int CHOLMOD(updown_mark)
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    Int *colmark,	/* Int array of size n. */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(updown_mask) (update, C, colmark, NULL, L, X, DeltaB,
	Common)) ;
}


/* ========================================================================== */
/* === cholmod_updown_mask ================================================== */
/* ========================================================================== */

int CHOLMOD(updown_mask)
(
    /* ---- input ---- */
    int update,		/* TRUE for update, FALSE for downdate */
    cholmod_sparse *C,	/* the incoming sparse update */
    Int *colmark,	/* Int array of size n.  See cholmod_updown.c */
    Int *mask,		/* size n */
    /* ---- in/out --- */
    cholmod_factor *L,	/* factor to modify */
    cholmod_dense *X,	/* solution to Lx=b (size n-by-1) */
    cholmod_dense *DeltaB,  /* change in b, zero on output */
    /* --------------- */
    cholmod_common *Common
)
{
    double xj, fl ;
    double *Lx, *W, *Xx, *Nx ;
    Int *Li, *Lp, *Lnz, *Cp, *Ci, *Cnz, *Head, *Flag, *Stack, *Lnext, *Iwork,
	*Set_ps1 [32], *Set_ps2 [32], *ps1, *ps2 ;
    size_t maxrank ;
    Path_type OrderedPath [32], Path [32] ;
    Int n, wdim, k1, k2, npaths, i, j, row, packed, ccol, p, cncol, do_solve,
	mark, jj, j2, kk, nextj, p1, p2, c, use_colmark, newlnz,
	k, newpath, path_order, w_order, scattered, path, newparent, pp1, pp2,
	smax, maxrow, row1, nsets, s, p3, newlnz1, Set [32], top, len, lnz, m,
	botrow ;
    size_t w ;
    int ok = TRUE ;
    DEBUG (Int oldparent) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (C, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_REAL, FALSE) ;
    RETURN_IF_XTYPE_INVALID (C, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
    n = L->n ;
    cncol = C->ncol ;
    if (!(C->sorted))
    {
	ERROR (CHOLMOD_INVALID, "C must have sorted columns") ;
	return (FALSE) ;
    }
    if (n != (Int) (C->nrow))
    {
	ERROR (CHOLMOD_INVALID, "C and L dimensions do not match") ;
	return (FALSE) ;
    }
    do_solve = (X != NULL) && (DeltaB != NULL) ;
    if (do_solve)
    {
	RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
	RETURN_IF_XTYPE_INVALID (DeltaB, CHOLMOD_REAL, CHOLMOD_REAL, FALSE) ;
	Xx = X->x ;
	Nx = DeltaB->x ;
	if (X->nrow != L->n || X->ncol != 1 || DeltaB->nrow != L->n ||
		DeltaB->ncol != 1 || Xx == NULL || Nx == NULL)
	{
	    ERROR (CHOLMOD_INVALID, "X and/or DeltaB invalid") ;
	    return (FALSE) ;
	}
    }
    else
    {
	Xx = NULL ;
	Nx = NULL ;
    }
    Common->status = CHOLMOD_OK ;
    Common->modfl = 0 ;

    fl = 0 ;
    use_colmark = (colmark != NULL) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* Note: cholmod_rowadd and cholmod_rowdel use the second n doubles in
     * Common->Xwork for Cx, and then perform a rank-1 update here, which uses
     * the first n doubles in Common->Xwork.   Both the rowadd and rowdel
     * routines allocate enough workspace so that Common->Xwork isn't destroyed
     * below.  Also, both cholmod_rowadd and cholmod_rowdel use the second n
     * ints in Common->Iwork for Ci.
     */

    /* make sure maxrank is in the proper range */
    maxrank = CHOLMOD(maxrank) (n, Common) ;
    k = MIN (cncol, (Int) maxrank) ;	/* maximum k is wdim */
    wdim = Power2 [k] ;		/* number of columns needed in W */
    ASSERT (wdim <= (Int) maxrank) ;
    PRINT1 (("updown wdim final "ID" k "ID"\n", wdim, k)) ;

    /* w = wdim * n */
    w = CHOLMOD(mult_size_t) (n, wdim, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (n, n, w, Common) ;
    if (Common->status < CHOLMOD_OK || maxrank == 0)
    {
	/* out of memory, L is returned unchanged */
	return (FALSE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* convert to simplicial numeric LDL' factor, if not already */
    /* ---------------------------------------------------------------------- */

    if (L->xtype == CHOLMOD_PATTERN || L->is_super || L->is_ll) 
    {
	/* can only update/downdate a simplicial LDL' factorization */
	CHOLMOD(change_factor) (CHOLMOD_REAL, FALSE, FALSE, FALSE, FALSE, L,
		Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory, L is returned unchanged */
	    return (FALSE) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* mark = CHOLMOD(clear_flag) (Common) ; */
    CHOLMOD_CLEAR_FLAG (Common) ;
    mark = Common->mark ;

    PRINT1 (("updown, rank %g update %d\n", (double) C->ncol, update)) ;
    DEBUG (CHOLMOD(dump_factor) (L, "input L for updown", Common)) ;
    ASSERT (CHOLMOD(dump_sparse) (C, "input C for updown", Common) >= 0) ;

    Ci = C->i ;
    Cp = C->p ;
    Cnz = C->nz ;
    packed = C->packed ;
    ASSERT (IMPLIES (!packed, Cnz != NULL)) ;

    /* ---------------------------------------------------------------------- */
    /* quick return */
    /* ---------------------------------------------------------------------- */

    if (cncol <= 0 || n == 0)
    {
	/* nothing to do */
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get L */
    /* ---------------------------------------------------------------------- */

    Li = L->i ;
    Lx = L->x ;
    Lp = L->p ;
    Lnz = L->nz ;
    Lnext = L->next ;
    ASSERT (Lnz != NULL) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	/* size n, Flag [i] <= mark must hold */
    Head = Common->Head ;	/* size n, Head [i] == EMPTY must hold */
    W = Common->Xwork ;		/* size n-by-wdim, zero on input and output*/

    /* note that Iwork [n .. 2*n-1] (i/i/l) may be in use in rowadd/rowdel: */
    Iwork = Common->Iwork ;
    Stack = Iwork ;		/* size n, uninitialized (i/i/l) */

    /* ---------------------------------------------------------------------- */
    /* entire rank-cncol update, done as a sequence of rank-k updates */
    /* ---------------------------------------------------------------------- */

    ps1 = NULL ;
    ps2 = NULL ;

    for (k1 = 0 ; k1 < cncol ; k1 += k)
    {

	/* ------------------------------------------------------------------ */
	/* get the next k columns of C for the update/downdate */
	/* ------------------------------------------------------------------ */

	/* the last update/downdate might be less than rank-k */
	if (k > cncol - k1)
	{
	    k = cncol - k1 ;
	    wdim = Power2 [k] ;
	}
	k2 = k1 + k - 1 ;

	/* workspaces are in the following state, on input and output */
	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, wdim, Common)) ;

	/* ------------------------------------------------------------------ */
	/* create a zero-length path for each column of W */
	/* ------------------------------------------------------------------ */

	nextj = n ;
	path = 0 ;
	for (ccol = k1 ; ccol <= k2 ; ccol++)
	{
	    PRINT1 (("Column ["ID"]: "ID"\n", path, ccol)) ;
	    ASSERT (ccol >= 0 && ccol <= cncol) ;
	    pp1 = Cp [ccol] ;
	    pp2 = (packed) ? (Cp [ccol+1]) : (pp1 + Cnz [ccol]) ;
	    /* get the row index j of the first entry in C (:,ccol) */
	    if (pp2 > pp1)
	    {
		/* Column ccol of C has at least one entry. */
		j = Ci [pp1] ;
	    }
	    else
	    {
		/* Column ccol of C is empty.  Pretend it has one entry in
		 * the last column with numerical value of zero. */
		j = n-1 ;
	    }
	    ASSERT (j >= 0 && j < n) ;

	    /* find first column to work on */
	    nextj = MIN (nextj, j) ;

	    Path [path].ccol = ccol ;	/* which column of C this path is for */
	    Path [path].start = EMPTY ;	/* paths for C have zero length */
	    Path [path].end = EMPTY ;
	    Path [path].parent = EMPTY ;    /* no parent yet */
	    Path [path].rank = 1 ;	    /* one column of W */
	    Path [path].c = EMPTY ;	    /* no child of this path (case A) */
	    Path [path].next = Head [j] ;   /* this path is pending at col j */
	    Path [path].pending = j ;	    /* this path is pending at col j */
	    Head [j] = path ;		    /* this path is pending at col j */
	    PRINT1(("Path "ID" starts: start "ID" end "ID" parent "ID" c "ID""
		    "j "ID" ccol "ID"\n", path, Path [path].start,
		    Path [path].end, Path [path].parent,
		    Path [path].c, j, ccol)) ;

	    /* initialize botrow for this path */
	    Path [path].botrow = (use_colmark) ? colmark [ccol] : n ;

	    path++ ;
	}

	/* we start with paths 0 to k-1.  Next one (now unused) is npaths */
	npaths = k ;

	j = nextj ;
	ASSERT (j < n) ;
	scattered = FALSE ;

	/* ------------------------------------------------------------------ */
	/* symbolic update of columns of L */
	/* ------------------------------------------------------------------ */

	while (j < n)
	{
	    ASSERT (j >= 0 && j < n && Lnz [j] > 0) ;

	    /* the old column, Li [p1..p2-1].  D (j,j) is stored in Lx [p1] */
	    p1 = Lp [j] ;
	    newlnz = Lnz [j] ;
	    p2 = p1 + newlnz  ;

#ifndef NDEBUG
	    PRINT1 (("\n=========Column j="ID" p1 "ID" p2 "ID" lnz "ID" \n",
			j, p1, p2, newlnz)) ;
	    dump_col ("Old", j, p1, p2, Li, Lx, n, Common) ;
	    oldparent = (Lnz [j] > 1) ? (Li [p1 + 1]) : EMPTY ;
	    ASSERT (CHOLMOD(dump_work) (TRUE, FALSE, 0, Common)) ;
	    ASSERT (!scattered) ;
	    PRINT1 (("Col "ID": Checking paths, npaths: "ID"\n", j, npaths)) ;
	    for (kk = 0 ; kk < npaths ; kk++)
	    {
		Int kk2, found, j3 = Path [kk].pending ;
		PRINT2 (("Path "ID" pending at "ID".\n", kk, j3)) ;
		if (j3 != EMPTY)
		{
		    /* Path kk must be somewhere in link list for column j3 */
		    ASSERT (Head [j3] != EMPTY) ;
		    PRINT3 (("    List at "ID": ", j3)) ;
		    found = FALSE ;
		    for (kk2 = Head [j3] ; kk2 != EMPTY ; kk2 = Path [kk2].next)
		    {
			PRINT3 ((""ID" ", kk2)) ;
			ASSERT (Path [kk2].pending == j3) ;
			found = found || (kk2 == kk) ;
		    }
		    PRINT3 (("\n")) ;
		    ASSERT (found) ;
		}
	    }
	    PRINT1 (("\nCol "ID": Paths at this column, head "ID"\n",
			j, Head [j]));
	    ASSERT (Head [j] != EMPTY) ;
	    for (kk = Head [j] ; kk != EMPTY ; kk = Path [kk].next)
	    {
		PRINT1 (("path "ID": (c="ID" j="ID") npaths "ID"\n",
			    kk, Path[kk].c, j, npaths)) ;
		ASSERT (kk >= 0 && kk < npaths) ;
		ASSERT (Path [kk].pending == j) ;
	    }
#endif

	    /* -------------------------------------------------------------- */
	    /* determine the path we're on */
	    /* -------------------------------------------------------------- */

	    /* get the first old path at column j */
	    path = Head [j] ;

	    /* -------------------------------------------------------------- */
	    /* update/downdate of forward solve, Lx=b */
	    /* -------------------------------------------------------------- */

	    if (do_solve)
	    {
		xj = Xx [j] ;
		if (IS_NONZERO (xj))
		{
		    xj = Xx [j] ;
		    /* This is first time column j has been seen for entire */
		    /* rank-k update/downdate. */

		    /* DeltaB += Lold (j:botrow-1,j) * X (j) */
		    Nx [j] += xj ;			/* diagonal of L */

		    /* find the botrow for this column */
		    botrow = (use_colmark) ? Path [path].botrow : n ;

		    for (p = p1 + 1 ; p < p2 ; p++)
		    {
			i = Li [p] ;
			if (i >= botrow)
			{
			    break ;
			}
			Nx [i] += Lx [p] * xj ;
		    }

		    /* clear X[j] to flag col j of Lold as having been seen.  If
		     * X (j) was initially zero, then the above code is never
		     * executed for column j.  This is safe, since if xj=0 the
		     * code above does not do anything anyway.  */
		    Xx [j] = 0.0 ;
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* start a new path at this column if two or more paths merge */
	    /* -------------------------------------------------------------- */

	    newpath =
		/* start a new path if paths have merged */
		(Path [path].next != EMPTY)
		/* or if j is the first node on a path (case A). */
		|| (Path [path].c == EMPTY) ;

	    if (newpath)
	    {
		/* get the botrow of the first path at column j */
		botrow = (use_colmark) ? Path [path].botrow : n ;

		path = npaths++ ;
		ASSERT (npaths <= 3*k) ;
		Path [path].ccol = EMPTY ; /* no single col of C for this path*/
		Path [path].start = j ;	   /* path starts at this column j */
		Path [path].end = EMPTY ;  /* don't know yet where it ends */
		Path [path].parent = EMPTY ;/* don't know parent path yet */
		Path [path].rank = 0 ;	/* rank is sum of child path ranks */
		PRINT1 (("Path "ID" starts: start "ID" end "ID" parent "ID"\n",
		path, Path [path].start, Path [path].end, Path [path].parent)) ;

		/* set the botrow of the new path */
		Path [path].botrow = (use_colmark) ? botrow : n ;
	    }

	    /* -------------------------------------------------------------- */
	    /* for each path kk pending at column j */
	    /* -------------------------------------------------------------- */

	    /* make a list of the sets that need to be merged into column j */
	    nsets = 0 ;

	    for (kk = Head [j] ; kk != EMPTY ; kk = Path [kk].next)
	    {

		/* ---------------------------------------------------------- */
		/* path kk is at (c,j) */
		/* ---------------------------------------------------------- */

		c = Path [kk].c ;
		ASSERT (c < j) ;
		PRINT1 (("TUPLE on path "ID" (c="ID" j="ID")\n", kk, c, j)) ;
		ASSERT (Path [kk].pending == j) ;

		if (newpath)
		{
		    /* finalize path kk and find rank of this path */
		    Path [kk].end = c ;	/* end of old path is previous node c */
		    Path [kk].parent = path ;	/* parent is this path */
		    Path [path].rank += Path [kk].rank ;    /* sum up ranks */
		    Path [kk].pending = EMPTY ;
		    PRINT1 (("Path "ID" done:start "ID" end "ID" parent "ID"\n",
		    kk, Path [kk].start, Path [kk].end, Path [kk].parent)) ;
		}

		if (c == EMPTY)
		{

		    /* ------------------------------------------------------ */
		    /* CASE A: first node in path */
		    /* ------------------------------------------------------ */

		    /* update:  add pattern of incoming column */

		    /* Column ccol of C is in Ci [pp1 ... pp2-1] */
		    ccol = Path [kk].ccol ;
		    pp1 = Cp [ccol] ;
		    pp2 = (packed) ? (Cp [ccol+1]) : (pp1 + Cnz [ccol]) ;
		    PRINT1 (("Case A, ccol = "ID" len "ID"\n", ccol, pp2-pp1)) ;
		    ASSERT (IMPLIES (pp2 > pp1, Ci [pp1] == j)) ;

		    if (!scattered)
		    {
			/* scatter the original pattern of column j of L */
			for (p = p1 ; p < p2 ; p++)
			{
			    Flag [Li [p]] = mark ;
			}
			scattered = TRUE ;
		    }

		    /* scatter column ccol of C (skip first entry, j) */
		    newlnz1 = newlnz ;
		    for (p = pp1 + 1 ; p < pp2 ; p++)
		    {
			row = Ci [p] ;
			if (Flag [row] < mark)
			{
			    /* this is a new entry in Lj' */
			    Flag [row] = mark ;
			    newlnz++ ;
			}
		    }
		    if (newlnz1 != newlnz)
		    {
			/* column ccol of C adds something to column j of L */
			Set [nsets++] = FLIP (ccol) ;
		    }

		}
		else if (Head [c] == 1)
		{

		    /* ------------------------------------------------------ */
		    /* CASE B: c is old, but changed, child of j */
		    /* CASE C: new child of j */
		    /* ------------------------------------------------------ */

		    /* Head [c] is 1 if col c of L has new entries,
		     * EMPTY otherwise */
		    Flag [c] = 0 ;
		    Head [c] = EMPTY ;

		    /* update: add Lc' */

		    /* column c of L is in Li [pp1 .. pp2-1] */
		    pp1 = Lp [c] ;
		    pp2 = pp1 + Lnz [c] ;
		    PRINT1 (("Case B/C: c = "ID"\n", c)) ;
		    DEBUG (dump_col ("Child", c, pp1, pp2, Li, Lx, n, Common)) ;
		    ASSERT (j == Li [pp1 + 1]) ; /* j is new parent of c */

		    if (!scattered)
		    {
			/* scatter the original pattern of column j of L */
			for (p = p1 ; p < p2 ; p++)
			{
			    Flag [Li [p]] = mark ;
			}
			scattered = TRUE ;
		    }

		    /* scatter column c of L (skip first two entries, c and j)*/
		    newlnz1 = newlnz ;
		    for (p = pp1 + 2 ; p < pp2 ; p++)
		    {
			row = Li [p] ;
			if (Flag [row] < mark)
			{
			    /* this is a new entry in Lj' */
			    Flag [row] = mark ;
			    newlnz++ ;
			}
		    }
		    PRINT2 (("\n")) ;

		    if (newlnz1 != newlnz)
		    {
			/* column c of L adds something to column j of L */
			Set [nsets++] = c ;
		    }
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* update the pattern of column j of L */
	    /* -------------------------------------------------------------- */

	    /* Column j of L will be in Li/Lx [p1 .. p3-1] */
	    p3 = p1 + newlnz ;
	    ASSERT (IMPLIES (nsets == 0, newlnz == Lnz [j])) ;
	    PRINT1 (("p1 "ID" p2 "ID" p3 "ID" nsets "ID"\n", p1, p2, p3,nsets));

	    /* -------------------------------------------------------------- */
	    /* ensure we have enough space for the longer column */
	    /* -------------------------------------------------------------- */

	    if (nsets > 0 && p3 > Lp [Lnext [j]])
	    {
		PRINT1 (("Col realloc: j "ID" newlnz "ID"\n", j, newlnz)) ;
		if (!CHOLMOD(reallocate_column) (j, newlnz, L, Common))
		{
		    /* out of memory, L is now simplicial symbolic */
		    CHOLMOD(clear_flag) (Common) ;
		    for (j = 0 ; j <= n ; j++)
		    {
			Head [j] = EMPTY ;
		    }
		    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, wdim, Common)) ;
		    return (FALSE) ;
		}
		/* L->i and L->x may have moved.  Column j has moved too */
		Li = L->i ;
		Lx = L->x ;
		p1 = Lp [j] ;
		p2 = p1 + Lnz [j] ;
		p3 = p1 + newlnz ;
	    }

	    /* -------------------------------------------------------------- */
	    /* create set pointers */
	    /* -------------------------------------------------------------- */

	    for (s = 0 ; s < nsets ; s++)
	    {
		/* Pattern of Set s is *(Set_ps1 [s] ... Set_ps2 [s]-1) */
		c = Set [s] ;
		if (c < EMPTY)
		{
		    /* column ccol of C, skip first entry (j) */
		    ccol = FLIP (c) ;
		    pp1 = Cp [ccol] ;
		    pp2 = (packed) ? (Cp [ccol+1]) : (pp1 + Cnz [ccol]) ;
		    ASSERT (pp2 - pp1 > 1) ;
		    Set_ps1 [s] = &(Ci [pp1 + 1]) ;
		    Set_ps2 [s] = &(Ci [pp2]) ;
		    PRINT1 (("set "ID" is ccol "ID"\n", s, ccol)) ;
		}
		else
		{
		    /* column c of L, skip first two entries (c and j)  */
		    pp1 = Lp [c] ;
		    pp2 = pp1 + Lnz [c]  ;
		    ASSERT (Lnz [c] > 2) ;
		    Set_ps1 [s] = &(Li [pp1 + 2]) ;
		    Set_ps2 [s] = &(Li [pp2]) ;
		    PRINT1 (("set "ID" is L "ID"\n", s, c)) ;
		}
		DEBUG (dump_set (s, Set_ps1, Set_ps2, j, n, Common)) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* multiset merge */
	    /* -------------------------------------------------------------- */

	    /* Merge the sets into a single sorted set, Lj'.  Before the merge
	     * starts, column j is located in Li/Lx [p1 ... p2-1] and the
	     * space Li/Lx [p2 ... p3-1] is empty.  p1 is Lp [j], p2 is
	     * Lp [j] + Lnz [j] (the old length of the column), and p3 is
	     * Lp [j] + newlnz (the new and longer length of the column).
	     *
	     * The sets 0 to nsets-1 are defined by the Set_ps1 and Set_ps2
	     * pointers.  Set s is located in *(Set_ps1 [s] ... Set_ps2 [s]-1).
	     * It may be a column of C, or a column of L.  All row indices i in
	     * the sets are in the range i > j and i < n.  All sets are sorted.
	     *
	     * The merge into column j of L is done in place.
	     *
	     * During the merge, p2 and p3 are updated.  Li/Lx [p1..p2-1]
	     * reflects the indices of the old column j of L that are yet to
	     * be merged into the new column.  Entries in their proper place in
	     * the new column j of L are located in Li/Lx [p3 ... p1+newlnz-1].
	     * The merge finishes when p2 == p3.
	     *
	     * During the merge, set s consumed as it is merged into column j of
	     * L.  Its unconsumed contents are *(Set_ps1 [s] ... Set_ps2 [s]-1).
	     * When a set is completely consumed, it is removed from the set of
	     * sets, and nsets is decremented.
	     *
	     * The multiset merge and 2-set merge finishes when p2 == p3.
	     */

	    PRINT1 (("Multiset merge p3 "ID" p2 "ID" nsets "ID"\n",
			p3, p2, nsets)) ;

	    while (p3 > p2 && nsets > 1)
	    {

#ifndef NDEBUG
		PRINT2 (("\nMultiset merge.  nsets = "ID"\n", nsets)) ;
		PRINT2 (("Source col p1 = "ID", p2 = "ID", p3= "ID"\n",
			    p1, p2, p3)) ;
		for (p = p1 + 1 ; p < p2 ; p++)
		{
		    PRINT2 (("    p: "ID" source row "ID" %g\n",
				p, Li[p], Lx[p])) ;
		    ASSERT (Li [p] > j && Li [p] < n) ;
		}
		PRINT2 (("---\n")) ;
		for (p = p3 ; p < p1 + newlnz ; p++)
		{
		    PRINT2 (("    p: "ID" target row "ID" %g\n",
				p, Li[p], Lx[p])) ;
		    ASSERT (Li [p] > j && Li [p] <  n) ;
		}
		for (s = 0 ; s < nsets ; s++)
		{
		    dump_set (s, Set_ps1, Set_ps2, j, n, Common) ;
		}
#endif

		/* get the entry at the tail end of source column Lj */
		row1 = Li [p2 - 1] ;
		ASSERT (row1 >= j && p2 >= p1) ;

		/* find the largest row in all the sets */
		maxrow = row1 ;
		smax = EMPTY ;
		for (s = nsets-1 ; s >= 0 ; s--)
		{
		    ASSERT (Set_ps1 [s] < Set_ps2 [s]) ;
		    row = *(Set_ps2 [s] - 1) ;
		    if (row == maxrow)
		    {
			/* skip past this entry in set s (it is a duplicate) */
			Set_ps2 [s]-- ;
			if (Set_ps1 [s] == Set_ps2 [s])
			{
			    /* nothing more in this set */
			    nsets-- ;
			    Set_ps1 [s] = Set_ps1 [nsets] ;
			    Set_ps2 [s] = Set_ps2 [nsets] ;
			    if (smax == nsets)
			    {
				/* Set smax redefined; it is now this set */
				smax = s ;
			    }
			}
		    }
		    else if (row > maxrow)
		    {
			maxrow = row ;
			smax = s ;
		    }
		}
		ASSERT (maxrow > j) ;

		/* move the row onto the stack of the target column */
		if (maxrow == row1)
		{
		    /* next entry is in Lj, move to the bottom of Lj' */
		    ASSERT (smax == EMPTY) ;
		    p2-- ;
		    p3-- ;
		    Li [p3] = maxrow ;
		    Lx [p3] = Lx [p2] ;
		}
		else
		{
		    /* new entry in Lj' */
		    ASSERT (smax >= 0 && smax < nsets) ;
		    Set_ps2 [smax]-- ;
		    p3-- ;
		    Li [p3] = maxrow ;
		    Lx [p3] = 0.0 ;
		    if (Set_ps1 [smax] == Set_ps2 [smax])
		    {
			/* nothing more in this set */
			nsets-- ;
			Set_ps1 [smax] = Set_ps1 [nsets] ;
			Set_ps2 [smax] = Set_ps2 [nsets] ;
			PRINT1 (("Set "ID" now empty\n", smax)) ;
		    }
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* 2-set merge: */
	    /* -------------------------------------------------------------- */

	    /* This the same as the multi-set merge, except there is only one
	     * set s = 0 left.  The source column j and the set 0 are being
	     * merged into the target column j. */

	    if (nsets > 0)
	    {
		ps1 = Set_ps1 [0] ;
		ps2 = Set_ps2 [0] ;
	    }

	    while (p3 > p2)
	    {

#ifndef NDEBUG
		PRINT2 (("\n2-set merge.\n")) ;
		ASSERT (nsets == 1) ;
		PRINT2 (("Source col p1 = "ID", p2 = "ID", p3= "ID"\n",
			    p1, p2, p3)) ;
		for (p = p1 + 1 ; p < p2 ; p++)
		{
		    PRINT2 (("    p: "ID" source row "ID" %g\n",
				p, Li[p], Lx[p])) ;
		    ASSERT (Li [p] > j && Li [p] < n) ;
		}
		PRINT2 (("---\n")) ;
		for (p = p3 ; p < p1 + newlnz ; p++)
		{
		    PRINT2 (("    p: "ID" target row "ID" %g\n",
				p, Li[p], Lx[p])) ;
		    ASSERT (Li [p] > j && Li [p] <  n) ;
		}
		dump_set (0, Set_ps1, Set_ps2, j, n, Common) ;
#endif

		if (p2 == p1 + 1)
		{
		    /* the top of Lj is empty; copy the set and quit */
		    while (p3 > p2)
		    {
			/* new entry in Lj' */
			row = *(--ps2) ;
			p3-- ;
			Li [p3] = row ;
			Lx [p3] = 0.0 ;
		    }
		}
		else
		{
		    /* get the entry at the tail end of Lj */
		    row1 = Li [p2 - 1] ;
		    ASSERT (row1 > j && row1 < n) ;
		    /* get the entry at the tail end of the incoming set */
		    ASSERT (ps1 < ps2) ;
		    row = *(ps2-1) ;
		    ASSERT (row > j && row1 < n) ;
		    /* move the larger of the two entries to the target set */
		    if (row1 >= row)
		    {
			/* next entry is in Lj, move to the bottom */
			if (row1 == row)
			{
			    /* skip past this entry in the set */
			    ps2-- ;
			}
			p2-- ;
			p3-- ;
			Li [p3] = row1 ;
			Lx [p3] = Lx [p2] ;
		    }
		    else
		    {
			/* new entry in Lj' */
			ps2-- ;
			p3-- ;
			Li [p3] = row ;
			Lx [p3] = 0.0 ;
		    }
		}
	    }

	    /* -------------------------------------------------------------- */
	    /* The new column j of L is now in Li/Lx [p1 ... p2-1] */
	    /* -------------------------------------------------------------- */

	    p2 = p1 + newlnz ;
	    DEBUG (dump_col ("After merge: ", j, p1, p2, Li, Lx, n, Common)) ;

	    fl += Path [path].rank * (6 + 4 * (double) newlnz) ;

	    /* -------------------------------------------------------------- */
	    /* clear Flag; original pattern of column j L no longer marked */
	    /* -------------------------------------------------------------- */

	    mark = CHOLMOD(clear_flag) (Common) ;
	    scattered = FALSE ;

	    /* -------------------------------------------------------------- */
	    /* find the new parent */
	    /* -------------------------------------------------------------- */

	    newparent = (newlnz > 1) ? (Li [p1 + 1]) : EMPTY ;
	    PRINT1 (("\nNew parent, Lnz: "ID": "ID" "ID"\n",
			j, newparent,newlnz));
	    ASSERT (oldparent == EMPTY || newparent <= oldparent) ;

	    /* -------------------------------------------------------------- */
	    /* go to the next node in the path */
	    /* -------------------------------------------------------------- */

	    /* path moves to (j,nextj) unless j is a root */
	    nextj = (newparent == EMPTY) ? n : newparent ;

	    /* place path at head of list for nextj, or terminate the path */
	    PRINT1 (("\n j = "ID" nextj = "ID"\n\n", j, nextj)) ;
	    Path [path].c = j ;
	    if (nextj < n)
	    {
		/* put path on link list of pending paths at column nextj */
		Path [path].next = Head [nextj] ;
		Path [path].pending = nextj ;
		Head [nextj] = path ;
		PRINT1 (("Path "ID" continues to ("ID","ID").  Rank "ID"\n",
		    path, Path [path].c, nextj, Path [path].rank)) ;
	    }
	    else
	    {
		/* path has ended here, at a root */
		Path [path].next = EMPTY ;
		Path [path].pending = EMPTY ;
		Path [path].end = j ;
		PRINT1 (("Path "ID" ends at root ("ID").  Rank "ID"\n",
		    path, Path [path].end, Path [path].rank)) ;
	    }

	    /* The link list Head [j] can now be emptied.  Set Head [j] to 1
	     * if column j has changed (it is no longer used as a link list). */
	    PRINT1 (("column "ID", oldlnz = "ID"\n", j, Lnz [j])) ;
	    Head [j] = (Lnz [j] != newlnz) ? 1 : EMPTY ;
	    Lnz [j] = newlnz ;
	    PRINT1 (("column "ID", newlnz = "ID"\n", j, newlnz)) ;
	    DEBUG (dump_col ("New", j, p1, p2, Li, Lx, n, Common)) ;

	    /* move to the next column */
	    if (k == Path [path].rank)
	    {
		/* only one path left */
		j = nextj ;
	    }
	    else
	    {
		/* The current path is moving from column j to column nextj
		 * (nextj is n if the path has ended).  However, there may be
		 * other paths pending in columns j+1 to nextj-1.  There are
		 * two methods for looking for the next column with a pending
		 * update.  The first one looks at all columns j+1 to nextj-1
		 * for a non-empty link list.  This can be costly if j and
		 * nextj differ by a large amount (it can be O(n), but this
		 * entire routine may take Omega(1) time).  The second method
		 * looks at all paths and finds the smallest column at which any
		 * path is pending.  It takes O(# of paths), which is bounded
		 * by 23: one for each column of C (up to 8), and then 15 for a
		 * balanced binary tree with 8 leaves.  However, if j and
		 * nextj differ by a tiny amount (nextj is often j+1 near
		 * the end of the matrix), looking at columns j+1 to nextj
		 * would be faster.  Both methods give the same answer. */

		if (nextj - j < npaths)
		{
		    /* there are fewer columns to search than paths */
		    PRINT1 (("check j="ID" to nextj="ID"\n", j, nextj)) ;
		    for (j2 = j + 1 ; j2 < nextj ; j2++)
		    {
			PRINT1 (("check j="ID" "ID"\n", j2, Head [j2])) ;
			if (Head [j2] != EMPTY)
			{
			    PRINT1 (("found, j="ID"\n", j2)) ;
			    ASSERT (Path [Head [j2]].pending == j2) ;
			    break ;
			}
		    }
		}
		else
		{
		    /* there are fewer paths than columns to search */
		    j2 = nextj ;
		    for (kk = 0 ; kk < npaths ; kk++)
		    {
			jj = Path [kk].pending ;
			PRINT2 (("Path "ID" pending at "ID"\n", kk, jj)) ;
			if (jj != EMPTY) j2 = MIN (j2, jj) ;
		    }
		}
		j = j2 ;
	    }
	}

	/* ensure workspaces are back to the values required on input */
	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, TRUE, Common)) ;

	/* ------------------------------------------------------------------ */
	/* depth-first-search of tree to order the paths */
	/* ------------------------------------------------------------------ */

	/* create lists of child paths */
	PRINT1 (("\n\nDFS search:\n\n")) ;
	for (path = 0 ; path < npaths ; path++)
	{
	    Path [path].c = EMPTY ;	    /* first child of path */
	    Path [path].next = EMPTY ;	    /* next sibling of path */
	    Path [path].order = EMPTY ;	    /* path is not ordered yet */
	    Path [path].wfirst = EMPTY ;    /* 1st column of W not found yet */

#ifndef NDEBUG
	    j = Path [path].start ;
	    PRINT1 (("Path "ID" : start "ID" end "ID" parent "ID" ccol "ID"\n", 
	    path, j, Path [path].end, Path [path].parent, Path [path].ccol)) ;
	    for ( ; ; )
	    {
		PRINT1 (("	column "ID"\n", j)) ;
		ASSERT (j == EMPTY || (j >= 0 && j < n)) ;
		if (j == Path [path].end)
		{
		    break ;
		}
		ASSERT (j >= 0 && j < n) ;
		j = (Lnz [j] > 1) ? (Li [Lp [j] + 1]) : EMPTY ;
	    }
#endif
	}

	for (path = 0 ; path < npaths ; path++)
	{
	    p = Path [path].parent ;	/* add path to child list of parent */
	    if (p != EMPTY)
	    {
		ASSERT (p < npaths) ;
		Path [path].next = Path [p].c ;
		Path [p].c = path ;
	    }
	}

	path_order = k ;
	w_order = 0 ;
	for (path = npaths-1 ; path >= 0 ; path--)
	{
	    if (Path [path].order == EMPTY)
	    {
		/* this path is the root of a subtree of Tbar */
		PRINT1 (("Root path "ID"\n", path)) ;
		ASSERT (path >= k) ;
		dfs (Path, k, path, &path_order, &w_order, 0, npaths) ;
	    }
	}
	ASSERT (path_order == npaths) ;
	ASSERT (w_order == k) ;

	/* reorder the paths */
	for (path = 0 ; path < npaths ; path++)
	{
	    /* old order is path, new order is Path [path].order */
	    OrderedPath [Path [path].order] = Path [path] ;
	}

#ifndef NDEBUG
	for (path = 0 ; path < npaths ; path++)
	{
	    PRINT1 (("Ordered Path "ID": start "ID" end "ID" wfirst "ID" rank "
		    ""ID" ccol "ID"\n", path, OrderedPath [path].start,
		    OrderedPath [path].end, OrderedPath [path].wfirst,
		    OrderedPath [path].rank, OrderedPath [path].ccol)) ;
	    if (path < k)
	    {
		ASSERT (OrderedPath [path].ccol >= 0) ;
	    }
	    else
	    {
		ASSERT (OrderedPath [path].ccol == EMPTY) ;
	    }
	}
#endif

	/* ------------------------------------------------------------------ */
	/* numeric update/downdate for all paths */
	/* ------------------------------------------------------------------ */

	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, wdim, Common)) ;

	switch (wdim)
	{
	    case 1:
		updown_1_r (update, C, k, L, W, OrderedPath, npaths, mask,
		    Common) ;
		break ;
	    case 2:
		updown_2_r (update, C, k, L, W, OrderedPath, npaths, mask,
		    Common) ;
		break ;
	    case 4:
		updown_4_r (update, C, k, L, W, OrderedPath, npaths, mask,
		    Common) ;
		break ;
	    case 8:
		updown_8_r (update, C, k, L, W, OrderedPath, npaths, mask,
		    Common) ;
		break ;
	}

	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, wdim, Common)) ;
    }

    /* ---------------------------------------------------------------------- */
    /* update/downdate the forward solve */
    /* ---------------------------------------------------------------------- */

    if (do_solve)
    {
	/* We now have DeltaB += Lold (:,j) * X (j) for all columns j in union
	 * of all paths seen during the entire rank-cncol update/downdate. For
	 * each j in path, do DeltaB -= Lnew (:,j)*DeltaB(j) 
	 * in topological order. */

#ifndef NDEBUG
	PRINT1 (("\ndo_solve, DeltaB + Lold(:,Path)*X(Path):\n")) ;
	for (i = 0 ; i < n ; i++)
	{
	    PRINT1 (("do_solve: "ID" %30.20e\n", i, Nx [i])) ;
	}
#endif

	/* Note that the downdate, if it deleted entries, would need to compute
	 * the Stack prior to doing any downdates. */

	/* find the union of all the paths in the new L */
	top = n ;	/* "top" is stack pointer, not a row or column index */
	for (ccol = 0 ; ccol < cncol ; ccol++)
	{

	    /* -------------------------------------------------------------- */
	    /* j = first row index of C (:,ccol) */
	    /* -------------------------------------------------------------- */

	    pp1 = Cp [ccol] ;
	    pp2 = (packed) ? (Cp [ccol+1]) : (pp1 + Cnz [ccol]) ;
	    if (pp2 > pp1)
	    {
		/* Column ccol of C has at least one entry. */
		j = Ci [pp1] ;
	    }
	    else
	    {
		/* Column ccol of C is empty */
		j = n-1 ;
	    }
	    PRINT1 (("\ndo_solve:      ccol= "ID"\n", ccol)) ;
	    ASSERT (j >= 0 && j < n) ;
	    len = 0 ;

	    /* -------------------------------------------------------------- */
	    /* find the new rowmark */
	    /* -------------------------------------------------------------- */

	    /* Each column of C can redefine the region of L that takes part in
	     * the update/downdate of the triangular solve Lx=b.  If
	     * i = colmark [ccol] for column C(:,ccol), then i = rowmark [j] is
	     * redefined for all columns along the path modified by C(:,ccol).
	     * If more than one column modifies any given column j of L, then
	     * the rowmark of j is determined by the colmark of the least-
	     * numbered column that affects column j.  That is, if both
	     * C(:,ccol1) and C(:,ccol2) affect column j of L, then
	     * rowmark [j] = colmark [MIN (ccol1, ccol2)].
	     *
	     * rowmark [j] is not modified if rowmark or colmark are NULL,
	     * or if colmark [ccol] is EMPTY.
	     */

	    botrow = (use_colmark) ? (colmark [ccol]) : EMPTY ;

	    /* -------------------------------------------------------------- */
	    /* traverse from j towards root, stopping if node already visited */
	    /* -------------------------------------------------------------- */

	    while (j != EMPTY && Flag [j] < mark)
	    {
		PRINT1 (("do_solve: subpath j= "ID"\n", j)) ;
		ASSERT (j >= 0 && j < n) ;
		Stack [len++] = j ;		/* place j on the stack */
		Flag [j] = mark ;		/* flag j as visited */

		/* if using colmark, mark column j with botrow */
		ASSERT (Li [Lp [j]] == j) ;	/* diagonal is always present */
		if (use_colmark)
		{
		    Li [Lp [j]] = botrow ;	/* use the space for botrow */
		}

		/* go up the tree, to the parent of j */
		j = (Lnz [j] > 1) ? (Li [Lp [j] + 1]) : EMPTY ;
	    }

	    /* -------------------------------------------------------------- */
	    /* move the path down to the bottom of the stack */
	    /* -------------------------------------------------------------- */

	    ASSERT (len <= top) ;
	    while (len > 0)
	    {
		Stack [--top] = Stack [--len] ;
	    }
	}

#ifndef NDEBUG
	/* Union of paths now in Stack [top..n-1] in topological order */
	PRINT1 (("\nTopological order:\n")) ;
	for (i = top ; i < n ; i++)
	{
	    PRINT1 (("column "ID" in full path\n", Stack [i])) ;
	}
#endif

	/* Do the forward solve for the full path part of L */
	for (m = top ; m < n ; m++)
	{
	    j = Stack [m] ;
	    ASSERT (j >= 0 && j < n) ;
	    PRINT1 (("do_solve: path j= "ID"\n", j)) ;
	    p1 = Lp [j] ;
	    lnz = Lnz [j] ;
	    p2 = p1 + lnz ;
	    xj = Nx [j] ;

	    /* copy new solution onto old one, for all cols in full path */
	    Xx [j] = xj ;
	    Nx [j] = 0. ;

	    /* DeltaB -= Lnew (j+1:botrow-1,j) * deltab(j) */
	    if (use_colmark)
	    {
		botrow = Li [p1] ;	/* get botrow */
		Li [p1] = j ;		/* restore diagonal entry */
		for (p = p1 + 1 ; p < p2 ; p++)
		{
		    i = Li [p] ;
		    if (i >= botrow) break ;
		    Nx [i] -= Lx [p] * xj ;
		}
	    }
	    else
	    {
		for (p = p1 + 1 ; p < p2 ; p++)
		{
		    Nx [Li [p]] -= Lx [p] * xj ;
		}
	    }
	}

	/* clear the Flag */
	mark = CHOLMOD(clear_flag) (Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* successful update/downdate */
    /* ---------------------------------------------------------------------- */

    Common->modfl = fl ;
    DEBUG (for (j = 0 ; j < n ; j++) ASSERT (IMPLIES (do_solve, Nx[j] == 0.))) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, TRUE, Common)) ;
    DEBUG (CHOLMOD(dump_factor) (L, "output L for updown", Common)) ;
    return (TRUE) ;
}
#endif
