/* ========================================================================== */
/* === Cholesky/cholmod_analyze ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2013, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Order and analyze a matrix (either simplicial or supernodal), in prepartion
 * for numerical factorization via cholmod_factorize or via the "expert"
 * routines cholmod_rowfac and cholmod_super_numeric.
 *
 * symmetric case:    A or A(p,p)
 * unsymmetric case:  AA', A(p,:)*A(p,:)', A(:,f)*A(:,f)', or A(p,f)*A(p,f)'
 *
 * For the symmetric case, only the upper or lower triangular part of A is
 * accessed (depending on the type of A).  LL'=A (or permuted A) is analzed.
 * For the unsymmetric case (LL'=AA' or permuted A).
 *
 * There can be no duplicate entries in p or f.  p is of length m if A is
 * m-by-n.  f can be length 0 to n.
 *
 * In both cases, the columns of A need not be sorted.  A can be in packed
 * or unpacked form.
 *
 * Ordering options include:
 *
 *	natural:    A is not permuted to reduce fill-in
 *	given:	    a permutation can be provided to this routine (UserPerm)
 *	AMD:	    approximate minumum degree (AMD for the symmetric case,
 *		    COLAMD for the AA' case).
 *	METIS:	    nested dissection with METIS_NodeND
 *	NESDIS:	    nested dissection using METIS_NodeComputeSeparator,
 *		    typically followed by a constrained minimum degree
 *		    (CAMD for the symmetric case, CCOLAMD for the AA' case).
 *
 * Multiple ordering options can be tried (up to 9 of them), and the best one
 * is selected (the one that gives the smallest number of nonzeros in the
 * simplicial factor L).  If one method fails, cholmod_analyze keeps going, and
 * picks the best among the methods that succeeded.  This routine fails (and
 * returns NULL) if either initial memory allocation fails, all ordering methods
 * fail, or the supernodal analysis (if requested) fails.  By default, the 9
 * methods available are:
 *
 *	1) given permutation (skipped if UserPerm is NULL)
 *	2) AMD (symmetric case) or COLAMD (unsymmetric case)
 *	3) METIS with default parameters
 *	4) NESDIS with default parameters (stopping the partitioning when
 *	    the graph is of size nd_small = 200 or less, remove nodes with
 *	    more than max (16, prune_dense * sqrt (n)) nodes where
 *	    prune_dense = 10, and follow partitioning with CCOLAMD, a
 *	    constrained minimum degree ordering).
 *	5) natural
 *	6) NESDIS, nd_small = 20000, prune_dense = 10
 *	7) NESDIS, nd_small =     4, prune_dense = 10, no min degree
 *	8) NESDIS, nd_small =   200, prune_dense = 0
 *	9) COLAMD for A*A' or AMD for A
 *
 * By default, the first two are tried, and METIS is tried if AMD reports a high
 * flop count and fill-in.  Let fl denote the flop count for the AMD, ordering,
 * nnz(L) the # of nonzeros in L, and nnz(tril(A)) (or A*A').  If
 * fl/nnz(L) >= 500 and nnz(L)/nnz(tril(A)) >= 5, then METIS is attempted.  The
 * best ordering is used (UserPerm if given, AMD, and METIS if attempted).  If
 * you do not have METIS, only the first two will be tried (user permutation,
 * if provided, and AMD/COLAMD).  This default behavior is obtained when
 * Common->nmethods is zero.  In this case, methods 0, 1, and 2 in
 * Common->method [..] are reset to User-provided, AMD, and METIS (or NESDIS
 * if Common->default_nesdis is set to the non-default value of TRUE),
 * respectively.
 *
 * You can modify these 9 methods and the number of methods tried by changing
 * parameters in the Common argument.  If you know the best ordering for your
 * matrix, set Common->nmethods to 1 and set Common->method[0].ordering to the
 * requested ordering method.  Parameters for each method can also be modified
 * (refer to cholmod.h for details).
 *
 * Note that it is possible for METIS to terminate your program if it runs out
 * of memory.  This is not the case for any CHOLMOD or minimum degree ordering
 * routine (AMD, COLAMD, CAMD, CCOLAMD, or CSYMAMD).  Since NESDIS relies on
 * METIS, it too can terminate your program.
 *
 * The factor L is returned as simplicial symbolic (L->is_super FALSE) if
 * Common->supernodal <= CHOLMOD_SIMPLICIAL (0) or as supernodal symbolic if
 * Common->supernodal >= CHOLMOD_SUPERNODAL (2).  If Common->supernodal is
 * equal to CHOLMOD_AUTO (1), then L is simplicial if the flop count per
 * nonzero in L is less than Common->supernodal_switch (default: 40), and
 * is returned as a supernodal factor otherwise.
 *
 * In both cases, L->xtype is CHOLMOD_PATTERN.
 * A subsequent call to cholmod_factorize will perform a
 * simplicial or supernodal factorization, depending on the type of L.
 *
 * For the simplicial case, L contains the fill-reducing permutation (L->Perm)
 * and the counts of nonzeros in each column of L (L->ColCount).  For the
 * supernodal case, L also contains the nonzero pattern of each supernode.
 *
 * workspace: Flag (nrow), Head (nrow+1)
 *	if symmetric:   Iwork (6*nrow)
 *	if unsymmetric: Iwork (6*nrow+ncol).
 *	calls various ordering routines, which typically allocate O(nnz(A))
 *	temporary workspace ((2 to 3)*nnz(A) * sizeof (Int) is typical, but it
 *	can be much higher if A*A' must be explicitly formed for METIS).  Also
 *	allocates up to 2 temporary (permuted/transpose) copies of the nonzero
 *	pattern of A, and up to 3*n*sizeof(Int) additional workspace.
 *
 * Supports any xtype (pattern, real, complex, or zomplex)
 */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"

#ifndef NSUPERNODAL
#include "cholmod_supernodal.h"
#endif

#ifndef NPARTITION
#include "cholmod_partition.h"
#endif


/* ========================================================================== */
/* === cholmod_analyze ====================================================== */
/* ========================================================================== */

/* Orders and analyzes A, AA', PAP', or PAA'P' and returns a symbolic factor
 * that can later be passed to cholmod_factorize. */

cholmod_factor *CHOLMOD(analyze)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order and analyze */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(analyze_p2) (TRUE, A, NULL, NULL, 0, Common)) ;
}


/* ========================================================================== */
/* === cholmod_analyze_p ==================================================== */
/* ========================================================================== */

/* Orders and analyzes A, AA', PAP', PAA'P', FF', or PFF'P and returns a
 * symbolic factor that can later be passed to cholmod_factorize, where
 * F = A(:,fset) if fset is not NULL and A->stype is zero.
 * UserPerm is tried if non-NULL.  */

cholmod_factor *CHOLMOD(analyze_p)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order and analyze */
    Int *UserPerm,	/* user-provided permutation, size A->nrow */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(analyze_p2) (TRUE, A, UserPerm, fset, fsize, Common)) ;
}


/* ========================================================================== */
/* === permute_matrices ===================================================== */
/* ========================================================================== */

/* Permute and transpose a matrix.  Allocates the A1 and A2 matrices, if needed,
 * or returns them as NULL if not needed.
 */

static int permute_matrices
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to permute */
    Int ordering,	/* ordering method used */
    Int *Perm,		/* fill-reducing permutation */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    Int do_rowcolcounts,/* if TRUE, compute both S and F.  If FALSE, only
			 * S is needed for the symmetric case, and only F for
			 * the unsymmetric case */
    /* ---- output --- */
    cholmod_sparse **A1_handle,	    /* see comments below for A1, A2, S, F */
    cholmod_sparse **A2_handle,
    cholmod_sparse **S_handle,
    cholmod_sparse **F_handle,
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *A1, *A2, *S, *F ;

    *A1_handle = NULL ;
    *A2_handle = NULL ;
    *S_handle = NULL ;
    *F_handle = NULL ;
    A1 = NULL ;
    A2 = NULL ;

    if (ordering == CHOLMOD_NATURAL)
    {

	/* ------------------------------------------------------------------ */
	/* natural ordering of A */
	/* ------------------------------------------------------------------ */

	if (A->stype < 0)
	{
	    /* symmetric lower case: A already in lower form, so S=A' */
	    /* workspace: Iwork (nrow) */
	    A2 = CHOLMOD(ptranspose) (A, 0, NULL, NULL, 0, Common) ;
	    F = A ;
	    S = A2 ;
	}
	else if (A->stype > 0)
	{
	    /* symmetric upper case: F = pattern of triu (A)', S = A */
	    /* workspace: Iwork (nrow) */
	    if (do_rowcolcounts)
	    {
		/* F not needed for symmetric case if do_rowcolcounts FALSE */
		A1 = CHOLMOD(ptranspose) (A, 0, NULL, fset, fsize, Common) ;
	    }
	    F = A1 ;
	    S = A ;
	}
	else
	{
	    /* unsymmetric case: F = pattern of A (:,f)',  S = A */
	    /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
	    A1 = CHOLMOD(ptranspose) (A, 0, NULL, fset, fsize, Common) ;
	    F = A1 ;
	    S = A ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* A is permuted */
	/* ------------------------------------------------------------------ */

	if (A->stype < 0)
	{
	    /* symmetric lower case: S = tril (A (p,p))' and F = S' */
	    /* workspace: Iwork (2*nrow) */
	    A2 = CHOLMOD(ptranspose) (A, 0, Perm, NULL, 0, Common) ;
	    S = A2 ;
	    /* workspace: Iwork (nrow) */
	    if (do_rowcolcounts)
	    {
		/* F not needed for symmetric case if do_rowcolcounts FALSE */
		A1 = CHOLMOD(ptranspose) (A2, 0, NULL, NULL, 0, Common) ;
	    }
	    F = A1 ;
	}
	else if (A->stype > 0)
	{
	    /* symmetric upper case: F = triu (A (p,p))' and S = F' */
	    /* workspace: Iwork (2*nrow) */
	    A1 = CHOLMOD(ptranspose) (A, 0, Perm, NULL, 0, Common) ;
	    F = A1 ;
	    /* workspace: Iwork (nrow) */
	    A2 = CHOLMOD(ptranspose) (A1, 0, NULL, NULL, 0, Common) ;
	    S = A2 ;
	}
	else
	{
	    /* unsymmetric case:     F = A (p,f)'         and S = F' */
	    /* workspace: Iwork (nrow if no fset, MAX(nrow,ncol) if fset) */
	    A1 = CHOLMOD(ptranspose) (A, 0, Perm, fset, fsize, Common) ;
	    F = A1 ;
	    if (do_rowcolcounts)
	    {
		/* S not needed for unsymmetric case if do_rowcolcounts FALSE */
		/* workspace: Iwork (nrow) */
		A2 = CHOLMOD(ptranspose) (A1, 0, NULL, NULL, 0, Common) ;
	    }
	    S = A2 ;
	}
    }

    /* If any cholmod_*transpose fails, one or more matrices will be NULL */
    *A1_handle = A1 ;
    *A2_handle = A2 ;
    *S_handle = S ;
    *F_handle = F ;
    return (Common->status == CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_analyze_ordering ============================================= */
/* ========================================================================== */

/* Given a matrix A and its fill-reducing permutation, compute the elimination
 * tree, its (non-weighted) postordering, and the number of nonzeros in each
 * column of L.  Also computes the flop count, the total nonzeros in L, and
 * the nonzeros in A (Common->fl, Common->lnz, and Common->anz).
 *
 * The column counts of L, flop count, and other statistics from
 * cholmod_rowcolcounts are not computed if ColCount is NULL.
 *
 * workspace: Iwork (2*nrow if symmetric, 2*nrow+ncol if unsymmetric),
 *	Flag (nrow), Head (nrow+1)
 */

int CHOLMOD(analyze_ordering)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    int ordering,	/* ordering method used */
    Int *Perm,		/* size n, fill-reducing permutation to analyze */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    Int *Parent,	/* size n, elimination tree */
    Int *Post,		/* size n, postordering of elimination tree */
    Int *ColCount,	/* size n, nnz in each column of L */
    /* ---- workspace  */
    Int *First,		/* size n workspace for cholmod_postorder */
    Int *Level,		/* size n workspace for cholmod_postorder */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *A1, *A2, *S, *F ;
    Int n, ok, do_rowcolcounts ;

    /* check inputs */
    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;

    n = A->nrow ;

    do_rowcolcounts = (ColCount != NULL) ;

    /* permute A according to Perm and fset */
    ok = permute_matrices (A, ordering, Perm, fset, fsize, do_rowcolcounts,
	    &A1, &A2, &S, &F, Common) ;

    /* find etree of S (symmetric upper/lower case) or F (unsym case) */
    /* workspace: symmmetric: Iwork (nrow), unsym: Iwork (nrow+ncol) */
    ok = ok && CHOLMOD(etree) (A->stype ? S:F, Parent, Common) ;

    /* postorder the etree (required by cholmod_rowcolcounts) */
    /* workspace: Iwork (2*nrow) */
    ok = ok && (CHOLMOD(postorder) (Parent, n, NULL, Post, Common) == n) ;

    /* cholmod_postorder doesn't set Common->status if it returns < n */
    Common->status = (!ok && Common->status == CHOLMOD_OK) ?
	CHOLMOD_INVALID : Common->status ;

    /* analyze LL'=S or SS' or S(:,f)*S(:,f)' */
    /* workspace:
     *	if symmetric:   Flag (nrow), Iwork (2*nrow)
     *	if unsymmetric: Flag (nrow), Iwork (2*nrow+ncol), Head (nrow+1)
     */
    if (do_rowcolcounts)
    {
	ok = ok && CHOLMOD(rowcolcounts) (A->stype ? F:S, fset, fsize, Parent,
	    Post, NULL, ColCount, First, Level, Common) ;
    }

    /* free temporary matrices and return result */
    CHOLMOD(free_sparse) (&A1, Common) ;
    CHOLMOD(free_sparse) (&A2, Common) ;
    return (ok) ;
}


/* ========================================================================== */
/* === Free workspace and return L ========================================== */
/* ========================================================================== */

#define FREE_WORKSPACE_AND_RETURN \
{ \
    Common->no_workspace_reallocate = FALSE ; \
    CHOLMOD(free) (n, sizeof (Int), Lparent,  Common) ; \
    CHOLMOD(free) (n, sizeof (Int), Perm,     Common) ; \
    CHOLMOD(free) (n, sizeof (Int), ColCount, Common) ; \
    if (Common->status < CHOLMOD_OK) \
    { \
	CHOLMOD(free_factor) (&L, Common) ; \
    } \
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ; \
    return (L) ; \
}


/* ========================================================================== */
/* === cholmod_analyze_p2 =================================================== */
/* ========================================================================== */

/* Ordering and analysis for sparse Cholesky or sparse QR.  CHOLMOD itself
 * always uses for_cholesky = TRUE.  The for_cholesky = FALSE option is
 * for SuiteSparseQR only. */

cholmod_factor *CHOLMOD(analyze_p2)
(
    /* ---- input ---- */
    int for_cholesky,   /* if TRUE, then analyze for Cholesky; else for QR */
    cholmod_sparse *A,	/* matrix to order and analyze */
    Int *UserPerm,	/* user-provided permutation, size A->nrow */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* --------------- */
    cholmod_common *Common
)
{
    double lnz_best ;
    Int *First, *Level, *Work4n, *Cmember, *CParent, *ColCount, *Lperm, *Parent,
	*Post, *Perm, *Lparent, *Lcolcount ;
    cholmod_factor *L ;
    Int k, n, ordering, method, nmethods, status, default_strategy, ncol, uncol,
	skip_analysis, skip_best ;
    Int amd_backup ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    Common->status = CHOLMOD_OK ;
    status = CHOLMOD_OK ;
    Common->selected = EMPTY ;
    Common->called_nd = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    ncol = A->ncol ;
    uncol = (A->stype == 0) ? (A->ncol) : 0 ;

    /* ---------------------------------------------------------------------- */
    /* set the default strategy */
    /* ---------------------------------------------------------------------- */

    lnz_best = (double) EMPTY ;
    skip_best = FALSE ;
    nmethods = MIN (Common->nmethods, CHOLMOD_MAXMETHODS) ;
    nmethods = MAX (0, nmethods) ;

#ifndef NDEBUG
    PRINT1 (("cholmod_analyze_p2 :: nmethods "ID"\n", nmethods)) ;
    for (method = 0 ; method < nmethods ; method++)
    {
        PRINT1 (("  "ID": ordering "ID"\n",     
            method, Common->method [method].ordering)) ;
    }
#endif

    default_strategy = (nmethods == 0) ;
    if (default_strategy)
    {
	/* default strategy: try UserPerm, if given.  Try AMD for A, or AMD
	 * to order A*A'.  Try METIS for the symmetric case only if AMD reports
         * a high degree of fill-in and flop count.  METIS is not tried if the
         * Partition Module isn't installed.   If Common->default_nesdis is
         * TRUE, then NESDIS is used as the 3rd ordering instead. */
	Common->method [0].ordering = CHOLMOD_GIVEN ;/* skip if UserPerm NULL */
	Common->method [1].ordering = CHOLMOD_AMD ;
	Common->method [2].ordering = 
	    (Common->default_nesdis ? CHOLMOD_NESDIS : CHOLMOD_METIS) ;
        amd_backup = FALSE ;
#ifndef NPARTITION
	nmethods = 3 ;
#else
	nmethods = 2 ;
#endif
    }
    else
    {
        /* If only METIS and NESDIS are selected, or if 2 or more methods are
         * being tried, then enable AMD backup */
        amd_backup = (nmethods > 1) || (nmethods == 1 &&
            (Common->method [0].ordering == CHOLMOD_METIS ||
             Common->method [0].ordering == CHOLMOD_NESDIS)) ;
    }

#ifdef NSUPERNODAL
    /* CHOLMOD Supernodal module not installed, just do simplicial analysis */
    Common->supernodal = CHOLMOD_SIMPLICIAL ;
#endif

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* Note: enough space needs to be allocated here so that routines called by
     * cholmod_analyze do not reallocate the space.
     */

    /* s = 6*n + uncol */
    s = CHOLMOD(mult_size_t) (n, 6, &ok) ;
    s = CHOLMOD(add_size_t) (s, uncol, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }

    CHOLMOD(allocate_work) (n, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ensure that subsequent routines, called by cholmod_analyze, do not
     * reallocate any workspace.  This is set back to FALSE in the
     * FREE_WORKSPACE_AND_RETURN macro, which is the only way this function
     * returns to its caller. */
    Common->no_workspace_reallocate = TRUE ;

    /* Use the last 4*n Int's in Iwork for Parent, First, Level, and Post, since
     * other CHOLMOD routines will use the first 2n+uncol space.  The ordering
     * routines (cholmod_amd, cholmod_colamd, cholmod_ccolamd, cholmod_metis)
     * are an exception.  They can use all 6n + ncol space, since the contents
     * of Parent, First, Level, and Post are not needed across calls to those
     * routines. */
    Work4n = Common->Iwork ;
    Work4n += 2*((size_t) n) + uncol ;
    Parent = Work4n ;
    First  = Work4n + n ;
    Level  = Work4n + 2*((size_t) n) ;
    Post   = Work4n + 3*((size_t) n) ;

    /* note that this assignment means that cholmod_nested_dissection,
     * cholmod_ccolamd, and cholmod_camd can use only the first 4n+uncol
     * space in Common->Iwork */
    Cmember = Post ;
    CParent = Level ;

    /* ---------------------------------------------------------------------- */
    /* allocate more workspace, and an empty simplicial symbolic factor */
    /* ---------------------------------------------------------------------- */

    L = CHOLMOD(allocate_factor) (n, Common) ;
    Lparent  = CHOLMOD(malloc) (n, sizeof (Int), Common) ;
    Perm     = CHOLMOD(malloc) (n, sizeof (Int), Common) ;
    ColCount = CHOLMOD(malloc) (n, sizeof (Int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	FREE_WORKSPACE_AND_RETURN ;
    }
    Lperm = L->Perm ;
    Lcolcount = L->ColCount ;
    Common->anz = EMPTY ;

    /* ---------------------------------------------------------------------- */
    /* try all the requested ordering options and backup to AMD if needed */
    /* ---------------------------------------------------------------------- */

    /* turn off error handling [ */
    Common->try_catch = TRUE ;

    for (method = 0 ; method <= nmethods ; method++)
    {

	/* ------------------------------------------------------------------ */
	/* determine the method to try */
	/* ------------------------------------------------------------------ */

	Common->fl = EMPTY ;
	Common->lnz = EMPTY ;
	skip_analysis = FALSE ;

	if (method == nmethods)
	{
	    /* All methods failed: backup to AMD */
	    if (Common->selected == EMPTY && amd_backup)
	    {
		PRINT1 (("All methods requested failed: backup to AMD\n")) ;
		ordering = CHOLMOD_AMD ;
	    }
	    else
	    {
		break ;
	    }
	}
	else
	{
	    ordering = Common->method [method].ordering ;
	}
	Common->current = method ;
	PRINT1 (("method "ID": Try method: "ID"\n", method, ordering)) ;

	/* ------------------------------------------------------------------ */
	/* find the fill-reducing permutation */
	/* ------------------------------------------------------------------ */

	if (ordering == CHOLMOD_NATURAL)
	{

	    /* -------------------------------------------------------------- */
	    /* natural ordering */
	    /* -------------------------------------------------------------- */

	    for (k = 0 ; k < n ; k++)
	    {
		Perm [k] = k ;
	    }

	}
	else if (ordering == CHOLMOD_GIVEN)
	{

	    /* -------------------------------------------------------------- */
	    /* use given ordering of A, if provided */
	    /* -------------------------------------------------------------- */

	    if (UserPerm == NULL)
	    {
		/* this is not an error condition */
		PRINT1 (("skip, no user perm given\n")) ;
		continue ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		/* UserPerm is checked in cholmod_ptranspose */
		Perm [k] = UserPerm [k] ;
	    }

	}
	else if (ordering == CHOLMOD_AMD)
	{

	    /* -------------------------------------------------------------- */
	    /* AMD ordering of A, A*A', or A(:,f)*A(:,f)' */
	    /* -------------------------------------------------------------- */

            amd_backup = FALSE ;    /* no need to try AMD twice ... */
	    CHOLMOD(amd) (A, fset, fsize, Perm, Common) ;
	    skip_analysis = TRUE ;

	}
	else if (ordering == CHOLMOD_COLAMD)
	{

	    /* -------------------------------------------------------------- */
	    /* AMD for symmetric case, COLAMD for A*A' or A(:,f)*A(:,f)' */
	    /* -------------------------------------------------------------- */

	    if (A->stype)
	    {
		CHOLMOD(amd) (A, fset, fsize, Perm, Common) ;
		skip_analysis = TRUE ;
	    }
	    else
	    {
		/* Alternative:
		CHOLMOD(ccolamd) (A, fset, fsize, NULL, Perm, Common) ;
		*/
		/* do not postorder, it is done later, below */
		/* workspace: Iwork (4*nrow+uncol), Flag (nrow), Head (nrow+1)*/
		CHOLMOD(colamd) (A, fset, fsize, FALSE, Perm, Common) ;
	    }

	}
	else if (ordering == CHOLMOD_METIS)
	{

	    /* -------------------------------------------------------------- */
	    /* use METIS_NodeND directly (via a CHOLMOD wrapper) */
	    /* -------------------------------------------------------------- */

#ifndef NPARTITION
	    /* postorder parameter is false, because it will be later, below */
	    /* workspace: Iwork (4*nrow+uncol), Flag (nrow), Head (nrow+1) */
	    Common->called_nd = TRUE ;
	    CHOLMOD(metis) (A, fset, fsize, FALSE, Perm, Common) ;
#else
	    Common->status = CHOLMOD_NOT_INSTALLED ;
#endif

	}
	else if (ordering == CHOLMOD_NESDIS)
	{

	    /* -------------------------------------------------------------- */
	    /* use CHOLMOD's nested dissection */
	    /* -------------------------------------------------------------- */

	    /* this method is based on METIS' node bissection routine
	     * (METIS_NodeComputeSeparator).  In contrast to METIS_NodeND,
	     * it calls CAMD or CCOLAMD on the whole graph, instead of MMD
	     * on just the leaves. */
#ifndef NPARTITION
	    /* workspace: Flag (nrow), Head (nrow+1), Iwork (2*nrow) */
	    Common->called_nd = TRUE ;
	    CHOLMOD(nested_dissection) (A, fset, fsize, Perm, CParent, Cmember,
		    Common) ;
#else
	    Common->status = CHOLMOD_NOT_INSTALLED ;
#endif

	}
	else
	{

	    /* -------------------------------------------------------------- */
	    /* invalid ordering method */
	    /* -------------------------------------------------------------- */

	    Common->status = CHOLMOD_INVALID ;
	    PRINT1 (("No such ordering: "ID"\n", ordering)) ;
	}

	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory, or method failed */
	    status = MIN (status, Common->status) ;
	    Common->status = CHOLMOD_OK ;
	    continue ;
	}

	/* ------------------------------------------------------------------ */
	/* analyze the ordering */
	/* ------------------------------------------------------------------ */

	if (!skip_analysis)
	{
	    if (!CHOLMOD(analyze_ordering) (A, ordering, Perm, fset, fsize,
		    Parent, Post, ColCount, First, Level, Common))
	    {
		/* ordering method failed; clear status and try next method */
		status = MIN (status, Common->status) ;
		Common->status = CHOLMOD_OK ;
		continue ;
	    }
	}

	ASSERT (Common->fl >= 0 && Common->lnz >= 0) ;
	Common->method [method].fl  = Common->fl ;
	Common->method [method].lnz = Common->lnz ;
	PRINT1 (("lnz %g fl %g\n", Common->lnz, Common->fl)) ;

	/* ------------------------------------------------------------------ */
	/* pick the best method */
	/* ------------------------------------------------------------------ */

	/* fl.pt. compare, but lnz can never be NaN */
	if (Common->selected == EMPTY || Common->lnz < lnz_best)
	{
	    Common->selected = method ;
	    PRINT1 (("this is best so far, method "ID"\n", method)) ;
	    L->ordering = ordering ;
	    lnz_best = Common->lnz ;
	    for (k = 0 ; k < n ; k++)
	    {
		Lperm [k] = Perm [k] ;
	    }
	    /* save the results of cholmod_analyze_ordering, if it was called */
	    skip_best = skip_analysis ;
	    if (!skip_analysis)
	    {
		/* save the column count; becomes permanent part of L */
		for (k = 0 ; k < n ; k++)
		{
		    Lcolcount [k] = ColCount [k] ;
		}
		/* Parent is needed for weighted postordering and for supernodal
		 * analysis.  Does not become a permanent part of L */
		for (k = 0 ; k < n ; k++)
		{
		    Lparent [k] = Parent [k] ;
		}
	    }
	}

	/* ------------------------------------------------------------------ */
	/* determine if METIS is to be skipped */
	/* ------------------------------------------------------------------ */

	if (default_strategy && ordering == CHOLMOD_AMD)
	{
	    if ((Common->fl < 500 * Common->lnz) ||
		(Common->lnz < 5 * Common->anz))
	    {
		/* AMD found an ordering with less than 500 flops per nonzero in
		 * L, or one with a fill-in ratio (nnz(L)/nnz(A)) of less than
		 * 5.  This is pretty good, and it's unlikely that METIS will do
		 * better (this heuristic is based on tests on all symmetric
		 * positive definite matrices in the UF sparse matrix
		 * collection, and it works well across a wide range of
		 * problems).  METIS can take much more time than AMD. */
		break ;
	    }
	}
    }

    /* turn error printing back on ] */
    Common->try_catch = FALSE ;

    /* ---------------------------------------------------------------------- */
    /* return if no ordering method succeeded */
    /* ---------------------------------------------------------------------- */

    if (Common->selected == EMPTY)
    {
	/* All methods failed.  
	 * If two or more methods failed, they may have failed for different
	 * reasons.  Both would clear Common->status and skip to the next
	 * method.  Common->status needs to be restored here to the worst error
	 * obtained in any of the methods.  CHOLMOD_INVALID is worse
	 * than CHOLMOD_OUT_OF_MEMORY, since the former implies something may
	 * be wrong with the user's input.  CHOLMOD_OUT_OF_MEMORY is simply an
	 * indication of lack of resources. */
        if (status >= CHOLMOD_OK)
        {
            /* this can occur if nmethods = 1, ordering = CHOLMOD_GIVEN,
               but UserPerm is NULL */
            status = CHOLMOD_INVALID ;
        }
	ERROR (status, "all methods failed") ;
	FREE_WORKSPACE_AND_RETURN ;
    }

    /* ---------------------------------------------------------------------- */
    /* do the analysis for AMD, if skipped */
    /* ---------------------------------------------------------------------- */

    Common->fl  = Common->method [Common->selected].fl  ;
    Common->lnz = Common->method [Common->selected].lnz ;
    ASSERT (Common->lnz >= 0) ;

    if (skip_best)
    {
	if (!CHOLMOD(analyze_ordering) (A, L->ordering, Lperm, fset, fsize,
		Lparent, Post, Lcolcount, First, Level, Common))
	{
	    /* out of memory, or method failed */
	    FREE_WORKSPACE_AND_RETURN ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* postorder the etree, weighted by the column counts */
    /* ---------------------------------------------------------------------- */

    if (Common->postorder)
    {
	/* combine the fill-reducing ordering with the weighted postorder */
	/* workspace: Iwork (2*nrow) */
	if (CHOLMOD(postorder) (Lparent, n, Lcolcount, Post, Common) == n)
	{
	    /* use First and Level as workspace [ */
	    Int *Wi = First, *InvPost = Level ;
	    Int newchild, oldchild, newparent, oldparent ;

	    for (k = 0 ; k < n ; k++)
	    {
		Wi [k] = Lperm [Post [k]] ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		Lperm [k] = Wi [k] ;
	    }

	    for (k = 0 ; k < n ; k++)
	    {
		Wi [k] = Lcolcount [Post [k]] ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		Lcolcount [k] = Wi [k] ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		InvPost [Post [k]] = k ;
	    }

	    /* updated Lparent needed only for supernodal case */
	    for (newchild = 0 ; newchild < n ; newchild++)
	    {
		oldchild = Post [newchild] ;
		oldparent = Lparent [oldchild] ;
		newparent = (oldparent == EMPTY) ? EMPTY : InvPost [oldparent] ;
		Wi [newchild] = newparent ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		Lparent [k] = Wi [k] ;
	    }
	    /* done using Iwork as workspace ] */

	    /* L is now postordered, no longer in natural ordering */
	    if (L->ordering == CHOLMOD_NATURAL)
	    {
		L->ordering = CHOLMOD_POSTORDERED ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* supernodal analysis, if requested or if selected automatically */
    /* ---------------------------------------------------------------------- */

#ifndef NSUPERNODAL
    if (Common->supernodal > CHOLMOD_AUTO
    || (Common->supernodal == CHOLMOD_AUTO &&
	Common->lnz > 0 &&
	(Common->fl / Common->lnz) >= Common->supernodal_switch))
    {
	cholmod_sparse *S, *F, *A2, *A1 ;

	permute_matrices (A, L->ordering, Lperm, fset, fsize, TRUE,
		&A1, &A2, &S, &F, Common) ;

	/* workspace: Flag (nrow), Head (nrow), Iwork (5*nrow) */
	CHOLMOD(super_symbolic2) (for_cholesky, S, F, Lparent, L, Common) ;
	PRINT1 (("status %d\n", Common->status)) ;

	CHOLMOD(free_sparse) (&A1, Common) ;
	CHOLMOD(free_sparse) (&A2, Common) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* free temporary matrices and workspace, and return result L */
    /* ---------------------------------------------------------------------- */

    FREE_WORKSPACE_AND_RETURN ;
}
#endif
