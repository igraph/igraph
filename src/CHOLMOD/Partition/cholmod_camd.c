/* ========================================================================== */
/* === Partition/cholmod_camd =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Partition Module.  Copyright (C) 2005-2013, Timothy A. Davis
 * The CHOLMOD/Partition Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* CHOLMOD interface to the CAMD ordering routine.  Orders A if the matrix is
 * symmetric.  On output, Perm [k] = i if row/column i of A is the kth
 * row/column of P*A*P'.  This corresponds to A(p,p) in MATLAB notation.
 *
 * If A is unsymmetric, cholmod_camd orders A*A'.  On output, Perm [k] = i if
 * row/column i of A*A' is the kth row/column of P*A*A'*P'.  This corresponds to
 * A(p,:)*A(p,:)' in MATLAB notation.  If f is present, A(p,f)*A(p,f)' is
 * ordered.
 *
 * Computes the flop count for a subsequent LL' factorization, the number
 * of nonzeros in L, and the number of nonzeros in the matrix ordered (A,
 * A*A' or A(:,f)*A(:,f)').
 *
 * workspace: Iwork (4*nrow). Head (nrow).
 *
 * Allocates a temporary copy of A+A' or A*A' (with
 * both upper and lower triangular parts) as input to CAMD.
 * Also allocates 3*(n+1) additional integer workspace (not in Common).
 *
 * Supports any xtype (pattern, real, complex, or zomplex)
 */

#ifndef NCAMD

#include "cholmod_internal.h"
#include "camd.h"
#include "cholmod_camd.h"

#if (CAMD_VERSION < CAMD_VERSION_CODE (2,0))
#error "CAMD v2.0 or later is required"
#endif

/* ========================================================================== */
/* === cholmod_camd ========================================================= */
/* ========================================================================== */

int CHOLMOD(camd)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    Int *Cmember,	/* size nrow.  see cholmod_ccolamd.c for description.*/
    /* ---- output ---- */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
)
{
    double Info [CAMD_INFO], Control2 [CAMD_CONTROL], *Control ;
    Int *Cp, *Len, *Nv, *Head, *Elen, *Degree, *Wi, *Next, *BucketSet,
	*Work3n, *p ;
    cholmod_sparse *C ;
    Int j, n, cnz ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    n = A->nrow ;

    /* s = 4*n */
    s = CHOLMOD(mult_size_t) (n, 4, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    RETURN_IF_NULL (Perm, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;
    if (n == 0)
    {
	/* nothing to do */
	Common->fl = 0 ;
	Common->lnz = 0 ;
	Common->anz = 0 ;
	return (TRUE) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    /* cholmod_analyze has allocated Cmember at Common->Iwork + 5*n+uncol, and
     * CParent at Common->Iwork + 4*n+uncol, where uncol is 0 if A is symmetric
     * or A->ncol otherwise.  Thus, only the first 4n integers in Common->Iwork
     * can be used here. */

    CHOLMOD(allocate_work) (n, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    p = Common->Iwork ;
    Degree = p ; p += n ;	/* size n */
    Elen   = p ; p += n ;	/* size n */
    Len    = p ; p += n ;	/* size n */
    Nv     = p ; p += n ;	/* size n */

    Work3n = CHOLMOD(malloc) (n+1, 3*sizeof (Int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    p = Work3n ;
    Next = p ; p += n ;		/* size n */
    Wi   = p ; p += (n+1) ;	/* size n+1 */
    BucketSet = p ;		/* size n */

    Head = Common->Head ;	/* size n+1 */

    /* ---------------------------------------------------------------------- */
    /* construct the input matrix for CAMD */
    /* ---------------------------------------------------------------------- */

    if (A->stype == 0)
    {
	/* C = A*A' or A(:,f)*A(:,f)', add extra space of nnz(C)/2+n to C */
	C = CHOLMOD(aat) (A, fset, fsize, -2, Common) ;
    }
    else
    {
	/* C = A+A', but use only the upper triangular part of A if A->stype = 1
	 * and only the lower part of A if A->stype = -1.  Add extra space of
	 * nnz(C)/2+n to C. */
	C = CHOLMOD(copy) (A, 0, -2, Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory, fset invalid, or other error */
	CHOLMOD(free) (n+1, 3*sizeof (Int), Work3n, Common) ;
	return (FALSE) ;
    }

    Cp = C->p ;
    for (j = 0 ; j < n ; j++)
    {
	Len [j] = Cp [j+1] - Cp [j] ;
    }

    /* C does not include the diagonal, and both upper and lower parts.
     * Common->anz includes the diagonal, and just the lower part of C */
    cnz = Cp [n] ;
    Common->anz = cnz / 2 + n ;

    /* ---------------------------------------------------------------------- */
    /* order C using CAMD */
    /* ---------------------------------------------------------------------- */

    /* get parameters */
    if (Common->current < 0 || Common->current >= CHOLMOD_MAXMETHODS)
    {
	/* use CAMD defaults */
	Control = NULL ;
    }
    else
    {
	Control = Control2 ;
	Control [CAMD_DENSE] = Common->method [Common->current].prune_dense ;
	Control [CAMD_AGGRESSIVE] = Common->method [Common->current].aggressive;
    }

    /* CAMD_2 does not use camd_malloc and camd_free, but set these pointers
     * just be safe. */
    camd_malloc = Common->malloc_memory ;
    camd_free = Common->free_memory ;
    camd_calloc = Common->calloc_memory ;
    camd_realloc = Common->realloc_memory ;

    /* CAMD_2 doesn't print anything either, but future versions might,
     * so set the camd_printf pointer too. */
    camd_printf = Common->print_function ;

#ifdef LONG
    /* DEBUG (camd_l_debug_init ("cholmod_l_camd")) ; */
    camd_l2 (n, C->p,  C->i, Len, C->nzmax, cnz, Nv, Next, Perm, Head, Elen,
	    Degree, Wi, Control, Info, Cmember, BucketSet) ;
#else
    /* DEBUG (camd_debug_init ("cholmod_camd")) ; */
    camd_2 (n, C->p,  C->i, Len, C->nzmax, cnz, Nv, Next, Perm, Head, Elen,
	    Degree, Wi, Control, Info, Cmember, BucketSet) ;
#endif

    /* LL' flop count.  Need to subtract n for LL' flop count.  Note that this
     * is a slight upper bound which is often exact (see CAMD/Source/camd_2.c
     * for details).  cholmod_analyze computes an exact flop count and
     * fill-in. */
    Common->fl = Info [CAMD_NDIV] + 2 * Info [CAMD_NMULTSUBS_LDL] + n ;

    /* Info [CAMD_LNZ] excludes the diagonal */
    Common->lnz = n + Info [CAMD_LNZ] ;

    /* ---------------------------------------------------------------------- */
    /* free the CAMD workspace and clear the persistent workspace in Common */
    /* ---------------------------------------------------------------------- */

    ASSERT (IMPLIES (Common->status == CHOLMOD_OK,
		CHOLMOD(dump_perm) (Perm, n, n, "CAMD2 perm", Common))) ;
    CHOLMOD(free_sparse) (&C, Common) ;
    for (j = 0 ; j <= n ; j++)
    {
	Head [j] = EMPTY ;
    }
    CHOLMOD(free) (n+1, 3*sizeof (Int), Work3n, Common) ;
    return (TRUE) ;
}
#endif
