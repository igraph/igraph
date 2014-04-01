/* ========================================================================== */
/* === Cholesky/cholmod_etree =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Compute the elimination tree of A or A'*A
 *
 * In the symmetric case, the upper triangular part of A is used.  Entries not
 * in this part of the matrix are ignored.  Computing the etree of a symmetric
 * matrix from just its lower triangular entries is not supported.
 *
 * In the unsymmetric case, all of A is used, and the etree of A'*A is computed.
 *
 * References:
 *
 * J. Liu, "A compact row storage scheme for Cholesky factors", ACM Trans.
 * Math. Software, vol 12, 1986, pp. 127-148.
 *
 * J. Liu, "The role of elimination trees in sparse factorization", SIAM J.
 * Matrix Analysis & Applic., vol 11, 1990, pp. 134-172.
 *
 * J. Gilbert, X. Li, E. Ng, B. Peyton, "Computing row and column counts for
 * sparse QR and LU factorization", BIT, vol 41, 2001, pp. 693-710.
 *
 * workspace: symmetric: Iwork (nrow), unsymmetric: Iwork (nrow+ncol)
 *
 * Supports any xtype (pattern, real, complex, or zomplex)
 */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"

/* ========================================================================== */
/* === update_etree ========================================================= */
/* ========================================================================== */

static void update_etree
(
    /* inputs, not modified */
    Int k,		/* process the edge (k,i) in the input graph */
    Int i,
    /* inputs, modified on output */
    Int Parent [ ],	/* Parent [t] = p if p is the parent of t */
    Int Ancestor [ ]	/* Ancestor [t] is the ancestor of node t in the
			   partially-constructed etree */
)
{
    Int a ;
    for ( ; ; )		/* traverse the path from k to the root of the tree */
    {
	a = Ancestor [k] ;
	if (a == i)
	{
	    /* final ancestor reached; no change to tree */
	    return ;
	}
	/* perform path compression */
	Ancestor [k] = i ;
	if (a == EMPTY)
	{
	    /* final ancestor undefined; this is a new edge in the tree */
	    Parent [k] = i ;
	    return ;
	}
	/* traverse up to the ancestor of k */
	k = a ;
    }
}

/* ========================================================================== */
/* === cholmod_etree ======================================================== */
/* ========================================================================== */

/* Find the elimination tree of A or A'*A */

int CHOLMOD(etree)
(
    /* ---- input ---- */
    cholmod_sparse *A,
    /* ---- output --- */
    Int *Parent,	/* size ncol.  Parent [j] = p if p is the parent of j */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Ap, *Ai, *Anz, *Ancestor, *Prev, *Iwork ;
    Int i, j, jprev, p, pend, nrow, ncol, packed, stype ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    stype = A->stype ;

    /* s = A->nrow + (stype ? 0 : A->ncol) */
    s = CHOLMOD(add_size_t) (A->nrow, (stype ? 0 : A->ncol), &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (0, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;	/* out of memory */
    }

    ASSERT (CHOLMOD(dump_sparse) (A, "etree", Common) >= 0) ;
    Iwork = Common->Iwork ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;	/* the number of columns of A */
    nrow = A->nrow ;	/* the number of rows of A */
    Ap = A->p ;		/* size ncol+1, column pointers for A */
    Ai = A->i ;		/* the row indices of A */
    Anz = A->nz ;	/* number of nonzeros in each column of A */
    packed = A->packed ;
    Ancestor = Iwork ;	/* size ncol (i/i/l) */

    for (j = 0 ; j < ncol ; j++)
    {
	Parent [j] = EMPTY ;
	Ancestor [j] = EMPTY ;
    }

    /* ---------------------------------------------------------------------- */
    /* compute the etree */
    /* ---------------------------------------------------------------------- */

    if (stype > 0)
    {

	/* ------------------------------------------------------------------ */
	/* symmetric (upper) case: compute etree (A) */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < ncol ; j++)
	{
	    /* for each row i in column j of triu(A), excluding the diagonal */
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i < j)
		{
		    update_etree (i, j, Parent, Ancestor) ;
		}
	    }
	}

    }
    else if (stype == 0)
    {

	/* ------------------------------------------------------------------ */
	/* unsymmetric case: compute etree (A'*A) */
	/* ------------------------------------------------------------------ */

	Prev = Iwork + ncol ;	/* size nrow (i/i/l) */
	for (i = 0 ; i < nrow ; i++)
	{
	    Prev [i] = EMPTY ;
	}
	for (j = 0 ; j < ncol ; j++)
	{
	    /* for each row i in column j of A */
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		/* a graph is constructed dynamically with one path per row
		 * of A.  If the ith row of A contains column indices
		 * (j1,j2,j3,j4) then the new graph has edges (j1,j2), (j2,j3),
		 * and (j3,j4).  When at node i of this path-graph, all edges
		 * (jprev,j) are considered, where jprev<j */
		i = Ai [p] ;
		jprev = Prev [i] ;
		if (jprev != EMPTY)
		{
		    update_etree (jprev, j, Parent, Ancestor) ;
		}
		Prev [i] = j ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* symmetric case with lower triangular part not supported */
	/* ------------------------------------------------------------------ */

	ERROR (CHOLMOD_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }

    ASSERT (CHOLMOD(dump_parent) (Parent, ncol, "Parent", Common)) ;
    return (TRUE) ;
}
#endif
