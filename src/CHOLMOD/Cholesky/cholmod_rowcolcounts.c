/* ========================================================================== */
/* === Cholesky/cholmod_rowcolcounts ======================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Compute the row and column counts of the Cholesky factor L of the matrix
 * A or A*A'.  The etree and its postordering must already be computed (see
 * cholmod_etree and cholmod_postorder) and given as inputs to this routine.
 *
 * For the symmetric case (LL'=A), A is accessed by column.  Only the lower
 * triangular part of A is used.  Entries not in this part of the matrix are
 * ignored.  This is the same as storing the upper triangular part of A by
 * rows, with entries in the lower triangular part being ignored.  NOTE: this
 * representation is the TRANSPOSE of the input to cholmod_etree.
 *
 * For the unsymmetric case (LL'=AA'), A is accessed by column.  Equivalently,
 * if A is viewed as a matrix in compressed-row form, this routine computes
 * the row and column counts for L where LL'=A'A.  If the input vector f is
 * present, then F*F' is analyzed instead, where F = A(:,f).
 *
 * The set f is held in fset and fsize.
 *	fset = NULL means ":" in MATLAB. fset is ignored.
 *	fset != NULL means f = fset [0..fset-1].
 *	fset != NULL and fsize = 0 means f is the empty set.
 *	Common->status is set to CHOLMOD_INVALID if fset is invalid.
 *
 * In both cases, the columns of A need not be sorted.
 * A can be packed or unpacked.
 *
 * References:
 * J. Gilbert, E. Ng, B. Peyton, "An efficient algorithm to compute row and
 * column counts for sparse Cholesky factorization", SIAM J. Matrix Analysis &
 * Applic., vol 15, 1994, pp. 1075-1091.
 *
 * J. Gilbert, X. Li, E. Ng, B. Peyton, "Computing row and column counts for
 * sparse QR and LU factorization", BIT, vol 41, 2001, pp. 693-710.
 *
 * workspace:
 *	if symmetric:   Flag (nrow), Iwork (2*nrow)
 *	if unsymmetric: Flag (nrow), Iwork (2*nrow+ncol), Head (nrow+1)
 *
 * Supports any xtype (pattern, real, complex, or zomplex).
 */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"

/* ========================================================================== */
/* === initialize_node ====================================================== */
/* ========================================================================== */

static int initialize_node  /* initial work for kth node in postordered etree */
(
    Int k,		/* at the kth step of the algorithm (and kth node) */
    Int Post [ ],	/* Post [k] = i, the kth node in postordered etree */
    Int Parent [ ],	/* Parent [i] is the parent of i in the etree */
    Int ColCount [ ],	/* ColCount [c] is the current weight of node c */
    Int PrevNbr [ ]	/* PrevNbr [u] = k if u was last considered at step k */
)
{
    Int p, parent ;
    /* determine p, the kth node in the postordered etree */
    p = Post [k] ;
    /* adjust the weight if p is not a root of the etree */
    parent = Parent [p] ;
    if (parent != EMPTY)
    {
	ColCount [parent]-- ;
    }
    /* flag node p to exclude self edges (p,p) */
    PrevNbr [p] = k ;
    return (p) ;
}


/* ========================================================================== */
/* === process_edge ========================================================= */
/* ========================================================================== */

/* edge (p,u) is being processed.  p < u is a descendant of its ancestor u in
 * the etree.  node p is the kth node in the postordered etree.  */

static void process_edge
(
    Int p,		/* process edge (p,u) of the matrix */
    Int u,
    Int k,		/* we are at the kth node in the postordered etree */
    Int First [ ],	/* First [i] = k if the postordering of first
			 * descendent of node i is k */
    Int PrevNbr [ ],	/* u was last considered at step k = PrevNbr [u] */
    Int ColCount [ ],	/* ColCount [c] is the current weight of node c */
    Int PrevLeaf [ ],	/* s = PrevLeaf [u] means that s was the last leaf
			 * seen in the subtree rooted at u.  */
    Int RowCount [ ],	/* RowCount [i] is # of nonzeros in row i of L,
			 * including the diagonal.  Not computed if NULL. */
    Int SetParent [ ],	/* the FIND/UNION data structure, which forms a set
			 * of trees.  A root i has i = SetParent [i].  Following
			 * a path from i to the root q of the subtree containing
			 * i means that q is the SetParent representative of i.
			 * All nodes in the tree could have their SetParent
			 * equal to the root q; the tree representation is used
			 * to save time.  When a path is traced from i to its
			 * root q, the path is re-traversed to set the SetParent
			 * of the whole path to be the root q. */
    Int Level [ ]	 /* Level [i] = length of path from node i to root */
)
{
    Int prevleaf, q, s, sparent ;
    if (First [p] > PrevNbr [u])
    {
	/* p is a leaf of the subtree of u */
	ColCount [p]++ ;
	prevleaf = PrevLeaf [u] ;
	if (prevleaf == EMPTY)
	{
	    /* p is the first leaf of subtree of u; RowCount will be incremented
	     * by the length of the path in the etree from p up to u. */
	    q = u ;
	}
	else
	{
	    /* q = FIND (prevleaf): find the root q of the
	     * SetParent tree containing prevleaf */
	    for (q = prevleaf ; q != SetParent [q] ; q = SetParent [q])
	    {
		;
	    }
	    /* the root q has been found; re-traverse the path and
	     * perform path compression */
	    s = prevleaf ;
	    for (s = prevleaf ; s != q ; s = sparent)
	    {
		sparent = SetParent [s] ;
		SetParent [s] = q ;
	    }
	    /* adjust the RowCount and ColCount; RowCount will be incremented by
	     * the length of the path from p to the SetParent root q, and
	     * decrement the ColCount of q by one. */
	    ColCount [q]-- ;
	}
	if (RowCount != NULL)
	{
	    /* if RowCount is being computed, increment it by the length of
	     * the path from p to q */
	    RowCount [u] += (Level [p] - Level [q]) ;
	}
	/* p is a leaf of the subtree of u, so mark PrevLeaf [u] to be p */
	PrevLeaf [u] = p ;
    }
    /* flag u has having been processed at step k */
    PrevNbr [u] = k ;
}


/* ========================================================================== */
/* === finalize_node ======================================================== */
/* ========================================================================== */

static void finalize_node    /* compute UNION (p, Parent [p]) */
(
    Int p,
    Int Parent [ ],	/* Parent [p] is the parent of p in the etree */
    Int SetParent [ ]	/* see process_edge, above */
)
{
    /* all nodes in the SetParent tree rooted at p now have as their final
     * root the node Parent [p].  This computes UNION (p, Parent [p]) */
    if (Parent [p] != EMPTY)
    {
	SetParent [p] = Parent [p] ;
    }
}


/* ========================================================================== */
/* === cholmod_rowcolcounts ================================================= */
/* ========================================================================== */

int CHOLMOD(rowcolcounts)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    Int *Parent,	/* size nrow.  Parent [i] = p if p is the parent of i */
    Int *Post,		/* size nrow.  Post [k] = i if i is the kth node in
			 * the postordered etree. */
    /* ---- output --- */
    Int *RowCount,	/* size nrow. RowCount [i] = # entries in the ith row of
			 * L, including the diagonal. */
    Int *ColCount,	/* size nrow. ColCount [i] = # entries in the ith
			 * column of L, including the diagonal. */
    Int *First,		/* size nrow.  First [i] = k is the least postordering
			 * of any descendant of i. */
    Int *Level,		/* size nrow.  Level [i] is the length of the path from
			 * i to the root, with Level [root] = 0. */
    /* --------------- */
    cholmod_common *Common
)
{
    double fl, ff ;
    Int *Ap, *Ai, *Anz, *PrevNbr, *SetParent, *Head, *PrevLeaf, *Anext, *Ipost,
	*Iwork ;
    Int i, j, r, k, len, s, p, pend, inew, stype, nf, anz, inode, parent,
	nrow, ncol, packed, use_fset, jj ;
    size_t w ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_NULL (Post, FALSE) ;
    RETURN_IF_NULL (ColCount, FALSE) ;
    RETURN_IF_NULL (First, FALSE) ;
    RETURN_IF_NULL (Level, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    stype = A->stype ;
    if (stype > 0)
    {
	/* symmetric with upper triangular part not supported */
	ERROR (CHOLMOD_INVALID, "symmetric upper not supported") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;	/* the number of rows of A */
    ncol = A->ncol ;	/* the number of columns of A */

    /* w = 2*nrow + (stype ? 0 : ncol) */
    w = CHOLMOD(mult_size_t) (nrow, 2, &ok) ;
    w = CHOLMOD(add_size_t) (w, (stype ? 0 : ncol), &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (nrow, w, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }

    ASSERT (CHOLMOD(dump_perm) (Post, nrow, nrow, "Post", Common)) ;
    ASSERT (CHOLMOD(dump_parent) (Parent, nrow, "Parent", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;	/* size ncol+1, column pointers for A */
    Ai = A->i ;	/* the row indices of A, of size nz=Ap[ncol+1] */
    Anz = A->nz ;
    packed = A->packed ;
    ASSERT (IMPLIES (!packed, Anz != NULL)) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    SetParent = Iwork ;		    /* size nrow (i/i/l) */
    PrevNbr   = Iwork + nrow ;	    /* size nrow (i/i/l) */
    Anext     = Iwork + 2*((size_t) nrow) ;    /* size ncol (i/i/l) (unsym only) */
    PrevLeaf  = Common->Flag ;	    /* size nrow */
    Head      = Common->Head ;	    /* size nrow+1 (unsym only)*/

    /* ---------------------------------------------------------------------- */
    /* find the first descendant and level of each node in the tree */
    /* ---------------------------------------------------------------------- */

    /* First [i] = k if the postordering of first descendent of node i is k */
    /* Level [i] = length of path from node i to the root (Level [root] = 0) */

    for (i = 0 ; i < nrow ; i++)
    {
	First [i] = EMPTY ;
    }

    /* postorder traversal of the etree */
    for (k = 0 ; k < nrow ; k++)
    {
	/* node i of the etree is the kth node in the postordered etree */
	i = Post [k] ;

	/* i is a leaf if First [i] is still EMPTY */
	/* ColCount [i] starts at 1 if i is a leaf, zero otherwise */
	ColCount [i] = (First [i] == EMPTY) ? 1 : 0 ;

	/* traverse the path from node i to the root, stopping if we find a
	 * node r whose First [r] is already defined. */
	len = 0 ;
	for (r = i ; (r != EMPTY) && (First [r] == EMPTY) ; r = Parent [r])
	{
	    First [r] = k ;
	    len++ ;
	}
	if (r == EMPTY)
	{
	    /* we hit a root node, the level of which is zero */
	    len-- ;
	}
	else
	{
	    /* we stopped at node r, where Level [r] is already defined */
	    len += Level [r] ;
	}
	/* re-traverse the path from node i to r; set the level of each node */
	for (s = i ; s != r ; s = Parent [s])
	{
	    Level [s] = len-- ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* AA' case: sort columns of A according to first postordered row index */
    /* ---------------------------------------------------------------------- */

    fl = 0.0 ;
    if (stype == 0)
    {
	/* [ use PrevNbr [0..nrow-1] as workspace for Ipost */
	Ipost = PrevNbr ;
	/* Ipost [i] = k if i is the kth node in the postordered etree. */
	for (k = 0 ; k < nrow ; k++)
	{
	    Ipost [Post [k]] = k ;
	}
	use_fset = (fset != NULL) ;
	if (use_fset)
	{
	    nf = fsize ;
	    /* clear Anext to check fset */
	    for (j = 0 ; j < ncol ; j++)
	    {
		Anext [j] = -2 ;
	    }
	    /* find the first postordered row in each column of A (post,f)
	     * and place the column in the corresponding link list */
	    for (jj = 0 ; jj < nf ; jj++)
	    {
		j = fset [jj] ;
		if (j < 0 || j > ncol || Anext [j] != -2)
		{
		    /* out-of-range or duplicate entry in fset */
		    ERROR (CHOLMOD_INVALID, "fset invalid") ;
		    return (FALSE) ;
		}
		/* flag column j as having been seen */
		Anext [j] = EMPTY ;
	    }
	    /* fset is now valid */
	    ASSERT (CHOLMOD(dump_perm) (fset, nf, ncol, "fset", Common)) ;
	}
	else
	{
	    nf = ncol ;
	}
	for (jj = 0 ; jj < nf ; jj++)
	{
	    j = (use_fset) ? (fset [jj]) : jj ;
	    /* column j is in the fset; find the smallest row (if any) */
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    ff = (double) MAX (0, pend - p) ;
	    fl += ff*ff + ff ;
	    if (pend > p)
	    {
		k = Ipost [Ai [p]] ;
		for ( ; p < pend ; p++)
		{
		    inew = Ipost [Ai [p]] ;
		    k = MIN (k, inew) ;
		}
		/* place column j in link list k */
		ASSERT (k >= 0 && k < nrow) ;
		Anext [j] = Head [k] ;
		Head [k] = j ;
	    }
	}
	/* Ipost no longer needed for inverse postordering ]
	 * Head [k] contains a link list of all columns whose first
	 * postordered row index is equal to k, for k = 0 to nrow-1. */
    }

    /* ---------------------------------------------------------------------- */
    /* compute the row counts and node weights */
    /* ---------------------------------------------------------------------- */

    if (RowCount != NULL)
    {
	for (i = 0 ; i < nrow ; i++)
	{
	    RowCount [i] = 1 ;
	}
    }
    for (i = 0 ; i < nrow ; i++)
    {
	PrevLeaf [i] = EMPTY ;
	PrevNbr [i] = EMPTY ;
	SetParent [i] = i ;	/* every node is in its own set, by itself */
    }

    if (stype != 0)
    {

	/* ------------------------------------------------------------------ */
	/* symmetric case: LL' = A */
	/* ------------------------------------------------------------------ */

	/* also determine the number of entries in triu(A) */
	anz = nrow ;
	for (k = 0 ; k < nrow ; k++)
	{
	    /* j is the kth node in the postordered etree */
	    j = initialize_node (k, Post, Parent, ColCount, PrevNbr) ;

	    /* for all nonzeros A(i,j) below the diagonal, in column j of A */
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i > j)
		{
		    /* j is a descendant of i in etree(A) */
		    anz++ ;
		    process_edge (j, i, k, First, PrevNbr, ColCount,
			    PrevLeaf, RowCount, SetParent, Level) ;
		}
	    }
	    /* update SetParent: UNION (j, Parent [j]) */
	    finalize_node (j, Parent, SetParent) ;
	}
	Common->anz = anz ;
    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* unsymmetric case: LL' = AA' */
	/* ------------------------------------------------------------------ */

	for (k = 0 ; k < nrow ; k++)
	{
	    /* inode is the kth node in the postordered etree */
	    inode = initialize_node (k, Post, Parent, ColCount, PrevNbr) ;

	    /* for all cols j whose first postordered row is k: */
	    for (j = Head [k] ; j != EMPTY ; j = Anext [j])
	    {
		/* k is the first postordered row in column j of A */
		/* for all rows i in column j: */
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    /* has i already been considered at this step k */
		    if (PrevNbr [i] < k)
		    {
			/* inode is a descendant of i in etree(AA') */
			/* process edge (inode,i) and set PrevNbr[i] to k */
			process_edge (inode, i, k, First, PrevNbr, ColCount,
				PrevLeaf, RowCount, SetParent, Level) ;
		    }
		}
	    }
	    /* clear link list k */
	    Head [k] = EMPTY ;
	    /* update SetParent: UNION (inode, Parent [inode]) */
	    finalize_node (inode, Parent, SetParent) ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* finish computing the column counts */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < nrow ; j++)
    {
	parent = Parent [j] ;
	if (parent != EMPTY)
	{
	    /* add the ColCount of j to its parent */
	    ColCount [parent] += ColCount [j] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace */
    /* ---------------------------------------------------------------------- */

    Common->mark = EMPTY ;
    /* CHOLMOD(clear_flag) (Common) ; */
    CHOLMOD_CLEAR_FLAG (Common) ;

    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* flop count and nnz(L) for subsequent LL' numerical factorization */
    /* ---------------------------------------------------------------------- */

    /* use double to avoid integer overflow.  lnz cannot be NaN. */
    Common->aatfl = fl ;
    Common->lnz = 0. ;
    fl = 0 ;
    for (j = 0 ; j < nrow ; j++)
    {
	ff = (double) (ColCount [j]) ;
	Common->lnz += ff ;
	fl += ff*ff ;
    }

    Common->fl = fl ;
    PRINT1 (("rowcol fl %g lnz %g\n", Common->fl, Common->lnz)) ;

    return (TRUE) ;
}
#endif
