/* ========================================================================== */
/* === Cholesky/cholmod_postorder =========================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Compute the postorder of a tree. */

#ifndef NCHOLESKY

#include "cholmod_internal.h"
#include "cholmod_cholesky.h"


/* ========================================================================== */
/* === dfs ================================================================== */
/* ========================================================================== */

/* The code below includes both a recursive and non-recursive depth-first-search
 * of a tree.  The recursive code is simpler, but can lead to stack overflow.
 * It is left here for reference, to understand what the non-recursive code
 * is computing.  To try the recursive version, uncomment the following
 * #define, or compile the code with -DRECURSIVE.  Be aware that stack
 * overflow may occur.
#define RECURSIVE
 */

#ifdef RECURSIVE

/* recursive version: a working code for reference only, not actual use */

static Int dfs			/* return the new value of k */
(
    Int p,		/* start a DFS at node p */
    Int k,		/* start the node numbering at k */
    Int Post [ ],	/* Post ordering, modified on output */
    Int Head [ ],	/* Head [p] = youngest child of p; EMPTY on output */
    Int Next [ ],	/* Next [j] = sibling of j; unmodified */
    Int Pstack [ ]	/* unused */
)
{
    Int j ;
    /* start a DFS at each child of node p */
    for (j = Head [p] ; j != EMPTY ; j = Next [j])
    {
	/* start a DFS at child node j */
	k = dfs (j, k, Post, Head, Next, Pstack) ;
    }
    Post [k++] = p ;	/* order node p as the kth node */
    Head [p] = EMPTY ;	/* link list p no longer needed */
    return (k) ;	/* the next node will be numbered k */
}

#else

/* non-recursive version for actual use */

static Int dfs		/* return the new value of k */
(
    Int p,		/* start the DFS at a root node p */
    Int k,		/* start the node numbering at k */
    Int Post [ ],	/* Post ordering, modified on output */
    Int Head [ ],	/* Head [p] = youngest child of p; EMPTY on output */
    Int Next [ ],	/* Next [j] = sibling of j; unmodified */
    Int Pstack [ ]	/* workspace of size n, undefined on input or output */
)
{
    Int j, phead ;

    /* put the root node on the stack */
    Pstack [0] = p ;
    phead = 0 ;

    /* while the stack is not empty, do: */
    while (phead >= 0)
    {
	/* grab the node p from top of the stack and get its youngest child j */
	p = Pstack [phead] ;
	j = Head [p] ;
	if (j == EMPTY)
	{
	    /* all children of p ordered.  remove p from stack and order it */
	    phead-- ;
	    Post [k++] = p ;	/* order node p as the kth node */
	}
	else
	{
	    /* leave p on the stack.  Start a DFS at child node j by putting
	     * j on the stack and removing j from the list of children of p. */
	    Head [p] = Next [j] ;
	    Pstack [++phead] = j ;
	}
    }
    return (k) ;	/* the next node will be numbered k */
}

#endif

/* ========================================================================== */
/* === cholmod_postorder ==================================================== */
/* ========================================================================== */

/* Postorder a tree.  The tree is either an elimination tree (the output from
 * from cholmod_etree) or a component tree (from cholmod_nested_dissection).
 *
 * An elimination tree is a complete tree of n nodes with Parent [j] > j or
 * Parent [j] = EMPTY if j is a root.  On output Post [0..n-1] is a complete
 * permutation vector.
 *
 * A component tree is a subset of 0..n-1.  Parent [j] = -2 if node j is not
 * in the component tree.  Parent [j] = EMPTY if j is a root of the component
 * tree, and Parent [j] is in the range 0 to n-1 if j is in the component
 * tree but not a root.  On output, Post [k] is defined only for nodes in
 * the component tree.  Post [k] = j if node j is the kth node in the
 * postordered component tree, where k is in the range 0 to the number of
 * components minus 1.
 *
 * Node j is ignored and not included in the postorder if Parent [j] < EMPTY.
 *
 * As a result, check_parent (Parent, n,...) may fail on input, since
 * cholmod_check_parent assumes Parent is an elimination tree.  Similarly,
 * cholmod_check_perm (Post, ...) may fail on output, since Post is a partial
 * permutation if Parent is a component tree.
 *
 * An optional node weight can be given.  When starting a postorder at node j,
 * the children of j are ordered in increasing order of their weight.
 * If no weights are given (Weight is NULL) then children are ordered in
 * increasing order of their node number.  The weight of a node must be in the
 * range 0 to n-1.  Weights outside that range are silently converted to that
 * range (weights < 0 are treated as zero, and weights >= n are treated as n-1).
 *
 *
 * workspace: Head (n), Iwork (2*n)
 */

SuiteSparse_long CHOLMOD(postorder)	/* return # of nodes postordered */
(
    /* ---- input ---- */
    Int *Parent,	/* size n. Parent [j] = p if p is the parent of j */
    size_t n,
    Int *Weight,	/* size n, optional. Weight [j] is weight of node j */
    /* ---- output --- */
    Int *Post,		/* size n. Post [k] = j is kth in postordered tree */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Head, *Next, *Pstack, *Iwork ;
    Int j, p, k, w, nextj ;
    size_t s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (Parent, EMPTY) ;
    RETURN_IF_NULL (Post, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = 2*n */
    s = CHOLMOD(mult_size_t) (n, 2, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (EMPTY) ;
    }

    CHOLMOD(allocate_work) (n, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (EMPTY) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Head  = Common->Head ;	/* size n+1, initially all EMPTY */
    Iwork = Common->Iwork ;
    Next  = Iwork ;		/* size n (i/i/l) */
    Pstack = Iwork + n ;	/* size n (i/i/l) */

    /* ---------------------------------------------------------------------- */
    /* construct a link list of children for each node */
    /* ---------------------------------------------------------------------- */

    if (Weight == NULL)
    {

	/* in reverse order so children are in ascending order in each list */
	for (j = n-1 ; j >= 0 ; j--)
	{
	    p = Parent [j] ;
	    if (p >= 0 && p < ((Int) n))
	    {
		/* add j to the list of children for node p */
		Next [j] = Head [p] ;
		Head [p] = j ;
	    }
	}

	/* Head [p] = j if j is the youngest (least-numbered) child of p */
	/* Next [j1] = j2 if j2 is the next-oldest sibling of j1 */

    }
    else
    {

	/* First, construct a set of link lists according to Weight.
	 *
	 * Whead [w] = j if node j is the first node in bucket w.
	 * Next [j1] = j2 if node j2 follows j1 in a link list.
	 */

	Int *Whead = Pstack ;	    /* use Pstack as workspace for Whead [ */

	for (w = 0 ; w < ((Int) n) ; w++)
	{
	    Whead [w] = EMPTY ;
	}
	/* do in forward order, so nodes that ties are ordered by node index */
	for (j = 0 ; j < ((Int) n) ; j++)
	{
	    p = Parent [j] ;
	    if (p >= 0 && p < ((Int) n))
	    {
		w = Weight [j] ;
		w = MAX (0, w) ;
		w = MIN (w, ((Int) n) - 1) ;
		/* place node j at the head of link list for weight w */
		Next [j] = Whead [w] ;
		Whead [w] = j ;
	    }
	}

	/* traverse weight buckets, placing each node in its parent's list */
	for (w = n-1 ; w >= 0 ; w--)
	{
	    for (j = Whead [w] ; j != EMPTY ; j = nextj)
	    {
		nextj = Next [j] ;
		/* put node j in the link list of its parent */
		p = Parent [j] ;
		ASSERT (p >= 0 && p < ((Int) n)) ;
		Next [j] = Head [p] ;
		Head [p] = j ;
	    }
	}

	/* Whead no longer needed ] */
	/* Head [p] = j if j is the lightest child of p */
	/* Next [j1] = j2 if j2 is the next-heaviest sibling of j1 */
    }

    /* ---------------------------------------------------------------------- */
    /* start a DFS at each root node of the etree */
    /* ---------------------------------------------------------------------- */

    k = 0 ;
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	if (Parent [j] == EMPTY)
	{
	    /* j is the root of a tree; start a DFS here */
	    k = dfs (j, k, Post, Head, Next, Pstack) ;
	}
    }

    /* this would normally be EMPTY already, unless Parent is invalid */
    for (j = 0 ; j < ((Int) n) ; j++)
    {
	Head [j] = EMPTY ;
    }

    PRINT1 (("postordered "ID" nodes\n", k)) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    return (k) ;
}
#endif
