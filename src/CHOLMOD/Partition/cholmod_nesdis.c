/* ========================================================================== */
/* === Partition/cholmod_nesdis ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Partition Module.
 * Copyright (C) 2005-2006, Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Partition Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* CHOLMOD nested dissection and graph partitioning.
 *
 * cholmod_bisect:
 *
 *	Finds a set of nodes that partitions the graph into two parts.
 *	Compresses the graph first.  Requires METIS.
 *
 * cholmod_nested_dissection:
 *
 *	Nested dissection, using its own compression and connected-commponents
 *	algorithms, an external graph partitioner (METIS), and a constrained
 *	minimum degree ordering algorithm (CCOLAMD or CSYMAMD).  Typically
 *	gives better orderings than METIS_NodeND (about 5% to 10% fewer
 *	nonzeros in L).
 *
 * cholmod_collapse_septree:
 *
 *	Prune the separator tree returned by cholmod_nested_dissection.
 *
 * This file contains several routines private to this file:
 *
 *	partition	compress and partition a graph
 *	clear_flag	clear Common->Flag, but do not modify negative entries
 *	find_components	find the connected components of a graph
 *
 * Supports any xtype (pattern, real, complex, or zomplex).
 */

#ifndef NPARTITION

#include "cholmod_internal.h"
#include "cholmod_partition.h"
#include "cholmod_cholesky.h"

/* ========================================================================== */
/* === partition ============================================================ */
/* ========================================================================== */

/* Find a set of nodes that partition a graph.  The graph must be symmetric
 * with no diagonal entries.  To compress the graph first, compress is TRUE
 * and on input Hash [j] holds the hash key for node j, which must be in the
 * range 0 to csize-1. The input graph (Cp, Ci) is destroyed.  Cew is all 1's
 * on input and output.  Cnw [j] > 0 is the initial weight of node j.  On
 * output, Cnw [i] = 0 if node i is absorbed into j and the original weight
 * Cnw [i] is added to Cnw [j].  If compress is FALSE, the graph is not
 * compressed and Cnw and Hash are unmodified.  The partition itself is held in
 * the output array Part of size n.  Part [j] is 0, 1, or 2, depending on
 * whether node j is in the left part of the graph, the right part, or the
 * separator, respectively.  Note that the input graph need not be connected,
 * and the output subgraphs (the three parts) may also be unconnected.
 *
 * Returns the size of the separator, in terms of the sum of the weights of
 * the nodes.  It is guaranteed to be between 1 and the total weight of all
 * the nodes.  If it is of size less than the total weight, then both the left
 * and right parts are guaranteed to be non-empty (this guarantee depends on
 * cholmod_metis_bisector).
 */

static SuiteSparse_long partition    /* size of separator or -1 if failure */
(
    /* inputs, not modified on output */
#ifndef NDEBUG
    Int csize,		/* upper bound on # of edges in the graph;
			 * csize >= MAX (n, nnz(C)) must hold. */
#endif
    int compress,	/* if TRUE the compress the graph first */

    /* input/output */
    Int Hash [ ],	/* Hash [i] = hash >= 0 is the hash function for node
			 * i on input.  On output, Hash [i] = FLIP (j) if node
			 * i is absorbed into j.  Hash [i] >= 0 if i has not
			 * been absorbed. */

    /* input graph, compressed graph of cn nodes on output */
    cholmod_sparse *C,

    /* input/output */
    Int Cnw [ ],	/* size n.  Cnw [j] > 0 is the weight of node j on
			 * input.  On output, if node i is absorbed into
			 * node j, then Cnw [i] = 0 and the original weight of
			 * node i is added to Cnw [j].  The sum of Cnw [0..n-1]
			 * is not modified. */

    /* workspace */
    Int Cew [ ],	/* size csize, all 1's on input and output */

    /* more workspace, undefined on input and output */
    Int Cmap [ ],	/* size n (i/i/l) */

    /* output */
    Int Part [ ],	/* size n, Part [j] = 0, 1, or 2. */

    cholmod_common *Common
)
{
    Int n, hash, head, i, j, k, p, pend, ilen, ilast, pi, piend,
	jlen, ok, cn, csep, pdest, nodes_pruned, nz, total_weight, jscattered ;
    Int *Cp, *Ci, *Next, *Hhead ;

#ifndef NDEBUG
    Int cnt, pruned ;
    double work = 0, goodwork = 0 ;
#endif

    /* ---------------------------------------------------------------------- */
    /* quick return for small or empty graphs */
    /* ---------------------------------------------------------------------- */

    n = C->nrow ;
    Cp = C->p ;
    Ci = C->i ;
    nz = Cp [n] ;

    PRINT2 (("Partition start, n "ID" nz "ID"\n", n, nz)) ;

    total_weight = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	ASSERT (Cnw [j] > 0) ;
	total_weight += Cnw [j] ;
    }

    if (n <= 2)
    {
	/* very small graph */
	for (j = 0 ; j < n ; j++)
	{
	    Part [j] = 2 ;
	}
	return (total_weight) ;
    }
    else if (nz <= 0)
    {
	/* no edges, this is easy */
	PRINT2 (("diagonal matrix\n")) ;
	k = n/2 ;
	for (j = 0 ; j < k ; j++)
	{
	    Part [j] = 0 ;
	}
	for ( ; j < n ; j++)
	{
	    Part [j] = 1 ;
	}
	/* ensure the separator is not empty (required by nested dissection) */
	Part [n-1] = 2 ;
	return (Cnw [n-1]) ;
    }

#ifndef NDEBUG
    ASSERT (n > 1 && nz > 0) ;
    PRINT2 (("original graph:\n")) ;
    for (j = 0 ; j < n ; j++)
    {
	PRINT2 ((""ID": ", j)) ;
	for (p = Cp [j] ; p < Cp [j+1] ; p++)
	{
	    i = Ci [p] ;
	    PRINT3 ((""ID" ", i)) ;
	    ASSERT (i >= 0 && i < n && i != j) ;
	}
	PRINT2 (("hash: "ID"\n", Hash [j])) ;
    }
    DEBUG (for (p = 0 ; p < csize ; p++) ASSERT (Cew [p] == 1)) ;
#endif

    nodes_pruned = 0 ;

    if (compress)
    {

	/* ------------------------------------------------------------------ */
	/* get workspace */
	/* ------------------------------------------------------------------ */

	Next = Part ;	/* use Part as workspace for Next [ */
	Hhead = Cew ;	/* use Cew as workspace for Hhead [ */

	/* ------------------------------------------------------------------ */
	/* create the hash buckets */
	/* ------------------------------------------------------------------ */

	for (j = 0 ; j < n ; j++)
	{
	    /* get the hash key for node j */
	    hash = Hash [j] ;
	    ASSERT (hash >= 0 && hash < csize) ;
	    head = Hhead [hash] ;
	    if (head > EMPTY)
	    {
		/* hash bucket for this hash key is empty. */
		head = EMPTY ;
	    }
	    else
	    {
		/* hash bucket for this hash key is not empty.  get old head */
		head = FLIP (head) ;
		ASSERT (head >= 0 && head < n) ;
	    }
	    /* node j becomes the new head of the hash bucket.  FLIP it so that
	     * we can tell the difference between an empty or non-empty hash
	     * bucket. */
	    Hhead [hash] = FLIP (j) ;
	    Next [j] = head ;
	    ASSERT (head >= EMPTY && head < n) ;
	}

#ifndef NDEBUG
	for (cnt = 0, k = 0 ; k < n ; k++)
	{
	    ASSERT (Hash [k] >= 0 && Hash [k] < csize) ;    /* k is alive */
	    hash = Hash [k] ;
	    ASSERT (hash >= 0 && hash < csize) ;
	    head = Hhead [hash] ;
	    ASSERT (head < EMPTY) ;	/* hash bucket not empty */
	    j = FLIP (head) ;
	    ASSERT (j >= 0 && j < n) ;
	    if (j == k)
	    {
		PRINT2 (("hash "ID": ", hash)) ;
		for ( ; j != EMPTY ; j = Next [j])
		{
		    PRINT3 ((" "ID"", j)) ;
		    ASSERT (j >= 0 && j < n) ;
		    ASSERT (Hash [j] == hash) ;
		    cnt++ ;
		    ASSERT (cnt <= n) ;
		}
		PRINT2 (("\n")) ;
	    }
	}
	ASSERT (cnt == n) ;
#endif

	/* ------------------------------------------------------------------ */
	/* scan the non-empty hash buckets for indistinguishable nodes */
	/* ------------------------------------------------------------------ */

	/* If there are no hash collisions and no compression occurs, this takes
	 * O(n) time.  If no hash collisions, but some nodes are removed, this
	 * takes time O(n+e) where e is the sum of the degress of the nodes
	 * that are removed.  Even with many hash collisions (a rare case),
	 * this algorithm has never been observed to perform more than nnz(A)
	 * useless work.
	 *
	 * Cmap is used as workspace to mark nodes of the graph, [
	 * for comparing the nonzero patterns of two nodes i and j.
	 */

#define Cmap_MARK(i)   Cmap [i] = j
#define Cmap_MARKED(i) (Cmap [i] == j)

	for (i = 0 ; i < n ; i++)
	{
	    Cmap [i] = EMPTY ;
	}

	for (k = 0 ; k < n ; k++)
	{
	    hash = Hash [k] ;
	    ASSERT (hash >= FLIP (n-1) && hash < csize) ;
	    if (hash < 0)
	    {
		/* node k has already been absorbed into some other node */
		ASSERT (FLIP (Hash [k]) >= 0 && FLIP (Hash [k] < n)) ;
		continue ;
	    }
	    head = Hhead [hash] ;
	    ASSERT (head < EMPTY || head == 1) ;
	    if (head == 1)
	    {
		/* hash bucket is already empty */
		continue ;
	    }
	    PRINT2 (("\n--------------------hash "ID":\n", hash)) ;
	    for (j = FLIP (head) ; j != EMPTY && Next[j] > EMPTY ; j = Next [j])
	    {
		/* compare j with all nodes i following it in hash bucket */
		ASSERT (j >= 0 && j < n && Hash [j] == hash) ;
		p = Cp [j] ;
		pend = Cp [j+1] ;
		jlen = pend - p ;
		jscattered = FALSE ;
		DEBUG (for (i = 0 ; i < n ; i++) ASSERT (!Cmap_MARKED (i))) ;
		DEBUG (pruned = FALSE) ;
		ilast = j ;
		for (i = Next [j] ; i != EMPTY ; i = Next [i])
		{
		    ASSERT (i >= 0 && i < n && Hash [i] == hash && i != j) ;
		    pi = Cp [i] ;
		    piend = Cp [i+1] ;
		    ilen = piend - pi ;
		    DEBUG (work++) ;
		    if (ilen != jlen)
		    {
			/* i and j have different degrees */
			ilast = i ;
			continue ;
		    }
		    /* scatter the pattern of node j, if not already */
		    if (!jscattered)
		    {
			Cmap_MARK (j) ;
			for ( ; p < pend ; p++)
			{
			    Cmap_MARK (Ci [p]) ;
			}
			jscattered = TRUE ;
			DEBUG (work += jlen) ;
		    }
		    for (ok = Cmap_MARKED (i) ; ok && pi < piend ; pi++)
		    {
			ok = Cmap_MARKED (Ci [pi]) ;
			DEBUG (work++) ;
		    }
		    if (ok)
		    {
			/* found it.  kill node i and merge it into j */
			PRINT2 (("found "ID" absorbed into "ID"\n", i, j)) ;
			Hash [i] = FLIP (j) ;
			Cnw [j] += Cnw [i] ;
			Cnw [i] = 0 ;
			ASSERT (ilast != i && ilast >= 0 && ilast < n) ;
			Next [ilast] = Next [i] ; /* delete i from bucket */
			nodes_pruned++ ;
			DEBUG (goodwork += (ilen+1)) ;
			DEBUG (pruned = TRUE) ;
		    }
		    else
		    {
			/* i and j are different */
			ilast = i ;
		    }
		}
		DEBUG (if (pruned) goodwork += jlen) ;
	    }
	    /* empty the hash bucket, restoring Cew */
	    Hhead [hash] = 1 ;
	}

	DEBUG (if (((work - goodwork) / (double) nz) > 0.20) PRINT0 ((
	    "work %12g good %12g nz %12g (wasted work/nz: %6.2f )\n",
	    work, goodwork, (double) nz, (work - goodwork) / ((double) nz)))) ;

	/* All hash buckets now empty.  Cmap no longer needed as workspace. ]
	 * Cew no longer needed as Hhead; Cew is now restored to all ones. ]
	 * Part no longer needed as workspace for Next. ] */
    }

    /* Edge weights are all one, node weights reflect node absorption */
    DEBUG (for (p = 0 ; p < csize ; p++) ASSERT (Cew [p] == 1)) ;
    DEBUG (for (cnt = 0, j = 0 ; j < n ; j++) cnt += Cnw [j]) ;
    ASSERT (cnt == total_weight) ;

    /* ---------------------------------------------------------------------- */
    /* compress and partition the graph */
    /* ---------------------------------------------------------------------- */

    if (nodes_pruned == 0)
    {

	/* ------------------------------------------------------------------ */
	/* no pruning done at all.  Do not create the compressed graph */
	/* ------------------------------------------------------------------ */

	/* FUTURE WORK: could call CHACO, SCOTCH, ... here too */
	csep = CHOLMOD(metis_bisector) (C, Cnw, Cew, Part, Common) ;

    }
    else if (nodes_pruned == n-1)
    {

	/* ------------------------------------------------------------------ */
	/* only one node left.  This is a dense graph */
	/* ------------------------------------------------------------------ */

	PRINT2 (("completely dense graph\n")) ;
	csep = total_weight ;
	for (j = 0 ; j < n ; j++)
	{
	    Part [j] = 2 ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* compress the graph and partition the compressed graph */
	/* ------------------------------------------------------------------ */

	/* ------------------------------------------------------------------ */
	/* create the map from the uncompressed graph to the compressed graph */
	/* ------------------------------------------------------------------ */

	/* Cmap [j] = k if node j is alive and the kth node of compressed graph.
	 * The mapping is done monotonically (that is, k <= j) to simplify the
	 * uncompression later on.  Cmap [j] = EMPTY if node j is dead. */

	for (j = 0 ; j < n ; j++)
	{
	    Cmap [j] = EMPTY ;
	}
	k = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    if (Cnw [j] > 0)
	    {
		ASSERT (k <= j) ;
		Cmap [j] = k++ ;
	    }
	}
	cn = k ;	    /* # of nodes in compressed graph */
	PRINT2 (("compressed graph from "ID" to "ID" nodes\n", n, cn)) ;
	ASSERT (cn > 1 && cn == n - nodes_pruned) ;

	/* ------------------------------------------------------------------ */
	/* create the compressed graph */
	/* ------------------------------------------------------------------ */

	k = 0 ;
	pdest = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    if (Cnw [j] > 0)
	    {
		/* node j in the full graph is node k in the compressed graph */
		ASSERT (k <= j && Cmap [j] == k) ;
		p = Cp [j] ;
		pend = Cp [j+1] ;
		Cp [k] = pdest ;
		Cnw [k] = Cnw [j] ;
		for ( ; p < pend ; p++)
		{
		    /* prune dead nodes, and remap to new node numbering */
		    i = Ci [p] ;
		    ASSERT (i >= 0 && i < n && i != j) ;
		    i = Cmap [i] ;
		    ASSERT (i >= EMPTY && i < cn && i != k) ;
		    if (i > EMPTY)
		    {
			ASSERT (pdest <= p) ;
			Ci [pdest++] = i ;
		    }
		}
		k++ ;
	    }
	}
	Cp [cn] = pdest ;
	C->nrow = cn ;
	C->ncol = cn ;	/* affects mem stats unless restored when C free'd */

#ifndef NDEBUG
	PRINT2 (("pruned graph ("ID"/"ID") nodes, ("ID"/"ID") edges\n",
		    cn, n, pdest, nz)) ;
	PRINT2 (("compressed graph:\n")) ;
	for (cnt = 0, j = 0 ; j < cn ; j++)
	{
	    PRINT2 ((""ID": ", j)) ;
	    for (p = Cp [j] ; p < Cp [j+1] ; p++)
	    {
		i = Ci [p] ;
		PRINT3 ((""ID" ", i)) ;
		ASSERT (i >= 0 && i < cn && i != j) ;
	    }
	    PRINT2 (("weight: "ID"\n", Cnw [j])) ;
	    ASSERT (Cnw [j] > 0) ;
	    cnt += Cnw [j] ;
	}
	ASSERT (cnt == total_weight) ;
	for (j = 0 ; j < n ; j++) PRINT2 (("Cmap ["ID"] = "ID"\n", j, Cmap[j]));
	ASSERT (k == cn) ;
#endif

	/* ------------------------------------------------------------------ */
	/* find the separator of the compressed graph */
	/* ------------------------------------------------------------------ */

	/* FUTURE WORK: could call CHACO, SCOTCH, ... here too */
	csep = CHOLMOD(metis_bisector) (C, Cnw, Cew, Part, Common) ;

	if (csep < 0)
	{
	    /* failed */
	    return (-1) ;
	}

	PRINT2 (("Part: ")) ;
	DEBUG (for (j = 0 ; j < cn ; j++) PRINT2 ((""ID" ", Part [j]))) ;
	PRINT2 (("\n")) ;

	/* Cp and Ci no longer needed */

	/* ------------------------------------------------------------------ */
	/* find the separator of the uncompressed graph */
	/* ------------------------------------------------------------------ */

	/* expand the separator to live nodes in the uncompressed graph */
	for (j = n-1 ; j >= 0 ; j--)
	{
	    /* do this in reverse order so that Cnw can be expanded in place */
	    k = Cmap [j] ;
	    ASSERT (k >= EMPTY && k < n) ;
	    if (k > EMPTY)
	    {
		/* node k in compressed graph and is node j in full graph */
		ASSERT (k <= j) ;
		ASSERT (Hash [j] >= EMPTY) ;
		Part [j] = Part [k] ;
		Cnw [j] = Cnw [k] ;
	    }
	    else
	    {
		/* node j is a dead node */
		Cnw [j] = 0 ;
		DEBUG (Part [j] = EMPTY) ;
		ASSERT (Hash [j] < EMPTY) ;
	    }
	}

	/* find the components for the dead nodes */
	for (i = 0 ; i < n ; i++)
	{
	    if (Hash [i] < EMPTY)
	    {
		/* node i has been absorbed into node j */
		j = FLIP (Hash [i]) ;
		ASSERT (Part [i] == EMPTY && j >= 0 && j < n && Cnw [i] == 0) ;
		Part [i] = Part [j] ;
	    }
	    ASSERT (Part [i] >= 0 && Part [i] <= 2) ;
	}

#ifndef NDEBUG
	PRINT2 (("Part: ")) ;
	for (cnt = 0, j = 0 ; j < n ; j++)
	{
	    ASSERT (Part [j] != EMPTY) ;
	    PRINT2 ((""ID" ", Part [j])) ;
	    if (Part [j] == 2) cnt += Cnw [j] ;
	}
	PRINT2 (("\n")) ;
	PRINT2 (("csep "ID" "ID"\n", cnt, csep)) ;
	ASSERT (cnt == csep) ;
	for (cnt = 0, j = 0 ; j < n ; j++) cnt += Cnw [j] ;
	ASSERT (cnt == total_weight) ;
#endif

    }

    /* ---------------------------------------------------------------------- */
    /* return the separator (or -1 if error) */
    /* ---------------------------------------------------------------------- */

    PRINT2 (("Partition done, n "ID" csep "ID"\n", n, csep)) ;
    return (csep) ;
}


/* ========================================================================== */
/* === clear_flag =========================================================== */
/* ========================================================================== */

/* A node j has been removed from the graph if Flag [j] < EMPTY.
 * If Flag [j] >= EMPTY && Flag [j] < mark, then node j is alive but unmarked.
 * Flag [j] == mark means that node j is alive and marked.  Incrementing mark
 * means that all nodes are either (still) dead, or live but unmarked.
 *
 * If Map is NULL, then on output, Common->mark < Common->Flag [i] for all i
 * from 0 to Common->nrow.  This is the same output condition as
 * cholmod_clear_flag, except that this routine maintains the Flag [i] < EMPTY
 * condition as well, if that condition was true on input.
 *
 * If Map is non-NULL, then on output, Common->mark < Common->Flag [i] for all
 * i in the set Map [0..cn-1].
 *
 * workspace: Flag (nrow)
 */

static SuiteSparse_long clear_flag (Int *Map, Int cn, cholmod_common *Common)
{
    Int nrow, i ;
    Int *Flag ;
    PRINT2 (("old mark %ld\n", Common->mark)) ;
    Common->mark++ ;
    PRINT2 (("new mark %ld\n", Common->mark)) ;
    if (Common->mark <= 0)
    {
	nrow = Common->nrow ;
	Flag = Common->Flag ;
        if (Map != NULL)
        {
            for (i = 0 ; i < cn ; i++)
            {
                /* if Flag [Map [i]] < EMPTY, leave it alone */
                if (Flag [Map [i]] >= EMPTY)
                {
                    Flag [Map [i]] = EMPTY ;
                }
            }
            /* now Flag [Map [i]] <= EMPTY for all i */
        }
        else
        {
            for (i = 0 ; i < nrow ; i++)
            {
                /* if Flag [i] < EMPTY, leave it alone */
                if (Flag [i] >= EMPTY)
                {
                    Flag [i] = EMPTY ;
                }
            }
            /* now Flag [i] <= EMPTY for all i */
        }
	Common->mark = 0 ;
    }
    return (Common->mark) ;
}


/* ========================================================================== */
/* === find_components ====================================================== */
/* ========================================================================== */

/* Find all connected components of the current subgraph C.  The subgraph C
 * consists of the nodes of B that appear in the set Map [0..cn-1].  If Map
 * is NULL, then it is assumed to be the identity mapping
 * (Map [0..cn-1] = 0..cn-1).
 *
 * A node j does not appear in B if it has been ordered (Flag [j] < EMPTY,
 * which means that j has been ordered and is "deleted" from B).
 *
 * If the size of a component is large, it is placed on the component stack,
 * Cstack.  Otherwise, its nodes are ordered and it is not placed on the Cstack.
 *
 * A component S is defined by a "representative node" (repnode for short)
 * called the snode, which is one of the nodes in the subgraph.  Likewise, the
 * subgraph C is defined by its repnode, called cnode.
 * 
 * If Part is not NULL on input, then Part [i] determines how the components
 * are placed on the stack.  Components containing nodes i with Part [i] == 0
 * are placed first, followed by components with Part [i] == 1. 
 *
 * The first node placed in each of the two parts is flipped when placed in
 * the Cstack.  This allows the components of the two parts to be found simply
 * by traversing the Cstack.
 *
 * workspace: Flag (nrow)
 */

static void find_components
(
    /* inputs, not modified on output */
    cholmod_sparse *B,
    Int Map [ ],	    /* size n, only Map [0..cn-1] used */
    Int cn,		    /* # of nodes in C */
    Int cnode,		    /* root node of component C, or EMPTY if C is the
			     * entire graph B */

    Int Part [ ],	    /* size cn, optional */

    /* input/output */
    Int Bnz [ ],	    /* size n.  Bnz [j] = # nonzeros in column j of B.
			     * Reduce since B is pruned of dead nodes. */

    Int CParent [ ],	    /* CParent [i] = j if component with repnode j is
			     * the parent of the component with repnode i.
			     * CParent [i] = EMPTY if the component with
			     * repnode i is a root of the separator tree.
			     * CParent [i] is -2 if i is not a repnode. */
    Int Cstack [ ],	    /* component stack for nested dissection */
    Int *top,		    /* Cstack [0..top] contains root nodes of the
			     * the components currently in the stack */

    /* workspace, undefined on input and output: */
    Int Queue [ ],	    /* size n, for breadth-first search */

    cholmod_common *Common
)
{
    Int n, mark, cj, j, sj, sn, p, i, snode, pstart, pdest, pend, nd_components,
	part, first, save_mark ;
    Int *Bp, *Bi, *Flag ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    PRINT2 (("find components: cn %d\n", cn)) ;
    Flag = Common->Flag ;	    /* size n */

    /* force initialization of Flag [Map [0..cn-1]] */
    save_mark = Common->mark ;      /* save the current mark */
    Common->mark = EMPTY ;

    /* clear Flag; preserve Flag [Map [i]] if Flag [Map [i]] already < EMPTY */
    /* this takes O(cn) time */
    mark = clear_flag (Map, cn, Common) ;

    Bp = B->p ;
    Bi = B->i ;
    n = B->nrow ;
    ASSERT (cnode >= EMPTY && cnode < n) ;
    ASSERT (IMPLIES (cnode >= 0, Flag [cnode] < EMPTY)) ;

    /* get ordering parameters */
    nd_components = Common->method [Common->current].nd_components ;

    /* ---------------------------------------------------------------------- */
    /* find the connected components of C via a breadth-first search */
    /* ---------------------------------------------------------------------- */

    part = (Part == NULL) ? 0 : 1 ;

    /* examine each part (part 1 and then part 0) */
    for (part = (Part == NULL) ? 0 : 1 ; part >= 0 ; part--)
    {

	/* first is TRUE for the first connected component in each part */
	first = TRUE ;

	/* find all connected components in the current part */
	for (cj = 0 ; cj < cn ; cj++)
	{
	    /* get node snode, which is node cj of C.  It might already be in
	     * the separator of C (and thus ordered, with Flag [snode] < EMPTY)
	     */
	    snode = (Map == NULL) ? (cj) : (Map [cj]) ;
	    ASSERT (snode >= 0 && snode < n) ;

	    if (Flag [snode] >= EMPTY && Flag [snode] < mark
		    && ((Part == NULL) || Part [cj] == part))
	    {

		/* ---------------------------------------------------------- */
		/* find new connected component S */
		/* ---------------------------------------------------------- */

		/* node snode is the repnode of a connected component S, the
		 * parent of which is cnode, the repnode of C.  If cnode is
		 * EMPTY then C is the original graph B. */
		PRINT2 (("----------:::snode "ID" cnode "ID"\n", snode, cnode));

		ASSERT (CParent [snode] == -2) ;
		if (first || nd_components)
		{
		    /* If this is the first node in this part, then it becomes
		     * the repnode of all components in this part, and all
		     * components in this part form a single node in the
		     * separator tree.  If nd_components is TRUE, then all
		     * connected components form their own node in the
		     * separator tree.
		     */
		    CParent [snode] = cnode ;
		}

		/* place j in the queue and mark it */
		Queue [0] = snode ;
		Flag [snode] = mark ;
		sn = 1 ;

		/* breadth-first traversal, starting at node j */
		for (sj = 0 ; sj < sn ; sj++)
		{
		    /* get node j from head of Queue and traverse its edges */
		    j = Queue [sj] ;
		    PRINT2 (("    j: "ID"\n", j)) ;
		    ASSERT (j >= 0 && j < n) ;
		    ASSERT (Flag [j] == mark) ;
		    pstart = Bp [j] ;
		    pdest = pstart ;
		    pend = pstart + Bnz [j] ;
		    for (p = pstart ; p < pend ; p++)
		    {
			i = Bi [p] ;
			if (i != j && Flag [i] >= EMPTY)
			{
			    /* node is still in the graph */
			    Bi [pdest++] = i ;
			    if (Flag [i] < mark)
			    {
				/* node i is in this component S, and unflagged
				 * (first time node i has been seen in this BFS)
				 * place node i in the queue and mark it */
				Queue [sn++] = i ;
				Flag [i] = mark ;
			    }
			}
		    }
		    /* edges to dead nodes have been removed */
		    Bnz [j] = pdest - pstart ;
		}

		/* ---------------------------------------------------------- */
		/* order S if it is small; place it on Cstack otherwise */
		/* ---------------------------------------------------------- */

		PRINT2 (("sn "ID"\n", sn)) ;

		/* place the new component on the Cstack.  Flip the node if
		 * is the first connected component of the current part,
		 * or if all components are treated as their own node in
		 * the separator tree. */
		Cstack [++(*top)] =
			(first || nd_components) ? FLIP (snode) : snode ;
		first = FALSE ;
	    }
	}
    }

    /* restore the flag (normally taking O(1) time except for Int overflow) */
    Common->mark = save_mark++ ;
    clear_flag (NULL, 0, Common) ;
    DEBUG (for (i = 0 ; i < n ; i++) ASSERT (Flag [i] < Common->mark)) ;
}


/* ========================================================================== */
/* === cholmod_bisect ======================================================= */
/* ========================================================================== */

/* Finds a node bisector of A, A*A', A(:,f)*A(:,f)'.
 *
 * workspace: Flag (nrow),
 *	Iwork (nrow if symmetric, max (nrow,ncol) if unsymmetric).
 *	Allocates a temporary matrix B=A*A' or B=A,
 *	and O(nnz(A)) temporary memory space.
 */

SuiteSparse_long CHOLMOD(bisect)	/* returns # of nodes in separator */
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to bisect */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int compress,	/* if TRUE, compress the graph first */
    /* ---- output --- */
    Int *Partition,	/* size A->nrow.  Node i is in the left graph if
			 * Partition [i] = 0, the right graph if 1, and in the
			 * separator if 2. */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Bp, *Bi, *Hash, *Cmap, *Bnw, *Bew, *Iwork ;
    cholmod_sparse *B ;
    unsigned Int hash ;
    Int j, n, bnz, sepsize, p, pend ;
    size_t csize, s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_NULL (Partition, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* quick return */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    if (n == 0)
    {
	return (0) ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = n + MAX (n, A->ncol) */
    s = CHOLMOD(add_size_t) (A->nrow, MAX (A->nrow, A->ncol), &ok) ;
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

    Iwork = Common->Iwork ;
    Hash = Iwork ;		/* size n, (i/l/l) */
    Cmap = Iwork + n ;		/* size n, (i/i/l) */

    /* ---------------------------------------------------------------------- */
    /* convert the matrix to adjacency list form */
    /* ---------------------------------------------------------------------- */

    /* The input graph to must be symmetric, with no diagonal entries
     * present.  The columns need not be sorted. */

    /* B = A, A*A', or A(:,f)*A(:,f)', upper and lower parts present */

    if (A->stype)
    {
	/* Add the upper/lower part to a symmetric lower/upper matrix by
	 * converting to unsymmetric mode */
	/* workspace: Iwork (nrow) */
	B = CHOLMOD(copy) (A, 0, -1, Common) ;
    }
    else
    {
	/* B = A*A' or A(:,f)*A(:,f)', no diagonal */
	/* workspace: Flag (nrow), Iwork (max (nrow,ncol)) */
	B = CHOLMOD(aat) (A, fset, fsize, -1, Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	return (EMPTY) ;
    }
    Bp = B->p ;
    Bi = B->i ;
    bnz = Bp [n] ;
    ASSERT ((Int) (B->nrow) == n && (Int) (B->ncol) == n) ;

    /* B does not include the diagonal, and both upper and lower parts.
     * Common->anz includes the diagonal, and just the lower part of B */
    Common->anz = bnz / 2 + ((double) n) ;

    /* Bew should be at least size n for the hash function to work well */
    /* this cannot cause overflow, because the matrix is already created */
    csize = MAX (((size_t) n) + 1, (size_t) bnz) ;

    /* create the graph using Flag as workspace for node weights [ */
    Bnw = Common->Flag ;    /* size n workspace */

    /* compute hash for each node if compression requested */
    if (compress)
    {
	for (j = 0 ; j < n ; j++)
	{
	    hash = j ;
	    pend = Bp [j+1] ;
	    for (p = Bp [j] ; p < pend ; p++)
	    {
		hash += Bi [p] ;
		ASSERT (Bi [p] != j) ;
	    }
	    /* finalize the hash key for node j */
	    hash %= csize ;
	    Hash [j] = (Int) hash ;
	    ASSERT (Hash [j] >= 0 && Hash [j] < csize) ;
	}
    }

    /* allocate edge weights */
    Bew = CHOLMOD(malloc) (csize, sizeof (Int), Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	CHOLMOD(free_sparse) (&B, Common) ;
	CHOLMOD(free) (csize, sizeof (Int), Bew, Common) ;
	return (EMPTY) ;
    }

    /* graph has unit node and edge weights */
    for (j = 0 ; j < n ; j++)
    {
	Bnw [j] = 1 ;
    }
    for (s = 0 ; s < csize ; s++)
    {
	Bew [s] = 1 ;
    }

    /* ---------------------------------------------------------------------- */
    /* compress and partition the graph */
    /* ---------------------------------------------------------------------- */

    sepsize = partition (
#ifndef NDEBUG
	    csize,
#endif
	    compress, Hash, B, Bnw, Bew, Cmap, Partition, Common) ;

    /* contents of Bp, Bi, Bnw, and Bew no longer needed ] */

    /* If partition fails, free the workspace below and return sepsize < 0 */

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    B->ncol = n ;   /* restore size for memory usage statistics */
    CHOLMOD(free_sparse) (&B, Common) ;
    Common->mark = EMPTY ;
    CHOLMOD_CLEAR_FLAG (Common) ;
    CHOLMOD(free) (csize, sizeof (Int), Bew, Common) ;
    return (sepsize) ;
}


/* ========================================================================== */
/* === cholmod_nested_dissection ============================================ */
/* ========================================================================== */

/* This method uses a node bisector, applied recursively (but using a
 * non-recursive algorithm).  Once the graph is partitioned, it calls a
 * constrained min degree code (CAMD or CSYMAMD for A+A', and CCOLAMD for A*A')
 * to order all the nodes in the graph - but obeying the constraints determined
 * by the separators.  This routine is similar to METIS_NodeND, except for how
 * it treats the leaf nodes.  METIS_NodeND orders the leaves of the separator
 * tree with MMD, ignoring the rest of the matrix when ordering a single leaf.
 * This routine orders the whole matrix with CSYMAMD or CCOLAMD, all at once,
 * when the graph partitioning is done.
 *
 * This function also returns a postorderd separator tree (CParent), and a
 * mapping of nodes in the graph to nodes in the separator tree (Cmember).
 *
 * workspace: Flag (nrow), Head (nrow+1), Iwork (4*nrow + (ncol if unsymmetric))
 *	Allocates a temporary matrix B=A*A' or B=A,
 *	and O(nnz(A)) temporary memory space.
 *	Allocates an additional 3*n*sizeof(Int) temporary workspace
 */

SuiteSparse_long CHOLMOD(nested_dissection)
    /* returns # of components, or -1 if error */
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    /* ---- output --- */
    Int *Perm,		/* size A->nrow, output permutation */
    Int *CParent,	/* size A->nrow.  On output, CParent [c] is the parent
			 * of component c, or EMPTY if c is a root, and where
			 * c is in the range 0 to # of components minus 1 */
    Int *Cmember,	/* size A->nrow.  Cmember [j] = c if node j of A is
			 * in component c */
    /* --------------- */
    cholmod_common *Common
)
{
    double prune_dense, nd_oksep ;
    Int *Bp, *Bi, *Bnz, *Cstack, *Imap, *Map, *Flag, *Head, *Next, *Bnw, *Iwork,
	*Ipost, *NewParent, *Hash, *Cmap, *Cp, *Ci, *Cew, *Cnw, *Part, *Post,
	*Work3n ;
    unsigned Int hash ;
    Int n, bnz, top, i, j, k, cnode, cdense, p, cj, cn, ci, cnz, mark, c, uncol,
	sepsize, parent, ncomponents, threshold, ndense, pstart, pdest, pend,
	nd_compress, nd_camd, csize, jnext, nd_small, total_weight,
	nchild, child = EMPTY ;
    cholmod_sparse *B, *C ;
    size_t s ;
    int ok = TRUE ;
    DEBUG (Int cnt) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_NULL (Perm, EMPTY) ;
    RETURN_IF_NULL (CParent, EMPTY) ;
    RETURN_IF_NULL (Cmember, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* quick return */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    if (n == 0)
    {
	return (1) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* get ordering parameters */
    prune_dense = Common->method [Common->current].prune_dense ;
    nd_compress = Common->method [Common->current].nd_compress ;
    nd_oksep = Common->method [Common->current].nd_oksep ;
    nd_oksep = MAX (0, nd_oksep) ;
    nd_oksep = MIN (1, nd_oksep) ;
    nd_camd = Common->method [Common->current].nd_camd ;
    nd_small = Common->method [Common->current].nd_small ;
    nd_small = MAX (4, nd_small) ;

    PRINT0 (("nd_components %d nd_small %d nd_oksep %g\n", 
	Common->method [Common->current].nd_components,
	nd_small, nd_oksep)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = 4*n + uncol */
    uncol = (A->stype == 0) ? A->ncol : 0 ;
    s = CHOLMOD(mult_size_t) (n, 4, &ok) ;
    s = CHOLMOD(add_size_t) (s, uncol, &ok) ;
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
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Flag = Common->Flag ;	/* size n */
    Head = Common->Head ;	/* size n+1, all equal to -1 */

    Iwork = Common->Iwork ;
    Imap = Iwork ;		/* size n, same as Queue in find_components */
    Map  = Iwork + n ;		/* size n */
    Bnz  = Iwork + 2*((size_t) n) ;	/* size n */
    Hash = Iwork + 3*((size_t) n) ;	/* size n */

    Work3n = CHOLMOD(malloc) (n, 3*sizeof (Int), Common) ;
    Part = Work3n ;		/* size n */
    Bnw  = Part + n ;		/* size n */
    Cnw  = Bnw + n ;		/* size n */

    Cstack = Perm ;		/* size n, use Perm as workspace for Cstack [ */
    Cmap = Cmember ;		/* size n, use Cmember as workspace [ */

    if (Common->status < CHOLMOD_OK)
    {
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* convert B to symmetric form with both upper/lower parts present */
    /* ---------------------------------------------------------------------- */

    /* B = A+A', A*A', or A(:,f)*A(:,f)', upper and lower parts present */

    if (A->stype)
    {
	/* Add the upper/lower part to a symmetric lower/upper matrix by
	 * converting to unsymmetric mode */
	/* workspace: Iwork (nrow) */
	B = CHOLMOD(copy) (A, 0, -1, Common) ;
    }
    else
    {
	/* B = A*A' or A(:,f)*A(:,f)', no diagonal */
	/* workspace: Flag (nrow), Iwork (max (nrow,ncol)) */
	B = CHOLMOD(aat) (A, fset, fsize, -1, Common) ;
    }

    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free) (3*n, sizeof (Int), Work3n, Common) ;
	return (EMPTY) ;
    }
    Bp = B->p ;
    Bi = B->i ;
    bnz = CHOLMOD(nnz) (B, Common) ;
    ASSERT ((Int) (B->nrow) == n && (Int) (B->ncol) == n) ;
    csize = MAX (n, bnz) ;
    ASSERT (CHOLMOD(dump_sparse) (B, "B for nd:", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* initializations */
    /* ---------------------------------------------------------------------- */

    /* all nodes start out unmarked and unordered (Type 4, see below) */
    Common->mark = EMPTY ;
    CHOLMOD_CLEAR_FLAG (Common) ;
    ASSERT (Flag == Common->Flag) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    for (j = 0 ; j < n ; j++)
    {
	CParent [j] = -2 ;
    }

    /* prune dense nodes from B */
    if (IS_NAN (prune_dense) || prune_dense < 0)
    {
	/* only remove completely dense nodes */
	threshold = n-2 ;
    }
    else
    {
	/* remove nodes with degree more than threshold */
	threshold = (Int) (MAX (16, prune_dense * sqrt ((double) (n)))) ;
	threshold = MIN (n, threshold) ;
    }
    ndense = 0 ;
    cnode = EMPTY ;
    cdense = EMPTY ;

    for (j = 0 ; j < n ; j++)
    {
	Bnz [j] = Bp [j+1] - Bp [j] ;
	if (Bnz [j] > threshold)
	{
	    /* node j is dense, prune it from B */
	    PRINT2 (("j is dense %d\n", j)) ;
	    ndense++ ;
	    if (cnode == EMPTY)
	    {
		/* first dense node found becomes root of this component,
		 * which contains all of the dense nodes found here */
		cdense = j ;
		cnode = j ;
		CParent [cnode] = EMPTY ;
	    }
	    Flag [j] = FLIP (cnode) ;
	}
    }
    B->packed = FALSE ;
    ASSERT (B->nz == NULL) ;

    if (ndense == n)
    {
	/* all nodes removed: Perm is identity, all nodes in component zero,
	 * and the separator tree has just one node. */
	PRINT2 (("all nodes are dense\n")) ;
	for (k = 0 ; k < n ; k++)
	{
	    Perm [k] = k ;
	    Cmember [k] = 0 ;
	}
	CParent [0] = EMPTY ;
	CHOLMOD(free_sparse) (&B, Common) ;
	CHOLMOD(free) (3*n, sizeof (Int), Work3n, Common) ;
	Common->mark = EMPTY ;
	CHOLMOD_CLEAR_FLAG (Common) ;
	return (1) ;
    }

    /* Cp and Ci are workspace to construct the subgraphs to partition */
    C = CHOLMOD(allocate_sparse) (n, n, csize, FALSE, TRUE, 0, CHOLMOD_PATTERN,
	    Common) ;
    Cew  = CHOLMOD(malloc) (csize, sizeof (Int), Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	CHOLMOD(free_sparse) (&C, Common) ;
	CHOLMOD(free_sparse) (&B, Common) ;
	CHOLMOD(free) (csize, sizeof (Int), Cew, Common) ;
	CHOLMOD(free) (3*n, sizeof (Int), Work3n, Common) ;
	Common->mark = EMPTY ;
	CHOLMOD_CLEAR_FLAG (Common) ;
	PRINT2 (("out of memory for C, etc\n")) ;
	return (EMPTY) ;
    }

    Cp = C->p ;
    Ci = C->i ;

    /* create initial unit node and edge weights */
    for (j = 0 ; j < n ; j++)
    {
	Bnw [j] = 1 ;
    }
    for (p = 0 ; p < csize ; p++)
    {
	Cew [p] = 1 ;
    }

    /* push the initial connnected components of B onto the Cstack */
    top = EMPTY ;	/* Cstack is empty */
    /* workspace: Flag (nrow), Iwork (nrow); use Imap as workspace for Queue [*/
    find_components (B, NULL, n, cnode, NULL,
	    Bnz, CParent, Cstack, &top, Imap, Common) ;
    /* done using Imap as workspace for Queue ] */

    /* Nodes can now be of Type 0, 1, 2, or 4 (see definition below) */

    /* ---------------------------------------------------------------------- */
    /* while Cstack is not empty, do: */
    /* ---------------------------------------------------------------------- */

    while (top >= 0)
    {

	/* clear the Flag array, but do not modify negative entries in Flag  */
	mark = clear_flag (NULL, 0, Common) ;

	DEBUG (for (i = 0 ; i < n ; i++) Imap [i] = EMPTY) ;

	/* ------------------------------------------------------------------ */
	/* get node(s) from the top of the Cstack */
	/* ------------------------------------------------------------------ */

	/* i is the repnode of its (unordered) connected component.  Get
	 * all repnodes for all connected components of a single part.  If
	 * each connected component is to be ordered separately (nd_components
	 * is TRUE), then this while loop iterates just once. */

	cnode = EMPTY ;
	cn = 0 ;
	while (cnode == EMPTY)
	{
	    i = Cstack [top--] ;

	    if (i < 0)
	    {
		/* this is the last node in this component */
		i = FLIP (i) ;
		cnode = i ;
	    }

	    ASSERT (i >= 0 && i < n && Flag [i] >= EMPTY) ;

	    /* place i in the queue and mark it */
	    Map [cn] = i ;
	    Flag [i] = mark ;
	    Imap [i] = cn ;
	    cn++ ;
	}

	ASSERT (cnode != EMPTY) ;

	/* During ordering, there are five kinds of nodes in the graph of B,
	 * based on Flag [j] and CParent [j] for nodes j = 0 to n-1:
	 *
	 * Type 0: If cnode is a repnode of an unordered component, then
	 * CParent [cnode] is in the range EMPTY to n-1 and
	 * Flag [cnode] >= EMPTY.  This is a "live" node.
	 *
	 * Type 1: If cnode is a repnode of an ordered separator component,
	 * then Flag [cnode] < EMPTY and FLAG [cnode] = FLIP (cnode).
	 * CParent [cnode] is in the range EMPTY to n-1.  cnode is a root of
	 * the separator tree if CParent [cnode] == EMPTY.  This node is dead.
	 *
	 * Type 2: If node j isn't a repnode, has not been absorbed via
	 * graph compression into another node, but is in an ordered separator
	 * component, then cnode = FLIP (Flag [j]) gives the repnode of the
	 * component that contains j and CParent [j]  is -2.  This node is dead.
	 * Note that Flag [j] < EMPTY.
	 *
	 * Type 3: If node i has been absorbed via graph compression into some
	 * other node j = FLIP (Flag [i]) where j is not a repnode.
	 * CParent [j] is -2.  Node i may or may not be in an ordered
	 * component.  This node is dead.  Note that Flag [j] < EMPTY.
	 *
	 * Type 4: If node j is "live" (not in an ordered component, and not
	 * absorbed into any other node), then Flag [j] >= EMPTY.
	 *
	 * Only "live" nodes (of type 0 or 4) are placed in a subgraph to be
	 * partitioned.  Node j is alive if Flag [j] >= EMPTY, and dead if
	 * Flag [j] < EMPTY.
	 */

	/* ------------------------------------------------------------------ */
	/* create the subgraph for this connected component C */
	/* ------------------------------------------------------------------ */

	/* Do a breadth-first search of the graph starting at cnode.
	 * use Map [0..cn-1] for nodes in the component C [
	 * use Cnw and Cew for node and edge weights of the resulting subgraph [
	 * use Cp and Ci for the resulting subgraph [
	 * use Imap [i] for all nodes i in B that are in the component C [
	 */

	cnz = 0 ;
	total_weight = 0 ;
	for (cj = 0 ; cj < cn ; cj++)
	{
	    /* get node j from the head of the queue; it is node cj of C */
	    j = Map [cj] ;
	    ASSERT (Flag [j] == mark) ;
	    Cp [cj] = cnz ;
	    Cnw [cj] = Bnw [j] ;
	    ASSERT (Cnw [cj] >= 0) ;
	    total_weight += Cnw [cj] ;
	    pstart = Bp [j] ;
	    pdest = pstart ;
	    pend = pstart + Bnz [j] ;
	    hash = cj ;
	    for (p = pstart ; p < pend ; p++)
	    {
		i = Bi [p] ;
		/* prune diagonal entries and dead edges from B */
		if (i != j && Flag [i] >= EMPTY)
		{
		    /* live node i is in the current component */
		    Bi [pdest++] = i ;
		    if (Flag [i] != mark)
		    {
			/* First time node i has been seen, it is a new node
			 * of C.  place node i in the queue and mark it */
			Map [cn] = i ;
			Flag [i] = mark ;
			Imap [i] = cn ;
			cn++ ;
		    }
		    /* place the edge (cj,ci) in the adjacency list of cj */
		    ci = Imap [i] ;
		    ASSERT (ci >= 0 && ci < cn && ci != cj && cnz < csize) ;
		    Ci [cnz++] = ci ;
		    hash += ci ;
		}
	    }
	    /* edges to dead nodes have been removed */
	    Bnz [j] = pdest - pstart ;
	    /* finalize the hash key for column j */
	    hash %= csize ;
	    Hash [cj] = (Int) hash ;
	    ASSERT (Hash [cj] >= 0 && Hash [cj] < csize) ;
	}
	Cp [cn] = cnz ;
	C->nrow = cn ;
	C->ncol = cn ;	/* affects mem stats unless restored when C free'd */

	/* contents of Imap no longer needed ] */

#ifndef NDEBUG
	for (cj = 0 ; cj < cn ; cj++)
	{
	    j = Map [cj] ;
	    PRINT2 (("----------------------------C column cj: "ID" j: "ID"\n",
		cj, j)) ;
	    ASSERT (j >= 0 && j < n) ;
	    ASSERT (Flag [j] >= EMPTY) ;
	    for (p = Cp [cj] ; p < Cp [cj+1] ; p++)
	    {
		ci = Ci [p] ;
		i = Map [ci] ;
		PRINT3 (("ci: "ID" i: "ID"\n", ci, i)) ;
		ASSERT (ci != cj && ci >= 0 && ci < cn) ;
		ASSERT (i != j && i >= 0 && i < n) ;
		ASSERT (Flag [i] >= EMPTY) ;
	    }
	}
#endif

	PRINT0 (("consider cn %d nd_small %d ", cn, nd_small)) ;
	if (cn < nd_small)  /* could be 'total_weight < nd_small' instead */
	{
	    /* place all nodes in the separator */
	    PRINT0 ((" too small\n")) ;
	    sepsize = total_weight ;
	}
	else
	{

	    /* Cp and Ci now contain the component, with cn nodes and cnz
	     * nonzeros.  The mapping of a node cj into node j the main graph
	     * B is given by Map [cj] = j */
	    PRINT0 ((" cut\n")) ;

	    /* -------------------------------------------------------------- */
	    /* compress and partition the graph C */
	    /* -------------------------------------------------------------- */

	    /* The edge weights Cew [0..csize-1] are all 1's on input to and
	     * output from the partition routine. */

	    sepsize = partition (
#ifndef NDEBUG
		    csize,
#endif
		    nd_compress, Hash, C, Cnw, Cew,
		    Cmap, Part, Common) ;

	    /* contents of Cp and Ci no longer needed ] */

	    if (sepsize < 0)
	    {
		/* failed */
		C->ncol = n ;   /* restore size for memory usage statistics */
		CHOLMOD(free_sparse) (&C, Common) ;
		CHOLMOD(free_sparse) (&B, Common) ;
		CHOLMOD(free) (csize, sizeof (Int), Cew, Common) ;
		CHOLMOD(free) (3*n, sizeof (Int), Work3n, Common) ;
		Common->mark = EMPTY ;
		CHOLMOD_CLEAR_FLAG (Common) ;
		return (EMPTY) ;
	    }

	    /* -------------------------------------------------------------- */
	    /* compress B based on how C was compressed */
	    /* -------------------------------------------------------------- */

	    for (ci = 0 ; ci < cn ; ci++)
	    {
		if (Hash [ci] < EMPTY)
		{
		    /* ci is dead in C, having been absorbed into cj */
		    cj = FLIP (Hash [ci]) ;
		    PRINT2 (("In C, "ID" absorbed into "ID" (wgt now "ID")\n",
			    ci, cj, Cnw [cj])) ;
		    /* i is dead in B, having been absorbed into j */
		    i = Map [ci] ;
		    j = Map [cj] ;
		    PRINT2 (("In B, "ID" (wgt "ID") => "ID" (wgt "ID")\n",
				i, Bnw [i], j, Bnw [j], Cnw [cj])) ;
		    /* more than one node may be absorbed into j.  This is
		     * accounted for in Cnw [cj].  Assign it here rather
		     * than += Bnw [i] */
		    Bnw [i] = 0 ;
		    Bnw [j] = Cnw [cj] ;
		    Flag [i] = FLIP (j) ;
		}
	    }

	    DEBUG (for (cnt = 0, j = 0 ; j < n ; j++) cnt += Bnw [j]) ;
	    ASSERT (cnt == n) ;
	}

	/* contents of Cnw [0..cn-1] no longer needed ] */

	/* ------------------------------------------------------------------ */
	/* order the separator, and stack the components when C is split */
	/* ------------------------------------------------------------------ */

	/* one more component has been found: either the separator of C,
	 * or all of C */

	ASSERT (sepsize >= 0 && sepsize <= total_weight) ;

	PRINT0 (("sepsize %d tot %d : %8.4f ", sepsize, total_weight,
	    ((double) sepsize) / ((double) total_weight))) ;

	if (sepsize == total_weight || sepsize == 0 ||
	    sepsize > nd_oksep * total_weight)
	{
	    /* Order the nodes in the component.  The separator is too large,
	     * or empty.  Note that the partition routine cannot return a
	     * sepsize of zero, but it can return a separator consisting of the
	     * whole graph.  The "sepsize == 0" test is kept, above, in case the
	     * partition routine changes.  In either case, this component
	     * remains unsplit, and becomes a leaf of the separator tree. */
	    PRINT2 (("cnode %d sepsize zero or all of graph: "ID"\n",
		cnode, sepsize)) ;
	    for (cj = 0 ; cj < cn ; cj++)
	    {
		j = Map [cj] ;
		Flag [j] = FLIP (cnode) ;
		PRINT2 (("      node cj: "ID" j: "ID" ordered\n", cj, j)) ;
	    }
	    ASSERT (Flag [cnode] == FLIP (cnode)) ;
	    ASSERT (cnode != EMPTY && Flag [cnode] < EMPTY) ;
	    PRINT0 (("discarded\n")) ;

	}
	else
	{

	    /* Order the nodes in the separator of C and find a new repnode
	     * cnode that is in the separator of C.  This requires the separator
	     * to be non-empty. */
	    PRINT0 (("sepsize not tiny: "ID"\n", sepsize)) ;
	    parent = CParent [cnode] ;
	    ASSERT (parent >= EMPTY && parent < n) ;
	    CParent [cnode] = -2 ;
	    cnode = EMPTY ;
	    for (cj = 0 ; cj < cn ; cj++)
	    {
		j = Map [cj] ;
		if (Part [cj] == 2)
		{
		    /* All nodes in the separator become part of a component
		     * whose repnode is cnode */
		    PRINT2 (("node cj: "ID" j: "ID" ordered\n", cj, j)) ;
		    if (cnode == EMPTY)
		    {
			PRINT2(("------------new cnode: cj "ID" j "ID"\n",
				    cj, j)) ;
			cnode = j ;
		    }
		    Flag [j] = FLIP (cnode) ;
		}
		else
		{
		    PRINT2 (("      node cj: "ID" j: "ID" not ordered\n",
				cj, j)) ;
		}
	    }
	    ASSERT (cnode != EMPTY && Flag [cnode] < EMPTY) ;
	    ASSERT (CParent [cnode] == -2) ;
	    CParent [cnode] = parent ;

	    /* find the connected components when C is split, and push
	     * them on the Cstack.  Use Imap as workspace for Queue. [ */
	    /* workspace: Flag (nrow) */
	    find_components (B, Map, cn, cnode, Part, Bnz,
		    CParent, Cstack, &top, Imap, Common) ;
	    /* done using Imap as workspace for Queue ] */
	}
	/* contents of Map [0..cn-1] no longer needed ] */
    }

    /* done using Cmember as workspace for Cmap ] */
    /* done using Perm as workspace for Cstack ] */

    /* ---------------------------------------------------------------------- */
    /* place nodes removed via compression into their proper component */
    /* ---------------------------------------------------------------------- */

    /* At this point, all nodes are of Type 1, 2, or 3, as defined above. */

    for (i = 0 ; i < n ; i++)
    {
	/* find the repnode cnode that contains node i */
	j = FLIP (Flag [i]) ;
	PRINT2 (("\nfind component for "ID", in: "ID"\n", i, j)) ;
	ASSERT (j >= 0 && j < n) ;
	DEBUG (cnt = 0) ;
	while (CParent [j] == -2)
	{
	    j = FLIP (Flag [j]) ;
	    PRINT2 (("    walk up to "ID" ", j)) ;
	    ASSERT (j >= 0 && j < n) ;
	    PRINT2 ((" CParent "ID"\n", CParent [j])) ;
	    ASSERT (cnt < n) ;
	    DEBUG (cnt++) ;
	}
	cnode = j ;
	ASSERT (cnode >= 0 && cnode < n) ;
	ASSERT (CParent [cnode] >= EMPTY && CParent [cnode] < n) ;
	PRINT2 (("i "ID" is in component with cnode "ID"\n", i, cnode)) ;
	ASSERT (Flag [cnode] == FLIP (cnode)) ;

	/* Mark all nodes along the path from i to cnode as being in the
	 * component whos repnode is cnode.  Perform path compression.  */
	j = FLIP (Flag [i]) ;
	Flag [i] = FLIP (cnode) ;
	DEBUG (cnt = 0) ;
	while (CParent [j] == -2)
	{
	    ASSERT (j >= 0 && j < n) ;
	    jnext = FLIP (Flag [j]) ;
	    PRINT2 (("    "ID" walk "ID" set cnode to "ID"\n", i, j, cnode)) ;
	    ASSERT (cnt < n) ;
	    DEBUG (cnt++) ;
	    Flag [j] = FLIP (cnode) ;
	    j = jnext ;
	}
    }

    /* At this point, all nodes fall into Types 1 or 2, as defined above. */

#ifndef NDEBUG
    for (j = 0 ; j < n ; j++)
    {
	PRINT2 (("j %d CParent %d  ", j, CParent [j])) ;
	if (CParent [j] >= EMPTY && CParent [j] < n)
	{
	    /* case 1: j is a repnode of a component */
	    cnode = j ;
	    PRINT2 ((" a repnode\n")) ;
	}
	else
	{
	    /* case 2: j is not a repnode of a component */
	    cnode = FLIP (Flag [j]) ;
	    PRINT2 ((" repnode is %d\n", cnode)) ;
	    ASSERT (cnode >= 0 && cnode < n) ;
	    ASSERT (CParent [cnode] >= EMPTY && CParent [cnode] < n) ;
	}
	ASSERT (Flag [cnode] == FLIP (cnode)) ;
	/* case 3 no longer holds */
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* free workspace */
    /* ---------------------------------------------------------------------- */

    C->ncol = n ;   /* restore size for memory usage statistics */
    CHOLMOD(free_sparse) (&C, Common) ;
    CHOLMOD(free_sparse) (&B, Common) ;
    CHOLMOD(free) (csize, sizeof (Int), Cew, Common) ;
    CHOLMOD(free) (3*n, sizeof (Int), Work3n, Common) ;

    /* ---------------------------------------------------------------------- */
    /* handle dense nodes */
    /* ---------------------------------------------------------------------- */

    /* The separator tree has nodes with either no children or two or more
     * children - with one exception.  There may exist a single root node with
     * exactly one child, which holds the dense rows/columns of the matrix.
     * Delete this node if it exists. */

    if (ndense > 0)
    {
	ASSERT (CParent [cdense] == EMPTY) ;	/* cdense has no parent */
	/* find the children of cdense */
	nchild = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    if (CParent [j] == cdense)
	    {
		nchild++ ;
		child = j ;
	    }
	}
	if (nchild == 1)
	{
	    /* the cdense node has just one child; merge the two nodes */
	    PRINT1 (("root has one child\n")) ;
	    CParent [cdense] = -2 ;		/* cdense is deleted */
	    CParent [child] = EMPTY ;		/* child becomes a root */
	    for (j = 0 ; j < n ; j++)
	    {
		if (Flag [j] == FLIP (cdense))
		{
		    /* j is a dense node */
		    PRINT1 (("dense %d\n", j)) ;
		    Flag [j] = FLIP (child) ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* postorder the components */
    /* ---------------------------------------------------------------------- */

    DEBUG (for (cnt = 0, j = 0 ; j < n ; j++) if (CParent [j] != -2) cnt++) ;

    /* use Cmember as workspace for Post [ */
    Post = Cmember ;

    /* cholmod_postorder uses Head and Iwork [0..2n].  It does not use Flag,
     * which here holds the mapping of nodes to repnodes.  It ignores all nodes
     * for which CParent [j] < -1, so it operates just on the repnodes. */
    /* workspace: Head (n), Iwork (2*n) */
    ncomponents = CHOLMOD(postorder) (CParent, n, NULL, Post, Common) ;
    ASSERT (cnt == ncomponents) ;

    /* use Iwork [0..n-1] as workspace for Ipost ( */
    Ipost = Iwork ;
    DEBUG (for (j = 0 ; j < n ; j++) Ipost [j] = EMPTY) ;

    /* compute inverse postorder */
    for (c = 0 ; c < ncomponents ; c++)
    {
	cnode = Post [c] ;
	ASSERT (cnode >= 0 && cnode < n) ;
	Ipost [cnode] = c ;
	ASSERT (Head [c] == EMPTY) ;
    }

    /* adjust the parent array */
    /* Iwork [n..2n-1] used for NewParent [ */
    NewParent = Iwork + n ;
    for (c = 0 ; c < ncomponents ; c++)
    {
	parent = CParent [Post [c]] ;
	NewParent [c] = (parent == EMPTY) ? EMPTY : (Ipost [parent]) ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	CParent [c] = NewParent [c] ;
    }
    ASSERT (CHOLMOD(dump_parent) (CParent, ncomponents, "CParent", Common)) ;

    /* Iwork [n..2n-1] no longer needed for NewParent ] */
    /* Cmember no longer needed for Post ] */

#ifndef NDEBUG
    /* count the number of children of each node */
    for (c = 0 ; c < ncomponents ; c++)
    {
	Cmember [c] = 0 ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	if (CParent [c] != EMPTY) Cmember [CParent [c]]++ ;
    }
    for (c = 0 ; c < ncomponents ; c++)
    {
	/* a node is either a leaf, or has 2 or more children */
	ASSERT (Cmember [c] == 0 || Cmember [c] >= 2) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* place each node in its component */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j < n ; j++)
    {
	/* node j is in the cth component, whose repnode is cnode */
	cnode = FLIP (Flag [j]) ;
	PRINT2 (("j "ID"  flag "ID" cnode "ID"\n",
		    j, Flag [j], FLIP (Flag [j]))) ;
	ASSERT (cnode >= 0 && cnode < n) ;
	c = Ipost [cnode] ;
	ASSERT (c >= 0 && c < ncomponents) ;
	Cmember [j] = c ;
    }

    /* Flag no longer needed for the node-to-component mapping */

    /* done using Iwork [0..n-1] as workspace for Ipost ) */

    /* ---------------------------------------------------------------------- */
    /* clear the Flag array */
    /* ---------------------------------------------------------------------- */

    Common->mark = EMPTY ;
    CHOLMOD_CLEAR_FLAG (Common) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* find the permutation */
    /* ---------------------------------------------------------------------- */

    PRINT1 (("nd_camd: %d A->stype %d\n", nd_camd, A->stype)) ;

    if (nd_camd)
    {

	/* ------------------------------------------------------------------ */
	/* apply camd, csymamd, or ccolamd using the Cmember constraints */
	/* ------------------------------------------------------------------ */

	if (A->stype != 0)
	{
	    /* ordering A+A', so fset and fsize are ignored.
	     * Add the upper/lower part to a symmetric lower/upper matrix by
	     * converting to unsymmetric mode
	     * workspace: Iwork (nrow) */
	    B = CHOLMOD(copy) (A, 0, -1, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		PRINT0 (("make symmetric failed\n")) ;
		return (EMPTY) ;
	    }
	    ASSERT ((Int) (B->nrow) == n && (Int) (B->ncol) == n) ;
	    PRINT2 (("nested dissection (2)\n")) ;
	    B->stype = -1 ;
	    if (nd_camd == 2)
	    {
		/* workspace:  Head (nrow+1), Iwork (nrow) if symmetric-upper */
		ok = CHOLMOD(csymamd) (B, Cmember, Perm, Common) ;
	    }
	    else
	    {
		/* workspace: Head (nrow), Iwork (4*nrow) */
		ok = CHOLMOD(camd) (B, NULL, 0, Cmember, Perm, Common) ;
	    }
	    CHOLMOD(free_sparse) (&B, Common) ;
	    if (!ok)
	    {
		/* failed */
		PRINT0 (("camd/csymamd failed\n")) ;
		return (EMPTY) ;
	    }
	}
	else
	{
	    /* ordering A*A' or A(:,f)*A(:,f)' */
	    /* workspace: Iwork (nrow if no fset; MAX(nrow,ncol) if fset) */
	    if (!CHOLMOD(ccolamd) (A, fset, fsize, Cmember, Perm, Common))
	    {
		/* ccolamd failed */
		PRINT2 (("ccolamd failed\n")) ;
		return (EMPTY) ;
	    }
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* natural ordering of each component */
	/* ------------------------------------------------------------------ */

	/* use Iwork [0..n-1] for Next [ */
	Next = Iwork  ;

	/* ------------------------------------------------------------------ */
	/* place the nodes in link lists, one list per component */
	/* ------------------------------------------------------------------ */

	/* do so in reverse order, to preserve original ordering */
	for (j = n-1 ; j >= 0 ; j--)
	{
	    /* node j is in the cth component */
	    c = Cmember [j] ;
	    ASSERT (c >= 0 && c < ncomponents) ;
	    /* place node j in link list for component c */
	    Next [j] = Head [c] ;
	    Head [c] = j ;
	}

	/* ------------------------------------------------------------------ */
	/* order each node in each component */
	/* ------------------------------------------------------------------ */

	k = 0 ;
	for (c = 0 ; c < ncomponents ; c++)
	{
	    for (j = Head [c] ; j != EMPTY ; j = Next [j])
	    {
		Perm [k++] = j ;
	    }
	    Head [c] = EMPTY ;
	}
	ASSERT (k == n) ;

	/* done using Iwork [0..n-1] for Next ] */
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace and return number of components */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    return (ncomponents) ;
}

/* ========================================================================== */
/* === cholmod_collapse_septree ============================================= */
/* ========================================================================== */

/* cholmod_nested_dissection returns the separator tree that was used in the
 * constrained minimum degree algorithm.  Parameter settings (nd_small,
 * nd_oksep, etc) that give a good fill-reducing ordering may give too fine of
 * a separator tree for other uses (parallelism, multi-level LPDASA, etc).  This
 * function takes as input the separator tree computed by
 * cholmod_nested_dissection, and collapses selected subtrees into single
 * nodes.  A subtree is collapsed if its root node (the separator) is large
 * compared to the total number of nodes in the subtree, or if the subtree is
 * small.  Note that the separator tree may actually be a forest.
 *
 * nd_oksep and nd_small act just like the ordering parameters in Common.
 * Returns the new number of nodes in the separator tree.
 */

SuiteSparse_long CHOLMOD(collapse_septree)
(
    /* ---- input ---- */
    size_t n,		/* # of nodes in the graph */
    size_t ncomponents,	/* # of nodes in the separator tree (must be <= n) */
    double nd_oksep,    /* collapse if #sep >= nd_oksep * #nodes in subtree */
    size_t nd_small,    /* collapse if #nodes in subtree < nd_small */
    /* ---- in/out --- */
    Int *CParent,	/* size ncomponents; from cholmod_nested_dissection */
    Int *Cmember,	/* size n; from cholmod_nested_dissection */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *First, *Count, *Csubtree, *W, *Map ;
    Int c, j, k, nc, sepsize, total_weight, parent, nc_new, first ;
    int collapse = FALSE, ok = TRUE ;
    size_t s ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (CParent, EMPTY) ;
    RETURN_IF_NULL (Cmember, EMPTY) ;
    if (n < ncomponents)
    {
	ERROR (CHOLMOD_INVALID, "invalid separator tree") ;
	return (EMPTY) ;
    }
    Common->status = CHOLMOD_OK ;
    nc = ncomponents ;
    if (n <= 1 || ncomponents <= 1)
    {
	/* no change; tree is one node already */
	return (nc) ;
    }

    nd_oksep = MAX (0, nd_oksep) ;
    nd_oksep = MIN (1, nd_oksep) ;
    nd_small = MAX (4, nd_small) ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    /* s = 3*ncomponents */
    s = CHOLMOD(mult_size_t) (ncomponents, 3, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (EMPTY) ;
    }
    CHOLMOD(allocate_work) (0, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (EMPTY) ;
    }
    W = Common->Iwork ;
    Count    = W ; W += ncomponents ;	    /* size ncomponents */
    Csubtree = W ; W += ncomponents ;	    /* size ncomponents */
    First    = W ; W += ncomponents ;	    /* size ncomponents */

    /* ---------------------------------------------------------------------- */
    /* find the first descendant of each node of the separator tree */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < nc ; c++)
    {
	First [c] = EMPTY ;
    }
    for (k = 0 ; k < nc ; k++)
    {
	for (c = k ; c != EMPTY && First [c] == -1 ; c = CParent [c])
	{
	    ASSERT (c >= 0 && c < nc) ;
	    First [c] = k ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* find the number of nodes of the graph in each node of the tree */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < nc ; c++)
    {
	Count [c] = 0 ;
    }
    for (j = 0 ; j < (Int) n ; j++)
    {
	ASSERT (Cmember [j] >= 0 && Cmember [j] < nc) ;
	Count [Cmember [j]]++ ;
    }

    /* ---------------------------------------------------------------------- */
    /* find the number of nodes in each subtree */
    /* ---------------------------------------------------------------------- */

    for (c = 0 ; c < nc ; c++)
    {
	/* each subtree includes its root */
	Csubtree [c] = Count [c] ;
	PRINT1 ((ID" size "ID" parent "ID" first "ID"\n",
	    c, Count [c], CParent [c], First [c])) ;
    }

    for (c = 0 ; c < nc ; c++)
    {
	/* add the subtree of the child, c, into the count of its parent */
	parent = CParent [c] ;
	ASSERT (parent >= EMPTY && parent < nc) ;
	if (parent != EMPTY)
	{
	    Csubtree [parent] += Csubtree [c] ;
	}
    }

#ifndef NDEBUG
    /* the sum of the roots should be n */
    j = 0 ;
    for (c = 0 ; c < nc ; c++) if (CParent [c] == EMPTY) j += Csubtree [c] ;
    ASSERT (j == (Int) n) ;
#endif

    /* ---------------------------------------------------------------------- */
    /* find subtrees to collapse */
    /* ---------------------------------------------------------------------- */

    /* consider all nodes in reverse post-order */
    for (c = nc-1 ; c >= 0 ; c--)
    {
	/* consider the subtree rooted at node c */
	sepsize = Count [c] ;
	total_weight = Csubtree [c] ;
	PRINT1 (("Node "ID" sepsize "ID" subtree "ID" ratio %g\n", c, sepsize,
	    total_weight, ((double) sepsize)/((double) total_weight))) ;
	first = First [c] ;
	if (first < c &&    /* c must not be a leaf */
	   (sepsize > nd_oksep * total_weight || total_weight < (int) nd_small))
	{
	    /* this separator is too large, or the subtree is too small.
	     * collapse the tree, by converting the entire subtree rooted at
	     * c into a single node.  The subtree consists of all nodes from
	     * First[c] to the root c.  Flag all nodes from First[c] to c-1
	     * as dead.
	     */
	    collapse = TRUE ;
	    for (k = first ; k < c ; k++)
	    {
		CParent [k] = -2 ;
		PRINT1 (("   collapse node "ID"\n", k)) ;
	    }
	    /* continue at the next node, first-1 */
	    c = first ;
	}
    }

    PRINT1 (("collapse: %d\n", collapse)) ;

    /* ---------------------------------------------------------------------- */
    /* compress the tree */
    /* ---------------------------------------------------------------------- */

    Map = Count ;	/* Count no longer needed */

    nc_new = nc ;
    if (collapse)
    {
	nc_new = 0 ;
	for (c = 0 ; c < nc ; c++)
	{
	    Map [c] = nc_new ;
	    if (CParent [c] >= EMPTY)
	    {
		/* node c is alive, and becomes node Map[c] in the new tree.
		 * Increment nc_new for the next node c. */
		nc_new++ ;
	    }
	}
	PRINT1 (("Collapse the tree from "ID" to "ID" nodes\n", nc, nc_new)) ;
	ASSERT (nc_new > 0) ;
	for (c = 0 ; c < nc ; c++)
	{
	    parent = CParent [c] ;
	    if (parent >= EMPTY)
	    {
		/* node c is alive */
		CParent [Map [c]] = (parent == EMPTY) ? EMPTY : Map [parent] ;
	    }
	}
	for (j = 0 ; j < (Int) n ; j++)
	{
	    PRINT1 (("j "ID" Cmember[j] "ID" Map[Cmember[j]] "ID"\n",
		j, Cmember [j], Map [Cmember [j]])) ;
	    Cmember [j] = Map [Cmember [j]] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return new size of separator tree */
    /* ---------------------------------------------------------------------- */

    return (nc_new) ;
}
#endif
