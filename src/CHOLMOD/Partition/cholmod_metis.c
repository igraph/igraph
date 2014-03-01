/* ========================================================================== */
/* === Partition/cholmod_metis ============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Partition Module.
 * Copyright (C) 2005-2006, Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Partition Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* CHOLMOD interface to the METIS package (Version 4.0.1):
 *
 * cholmod_metis_bisector:
 *
 *	Wrapper for METIS_NodeComputeSeparator.  Finds a set of nodes that
 *	partitions the graph into two parts.
 *
 * cholmod_metis:
 *
 *	Wrapper for METIS_NodeND, METIS's own nested dissection algorithm.
 *	Typically faster than cholmod_nested_dissection, mostly because it
 *	uses minimum degree on just the leaves of the separator tree, rather
 *	than the whole matrix.
 *
 * Note that METIS does not return an error if it runs out of memory.  Instead,
 * it terminates the program.  This interface attempts to avoid that problem
 * by preallocating space that should be large enough for any memory allocations
 * within METIS, and then freeing that space, just before the call to METIS.
 * While this is not guaranteed to work, it is very unlikely to fail.  If you
 * encounter this problem, increase Common->metis_memory.  If you don't mind
 * having your program terminated, set Common->metis_memory to zero (a value of
 * 2.0 is usually safe).  Several other METIS workarounds are made in the
 * routines in this file.  See the description of metis_memory_ok, just below,
 * for more details.
 *
 * FUTURE WORK: interfaces to other partitioners (CHACO, SCOTCH, JOSTLE, ... )
 *
 * workspace: several size-nz and size-n temporary arrays.  Uses no workspace
 * in Common.
 *
 * Supports any xtype (pattern, real, complex, or zomplex).
 */

#ifndef NPARTITION

#include "cholmod_internal.h"
#undef ASSERT

#include "metis.h"
/* METIS has its own ASSERT that it reveals to the user, so remove it here: */
#undef ASSERT

/* and redefine it back again */
#ifndef NDEBUG
#define ASSERT(expression) (assert (expression))
#else
#define ASSERT(expression)
#endif

#include "cholmod_partition.h"
#include "cholmod_cholesky.h"


/* ========================================================================== */
/* === dumpgraph ============================================================ */
/* ========================================================================== */

/* For dumping the input graph to METIS_NodeND, to check with METIS's onmetis
 * and graphchk programs.  For debugging only.  To use, uncomment this #define:
#define DUMP_GRAPH
 */

#ifdef DUMP_GRAPH
#include <stdio.h>
/* After dumping the graph with this routine, run "onmetis metisgraph" */
static void dumpgraph (idxtype *Mp, idxtype *Mi, SuiteSparse_long n,
    cholmod_common *Common)
{
    SuiteSparse_long i, j, p, nz ;
    FILE *f ;
    nz = Mp [n] ;
    printf ("Dumping METIS graph n %ld nz %ld\n", n, nz) ;    /* DUMP_GRAPH */
    f = fopen ("metisgraph", "w") ;
    if (f == NULL)
    {
	ERROR (-99, "cannot open metisgraph") ;
	return ;
    }
    fprintf (f, "%ld %ld\n", n, nz/2) ;			    /* DUMP_GRAPH */
    for (j = 0 ; j < n ; j++)
    {
	for (p = Mp [j] ; p < Mp [j+1] ; p++)
	{
	    i = Mi [p] ;
	    fprintf (f, " %ld", i+1) ;			    /* DUMP_GRAPH */
	}
	fprintf (f, "\n") ;				    /* DUMP_GRAPH */
    }
    fclose (f) ;
}
#endif


/* ========================================================================== */
/* === metis_memory_ok ====================================================== */
/* ========================================================================== */

/* METIS_NodeND and METIS_NodeComputeSeparator will terminate your program it
 * they run out of memory.  In an attempt to workaround METIS' behavior, this
 * routine allocates a single block of memory of size equal to an observed
 * upper bound on METIS' memory usage.  It then immediately deallocates the
 * block.  If the allocation fails, METIS is not called.
 *
 * Median memory usage for a graph with n nodes and nz edges (counting each
 * edge twice, or both upper and lower triangular parts of a matrix) is
 * 4*nz + 40*n + 4096 integers.  A "typical" upper bound is 10*nz + 50*n + 4096
 * integers.  Nearly all matrices tested fit within that upper bound, with the
 * exception two in the UF sparse matrix collection: Schenk_IBMNA/c-64 and
 * Gupta/gupta2.  The latter exceeds the "upper bound" by a factor of just less
 * than 2.
 *
 * If you do not mind having your program terminated if it runs out of memory,
 * set Common->metis_memory to zero.  Its default value is 2, which allows for
 * some memory fragmentation, and also accounts for the Gupta/gupta2 matrix.
 *
 * An alternative, if CHOLMOD is used in MATLAB, is to use a version of METIS
 * (4.0.2, perhaps) proposed to George Karypis.  This version uses function
 * pointer for malloc and free.  They can be set to mxMalloc and mxFree
 * (see sputil_config in MATLAB/sputil.c).  On Linux, with gcc, you must also
 * compile CHOLMOD, METIS, AMD, COLAMD, and CCOLAMD with the -fexceptions
 * compiler flag.  With this configuration, mxMalloc safely aborts the
 * mexFunction, frees all memory allocted by METIS, and safely returns to
 * MATLAB.  You may then set Common->metis_memory = 0.
 */

#define GUESS(nz,n) (10 * (nz) + 50 * (n) + 4096)

static int metis_memory_ok
(
    Int n,
    Int nz,
    cholmod_common *Common
)
{
    double s ;
    void *p ;
    size_t metis_guard ;

    if (Common->metis_memory <= 0)
    {
	/* do not prevent METIS from running out of memory */
	return (TRUE) ;
    }

    n  = MAX (1, n) ;
    nz = MAX (0, nz) ;

    /* compute in double, to avoid integer overflow */
    s = GUESS ((double) nz, (double) n) ;
    s *= Common->metis_memory ;

    if (s * sizeof (idxtype) >= ((double) Size_max))
    {
	/* don't even attempt to malloc such a large block */
	return (FALSE) ;
    }

    /* recompute in size_t */
    metis_guard = GUESS ((size_t) nz, (size_t) n) ;
    metis_guard *= Common->metis_memory ;

    /* attempt to malloc the block */
    p = CHOLMOD(malloc) (metis_guard, sizeof (idxtype), Common) ;
    if (p == NULL)
    {
	/* failure - return out-of-memory condition */
	return (FALSE) ;
    }

    /* success - free the block */
    CHOLMOD(free) (metis_guard, sizeof (idxtype), p, Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_metis_bisector =============================================== */
/* ========================================================================== */

/* Finds a set of nodes that bisects the graph of A or AA' (direct interface
 * to METIS_NodeComputeSeparator).
 *
 * The input matrix A must be square, symmetric (with both upper and lower
 * parts present) and with no diagonal entries.  These conditions are NOT
 * checked.
 */

SuiteSparse_long CHOLMOD(metis_bisector)	/* returns separator size */
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to bisect */
    Int *Anw,		/* size A->nrow, node weights */
    Int *Aew,		/* size nz, edge weights */
    /* ---- output --- */
    Int *Partition,	/* size A->nrow */
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Ap, *Ai ;
    idxtype *Mp, *Mi, *Mnw, *Mew, *Mpart ;
    Int n, nleft, nright, j, p, csep, total_weight, lightest, nz ;
    int Opt [8], nn, csp ;
    size_t n1 ;
    DEBUG (Int nsep) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_NULL (Anw, EMPTY) ;
    RETURN_IF_NULL (Aew, EMPTY) ;
    RETURN_IF_NULL (Partition, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    if (A->stype || A->nrow != A->ncol)
    {
	/* A must be square, with both upper and lower parts present */
	ERROR (CHOLMOD_INVALID, "matrix must be square, symmetric,"
		" and with both upper/lower parts present") ;
	return (EMPTY) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* quick return */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    if (n == 0)
    {
	return (0) ;
    }
    n1 = ((size_t) n) + 1 ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    nz = Ap [n] ;

    /* ---------------------------------------------------------------------- */
    /* METIS does not have a 64-bit integer version */
    /* ---------------------------------------------------------------------- */

#ifdef LONG
    if (sizeof (Int) > sizeof (idxtype) && MAX (n,nz) > INT_MAX / sizeof (int))
    {
	/* CHOLMOD's matrix is too large for METIS */
	return (EMPTY) ;
    }
#endif

    /* ---------------------------------------------------------------------- */
    /* set default options */
    /* ---------------------------------------------------------------------- */

    Opt [0] = 0 ;	/* use defaults */
    Opt [1] = 3 ;	/* matching type */
    Opt [2] = 1 ;	/* init. partitioning algo*/
    Opt [3] = 2 ;	/* refinement algorithm */
    Opt [4] = 0 ;	/* no debug */
    Opt [5] = 0 ;	/* unused */
    Opt [6] = 0 ;	/* unused */
    Opt [7] = -1 ;	/* random seed */

    DEBUG (for (j = 0 ; j < n ; j++) ASSERT (Anw [j] > 0)) ;

    /* ---------------------------------------------------------------------- */
    /* copy Int to METIS idxtype, if necessary */
    /* ---------------------------------------------------------------------- */

    DEBUG (for (j = 0 ; j < nz ; j++) ASSERT (Aew [j] > 0)) ;
    if (sizeof (Int) == sizeof (idxtype))
    {
	/* this is the typical case */
	Mi    = (idxtype *) Ai ;
	Mew   = (idxtype *) Aew ;
	Mp    = (idxtype *) Ap ;
	Mnw   = (idxtype *) Anw ;
	Mpart = (idxtype *) Partition ;
    }
    else
    {
	/* idxtype and Int differ; copy the graph into the METIS idxtype */
	Mi    = CHOLMOD(malloc) (nz, sizeof (idxtype), Common) ;
	Mew   = CHOLMOD(malloc) (nz, sizeof (idxtype), Common) ;
	Mp    = CHOLMOD(malloc) (n1, sizeof (idxtype), Common) ;
	Mnw   = CHOLMOD(malloc) (n,  sizeof (idxtype), Common) ;
	Mpart = CHOLMOD(malloc) (n,  sizeof (idxtype), Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    CHOLMOD(free) (nz, sizeof (idxtype), Mi,    Common) ;
	    CHOLMOD(free) (nz, sizeof (idxtype), Mew,   Common) ;
	    CHOLMOD(free) (n1, sizeof (idxtype), Mp,    Common) ;
	    CHOLMOD(free) (n,  sizeof (idxtype), Mnw,   Common) ;
	    CHOLMOD(free) (n,  sizeof (idxtype), Mpart, Common) ;
	    return (EMPTY) ;
	}
	for (p = 0 ; p < nz ; p++)
	{
	    Mi [p] = Ai [p] ;
	}
	for (p = 0 ; p < nz ; p++)
	{
	    Mew [p] = Aew [p] ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Mp [j] = Ap [j] ;
	}
	for (j = 0 ; j <  n ; j++)
	{
	    Mnw [j] = Anw [j] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* METIS workaround: try to ensure METIS doesn't run out of memory */
    /* ---------------------------------------------------------------------- */

    if (!metis_memory_ok (n, nz, Common))
    {
	/* METIS might ask for too much memory and thus terminate the program */
	if (sizeof (Int) != sizeof (idxtype))
	{
	    CHOLMOD(free) (nz, sizeof (idxtype), Mi,    Common) ;
	    CHOLMOD(free) (nz, sizeof (idxtype), Mew,   Common) ;
	    CHOLMOD(free) (n1, sizeof (idxtype), Mp,    Common) ;
	    CHOLMOD(free) (n,  sizeof (idxtype), Mnw,   Common) ;
	    CHOLMOD(free) (n,  sizeof (idxtype), Mpart, Common) ;
	}
	return (EMPTY) ;
    }

    /* ---------------------------------------------------------------------- */
    /* partition the graph */
    /* ---------------------------------------------------------------------- */

#ifndef NDEBUG
    PRINT1 (("Metis graph, n = "ID"\n", n)) ;
    for (j = 0 ; j < n ; j++)
    {
	Int ppp ;
	PRINT2 (("M(:,"ID") node weight "ID"\n", j, (Int) Mnw [j])) ;
	ASSERT (Mnw [j] > 0) ;
	for (ppp = Mp [j] ; ppp < Mp [j+1] ; ppp++)
	{
	    PRINT3 ((" "ID" : "ID"\n", (Int) Mi [ppp], (Int) Mew [ppp])) ;
	    ASSERT (Mi [ppp] != j) ;
	    ASSERT (Mew [ppp] > 0) ;
	}
    }
#endif

    nn = n ;
    METIS_NodeComputeSeparator (&nn, Mp, Mi, Mnw, Mew, Opt, &csp, Mpart) ;
    n = nn ;
    csep = csp ;

    PRINT1 (("METIS csep "ID"\n", csep)) ;

    /* ---------------------------------------------------------------------- */
    /* copy the results back from idxtype, if required */
    /* ---------------------------------------------------------------------- */

    if (sizeof (Int) != sizeof (idxtype))
    {
	for (j = 0 ; j < n ; j++)
	{
	    Partition [j] = Mpart [j] ;
	}
	CHOLMOD(free) (nz, sizeof (idxtype), Mi,    Common) ;
	CHOLMOD(free) (nz, sizeof (idxtype), Mew,   Common) ;
	CHOLMOD(free) (n1, sizeof (idxtype), Mp,    Common) ;
	CHOLMOD(free) (n,  sizeof (idxtype), Mnw,   Common) ;
	CHOLMOD(free) (n,  sizeof (idxtype), Mpart, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* ensure a reasonable separator */
    /* ---------------------------------------------------------------------- */

    /* METIS can return a valid separator with no nodes in (for example) the
     * left part.  In this case, there really is no separator.  CHOLMOD
     * prefers, in this case, for all nodes to be in the separator (and both
     * left and right parts to be empty).  Also, if the graph is unconnected,
     * METIS can return a valid empty separator.  CHOLMOD prefers at least one
     * node in the separator.  Note that cholmod_nested_dissection only calls
     * this routine on connected components, but cholmod_bisect can call this
     * routine for any graph. */

    if (csep == 0)
    {
	/* The separator is empty, select lightest node as separator.  If
	 * ties, select the highest numbered node. */
	lightest = 0 ;
	for (j = 0 ; j < n ; j++)
	{
	    if (Anw [j] <= Anw [lightest])
	    {
		lightest = j ;
	    }
	}
	PRINT1 (("Force "ID" as sep\n", lightest)) ;
	Partition [lightest] = 2 ;
	csep = Anw [lightest] ;
    }

    /* determine the node weights in the left and right part of the graph */
    nleft = 0 ;
    nright = 0 ;
    DEBUG (nsep = 0) ;
    for (j = 0 ; j < n ; j++)
    {
	PRINT1 (("Partition ["ID"] = "ID"\n", j, Partition [j])) ;
	if (Partition [j] == 0)
	{
	    nleft += Anw [j] ;
	}
	else if (Partition [j] == 1)
	{
	    nright += Anw [j] ;
	}
#ifndef NDEBUG
	else
	{
	    ASSERT (Partition [j] == 2) ;
	    nsep += Anw [j] ;
	}
#endif
    }
    ASSERT (csep == nsep) ;

    total_weight = nleft + nright + csep ;

    if (csep < total_weight)
    {
	/* The separator is less than the whole graph.  Make sure the left and
	 * right parts are either both empty or both non-empty. */
	PRINT1 (("nleft "ID" nright "ID" csep "ID" tot "ID"\n",
		nleft, nright, csep, total_weight)) ;
	ASSERT (nleft + nright + csep == total_weight) ;
	ASSERT (nleft > 0 || nright > 0) ;
	if ((nleft == 0 && nright > 0) || (nleft > 0 && nright == 0))
	{
	    /* left or right is empty; put all nodes in the separator */
	    PRINT1 (("Force all in sep\n")) ;
	    csep = total_weight ;
	    for (j = 0 ; j < n ; j++)
	    {
		Partition [j] = 2 ;
	    }
	}
    }

    ASSERT (CHOLMOD(dump_partition) (n, Ap, Ai, Anw, Partition, csep, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* return the sum of the weights of nodes in the separator */
    /* ---------------------------------------------------------------------- */

    return (csep) ;
}


/* ========================================================================== */
/* === cholmod_metis ======================================================== */
/* ========================================================================== */

/* CHOLMOD wrapper for the METIS_NodeND ordering routine.  Creates A+A',
 * A*A' or A(:,f)*A(:,f)' and then calls METIS_NodeND on the resulting graph.
 * This routine is comparable to cholmod_nested_dissection, except that it
 * calls METIS_NodeND directly, and it does not return the separator tree.
 *
 * workspace:  Flag (nrow), Iwork (4*n+uncol)
 *	Allocates a temporary matrix B=A*A' or B=A.
 */

int CHOLMOD(metis)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to order */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int postorder,	/* if TRUE, follow with etree or coletree postorder */
    /* ---- output --- */
    Int *Perm,		/* size A->nrow, output permutation */
    /* --------------- */
    cholmod_common *Common
)
{
    double d ;
    Int *Iperm, *Iwork, *Bp, *Bi ;
    idxtype *Mp, *Mi, *Mperm, *Miperm ;
    cholmod_sparse *B ;
    Int i, j, n, nz, p, identity, uncol ;
    int Opt [8], nn, zero = 0 ;
    size_t n1, s ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (Perm, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* quick return */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;
    if (n == 0)
    {
	return (TRUE) ;
    }
    n1 = ((size_t) n) + 1 ;

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
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (n, s, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* convert the matrix to adjacency list form */
    /* ---------------------------------------------------------------------- */

    /* The input graph for METIS must be symmetric, with both upper and lower
     * parts present, and with no diagonal entries.  The columns need not be
     * sorted.
     * B = A+A', A*A', or A(:,f)*A(:,f)', upper and lower parts present */
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
    ASSERT (CHOLMOD(dump_sparse) (B, "B for NodeND", Common) >= 0) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (FALSE) ;
    }
    ASSERT (B->nrow == A->nrow) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Iperm = Iwork ;		/* size n (i/i/l) */

    Bp = B->p ;
    Bi = B->i ;
    nz = Bp [n] ;

    /* ---------------------------------------------------------------------- */
    /* METIS does not have a SuiteSparse_long integer version */
    /* ---------------------------------------------------------------------- */

#ifdef LONG
    if (sizeof (Int) > sizeof (idxtype) && MAX (n,nz) > INT_MAX / sizeof (int))
    {
	/* CHOLMOD's matrix is too large for METIS */
	CHOLMOD(free_sparse) (&B, Common) ;
	return (FALSE) ;
    }
#endif

    /* B does not include the diagonal, and both upper and lower parts.
     * Common->anz includes the diagonal, and just the lower part of B */
    Common->anz = nz / 2 + n ;

    /* ---------------------------------------------------------------------- */
    /* set control parameters for METIS_NodeND */
    /* ---------------------------------------------------------------------- */

    Opt [0] = 0 ;	/* use defaults */
    Opt [1] = 3 ;	/* matching type */
    Opt [2] = 1 ;	/* init. partitioning algo*/
    Opt [3] = 2 ;	/* refinement algorithm */
    Opt [4] = 0 ;	/* no debug */
    Opt [5] = 1 ;	/* initial compression */
    Opt [6] = 0 ;	/* no dense node removal */
    Opt [7] = 1 ;	/* number of separators @ each step */

    /* ---------------------------------------------------------------------- */
    /* allocate the METIS input arrays, if needed */
    /* ---------------------------------------------------------------------- */

    if (sizeof (Int) == sizeof (idxtype))
    {
	/* This is the typical case. */
	Miperm = (idxtype *) Iperm ;
	Mperm  = (idxtype *) Perm ;
	Mp     = (idxtype *) Bp ;
	Mi     = (idxtype *) Bi ;
    }
    else
    {
	/* allocate graph for METIS only if Int and idxtype differ */
	Miperm = CHOLMOD(malloc) (n,  sizeof (idxtype), Common) ;
	Mperm  = CHOLMOD(malloc) (n,  sizeof (idxtype), Common) ;
	Mp     = CHOLMOD(malloc) (n1, sizeof (idxtype), Common) ;
	Mi     = CHOLMOD(malloc) (nz, sizeof (idxtype), Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    CHOLMOD(free_sparse) (&B, Common) ;
	    CHOLMOD(free) (n,  sizeof (idxtype), Miperm, Common) ;
	    CHOLMOD(free) (n,  sizeof (idxtype), Mperm, Common) ;
	    CHOLMOD(free) (n1, sizeof (idxtype), Mp, Common) ;
	    CHOLMOD(free) (nz, sizeof (idxtype), Mi, Common) ;
	    return (FALSE) ;
	}
	for (j = 0 ; j <= n ; j++)
	{
	    Mp [j] = Bp [j] ;
	}
	for (p = 0 ; p < nz ; p++)
	{
	    Mi [p] = Bi [p] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* METIS workarounds */
    /* ---------------------------------------------------------------------- */

    identity = FALSE ;
    if (nz == 0)
    {
	/* The matrix has no off-diagonal entries.  METIS_NodeND fails in this
	 * case, so avoid using it.  The best permutation is identity anyway,
	 * so this is an easy fix. */
	identity = TRUE ;
	PRINT1 (("METIS:: no nz\n")) ;
    }
    else if (Common->metis_nswitch > 0)
    {
	/* METIS_NodeND in METIS 4.0.1 gives a seg fault with one matrix of
	 * order n = 3005 and nz = 6,036,025, including the diagonal entries.
	 * The workaround is to return the identity permutation instead of using
	 * METIS for matrices of dimension 3000 or more and with density of 66%
	 * or more - admittedly an uncertain fix, but such matrices are so dense
	 * that any reasonable ordering will do, even identity (n^2 is only 50%
	 * higher than nz in this case).  CHOLMOD's nested dissection method
	 * (cholmod_nested_dissection) has no problems with the same matrix,
	 * even though it too uses METIS_NodeComputeSeparator.  The matrix is
	 * derived from LPnetlib/lpi_cplex1 in the UF sparse matrix collection.
	 * If C is the lpi_cplex matrix (of order 3005-by-5224), A = (C*C')^2
	 * results in the seg fault.  The seg fault also occurs in the stand-
	 * alone onmetis program that comes with METIS.  If a future version of
	 * METIS fixes this problem, then set Common->metis_nswitch to zero.
	 */
	d = ((double) nz) / (((double) n) * ((double) n)) ;
	if (n > (Int) (Common->metis_nswitch) && d > Common->metis_dswitch)
	{
	    identity = TRUE ;
	    PRINT1 (("METIS:: nswitch/dswitch activated\n")) ;
	}
    }

    if (!identity && !metis_memory_ok (n, nz, Common))
    {
	/* METIS might ask for too much memory and thus terminate the program */
	identity = TRUE ;
    }

    /* ---------------------------------------------------------------------- */
    /* find the permutation */
    /* ---------------------------------------------------------------------- */

    if (identity)
    {
	/* no need to do the postorder */
	postorder = FALSE ;
	for (i = 0 ; i < n ; i++)
	{
	    Mperm [i] = i ;
	}
    }
    else
    {
#ifdef DUMP_GRAPH
	/* DUMP_GRAPH */ printf ("Calling METIS_NodeND n "ID" nz "ID""
	"density %g\n", n, nz, ((double) nz) / (((double) n) * ((double) n)));
	dumpgraph (Mp, Mi, n, Common) ;
#endif

	nn = n ;
	METIS_NodeND (&nn, Mp, Mi, &zero, Opt, Mperm, Miperm) ;
	n = nn ;

	PRINT0 (("METIS_NodeND done\n")) ;
    }

    /* ---------------------------------------------------------------------- */
    /* free the METIS input arrays */
    /* ---------------------------------------------------------------------- */

    if (sizeof (Int) != sizeof (idxtype))
    {
	for (i = 0 ; i < n ; i++)
	{
	    Perm [i] = (Int) (Mperm [i]) ;
	}
	CHOLMOD(free) (n,   sizeof (idxtype), Miperm, Common) ;
	CHOLMOD(free) (n,   sizeof (idxtype), Mperm, Common) ;
	CHOLMOD(free) (n+1, sizeof (idxtype), Mp, Common) ;
	CHOLMOD(free) (nz,  sizeof (idxtype), Mi, Common) ;
    }

    CHOLMOD(free_sparse) (&B, Common) ;

    /* ---------------------------------------------------------------------- */
    /* etree or column-etree postordering, using the Cholesky Module */
    /* ---------------------------------------------------------------------- */

    if (postorder)
    {
	Int *Parent, *Post, *NewPerm ;
	Int k ;

	Parent = Iwork + 2*((size_t) n) + uncol ;   /* size n = nrow */
	Post   = Parent + n ;			    /* size n */

	/* workspace: Iwork (2*nrow+uncol), Flag (nrow), Head (nrow+1) */
	CHOLMOD(analyze_ordering) (A, CHOLMOD_METIS, Perm, fset, fsize,
		Parent, Post, NULL, NULL, NULL, Common) ;
	if (Common->status == CHOLMOD_OK)
	{
	    /* combine the METIS permutation with its postordering */
	    NewPerm = Parent ;	    /* use Parent as workspace */
	    for (k = 0 ; k < n ; k++)
	    {
		NewPerm [k] = Perm [Post [k]] ;
	    }
	    for (k = 0 ; k < n ; k++)
	    {
		Perm [k] = NewPerm [k] ;
	    }
	}
    }

    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;
    PRINT1 (("cholmod_metis done\n")) ;
    return (Common->status == CHOLMOD_OK) ;
}
#endif
