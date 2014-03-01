/* ========================================================================== */
/* === Supernodal/cholmod_super_symbolic ==================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Supernodal Module. Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Supernodal Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* Supernodal symbolic analysis of the LL' factorization of A, A*A',
 * A(:,f)*A(:,f)'.
 *
 * This routine must be preceded by a simplicial symbolic analysis
 * (cholmod_rowcolcounts).  See cholmod_analyze.c for an example of how to use
 * this routine.
 *
 * The user need not call this directly; cholmod_analyze is a "simple" wrapper
 * for this routine.
 *
 * Symmetric case:
 *
 *	A is stored in column form, with entries stored in the upper triangular
 *	part.  Entries in the lower triangular part are ignored.
 *
 * Unsymmetric case:
 *
 *	A is stored in column form.  If F is equal to the transpose of A, then
 *	A*A' is analyzed.  F can include a subset of the columns of A
 *	(F=A(:,f)'), in which case F*F' is analyzed.
 *
 * Requires Parent and L->ColCount to be defined on input; these are the
 * simplicial Parent and ColCount arrays as computed by cholmod_rowcolcounts.
 * Does not use L->Perm; the input matrices A and F must already be properly
 * permuted.  Allocates and computes the supernodal pattern of L (L->super,
 * L->pi, L->px, and L->s).  Does not allocate the real part (L->x).
 *
 * Supports any xtype (pattern, real, complex, or zomplex).
 */

#ifndef NSUPERNODAL

#include "cholmod_internal.h"
#include "cholmod_supernodal.h"


/* ========================================================================== */
/* === subtree ============================================================== */
/* ========================================================================== */

/* In the symmetric case, traverse the kth row subtree from the nonzeros in
 * A (0:k1-1,k) and add the new entries found to the pattern of the kth row
 * of L.  The current supernode s contains the diagonal block k1:k2-1, so it
 * can be skipped.
 *
 * In the unsymmetric case, the nonzero pattern of A*F is computed one column
 * at a time (thus, the total time spent in this function is bounded below by
 * the time taken to multiply A*F, which can be high if A is tall and thin).
 * The kth column is A*F(:,k), or the set union of all columns A(:,j) for which
 * F(j,k) is nonzero.  This routine is called once for each entry j.  Only the
 * upper triangular part is needed, so only A (0:k1-1,j) is accessed, where
 * k1:k2-1 are the columns of the current supernode s (k is in the range k1 to
 * k2-1).
 *
 * If A is sorted, then the total time taken by this function is proportional
 * to the number of nonzeros in the strictly block upper triangular part of A,
 * plus the number of entries in the strictly block lower triangular part of
 * the supernodal part of L.  This excludes entries in the diagonal blocks
 * corresponding to the columns in each supernode.  That is, if k1:k2-1 are
 * in a single supernode, then only A (0:k1-1,k1:k2-1) are accessed.
 *
 * For the unsymmetric case, only the strictly block upper triangular part
 * of A*F is constructed.
 *
 * Only adds column indices corresponding to the leading columns of each
 * relaxed supernode.
 */

static void subtree
(
    /* inputs, not modified: */
    Int j,		/* j = k for symmetric case */
    Int k,
    Int Ap [ ],
    Int Ai [ ],
    Int Anz [ ],
    Int SuperMap [ ],
    Int Sparent [ ],
    Int mark,
    Int sorted,         /* true if the columns of A are sorted */
    Int k1,             /* only consider A (0:k1-1,k) */

    /* input/output: */
    Int Flag [ ],
    Int Ls [ ],
    Int Lpi2 [ ]
)
{
    Int p, pend, i, si ;
    p = Ap [j] ;
    pend = (Anz == NULL) ? (Ap [j+1]) : (p + Anz [j]) ;

    for ( ; p < pend ; p++)
    {
	i = Ai [p] ;
	if (i < k1)
	{
	    /* (i,k) is an entry in the upper triangular part of A or A*F'.
	     * symmetric case:   A(i,k) is nonzero (j=k).
	     * unsymmetric case: A(i,j) and F(j,k) are both nonzero.
	     *
	     * Column i is in supernode si = SuperMap [i].  Follow path from si
	     * to root of supernodal etree, stopping at the first flagged
	     * supernode.  The root of the row subtree is supernode SuperMap[k],
	     * which is flagged already. This traversal will stop there, or it
	     * might stop earlier if supernodes have been flagged by previous
	     * calls to this routine for the same k. */
	    for (si = SuperMap [i] ; Flag [si] < mark ; si = Sparent [si])
	    {
		ASSERT (si <= SuperMap [k]) ;
		Ls [Lpi2 [si]++] = k ;
		Flag [si] = mark ;
	    }
	}
        else if (sorted)
        {
            break ;
        }
    }
}


/* clear workspace used by cholmod_super_symbolic */
#define FREE_WORKSPACE \
{ \
    /* CHOLMOD(clear_flag) (Common) ; */ \
    CHOLMOD_CLEAR_FLAG (Common) ; \
    for (k = 0 ; k <= nfsuper ; k++) \
    { \
	Head [k] = EMPTY ; \
    } \
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ; \
} \


/* ========================================================================== */
/* === cholmod_super_symbolic2 ============================================== */
/* ========================================================================== */

/* Analyze for supernodal Cholesky or multifrontal QR.  CHOLMOD itself always
 * analyzes for supernodal Cholesky, of course.  The "for_cholesky = TRUE"
 * option is used by SuiteSparseQR only. */

int CHOLMOD(super_symbolic2)
(
    /* ---- input ---- */
    int for_cholesky,   /* Cholesky if TRUE, QR if FALSE */
    cholmod_sparse *A,	/* matrix to analyze */
    cholmod_sparse *F,	/* F = A' or A(:,f)' */
    Int *Parent,	/* elimination tree */
    /* ---- in/out --- */
    cholmod_factor *L,	/* simplicial symbolic on input,
			 * supernodal symbolic on output */
    /* --------------- */
    cholmod_common *Common
)
{
    double zrelax0, zrelax1, zrelax2, xxsize ;
    Int *Wi, *Wj, *Super, *Snz, *Ap, *Ai, *Flag, *Head, *Ls, *Lpi, *Lpx, *Fnz,
	*Sparent, *Anz, *SuperMap, *Merged, *Nscol, *Zeros, *Fp, *Fj,
	*ColCount, *Lpi2, *Lsuper, *Iwork ;
    Int nsuper, d, n, j, k, s, mark, parent, p, pend, k1, k2, packed, nscol,
	nsrow, ndrow1, ndrow2, stype, ssize, xsize, sparent, plast, slast,
	csize, maxcsize, ss, nscol0, nscol1, ns, nfsuper, newzeros, totzeros,
	merge, snext, esize, maxesize, nrelax0, nrelax1, nrelax2, Asorted ;
    size_t w ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_NULL (L, FALSE) ;
    RETURN_IF_NULL (Parent, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (L, CHOLMOD_PATTERN, CHOLMOD_PATTERN, FALSE) ;
    stype = A->stype ;
    if (stype < 0)
    {
	/* invalid symmetry; symmetric lower form not supported */
	ERROR (CHOLMOD_INVALID, "symmetric lower not supported") ;
	return (FALSE) ;
    }
    if (stype == 0)
    {
	/* F must be present in the unsymmetric case */
	RETURN_IF_NULL (F, FALSE) ;
    }
    if (L->is_super)
    {
	/* L must be a simplicial symbolic factor */
	ERROR (CHOLMOD_INVALID, "L must be symbolic on input") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    n = A->nrow ;

    /* w = 5*n */
    w = CHOLMOD(mult_size_t) (n, 5, &ok) ;
    if (!ok)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (FALSE) ;
    }

    CHOLMOD(allocate_work) (n, w, 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (FALSE) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* A is now either A or triu(A(p,p)) for the symmetric case.  It is either
     * A or A(p,f) for the unsymmetric case (both in column form).  It can be
     * either packed or unpacked, and either sorted or unsorted.  Entries in
     * the lower triangular part may be present if A is symmetric, but these
     * are ignored. */

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;

    if (stype != 0)
    {
	/* F not accessed */
	Fp = NULL ;
	Fj = NULL ;
	Fnz = NULL ;
	packed = TRUE ;
    }
    else
    {
	/* F = A(:,f) or A(p,f) in packed row form, either sorted or unsorted */
	Fp = F->p ;
	Fj = F->i ;
	Fnz = F->nz ;
	packed = F->packed ;
    }

    ColCount = L->ColCount ;

    nrelax0 = Common->nrelax [0] ;
    nrelax1 = Common->nrelax [1] ;
    nrelax2 = Common->nrelax [2] ;

    zrelax0 = Common->zrelax [0] ;
    zrelax1 = Common->zrelax [1] ;
    zrelax2 = Common->zrelax [2] ;

    zrelax0 = IS_NAN (zrelax0) ? 0 : zrelax0 ;
    zrelax1 = IS_NAN (zrelax1) ? 0 : zrelax1 ;
    zrelax2 = IS_NAN (zrelax2) ? 0 : zrelax2 ;

    ASSERT (CHOLMOD(dump_parent) (Parent, n, "Parent", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    /* Sparent, Snz, and Merged could be allocated later, of size nfsuper */

    Iwork = Common->Iwork ;
    Wi      = Iwork ;	    /* size n (i/l/l).  Lpi2 is i/l/l */
    Wj      = Iwork + n ;   /* size n (i/l/l).  Zeros is i/l/l */
    Sparent = Iwork + 2*((size_t) n) ; /* size nfsuper <= n [ */
    Snz     = Iwork + 3*((size_t) n) ; /* size nfsuper <= n [ */
    Merged  = Iwork + 4*((size_t) n) ; /* size nfsuper <= n [ */

    Flag = Common->Flag ;   /* size n */
    Head = Common->Head ;   /* size n+1 */

    /* ---------------------------------------------------------------------- */
    /* find the fundamental supernodes */
    /* ---------------------------------------------------------------------- */

    /* count the number of children of each node, using Wi [ */
    for (j = 0 ; j < n ; j++)
    {
	Wi [j] = 0 ;
    }
    for (j = 0 ; j < n ; j++)
    {
	parent = Parent [j] ;
	if (parent != EMPTY)
	{
	    Wi [parent]++ ;
	}
    }

    Super = Head ;  /* use Head [0..nfsuper] as workspace for Super list ( */

    /* column 0 always starts a new supernode */
    nfsuper = (n == 0) ? 0 : 1 ;	/* number of fundamental supernodes */
    Super [0] = 0 ;

    for (j = 1 ; j < n ; j++)
    {
	/* check if j starts new supernode, or in the same supernode as j-1 */
	if (Parent [j-1] != j	    /* parent of j-1 is not j */
	    || (ColCount [j-1] != ColCount [j] + 1) /* j-1 not subset of j*/
	    || Wi [j] > 1)	    /* j has more than one child */
	{
	    /* j is the leading node of a supernode */
	    Super [nfsuper++] = j ;
	}
    }
    Super [nfsuper] = n ;

    /* contents of Wi no longer needed for child count ] */

    Nscol = Wi ; /* use Wi as size-nfsuper workspace for Nscol [ */

    /* ---------------------------------------------------------------------- */
    /* find the mapping of fundamental nodes to supernodes */
    /* ---------------------------------------------------------------------- */

    SuperMap = Wj ;	/* use Wj as workspace for SuperMap [ */

    /* SuperMap [k] = s if column k is contained in supernode s */
    for (s = 0 ; s < nfsuper ; s++)
    {
	for (k = Super [s] ; k < Super [s+1] ; k++)
	{
	    SuperMap [k] = s ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct the fundamental supernodal etree */
    /* ---------------------------------------------------------------------- */

    for (s = 0 ; s < nfsuper ; s++)
    {
	j = Super [s+1] - 1 ;	/* last node in supernode s */
	parent = Parent [j] ;	/* parent of last node */
	Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
	PRINT1 (("Sparent ["ID"] = "ID"\n", s, Sparent [s])) ;
    }

    /* contents of Wj no longer needed as workspace for SuperMap ]
     * SuperMap will be recomputed below, for the relaxed supernodes. */

    Zeros = Wj ;   /* use Wj for Zeros, workspace of size nfsuper [ */

    /* ---------------------------------------------------------------------- */
    /* relaxed amalgamation */
    /* ---------------------------------------------------------------------- */

    for (s = 0 ; s < nfsuper ; s++)
    {
	Merged [s] = EMPTY ;			/* s not merged into another */
	Nscol [s] = Super [s+1] - Super [s] ;	/* # of columns in s */
	Zeros [s] = 0 ;				/* # of zero entries in s */
	ASSERT (s <= Super [s]) ;
	Snz [s] = ColCount [Super [s]] ;  /* # of entries in leading col of s */
	PRINT2 (("lnz ["ID"] "ID"\n", s, Snz [s])) ;
    }

    for (s = nfsuper-2 ; s >= 0 ; s--)
    {
	/* should supernodes s and s+1 merge into a new node s? */
	PRINT1 (("\n========= Check relax of s "ID" and s+1 "ID"\n", s, s+1)) ;

	ss = Sparent [s] ;
	if (ss == EMPTY)
	{
	    PRINT1 (("s "ID" is a root, no merge with s+1 = "ID"\n", s, s+1)) ;
	    continue ;
	}

	/* find the current parent of s (perform path compression as needed) */
	for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = Merged [ss]) ;
	sparent = ss ;
	PRINT2 (("Current sparent of s "ID" is "ID"\n", s, sparent)) ;

	/* ss is the current parent of s */
	for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = snext)
	{
	    snext = Merged [ss] ;
	    PRINT2 (("ss "ID" is dead, merged into snext "ID"\n", ss, snext)) ;
	    Merged [ss] = sparent ;
	}

	/* if s+1 is not the current parent of s, do not merge */
	if (sparent != s+1)
	{
	    continue ;
	}

	nscol0 = Nscol [s] ;	/* # of columns in s */
	nscol1 = Nscol [s+1] ;	/* # of columns in s+1 */
	ns = nscol0 + nscol1 ;
	PRINT2 (("ns "ID" nscol0 "ID" nscol1 "ID"\n", ns, nscol0, nscol1)) ;

	totzeros = Zeros [s+1] ;	/* current # of zeros in s+1 */

	/* determine if supernodes s and s+1 should merge */
	if (ns <= nrelax0)
	{
	    PRINT2 (("ns is tiny ("ID"), so go ahead and merge\n", ns)) ;
	    merge = TRUE ;
	}
	else
	{
	    /* use double to avoid integer overflow */
	    double lnz0 = Snz [s] ;	/* # entries in leading column of s */
	    double lnz1 = Snz [s+1] ;	/* # entries in leading column of s+1 */
	    double xnewzeros = nscol0 * (lnz1 + nscol0 - lnz0) ;

	    /* use Int for the final update of Zeros [s] below */
	    newzeros = nscol0 * (Snz [s+1] + nscol0 - Snz [s]) ;
	    ASSERT (newzeros == xnewzeros) ;

	    PRINT2 (("lnz0 %g lnz1 %g xnewzeros %g\n", lnz0, lnz1, xnewzeros)) ;
	    if (xnewzeros == 0)
	    {
		/* no new zeros, so go ahead and merge */
		PRINT2 (("no new fillin, so go ahead and merge\n")) ;
		merge = TRUE ;
	    }
	    else
	    {
		/* # of zeros if merged */
		double xtotzeros = ((double) totzeros) + xnewzeros ;

		/* xtotsize: total size of merged supernode, if merged: */
		double xns = (double) ns ;
		double xtotsize  = (xns * (xns+1) / 2) + xns * (lnz1 - nscol1) ;
		double z = xtotzeros / xtotsize ;

		Int totsize ;
		totsize  = (ns * (ns+1) / 2) + ns * (Snz [s+1] - nscol1) ;

		PRINT2 (("oldzeros "ID" newzeros "ID" xtotsize %g z %g\n",
			    Zeros [s+1], newzeros, xtotsize, z)) ;

		/* use Int for the final update of Zeros [s] below */
		totzeros += newzeros ;

		/* do not merge if supernode would become too big
		 * (Int overflow).  Continue computing; not (yet) an error. */
		/* fl.pt. compare, but no NaN's can occur here */
		merge = ((ns <= nrelax1 && z < zrelax0) ||
			 (ns <= nrelax2 && z < zrelax1) ||
					  (z < zrelax2)) &&
			(xtotsize < Int_max / sizeof (double)) ;

	    }
	}

	if (merge)
	{
	    PRINT1 (("Merge node s ("ID") and s+1 ("ID")\n", s, s+1)) ;
	    Zeros [s] = totzeros ;
	    Merged [s+1] = s ;
	    Snz [s] = nscol0 + Snz [s+1] ;
	    Nscol [s] += Nscol [s+1] ;
	}
    }

    /* contents of Wj no longer needed for Zeros ] */
    /* contents of Wi no longer needed for Nscol ] */
    /* contents of Sparent no longer needed (recomputed below) */

    /* ---------------------------------------------------------------------- */
    /* construct the relaxed supernode list */
    /* ---------------------------------------------------------------------- */

    nsuper = 0 ;
    for (s = 0 ; s < nfsuper ; s++)
    {
	if (Merged [s] == EMPTY)
	{
	    PRINT1 (("live supernode: "ID" snz "ID"\n", s, Snz [s])) ;
	    Super [nsuper] = Super [s] ;
	    Snz [nsuper] = Snz [s] ;
	    nsuper++ ;
	}
    }
    Super [nsuper] = n ;
    PRINT1 (("Fundamental supernodes: "ID"  relaxed "ID"\n", nfsuper, nsuper)) ;

    /* Merged no longer needed ] */

    /* ---------------------------------------------------------------------- */
    /* find the mapping of relaxed nodes to supernodes */
    /* ---------------------------------------------------------------------- */

    /* use Wj as workspace for SuperMap { */

    /* SuperMap [k] = s if column k is contained in supernode s */
    for (s = 0 ; s < nsuper ; s++)
    {
	for (k = Super [s] ; k < Super [s+1] ; k++)
	{
	    SuperMap [k] = s ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* construct the relaxed supernodal etree */
    /* ---------------------------------------------------------------------- */

    for (s = 0 ; s < nsuper ; s++)
    {
	j = Super [s+1] - 1 ;	/* last node in supernode s */
	parent = Parent [j] ;	/* parent of last node */
	Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
	PRINT1 (("new Sparent ["ID"] = "ID"\n", s, Sparent [s])) ;
    }

    /* ---------------------------------------------------------------------- */
    /* determine the size of L->s and L->x */
    /* ---------------------------------------------------------------------- */

    ssize = 0 ;
    xsize = 0 ;
    xxsize = 0 ;
    for (s = 0 ; s < nsuper ; s++)
    {
	nscol = Super [s+1] - Super [s] ;
	nsrow = Snz [s] ;
	ASSERT (nscol > 0) ;
	ssize += nsrow ;
        if (for_cholesky)
        {
            xsize += nscol * nsrow ;
            /* also compute xsize in double to guard against Int overflow */
            xxsize += ((double) nscol) * ((double) nsrow) ;
        }
	if (ssize < 0 || (for_cholesky && xxsize > Int_max))
	{
	    /* Int overflow, clear workspace and return.
               QR factorization will not use xxsize, so that error is ignored.
               For Cholesky factorization, however, memory of space xxsize
               will be allocated, so this is a failure.  Both QR and Cholesky
               fail if ssize overflows. */
	    ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	    FREE_WORKSPACE ;
	    return (FALSE) ;
	}
	ASSERT (ssize > 0) ;
        ASSERT (IMPLIES (for_cholesky, xsize > 0)) ;
    }
    xsize = MAX (1, xsize) ;
    ssize = MAX (1, ssize) ;
    PRINT1 (("ix sizes: "ID" "ID" nsuper "ID"\n", ssize, xsize, nsuper)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate L (all except real part L->x) */
    /* ---------------------------------------------------------------------- */

    L->ssize = ssize ;
    L->xsize = xsize ;
    L->nsuper = nsuper ;

    CHOLMOD(change_factor) (CHOLMOD_PATTERN, TRUE, TRUE, TRUE, TRUE, L, Common);

    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory; L is still a valid simplicial symbolic factor */
	FREE_WORKSPACE ;
	return (FALSE) ;
    }

    DEBUG (CHOLMOD(dump_factor) (L, "L to symbolic super", Common)) ;
    ASSERT (L->is_ll && L->xtype == CHOLMOD_PATTERN && L->is_super) ;

    Lpi = L->pi ;
    Lpx = L->px ;
    Ls = L->s ;
    Ls [0] = 0 ;    /* flag for cholmod_check_factor; supernodes are defined */
    Lpx [0] = for_cholesky ? 0 : 123456 ;   /* magic number for sparse QR */
    Lsuper = L->super ;

    /* copy the list of relaxed supernodes into the final list in L */
    for (s = 0 ; s <= nsuper ; s++)
    {
	Lsuper [s] = Super [s] ;
    }

    /* Head no longer needed as workspace for fundamental Super list ) */

    Super = Lsuper ;	    /* Super is now the list of relaxed supernodes */

    /* ---------------------------------------------------------------------- */
    /* construct column pointers of relaxed supernodal pattern (L->pi) */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (s = 0 ; s < nsuper ; s++)
    {
	Lpi [s] = p ;
	p += Snz [s] ;
	PRINT1 (("Snz ["ID"] = "ID", Super ["ID"] = "ID"\n",
		    s, Snz [s], s, Super[s])) ;
    }
    Lpi [nsuper] = p ;
    ASSERT ((Int) (L->ssize) == MAX (1,p)) ;

    /* ---------------------------------------------------------------------- */
    /* construct pointers for supernodal values (L->px) */
    /* ---------------------------------------------------------------------- */

    if (for_cholesky)
    {
        /* L->px is not needed for QR factorization (it may lead to Int
           overflow, anyway, if xsize caused Int overflow above) */
        p = 0 ;
        for (s = 0 ; s < nsuper ; s++)
        {
            nscol = Super [s+1] - Super [s] ;   /* number of columns in s */
            nsrow = Snz [s] ;           /* # of rows, incl triangular part*/
            Lpx [s] = p ;               /* pointer to numerical part of s */
            p += nscol * nsrow ;
        }
        Lpx [s] = p ;
        ASSERT ((Int) (L->xsize) == MAX (1,p)) ;
    }

    /* Snz no longer needed ] */

    /* ---------------------------------------------------------------------- */
    /* symbolic analysis to construct the relaxed supernodal pattern (L->s) */
    /* ---------------------------------------------------------------------- */

    Lpi2 = Wi ;	    /* copy Lpi into Lpi2, using Wi as workspace for Lpi2 [ */
    for (s = 0 ; s < nsuper ; s++)
    {
	Lpi2 [s] = Lpi [s] ;
    }

    Asorted = A->sorted ;

    for (s = 0 ; s < nsuper ; s++)
    {
	/* sth supernode is in columns k1 to k2-1.
	 * compute nonzero pattern of L (k1:k2-1,:). */

	/* place rows k1 to k2-1 in leading column of supernode s */
	k1 = Super [s] ;
	k2 = Super [s+1] ;
	PRINT1 (("=========>>> Supernode "ID" k1 "ID" k2-1 "ID"\n",
		    s, k1, k2-1)) ;
	for (k = k1 ; k < k2 ; k++)
	{
	    Ls [Lpi2 [s]++] = k ;
	}

	/* compute nonzero pattern each row k1 to k2-1 */
	for (k = k1 ; k < k2 ; k++)
	{
	    /* compute row k of L.  In the symmetric case, the pattern of L(k,:)
	     * is the set of nodes reachable in the supernodal etree from any
	     * row i in the nonzero pattern of A(0:k,k).  In the unsymmetric
	     * case, the pattern of the kth column of A*A' is the set union
	     * of all columns A(0:k,j) for each nonzero F(j,k). */

	    /* clear the Flag array and mark the current supernode */
	    /* mark = CHOLMOD(clear_flag) (Common) ; */
	    CHOLMOD_CLEAR_FLAG (Common) ;
	    mark = Common->mark ;
	    Flag [s] = mark ;
	    ASSERT (s == SuperMap [k]) ;

	    /* traverse the row subtree for each nonzero in A or AA' */
	    if (stype != 0)
	    {
		subtree (k, k, Ap, Ai, Anz, SuperMap, Sparent, mark,
                        Asorted, k1, Flag, Ls, Lpi2) ;
	    }
	    else
	    {
		/* for each j nonzero in F (:,k) do */
		p = Fp [k] ;
		pend = (packed) ? (Fp [k+1]) : (p + Fnz [k]) ;
		for ( ; p < pend ; p++)
		{
		    subtree (Fj [p], k, Ap, Ai, Anz, SuperMap, Sparent, mark,
			    Asorted, k1, Flag, Ls, Lpi2) ;
		}
	    }
	}
    }
#ifndef NDEBUG
    for (s = 0 ; s < nsuper ; s++)
    {
	PRINT1 (("Lpi2[s] "ID" Lpi[s+1] "ID"\n", Lpi2 [s], Lpi [s+1])) ;
	ASSERT (Lpi2 [s] == Lpi [s+1]) ;
	CHOLMOD(dump_super) (s, Super, Lpi, Ls, NULL, NULL, 0, Common) ;
    }
#endif

    /* contents of Wi no longer needed for Lpi2 ] */
    /* Sparent no longer needed ] */

    /* ---------------------------------------------------------------------- */
    /* determine the largest update matrix (L->maxcsize) */
    /* ---------------------------------------------------------------------- */

    /* maxcsize could be determined before L->s is allocated and defined, which
     * would mean that all memory requirements for both the symbolic and numeric
     * factorizations could be computed using O(nnz(A)+O(n)) space.  However, it
     * would require a lot of extra work.  The analysis phase, above, would need
     * to be duplicated, but with Ls not kept; instead, the algorithm would keep
     * track of the current s and slast for each supernode d, and update them
     * when a new row index appears in supernode d.  An alternative would be to
     * do this computation only if the allocation of L->s failed, in which case
     * the following code would be skipped.
     *
     * The csize for a supernode is the size of its largest contribution to
     * a subsequent ancestor supernode.  For example, suppose the rows of #'s
     * in the figure below correspond to the columns of a subsequent supernode,
     * and the dots are the entries in that ancestore.
     *
     *	    c
     *	    c c
     *	    c c c
     *	    x x x
     *	    x x x
     *	    # # #   .
     *	    # # #   . .
     *	    * * *   . .
     *	    * * *   . .
     *	    * * *   . .
     *	            . .
     *
     * Then for this update, the csize is 3-by-2, or 6, because there are 3
     * rows of *'s which is the number of rows in the update, and there are
     * 2 rows of #'s, which is the number columns in the update.  The csize
     * of a supernode is the largest such contribution for any ancestor
     * supernode.  maxcsize, for the whole matrix, has a rough upper bound of
     * the maximum size of any supernode.  This bound is loose, because the
     * the contribution must be less than the size of the ancestor supernodal
     * that it's updating.  maxcsize of a completely dense matrix, with one
     * supernode, is zero.
     *
     * maxesize is the column dimension for the workspace E needed for the
     * solve.  E is of size nrhs-by-maxesize, where the nrhs is the number of
     * columns in the right-hand-side.  The maxesize is the largest esize of
     * any supernode.  The esize of a supernode is the number of row indices
     * it contains, excluding the column indices of the supernode itself.
     * For the following example, esize is 4:
     *
     *	    c
     *	    c c
     *	    c c c
     *	    x x x
     *	    x x x
     *	    x x x
     *	    x x x
     *
     * maxesize can be no bigger than n.
     */

    maxcsize = 1 ;
    maxesize = 1 ;

    /* Do not need to guard csize against Int overflow since xsize is OK. */

    if (for_cholesky)
    {
        /* this is not needed for QR factorization */
        for (d = 0 ; d < nsuper ; d++)
        {
            nscol = Super [d+1] - Super [d] ;
            p = Lpi [d] + nscol ;
            plast = p ;
            pend = Lpi [d+1] ;
            esize = pend - p ;
            maxesize = MAX (maxesize, esize) ;
            slast = (p == pend) ? (EMPTY) : (SuperMap [Ls [p]]) ;
            for ( ; p <= pend ; p++)
            {
                s = (p == pend) ? (EMPTY) : (SuperMap [Ls [p]]) ;
                if (s != slast)
                {
                    /* row i is the start of a new supernode */
                    ndrow1 = p - plast ;
                    ndrow2 = pend - plast ;
                    csize = ndrow2 * ndrow1 ;
                    PRINT1 (("Supernode "ID" ancestor "ID" C: "ID"-by-"ID
                        "  csize "ID"\n", d, slast, ndrow1, ndrow2, csize)) ;
                    maxcsize = MAX (maxcsize, csize) ;
                    plast = p ;
                    slast = s ;
                }
            }
        }
        PRINT1 (("max csize "ID"\n", maxcsize)) ;
    }

    /* Wj no longer needed for SuperMap } */

    L->maxcsize = maxcsize ;
    L->maxesize = maxesize ;
    L->is_super = TRUE ;
    ASSERT (L->xtype == CHOLMOD_PATTERN && L->is_ll) ;

    /* ---------------------------------------------------------------------- */
    /* supernodal symbolic factorization is complete */
    /* ---------------------------------------------------------------------- */

    FREE_WORKSPACE ;
    return (TRUE) ;
}

/* ========================================================================== */
/* === cholmod_super_symbolic =============================================== */
/* ========================================================================== */

/* Analyzes A, AA', or A(:,f)*A(:,f)' in preparation for a supernodal numeric
 * factorization.  The user need not call this directly; cholmod_analyze is
 * a "simple" wrapper for this routine.
 * 
 * This function does all the analysis for a supernodal Cholesky factorization.
 *
 * workspace: Flag (nrow), Head (nrow), Iwork (2*nrow),
 * and temporary space of size 3*nfsuper*sizeof(Int), where nfsuper <= n
 * is the number of fundamental supernodes.
 */

int CHOLMOD(super_symbolic)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to analyze */
    cholmod_sparse *F,	/* F = A' or A(:,f)' */
    Int *Parent,	/* elimination tree */
    /* ---- in/out --- */
    cholmod_factor *L,	/* simplicial symbolic on input,
			 * supernodal symbolic on output */
    /* --------------- */
    cholmod_common *Common
)
{
    return (CHOLMOD(super_symbolic2) (TRUE, A, F, Parent, L, Common)) ;
}
#endif
