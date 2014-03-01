/* ========================================================================== */
/* === MatrixOps/cholmod_ssmult ============================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/MatrixOps Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/MatrixOps Module is licensed under Version 2.0 of the GNU
 * General Public License.  See gpl.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * http://www.suitesparse.com
 * -------------------------------------------------------------------------- */

/* C = A*B.  Multiply two sparse matrices.
 *
 * A and B can be packed or unpacked, sorted or unsorted, and of any stype.
 * If A or B are symmetric, an internal unsymmetric copy is made first, however.
 * C is computed as if A and B are unsymmetric, and then if the stype input
 * parameter requests a symmetric form (upper or lower) the matrix is converted
 * into that form.
 *
 * C is returned as packed, and either unsorted or sorted, depending on the
 * "sorted" input parameter.  If C is returned sorted, then either C = (B'*A')'
 * or C = (A*B)'' is computed, depending on the number of nonzeros in A, B, and
 * C.
 *
 * workspace:
 *	if C unsorted: Flag (A->nrow), W (A->nrow) if values
 *	if C sorted:   Flag (B->ncol), W (B->ncol) if values
 *	Iwork (max (A->ncol, A->nrow, B->nrow, B->ncol))
 *	allocates temporary copies for A, B, and C, if required.
 *
 * Only pattern and real matrices are supported.  Complex and zomplex matrices
 * are supported only when the numerical values are not computed ("values"
 * is FALSE).
 */

#ifndef NMATRIXOPS

#include "cholmod_internal.h"
#include "cholmod_matrixops.h"


/* ========================================================================== */
/* === cholmod_ssmult ======================================================= */
/* ========================================================================== */

cholmod_sparse *CHOLMOD(ssmult)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* left matrix to multiply */
    cholmod_sparse *B,	/* right matrix to multiply */
    int stype,		/* requested stype of C */
    int values,		/* TRUE: do numerical values, FALSE: pattern only */
    int sorted,		/* if TRUE then return C with sorted columns */
    /* --------------- */
    cholmod_common *Common
)
{
    double bjt ;
    double *Ax, *Bx, *Cx, *W ;
    Int *Ap, *Anz, *Ai, *Bp, *Bnz, *Bi, *Cp, *Ci, *Flag ;
    cholmod_sparse *C, *A2, *B2, *A3, *B3, *C2 ;
    Int apacked, bpacked, j, i, pa, paend, pb, pbend, ncol, mark, cnz, t, p,
	nrow, anz, bnz, do_swap_and_transpose, n1, n2 ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_NULL (B, NULL) ;
    values = values &&
	(A->xtype != CHOLMOD_PATTERN) && (B->xtype != CHOLMOD_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, 
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    RETURN_IF_XTYPE_INVALID (B, CHOLMOD_PATTERN, 
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    if (A->ncol != B->nrow)
    {
	/* inner dimensions must agree */
	ERROR (CHOLMOD_INVALID, "A and B inner dimensions must match") ;
	return (NULL) ;
    }
    /* A and B must have the same numerical type if values is TRUE (both must
     * be CHOLMOD_REAL, this is implicitly checked above) */
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    if (A->nrow <= 1)
    {
	/* C will be implicitly sorted, so no need to sort it here */
	sorted = FALSE ;
    }
    if (sorted)
    {
	n1 = MAX (A->nrow, B->ncol) ;
    }
    else
    {
	n1 = A->nrow ;
    }
    n2 = MAX4 (A->ncol, A->nrow, B->nrow, B->ncol) ;
    CHOLMOD(allocate_work) (n1, n2, values ? n1 : 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (NULL) ;
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1 : 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    /* convert A to unsymmetric, if necessary */
    A2 = NULL ;
    B2 = NULL ;
    if (A->stype)
    {
	/* workspace: Iwork (max (A->nrow,A->ncol)) */
	A2 = CHOLMOD(copy) (A, 0, values, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common)) ;
	    return (NULL) ;
	}
	A = A2 ;
    }

    /* convert B to unsymmetric, if necessary */
    if (B->stype)
    {
	/* workspace: Iwork (max (B->nrow,B->ncol)) */
	B2 = CHOLMOD(copy) (B, 0, values, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    CHOLMOD(free_sparse) (&A2, Common) ;
	    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common)) ;
	    return (NULL) ;
	}
	B = B2 ;
    }

    ASSERT (CHOLMOD(dump_sparse) (A, "A", Common) >= 0) ;
    ASSERT (CHOLMOD(dump_sparse) (B, "B", Common) >= 0) ;

    /* get the A matrix */
    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    apacked = A->packed ;

    /* get the B matrix */
    Bp  = B->p ;
    Bnz = B->nz ;
    Bi  = B->i ;
    Bx  = B->x ;
    bpacked = B->packed ;

    /* get the size of C */
    nrow = A->nrow ;
    ncol = B->ncol ;

    /* get workspace */
    W = Common->Xwork ;		/* size nrow, unused if values is FALSE */
    Flag = Common->Flag ;	/* size nrow, Flag [0..nrow-1] < mark on input*/

    /* ---------------------------------------------------------------------- */
    /* count the number of entries in the result C */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	/* clear the Flag array */
	/* mark = CHOLMOD(clear_flag) (Common) ; */
	CHOLMOD_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	/* for each nonzero B(t,j) in column j, do: */
	pb = Bp [j] ;
	pbend = (bpacked) ? (Bp [j+1]) : (pb + Bnz [j]) ;
	for ( ; pb < pbend ; pb++)
	{
	    /* B(t,j) is nonzero */
	    t = Bi [pb] ;

	    /* add the nonzero pattern of A(:,t) to the pattern of C(:,j) */
	    pa = Ap [t] ;
	    paend = (apacked) ? (Ap [t+1]) : (pa + Anz [t]) ;
	    for ( ; pa < paend ; pa++)
	    {
		i = Ai [pa] ;
		if (Flag [i] != mark)
		{
		    Flag [i] = mark ;
		    cnz++ ;
		}
	    }
	}
	if (cnz < 0)
	{
	    break ;	    /* integer overflow case */
	}
    }

    /* mark = CHOLMOD(clear_flag) (Common) ; */
    CHOLMOD_CLEAR_FLAG (Common) ;
    mark = Common->mark ;

    /* ---------------------------------------------------------------------- */
    /* check for integer overflow */
    /* ---------------------------------------------------------------------- */

    if (cnz < 0)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	CHOLMOD(free_sparse) (&A2, Common) ;
	CHOLMOD(free_sparse) (&B2, Common) ;
	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common)) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* Determine how to return C sorted (if requested) */
    /* ---------------------------------------------------------------------- */

    do_swap_and_transpose = FALSE ;

    if (sorted)
    {
	/* Determine the best way to return C with sorted columns.  Computing
	 * C = (B'*A')' takes cnz + anz + bnz time (ignoring O(n) terms).
	 * Sorting C when done, C = (A*B)'', takes 2*cnz time.  Pick the one
	 * with the least amount of work. */

	anz = CHOLMOD(nnz) (A, Common) ;
	bnz = CHOLMOD(nnz) (B, Common) ;

	do_swap_and_transpose = (anz + bnz < cnz) ;

	if (do_swap_and_transpose)
	{

	    /* -------------------------------------------------------------- */
	    /* C = (B'*A')' */
	    /* -------------------------------------------------------------- */

	    /* workspace: Iwork (A->nrow) */
	    A3 = CHOLMOD(ptranspose) (A, values, NULL, NULL, 0, Common) ;
	    CHOLMOD(free_sparse) (&A2, Common) ;
	    A2 = A3 ;
	    if (Common->status < CHOLMOD_OK)
	    {
		/* out of memory */
		CHOLMOD(free_sparse) (&A2, Common) ;
		CHOLMOD(free_sparse) (&B2, Common) ;
		ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common));
		return (NULL) ;
	    }
	    /* workspace: Iwork (B->nrow) */
	    B3 = CHOLMOD(ptranspose) (B, values, NULL, NULL, 0, Common) ;
	    CHOLMOD(free_sparse) (&B2, Common) ;
	    B2 = B3 ;
	    if (Common->status < CHOLMOD_OK)
	    {
		/* out of memory */
		CHOLMOD(free_sparse) (&A2, Common) ;
		CHOLMOD(free_sparse) (&B2, Common) ;
		ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common));
		return (NULL) ;
	    }
	    A = B2 ;
	    B = A2 ;

	    /* get the new A matrix */
	    Ap  = A->p ;
	    Anz = A->nz ;
	    Ai  = A->i ;
	    Ax  = A->x ;
	    apacked = A->packed ;

	    /* get the new B matrix */
	    Bp  = B->p ;
	    Bnz = B->nz ;
	    Bi  = B->i ;
	    Bx  = B->x ;
	    bpacked = B->packed ;

	    /* get the size of C' */
	    nrow = A->nrow ;
	    ncol = B->ncol ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* allocate C */
    /* ---------------------------------------------------------------------- */

    C = CHOLMOD(allocate_sparse) (nrow, ncol, cnz, FALSE, TRUE, 0,
	    values ? A->xtype : CHOLMOD_PATTERN, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	CHOLMOD(free_sparse) (&A2, Common) ;
	CHOLMOD(free_sparse) (&B2, Common) ;
	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common)) ;
	return (NULL) ;
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* C = A*B */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;

    if (values)
    {

	/* pattern and values */
	for (j = 0 ; j < ncol ; j++)
	{
	    /* clear the Flag array */
	    /* mark = CHOLMOD(clear_flag (Common)) ; */
	    CHOLMOD_CLEAR_FLAG (Common) ;
	    mark = Common->mark ;

	    /* start column j of C */
	    Cp [j] = cnz ;

	    /* for each nonzero B(t,j) in column j, do: */
	    pb = Bp [j] ;
	    pbend = (bpacked) ? (Bp [j+1]) : (pb + Bnz [j]) ;
	    for ( ; pb < pbend ; pb++)
	    {
		/* B(t,j) is nonzero */
		t = Bi [pb] ;
		bjt = Bx [pb] ;

		/* add the nonzero pattern of A(:,t) to the pattern of C(:,j)
		 * and scatter the values into W */
		pa = Ap [t] ;
		paend = (apacked) ? (Ap [t+1]) : (pa + Anz [t]) ;
		for ( ; pa < paend ; pa++)
		{
		    i = Ai [pa] ;
		    if (Flag [i] != mark)
		    {
			Flag [i] = mark ;
			Ci [cnz++] = i ;
		    }
		    W [i] += Ax [pa] * bjt ;
		}
	    }

	    /* gather the values into C(:,j) */
	    for (p = Cp [j] ; p < cnz ; p++)
	    {
		i = Ci [p] ;
		Cx [p] = W [i] ;
		W [i] = 0 ;
	    }
	}

    }
    else
    {

	/* pattern only */
	for (j = 0 ; j < ncol ; j++)
	{
	    /* clear the Flag array */
	    /* mark = CHOLMOD(clear_flag) (Common) ; */
	    CHOLMOD_CLEAR_FLAG (Common) ;
	    mark = Common->mark ;

	    /* start column j of C */
	    Cp [j] = cnz ;

	    /* for each nonzero B(t,j) in column j, do: */
	    pb = Bp [j] ;
	    pbend = (bpacked) ? (Bp [j+1]) : (pb + Bnz [j]) ;
	    for ( ; pb < pbend ; pb++)
	    {
		/* B(t,j) is nonzero */
		t = Bi [pb] ;

		/* add the nonzero pattern of A(:,t) to the pattern of C(:,j) */
		pa = Ap [t] ;
		paend = (apacked) ? (Ap [t+1]) : (pa + Anz [t]) ;
		for ( ; pa < paend ; pa++)
		{
		    i = Ai [pa] ;
		    if (Flag [i] != mark)
		    {
			Flag [i] = mark ;
			Ci [cnz++] = i ;
		    }
		}
	    }
	}
    }

    Cp [ncol] = cnz ;
    ASSERT (MAX (1,cnz) == C->nzmax) ;

    /* ---------------------------------------------------------------------- */
    /* clear workspace and free temporary matrices */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(free_sparse) (&A2, Common) ;
    CHOLMOD(free_sparse) (&B2, Common) ;
    /* CHOLMOD(clear_flag) (Common) ; */
    CHOLMOD_CLEAR_FLAG (Common) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* convert C to a symmetric upper/lower matrix if requested */
    /* ---------------------------------------------------------------------- */

    /* convert C in place, which cannot fail since no memory is allocated */
    if (stype > 0)
    {
	/* C = triu (C), in place */
	(void) CHOLMOD(band_inplace) (0, ncol, values, C, Common) ;
	C->stype = 1 ;
    }
    else if (stype < 0)
    {
	/* C = tril (C), in place */
	(void) CHOLMOD(band_inplace) (-nrow, 0, values, C, Common) ;
	C->stype = -1 ;
    }
    ASSERT (Common->status >= CHOLMOD_OK) ;

    /* ---------------------------------------------------------------------- */
    /* sort C, if requested */
    /* ---------------------------------------------------------------------- */

    if (sorted)
    {
	if (do_swap_and_transpose)
	{
	    /* workspace: Iwork (C->ncol), which is A->nrow since C=(B'*A') */
	    C2 = CHOLMOD(ptranspose) (C, values, NULL, NULL, 0, Common) ;
	    CHOLMOD(free_sparse) (&C, Common) ;
	    if (Common->status < CHOLMOD_OK)
	    {
		/* out of memory */
		ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common));
		return (NULL) ;
	    }
	    C = C2 ;
	}
	else
	{
	    /* workspace: Iwork (max (C->nrow,C->ncol)) */
	    if (!CHOLMOD(sort) (C, Common))
	    {
		/* out of memory */
		CHOLMOD(free_sparse) (&C, Common) ;
		ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common));
		return (NULL) ;
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    DEBUG (CHOLMOD(dump_sparse) (C, "ssmult", Common) >= 0) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n1:0, Common)) ;
    return (C) ;
}
#endif
