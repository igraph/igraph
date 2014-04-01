/* ========================================================================== */
/* === Core/cholmod_aat ===================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* C = A*A' or C = A(:,f)*A(:,f)'
 *
 * A can be packed or unpacked, sorted or unsorted, but must be stored with
 * both upper and lower parts (A->stype of zero).  C is returned as packed,
 * C->stype of zero (both upper and lower parts present), and unsorted.  See
 * cholmod_ssmult in the MatrixOps Module for a more general matrix-matrix
 * multiply.
 *
 * You can trivially convert C into a symmetric upper/lower matrix by
 * changing C->stype = 1 or -1 after calling this routine.
 *
 * workspace:
 *	Flag (A->nrow),
 *	Iwork (max (A->nrow, A->ncol)) if fset present,
 *	Iwork (A->nrow) if no fset,
 *	W (A->nrow) if mode > 0,
 *	allocates temporary copy for A'.
 *
 * A can be pattern or real.  Complex or zomplex cases are supported only
 * if the mode is <= 0 (in which case the numerical values are ignored).
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"

cholmod_sparse *CHOLMOD(aat)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* input matrix; C=A*A' is constructed */
    Int *fset,		/* subset of 0:(A->ncol)-1 */
    size_t fsize,	/* size of fset */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag)
			 * -2: pattern only, no diagonal, add 50% + n extra
			 * space to C */
    /* --------------- */
    cholmod_common *Common
)
{
    double fjt ;
    double *Ax, *Fx, *Cx, *W ;
    Int *Ap, *Anz, *Ai, *Fp, *Fi, *Cp, *Ci, *Flag ;
    cholmod_sparse *C, *F ;
    Int packed, j, i, pa, paend, pf, pfend, n, mark, cnz, t, p, values, diag,
	extra ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    values = (mode > 0) && (A->xtype != CHOLMOD_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN,
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    if (A->stype)
    {
	ERROR (CHOLMOD_INVALID, "matrix cannot be symmetric") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    diag = (mode >= 0) ;
    n = A->nrow ;
    CHOLMOD(allocate_work) (n, MAX (A->ncol, A->nrow), values ? n : 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n : 0, Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (A, "A", Common) >= 0) ;

    /* get the A matrix */
    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;

    /* get workspace */
    W = Common->Xwork ;		/* size n, unused if values is FALSE */
    Flag = Common->Flag ;	/* size n, Flag [0..n-1] < mark on input*/

    /* ---------------------------------------------------------------------- */
    /* F = A' or A(:,f)' */
    /* ---------------------------------------------------------------------- */

    /* workspace: Iwork (nrow if no fset; MAX (nrow,ncol) if fset)*/
    F = CHOLMOD(ptranspose) (A, values, NULL, fset, fsize, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    Fp = F->p ;
    Fi = F->i ;
    Fx = F->x ;

    /* ---------------------------------------------------------------------- */
    /* count the number of entries in the result C */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	/* clear the Flag array */
	/* mark = CHOLMOD(clear_flag) (Common) ; */
	CHOLMOD_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	/* exclude the diagonal, if requested */
	if (!diag)
	{
	    Flag [j] = mark ;
	}

	/* for each nonzero F(t,j) in column j, do: */
	pfend = Fp [j+1] ;
	for (pf = Fp [j] ; pf < pfend ; pf++)
	{
	    /* F(t,j) is nonzero */
	    t = Fi [pf] ;

	    /* add the nonzero pattern of A(:,t) to the pattern of C(:,j) */
	    pa = Ap [t] ;
	    paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
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

    extra = (mode == -2) ? (cnz/2 + n) : 0 ;

    mark = CHOLMOD(clear_flag) (Common) ;

    /* ---------------------------------------------------------------------- */
    /* check for integer overflow */
    /* ---------------------------------------------------------------------- */

    if (cnz < 0 || (cnz + extra) < 0)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	CHOLMOD(clear_flag) (Common) ;
	CHOLMOD(free_sparse) (&F, Common) ;
	return (NULL) ;	    /* problem too large */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate C */
    /* ---------------------------------------------------------------------- */

    C = CHOLMOD(allocate_sparse) (n, n, cnz + extra, FALSE, TRUE, 0,
	    values ? A->xtype : CHOLMOD_PATTERN, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free_sparse) (&F, Common) ;
	return (NULL) ;	    /* out of memory */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* C = A*A' */
    /* ---------------------------------------------------------------------- */

    cnz = 0 ;

    if (values)
    {

	/* pattern and values */
	for (j = 0 ; j < n ; j++)
	{
	    /* clear the Flag array */
	    mark = CHOLMOD(clear_flag) (Common) ;

	    /* start column j of C */
	    Cp [j] = cnz ;

	    /* for each nonzero F(t,j) in column j, do: */
	    pfend = Fp [j+1] ;
	    for (pf = Fp [j] ; pf < pfend ; pf++)
	    {
		/* F(t,j) is nonzero */
		t = Fi [pf] ;
		fjt = Fx [pf] ;

		/* add the nonzero pattern of A(:,t) to the pattern of C(:,j)
		 * and scatter the values into W */
		pa = Ap [t] ;
		paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
		for ( ; pa < paend ; pa++)
		{
		    i = Ai [pa] ;
		    if (Flag [i] != mark)
		    {
			Flag [i] = mark ;
			Ci [cnz++] = i ;
		    }
		    W [i] += Ax [pa] * fjt ;
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
	for (j = 0 ; j < n ; j++)
	{
	    /* clear the Flag array */
	    mark = CHOLMOD(clear_flag) (Common) ;

	    /* exclude the diagonal, if requested */
	    if (!diag)
	    {
		Flag [j] = mark ;
	    }

	    /* start column j of C */
	    Cp [j] = cnz ;

	    /* for each nonzero F(t,j) in column j, do: */
	    pfend = Fp [j+1] ;
	    for (pf = Fp [j] ; pf < pfend ; pf++)
	    {
		/* F(t,j) is nonzero */
		t = Fi [pf] ;

		/* add the nonzero pattern of A(:,t) to the pattern of C(:,j) */
		pa = Ap [t] ;
		paend = (packed) ? (Ap [t+1]) : (pa + Anz [t]) ;
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

    Cp [n] = cnz ;
    ASSERT (IMPLIES (mode != -2, MAX (1,cnz) == C->nzmax)) ;

    /* ---------------------------------------------------------------------- */
    /* clear workspace and free temporary matrices and return result */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(free_sparse) (&F, Common) ;
    CHOLMOD(clear_flag) (Common) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, values ? n : 0, Common)) ;
    DEBUG (i = CHOLMOD(dump_sparse) (C, "aat", Common)) ;
    ASSERT (IMPLIES (mode < 0, i == 0)) ;
    return (C) ;
}
