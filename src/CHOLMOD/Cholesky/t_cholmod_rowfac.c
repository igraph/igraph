/* ========================================================================== */
/* === Cholesky/t_cholmod_rowfac ============================================ */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Cholesky Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * The CHOLMOD/Cholesky Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Template routine for cholmod_rowfac.  Supports any numeric xtype
 * (real, complex, or zomplex).
 *
 * workspace: Iwork (n), Flag (n), Xwork (n if real, 2*n if complex)
 */

#include "cholmod_template.h"

#ifdef MASK
static int TEMPLATE (cholmod_rowfac_mask)
#else
static int TEMPLATE (cholmod_rowfac)
#endif
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to factorize */
    cholmod_sparse *F,	/* used for A*A' case only. F=A' or A(:,f)' */
    double beta [2],	/* factorize beta*I+A or beta*I+AA' (beta [0] only) */
    size_t kstart,	/* first row to factorize */
    size_t kend,	/* last row to factorize is kend-1 */
#ifdef MASK
    /* These inputs are used for cholmod_rowfac_mask only */
    Int *mask,		/* size A->nrow. if mask[i] then W(i) is set to zero */
    Int *RLinkUp,	/* size A->nrow. link list of rows to compute */
#endif
    /* ---- in/out --- */
    cholmod_factor *L,
    /* --------------- */
    cholmod_common *Common
)
{
    double yx [2], lx [2], fx [2], dk [1], di [1], fl = 0 ;
#ifdef ZOMPLEX
    double yz [1], lz [1], fz [1] ;
#endif
    double *Ax, *Az, *Lx, *Lz, *Wx, *Wz, *Fx, *Fz ;
    Int *Ap, *Anz, *Ai, *Lp, *Lnz, *Li, *Lnext, *Flag, *Stack, *Fp, *Fi, *Fnz,
	*Iwork ;
    Int i, p, k, t, pf, pfend, top, s, mark, pend, n, lnz, is_ll, multadds,
	use_dbound, packed, stype, Fpacked, sorted, nzmax, len, parent ;
#ifndef REAL
    Int dk_imaginary ;
#endif

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    PRINT1 (("\nin cholmod_rowfac, kstart %d kend %d stype %d\n",
		kstart, kend, A->stype)) ;
    DEBUG (CHOLMOD(dump_factor) (L, "Initial L", Common)) ;

    n = A->nrow ;
    stype = A->stype ;

    if (stype > 0)
    {
	/* symmetric upper case: F is not needed.  It may be NULL */
	Fp = NULL ;
	Fi = NULL ;
	Fx = NULL ;
	Fz = NULL ;
	Fnz = NULL ;
	Fpacked = TRUE ;
    }
    else
    {
	/* unsymmetric case: F is required. */
	Fp = F->p ;
	Fi = F->i ;
	Fx = F->x ;
	Fz = F->z ;
	Fnz = F->nz ;
	Fpacked = F->packed ;
    }

    Ap = A->p ;		/* size A->ncol+1, column pointers of A */
    Ai = A->i ;		/* size nz = Ap [A->ncol], row indices of A */
    Ax = A->x ;		/* size nz, numeric values of A */
    Az = A->z ;
    Anz = A->nz ;
    packed = A->packed ;
    sorted = A->sorted ;

    use_dbound = IS_GT_ZERO (Common->dbound) ;

    /* get the current factors L (and D for LDL'); allocate space if needed */
    is_ll = L->is_ll ;
    if (L->xtype == CHOLMOD_PATTERN)
    {
	/* ------------------------------------------------------------------ */
	/* L is symbolic only; allocate and initialize L (and D for LDL') */
	/* ------------------------------------------------------------------ */

	/* workspace: none */
	CHOLMOD(change_factor) (A->xtype, is_ll, FALSE, FALSE, TRUE, L, Common);
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    return (FALSE) ;
	}
	ASSERT (L->minor == (size_t) n) ;
    }
    else if (kstart == 0 && kend == (size_t) n)
    {
	/* ------------------------------------------------------------------ */
	/* refactorization; reset L->nz and L->minor to restart factorization */
	/* ------------------------------------------------------------------ */

	L->minor = n ;
	Lnz = L->nz ;
	for (k = 0 ; k < n ; k++)
	{
	    Lnz [k] = 1 ;
	}
    }

    ASSERT (is_ll == L->is_ll) ;
    ASSERT (L->xtype != CHOLMOD_PATTERN) ;
    DEBUG (CHOLMOD(dump_factor) (L, "L ready", Common)) ;
    DEBUG (CHOLMOD(dump_sparse) (A, "A ready", Common)) ;
    DEBUG (if (stype == 0) CHOLMOD(dump_sparse) (F, "F ready", Common)) ;

    /* inputs, can be modified on output: */
    Lp = L->p ;		/* size n+1 */
    ASSERT (Lp != NULL) ;

    /* outputs, contents defined on input for incremental case only: */
    Lnz = L->nz ;	/* size n */
    Lnext = L->next ;	/* size n+2 */
    Li = L->i ;		/* size L->nzmax, can change in size */
    Lx = L->x ;		/* size L->nzmax or 2*L->nzmax, can change in size */
    Lz = L->z ;		/* size L->nzmax for zomplex case, can change in size */
    nzmax = L->nzmax ;
    ASSERT (Lnz != NULL && Li != NULL && Lx != NULL) ;

    /* ---------------------------------------------------------------------- */
    /* get workspace */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Stack = Iwork ;		/* size n (i/i/l) */
    Flag = Common->Flag ;	/* size n, Flag [i] < mark must hold */
    Wx = Common->Xwork ;	/* size n if real, 2*n if complex or 
				 * zomplex.  Xwork [i] == 0 must hold. */
    Wz = Wx + n ;		/* size n for zomplex case only */
    mark = Common->mark ;
    ASSERT ((Int) Common->xworksize >= (L->xtype == CHOLMOD_REAL ? 1:2)*n) ;

    /* ---------------------------------------------------------------------- */
    /* compute LDL' or LL' factorization by rows */
    /* ---------------------------------------------------------------------- */

#ifdef MASK
#define NEXT(k) k = RLinkUp [k]
#else
#define NEXT(k) k++
#endif

    for (k = kstart ; k < ((Int) kend) ; NEXT(k))
    {
	PRINT1 (("\n===============K "ID" Lnz [k] "ID"\n", k, Lnz [k])) ;

	/* ------------------------------------------------------------------ */
	/* compute pattern of kth row of L and scatter kth input column */
	/* ------------------------------------------------------------------ */

	/* column k of L is currently empty */
	ASSERT (Lnz [k] == 1) ;
	ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, 2*n, Common)) ;

	top = n ;		/* Stack is empty */
	Flag [k] = mark ;	/* do not include diagonal entry in Stack */

	/* use Li [Lp [i]+1] for etree */
#define PARENT(i) (Lnz [i] > 1) ? (Li [Lp [i] + 1]) : EMPTY

	if (stype > 0)
	{
	    /* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	    p = Ap [k] ;
	    pend = (packed) ? (Ap [k+1]) : (p + Anz [k]) ;
	    /* W [i] = Ax [i] ; scatter column of A */
#define SCATTER ASSIGN(Wx,Wz,i, Ax,Az,p)
	    SUBTREE ;
#undef SCATTER
	}
	else
	{
	    /* scatter kth col of triu (beta*I+AA'), get pattern L(k,:) */
	    pf = Fp [k] ;
	    pfend = (Fpacked) ? (Fp [k+1]) : (pf + Fnz [k]) ;
	    for ( ; pf < pfend ; pf++)
	    {
		/* get nonzero entry F (t,k) */
		t = Fi [pf] ;
		/* fk = Fx [pf] */
		ASSIGN (fx, fz, 0, Fx, Fz, pf) ;
		p = Ap [t] ;
		pend = (packed) ? (Ap [t+1]) : (p + Anz [t]) ;
		multadds = 0 ;
		/* W [i] += Ax [p] * fx ; scatter column of A*A' */
#define SCATTER MULTADD (Wx,Wz,i, Ax,Az,p, fx,fz,0) ; multadds++  ;
		SUBTREE ;
#undef SCATTER
#ifdef REAL
		fl += 2 * ((double) multadds) ;
#else
		fl += 8 * ((double) multadds) ;
#endif
	    }
	}

#undef PARENT

	/* ------------------------------------------------------------------ */
	/* if mask is present, set the corresponding entries in W to zero */
	/* ------------------------------------------------------------------ */

#ifdef MASK
        /* remove the dead element of Wx */
        if (mask != NULL)
        {

#if 0
	    /* older version */
            for (p = n; p > top;)
            {
                i = Stack [--p] ;
                if ( mask [i] >= 0 )
		{
		    CLEAR (Wx,Wz,i) ;	/* set W(i) to zero */
		}
            }
#endif

            for (s = top ; s < n ; s++)
            {
                i = Stack [s] ;
                if (mask [i] >= 0)
		{
		    CLEAR (Wx,Wz,i) ;	/* set W(i) to zero */
		}
            }

        }
#endif

	/* nonzero pattern of kth row of L is now in Stack [top..n-1].
	 * Flag [Stack [top..n-1]] is equal to mark, but no longer needed */

	/* mark = CHOLMOD(clear_flag) (Common) ; */
	CHOLMOD_CLEAR_FLAG (Common) ;
	mark = Common->mark ;

	/* ------------------------------------------------------------------ */
	/* compute kth row of L and store in column form */
	/* ------------------------------------------------------------------ */

	/* Solve L (0:k-1, 0:k-1) * y (0:k-1) = b (0:k-1) where
	 * b (0:k) = A (0:k,k) or A(0:k,:) * F(:,k) is in W and Stack.
	 *
	 * For LDL' factorization:
	 * L (k, 0:k-1) = y (0:k-1) ./ D (0:k-1)
	 * D (k) = b (k) - L (k, 0:k-1) * y (0:k-1)
	 *
	 * For LL' factorization:
	 * L (k, 0:k-1) = y (0:k-1)
	 * L (k,k) = sqrt (b (k) - L (k, 0:k-1) * L (0:k-1, k))
	 */

	/* dk = W [k] + beta */
	ADD_REAL (dk,0, Wx,k, beta,0) ;

#ifndef REAL
	/* In the unsymmetric case, the imaginary part of W[k] must be real,
	 * since F is assumed to be the complex conjugate transpose of A.  In
	 * the symmetric case, W[k] is the diagonal of A.  If the imaginary part
	 * of W[k] is nonzero, then the Cholesky factorization cannot be
	 * computed; A is not positive definite */
	dk_imaginary = (stype > 0) ? (IMAG_IS_NONZERO (Wx,Wz,k)) : FALSE ;
#endif

	/* W [k] = 0.0 ; */
	CLEAR (Wx,Wz,k) ;

	for (s = top ; s < n ; s++)
	{
	    /* get i for each nonzero entry L(k,i) */
	    i = Stack [s] ;

	    /* y = W [i] ; */
	    ASSIGN (yx,yz,0, Wx,Wz,i) ;

	    /* W [i] = 0.0 ; */
	    CLEAR (Wx,Wz,i) ;

	    lnz = Lnz [i] ;
	    p = Lp [i] ;
	    ASSERT (lnz > 0 && Li [p] == i) ;
	    pend = p + lnz ;

	    /* di = Lx [p] ; the diagonal entry L or D(i,i), which is real */
	    ASSIGN_REAL (di,0, Lx,p) ;

	    if (i >= (Int) L->minor || IS_ZERO (di [0]))
	    {
		/* For the LL' factorization, L(i,i) is zero.  For the LDL',
		 * D(i,i) is zero.  Skip column i of L, and set L(k,i) = 0. */
		CLEAR (lx,lz,0) ;
		p = pend ;
	    }
	    else if (is_ll)
	    {
#ifdef REAL
		fl += 2 * ((double) (pend - p - 1)) + 3 ;
#else
		fl += 8 * ((double) (pend - p - 1)) + 6 ;
#endif
		/* forward solve using L (i:(k-1),i) */
		/* divide by L(i,i), which must be real and nonzero */
		/* y /= di [0] */
		DIV_REAL (yx,yz,0, yx,yz,0, di,0) ;
		for (p++ ; p < pend ; p++)
		{
		    /* W [Li [p]] -= Lx [p] * y ; */
		    MULTSUB (Wx,Wz,Li[p], Lx,Lz,p, yx,yz,0) ;
		}
		/* do not scale L; compute dot product for L(k,k) */
		/* L(k,i) = conj(y) ; */
		ASSIGN_CONJ (lx,lz,0, yx,yz,0) ;
		/* d -= conj(y) * y ; */
		LLDOT (dk,0, yx,yz,0) ;
	    }
	    else
	    {
#ifdef REAL
		fl += 2 * ((double) (pend - p - 1)) + 3 ;
#else
		fl += 8 * ((double) (pend - p - 1)) + 6 ;
#endif
		/* forward solve using D (i,i) and L ((i+1):(k-1),i) */
		for (p++ ; p < pend ; p++)
		{
		    /* W [Li [p]] -= Lx [p] * y ; */
		    MULTSUB (Wx,Wz,Li[p], Lx,Lz,p, yx,yz,0) ;
		}
		/* Scale L (k,0:k-1) for LDL' factorization, compute D (k,k)*/
#ifdef REAL
		/* L(k,i) = y/d */
		lx [0] = yx [0] / di [0] ;
		/* d -= L(k,i) * y */
		dk [0] -= lx [0] * yx [0] ;
#else
		/* L(k,i) = conj(y) ; */
		ASSIGN_CONJ (lx,lz,0, yx,yz,0) ;
		/* L(k,i) /= di ; */
		DIV_REAL (lx,lz,0, lx,lz,0, di,0) ;
		/* d -= conj(y) * y / di */
		LDLDOT (dk,0, yx,yz,0, di,0) ;
#endif
	    }

	    /* determine if column i of L can hold the new L(k,i) entry */
	    if (p >= Lp [Lnext [i]])
	    {
		/* column i needs to grow */
		PRINT1 (("Factor Colrealloc "ID", old Lnz "ID"\n", i, Lnz [i]));
		if (!CHOLMOD(reallocate_column) (i, lnz + 1, L, Common))
		{
		    /* out of memory, L is now simplicial symbolic */
		    for (i = 0 ; i < n ; i++)
		    {
			/* W [i] = 0 ; */
			CLEAR (Wx,Wz,i) ;
		    }
		    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, n, Common)) ;
		    return (FALSE) ;
		}
		Li = L->i ;		/* L->i, L->x, L->z may have moved */
		Lx = L->x ;
		Lz = L->z ;
		p = Lp [i] + lnz ;	/* contents of L->p changed */
		ASSERT (p < Lp [Lnext [i]]) ;
	    }

	    /* store L (k,i) in the column form matrix of L */
	    Li [p] = k ;
	    /* Lx [p] = L(k,i) ; */
	    ASSIGN (Lx,Lz,p, lx,lz,0) ;
	    Lnz [i]++ ;
	}

	/* ------------------------------------------------------------------ */
	/* ensure abs (d) >= dbound if dbound is given, and store it in L */
	/* ------------------------------------------------------------------ */

	p = Lp [k] ;
	Li [p] = k ;

	if (k >= (Int) L->minor)
	{
	    /* the matrix is already not positive definite */
	    dk [0] = 0 ;
	}
	else if (use_dbound)
	{
	    /* modify the diagonal to force LL' or LDL' to exist */
	    dk [0] = CHOLMOD(dbound) (is_ll ? fabs (dk [0]) : dk [0], Common) ;
	}
	else if ((is_ll ? (IS_LE_ZERO (dk [0])) : (IS_ZERO (dk [0])))
#ifndef REAL
		|| dk_imaginary
#endif
		)
	{
	    /* the matrix has just been found to be not positive definite */
	    dk [0] = 0 ;
	    L->minor = k ;
	    ERROR (CHOLMOD_NOT_POSDEF, "not positive definite") ;
	}

	if (is_ll)
	{
	    /* this is counted as one flop, below */
	    dk [0] = sqrt (dk [0]) ;
	}

	/* Lx [p] = D(k,k) = d ; real part only */
	ASSIGN_REAL (Lx,p, dk,0) ;
	CLEAR_IMAG (Lx,Lz,p) ;
    }

#undef NEXT

    if (is_ll) fl += MAX ((Int) kend - (Int) kstart, 0) ;   /* count sqrt's */
    Common->rowfacfl = fl ;

    DEBUG (CHOLMOD(dump_factor) (L, "final cholmod_rowfac", Common)) ;
    ASSERT (CHOLMOD(dump_work) (TRUE, TRUE, n, Common)) ;
    return (TRUE) ;
}
#undef PATTERN
#undef REAL
#undef COMPLEX
#undef ZOMPLEX
