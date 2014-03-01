/* ========================================================================== */
/* === Core/cholmod_copy ==================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* C = A, which allocates C and copies A into C, with possible change of
 * stype.  The diagonal can optionally be removed.  The numerical entries
 * can optionally be copied.  This routine differs from cholmod_copy_sparse,
 * which makes an exact copy of a sparse matrix.
 *
 * A can be of any type (packed/unpacked, upper/lower/unsymmetric).  C is
 * packed and can be of any stype (upper/lower/unsymmetric), except that if
 * A is rectangular C can only be unsymmetric.  If the stype of A and C
 * differ, then the appropriate conversion is made.
 *
 * Symmetry of A (A->stype):
 * <0: lower: assume A is symmetric with just tril(A); the rest of A is ignored
 *  0  unsym: assume A is unsymmetric; consider all entries in A
 * >0  upper: assume A is symmetric with just triu(A); the rest of A is ignored
 *
 * Symmetry of C (stype parameter):
 * <0  lower: return just tril(C)
 *  0  unsym: return all of C
 * >0  upper: return just triu(C)
 *
 * In MATLAB:					    Using cholmod_copy:
 * ----------					    ----------------------------
 * C = A ;					    A unsymmetric, C unsymmetric
 * C = tril (A) ;				    A unsymmetric, C lower
 * C = triu (A) ;				    A unsymmetric, C upper
 * U = triu (A) ; L = tril (U',-1) ; C = L+U ;	    A upper, C unsymmetric
 * C = triu (A)' ;				    A upper, C lower
 * C = triu (A) ;				    A upper, C upper
 * L = tril (A) ; U = triu (L',1) ; C = L+U ;	    A lower, C unsymmetric
 * C = tril (A) ;				    A lower, C lower
 * C = tril (A)' ;				    A lower, C upper
 *
 * workspace: Iwork (max (nrow,ncol))
 *
 * A can have an xtype of pattern or real.  Complex and zomplex cases only
 * supported when mode <= 0 (in which case the numerical values are ignored).
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"


/* ========================================================================== */
/* === copy_sym_to_unsym ==================================================== */
/* ========================================================================== */

/* Construct an unsymmetric copy of a symmetric sparse matrix.  This does the
 * work for as C = cholmod_copy (A, 0, mode, Common) when A is symmetric.
 * In this case, extra space can be added to C.
 */

static cholmod_sparse *copy_sym_to_unsym
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to copy */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag)
			 * -2: pattern only, no diagonal, add 50% + n extra
			 * space to C */
    /* --------------- */
    cholmod_common *Common
)
{
    double aij ;
    double *Ax, *Cx ;
    Int *Ap, *Ai, *Anz, *Cp, *Ci, *Wj, *Iwork ;
    cholmod_sparse *C ;
    Int nrow, ncol, nz, packed, j, p, pend, i, pc, up, lo, values, diag,
	astype, extra ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = A->nrow ;
    ncol = A->ncol ;
    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    packed = A->packed ;
    values = (mode > 0) && (A->xtype != CHOLMOD_PATTERN) ;
    diag = (mode >= 0) ;

    astype = SIGN (A->stype) ;
    up = (astype > 0) ;
    lo = (astype < 0) ;
    ASSERT (astype != 0) ;

    /* ---------------------------------------------------------------------- */
    /* create an unsymmetric copy of a symmetric matrix */
    /* ---------------------------------------------------------------------- */

    Iwork = Common->Iwork ;
    Wj = Iwork ;		    /* size ncol (i/i/l) */

    /* In MATLAB notation, for converting a symmetric/upper matrix:
     *	U = triu (A) ;
     *	L = tril (U',-1) ;
     *	C = L + U ;
     *
     * For converting a symmetric/lower matrix to unsymmetric:
     *	L = tril (A) ;
     *	U = triu (L',1) ;
     *	C = L + U ;
     */
    ASSERT (up || lo) ;
    PRINT1 (("copy: convert symmetric to unsym\n")) ;

    /* count the number of entries in each column of C */
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = 0 ;
    }
    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if (i == j)
	    {
		/* the diagonal entry A(i,i) will appear just once
		 * (unless it is excluded with mode < 0) */
		if (diag)
		{
		    Wj [j]++ ;
		}
	    }
	    else if ((up && i < j) || (lo && i > j))
	    {
		/* upper case:  A(i,j) is in the strictly upper part;
		 * A(j,i) will be added to the strictly lower part of C.
		 * lower case is the opposite. */
		Wj [j]++ ;
		Wj [i]++ ;
	    }
	}
    }
    nz = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	nz += Wj [j] ;
    }

    extra = (mode == -2) ? (nz/2 + ncol) : 0 ;

    /* allocate C.  C is sorted if and only if A is sorted */
    C = CHOLMOD(allocate_sparse) (nrow, ncol, nz + extra, A->sorted, TRUE, 0,
	    values ? A->xtype : CHOLMOD_PATTERN, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* construct the column pointers for C */
    p = 0 ;
    for (j = 0 ; j < ncol ; j++)
    {
	Cp [j] = p ;
	p += Wj [j] ;
    }
    Cp [ncol] = p ;
    for (j = 0 ; j < ncol ; j++)
    {
	Wj [j] = Cp [j] ;
    }

    /* construct C */
    if (values)
    {

	/* pattern and values */
	ASSERT (diag) ;
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		aij = Ax [p] ;
		if (i == j)
		{
		    /* add diagonal entry A(i,i) to column i */
		    pc = Wj [i]++ ;
		    Ci [pc] = i ;
		    Cx [pc] = aij ;
		}
		else if ((up && i < j) || (lo && i > j))
		{
		    /* add A(i,j) to column j */
		    pc = Wj [j]++ ;
		    Ci [pc] = i ;
		    Cx [pc] = aij ;
		    /* add A(j,i) to column i */
		    pc = Wj [i]++ ;
		    Ci [pc] = j ;
		    Cx [pc] = aij ;
		}
	    }
	}

    }
    else
    {

	/* pattern only, possibly excluding the diagonal */
	for (j = 0 ; j < ncol ; j++)
	{
	    p = Ap [j] ;
	    pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	    for ( ; p < pend ; p++)
	    {
		i = Ai [p] ;
		if (i == j)
		{
		    /* add diagonal entry A(i,i) to column i
		     * (unless it is excluded with mode < 0) */
		    if (diag)
		    {
			Ci [Wj [i]++] = i ;
		    }
		}
		else if ((up && i < j) || (lo && i > j))
		{
		    /* add A(i,j) to column j */
		    Ci [Wj [j]++] = i ;
		    /* add A(j,i) to column i */
		    Ci [Wj [i]++] = j ;
		}
	    }
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    DEBUG (i = CHOLMOD(dump_sparse) (C, "copy_sym_to_unsym", Common)) ;
    PRINT1 (("mode %d nnzdiag "ID"\n", mode, i)) ;
    ASSERT (IMPLIES (mode < 0, i == 0)) ;
    return (C) ;
}


/* ========================================================================== */
/* === cholmod_copy ========================================================= */
/* ========================================================================== */

cholmod_sparse *CHOLMOD(copy)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to copy */
    int stype,		/* requested stype of C */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *C ;
    Int nrow, ncol, up, lo, values, diag, astype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    values = (mode > 0) && (A->xtype != CHOLMOD_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN,
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if ((stype || A->stype) && nrow != ncol)
    {
	/* inputs invalid */
	ERROR (CHOLMOD_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(allocate_work) (0, MAX (nrow,ncol), 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    diag = (mode >= 0) ;
    astype = SIGN (A->stype) ;
    stype = SIGN (stype) ;
    up = (astype > 0) ;
    lo = (astype < 0) ;

    /* ---------------------------------------------------------------------- */
    /* copy the matrix */
    /* ---------------------------------------------------------------------- */

    if (astype == stype)
    {

	/* ------------------------------------------------------------------ */
	/* symmetry of A and C are the same */
	/* ------------------------------------------------------------------ */

	/* copy A into C, keeping the same symmetry.  If A is symmetric
	 * entries in the ignored part of A are not copied into C */
	C = CHOLMOD(band) (A, -nrow, ncol, mode, Common) ;

    }
    else if (!astype)
    {

	/* ------------------------------------------------------------------ */
	/* convert unsymmetric matrix A into a symmetric matrix C */
	/* ------------------------------------------------------------------ */

	if (stype > 0)
	{
	    /* C = triu (A) */
	    C = CHOLMOD(band) (A, 0, ncol, mode, Common) ;
	}
	else
	{
	    /* C = tril (A) */
	    C = CHOLMOD(band) (A, -nrow, 0, mode, Common) ;
	}
	if (Common->status < CHOLMOD_OK)
	{
	    /* out of memory */
	    return (NULL) ;
	}
	C->stype = stype ;

    }
    else if (astype == -stype)
    {

	/* ------------------------------------------------------------------ */
	/* transpose a symmetric matrix */
	/* ------------------------------------------------------------------ */

	/* converting upper to lower or lower to upper */
	/* workspace: Iwork (nrow) */
	C = CHOLMOD(transpose) (A, values, Common) ;
	if (!diag)
	{
	    /* remove diagonal, if requested */
	    CHOLMOD(band_inplace) (-nrow, ncol, -1, C, Common) ;
	}

    }
    else
    {

	/* ------------------------------------------------------------------ */
	/* create an unsymmetric copy of a symmetric matrix */
	/* ------------------------------------------------------------------ */

	C = copy_sym_to_unsym (A, mode, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return if error */
    /* ---------------------------------------------------------------------- */

    if (Common->status < CHOLMOD_OK)
    {
	/* out of memory */
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    DEBUG (diag = CHOLMOD(dump_sparse) (C, "copy", Common)) ;
    PRINT1 (("mode %d nnzdiag "ID"\n", mode, diag)) ;
    ASSERT (IMPLIES (mode < 0, diag == 0)) ;
    return (C) ;
}
