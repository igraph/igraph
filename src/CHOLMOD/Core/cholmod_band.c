/* ========================================================================== */
/* === Core/cholmod_band ==================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* C = tril (triu (A,k1), k2)
 *
 * C is a matrix consisting of the diagonals of A from k1 to k2.
 *
 * k=0 is the main diagonal of A, k=1 is the superdiagonal, k=-1 is the
 * subdiagonal, and so on.  If A is m-by-n, then:
 *
 *	k1=-m		    C = tril (A,k2)
 *	k2=n		    C = triu (A,k1)
 *	k1=0 and k2=0	    C = diag(A), except C is a matrix, not a vector
 *
 * Values of k1 and k2 less than -m are treated as -m, and values greater
 * than n are treated as n.
 *
 * A can be of any symmetry (upper, lower, or unsymmetric); C is returned in
 * the same form, and packed.  If A->stype > 0, entries in the lower
 * triangular part of A are ignored, and the opposite is true if
 * A->stype < 0.  If A has sorted columns, then so does C.
 * C has the same size as A.
 *
 * If inplace is TRUE, then the matrix A is modified in place.
 * Only packed matrices can be converted in place.
 *
 * C can be returned as a numerical valued matrix (if A has numerical values
 * and mode > 0), as a pattern-only (mode == 0), or as a pattern-only but with
 * the diagonal entries removed (mode < 0).
 *
 * workspace: none
 *
 * A can have an xtype of pattern or real.  Complex and zomplex cases supported
 * only if mode <= 0 (in which case the numerical values are ignored).
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"

static cholmod_sparse *band		/* returns C, or NULL if failure */
(
    /* ---- input or in/out if inplace is TRUE --- */
    cholmod_sparse *A,
    /* ---- input ---- */
    SuiteSparse_long k1,    /* ignore entries below the k1-st diagonal */
    SuiteSparse_long k2,    /* ignore entries above the k2-nd diagonal */
    int mode,	    /* >0: numerical, 0: pattern, <0: pattern (no diagonal) */
    int inplace,    /* if TRUE, then convert A in place */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Cx ;
    Int packed, nz, j, p, pend, i, ncol, nrow, jlo, jhi, ilo, ihi, sorted,
	values, diag ;
    Int *Ap, *Anz, *Ai, *Cp, *Ci ;
    cholmod_sparse *C ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    values = (mode > 0) && (A->xtype != CHOLMOD_PATTERN) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN,
	    values ? CHOLMOD_REAL : CHOLMOD_ZOMPLEX, NULL) ;
    packed = A->packed ;
    diag = (mode >= 0) ;
    if (inplace && !packed)
    {
	/* cannot operate on an unpacked matrix in place */
	ERROR (CHOLMOD_INVALID, "cannot operate on unpacked matrix in-place") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */


    PRINT1 (("k1 %ld k2 %ld\n", k1, k2)) ;
    Ap  = A->p ;
    Anz = A->nz ;
    Ai  = A->i ;
    Ax  = A->x ;
    sorted = A->sorted ;


    if (A->stype > 0)
    {
	/* ignore any entries in strictly lower triangular part of A */
	k1 = MAX (k1, 0) ;
    }
    if (A->stype < 0)
    {
	/* ignore any entries in strictly upper triangular part of A */
	k2 = MIN (k2, 0) ;
    }
    ncol = A->ncol ;
    nrow = A->nrow ;

    /* ensure k1 and k2 are in the range -nrow to +ncol to
     * avoid possible integer overflow if k1 and k2 are huge */
    k1 = MAX (-nrow, k1) ;
    k1 = MIN (k1, ncol) ;
    k2 = MAX (-nrow, k2) ;
    k2 = MIN (k2, ncol) ;

    /* consider columns jlo to jhi.  columns outside this range are empty */
    jlo = MAX (k1, 0) ;
    jhi = MIN (k2+nrow, ncol) ;

    if (k1 > k2)
    {
	/* nothing to do */
	jlo = ncol ;
	jhi = ncol ;
    }

    /* ---------------------------------------------------------------------- */
    /* allocate C, or operate on A in place */
    /* ---------------------------------------------------------------------- */

    if (inplace)
    {
	/* convert A in place */
	C = A ;
    }
    else
    {
	/* count the number of entries in the result C */
	nz = 0 ;
	if (sorted)
	{
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i > ihi)
		    {
			break ;
		    }
		    if (i >= ilo && (diag || i != j))
		    {
			nz++ ;
		    }
		}
	    }
	}
	else
	{
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= ilo && i <= ihi && (diag || i != j))
		    {
			nz++ ;
		    }
		}
	    }
	}
	/* allocate C; A will not be modified.  C is sorted if A is sorted */
	C = CHOLMOD(allocate_sparse) (A->nrow, ncol, nz, sorted, TRUE,
		A->stype, values ? A->xtype : CHOLMOD_PATTERN, Common) ;
	if (Common->status < CHOLMOD_OK)
	{
	    return (NULL) ;	/* out of memory */
	}
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;

    /* ---------------------------------------------------------------------- */
    /* construct C */
    /* ---------------------------------------------------------------------- */

    /* columns 0 to jlo-1 are empty */
    for (j = 0 ; j < jlo ; j++)
    {
	Cp [j] = 0 ;
    }

    nz = 0 ;
    if (sorted)
    {
	if (values)
	{
	    /* pattern and values */
	    ASSERT (diag) ;
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i > ihi)
		    {
			break ;
		    }
		    if (i >= ilo)
		    {
			Ci [nz] = i ;
			Cx [nz] = Ax [p] ;
			nz++ ;
		    }
		}
	    }
	}
	else
	{
	    /* pattern only, perhaps with no diagonal */
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i > ihi)
		    {
			break ;
		    }
		    if (i >= ilo && (diag || i != j))
		    {
			Ci [nz++] = i ;
		    }
		}
	    }
	}
    }
    else
    {
	if (values)
	{
	    /* pattern and values */
	    ASSERT (diag) ;
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= ilo && i <= ihi)
		    {
			Ci [nz] = i ;
			Cx [nz] = Ax [p] ;
			nz++ ;
		    }
		}
	    }
	}
	else
	{
	    /* pattern only, perhaps with no diagonal */
	    for (j = jlo ; j < jhi ; j++)
	    {
		ilo = j-k2 ;
		ihi = j-k1 ;
		p = Ap [j] ;
		pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
		Cp [j] = nz ;
		for ( ; p < pend ; p++)
		{
		    i = Ai [p] ;
		    if (i >= ilo && i <= ihi && (diag || i != j))
		    {
			Ci [nz++] = i ;
		    }
		}
	    }
	}
    }

    /* columns jhi to ncol-1 are empty */
    for (j = jhi ; j <= ncol ; j++)
    {
	Cp [j] = nz ;
    }

    /* ---------------------------------------------------------------------- */
    /* reduce A in size if done in place */
    /* ---------------------------------------------------------------------- */

    if (inplace)
    {
	/* free the unused parts of A, and reduce A->i and A->x in size */
	ASSERT (MAX (1,nz) <= A->nzmax) ;
	CHOLMOD(reallocate_sparse) (nz, A, Common) ;
	ASSERT (Common->status >= CHOLMOD_OK) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return the result C */
    /* ---------------------------------------------------------------------- */

    DEBUG (i = CHOLMOD(dump_sparse) (C, "band", Common)) ;
    ASSERT (IMPLIES (mode < 0, i == 0)) ;
    return (C) ;
}


/* ========================================================================== */
/* === cholmod_band ========================================================= */
/* ========================================================================== */

cholmod_sparse *CHOLMOD(band)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to extract band matrix from */
    SuiteSparse_long k1,    /* ignore entries below the k1-st diagonal */
    SuiteSparse_long k2,    /* ignore entries above the k2-nd diagonal */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    /* --------------- */
    cholmod_common *Common
)
{
    return (band (A, k1, k2, mode, FALSE, Common)) ;
}


/* ========================================================================== */
/* === cholmod_band_inplace ================================================= */
/* ========================================================================== */

int CHOLMOD(band_inplace)
(
    /* ---- input ---- */
    SuiteSparse_long k1,    /* ignore entries below the k1-st diagonal */
    SuiteSparse_long k2,    /* ignore entries above the k2-nd diagonal */
    int mode,		/* >0: numerical, 0: pattern, <0: pattern (no diag) */
    /* ---- in/out --- */
    cholmod_sparse *A,	/* matrix from which entries not in band are removed */
    /* --------------- */
    cholmod_common *Common
)
{
    return (band (A, k1, k2, mode, TRUE, Common) != NULL) ;
}
