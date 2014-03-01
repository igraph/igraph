/* ========================================================================== */
/* === Core/cholmod_sparse ================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Core utility routines for the cholmod_sparse object:
 *
 * A sparse matrix is held in compressed column form.  In the basic type
 * ("packed", which corresponds to a MATLAB sparse matrix), an n-by-n matrix
 * with nz entries is held in three arrays: p of size n+1, i of size nz, and x
 * of size nz.  Row indices of column j are held in i [p [j] ... p [j+1]-1] and
 * in the same locations in x.  There may be no duplicate entries in a column.
 * Row indices in each column may be sorted or unsorted (CHOLMOD keeps track).
 *
 * Primary routines:
 * -----------------
 * cholmod_allocate_sparse	allocate a sparse matrix
 * cholmod_free_sparse		free a sparse matrix
 *
 * Secondary routines:
 * -------------------
 * cholmod_reallocate_sparse	change the size (# entries) of sparse matrix
 * cholmod_nnz			number of nonzeros in a sparse matrix
 * cholmod_speye		sparse identity matrix
 * cholmod_spzeros		sparse zero matrix
 * cholmod_copy_sparse		create a copy of a sparse matrix
 *
 * All xtypes are supported (pattern, real, complex, and zomplex)
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"


/* ========================================================================== */
/* === cholmod_allocate_sparse ============================================== */
/* ========================================================================== */

/* Allocate space for a matrix.  A->i and A->x are not initialized.  A->p
 * (and A->nz if A is not packed) are set to zero, so a matrix containing no
 * entries (all zero) is returned.  See also cholmod_spzeros.
 *
 * workspace: none
 */

cholmod_sparse *CHOLMOD(allocate_sparse)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of A */
    size_t ncol,	/* # of columns of A */
    size_t nzmax,	/* max # of nonzeros of A */
    int sorted,		/* TRUE if columns of A sorted, FALSE otherwise */
    int packed,		/* TRUE if A will be packed, FALSE otherwise */
    int stype,		/* stype of A */
    int xtype,		/* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *A ;
    Int *Ap, *Anz ;
    size_t nzmax0 ;
    Int j ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    if (stype != 0 && nrow != ncol)
    {
	ERROR (CHOLMOD_INVALID, "rectangular matrix with stype != 0 invalid") ;
	return (NULL) ;
    }
    if (xtype < CHOLMOD_PATTERN || xtype > CHOLMOD_ZOMPLEX)
    {
	ERROR (CHOLMOD_INVALID, "xtype invalid") ;
	return (NULL) ;
    }
    /* ensure the dimensions do not cause integer overflow */
    (void) CHOLMOD(add_size_t) (ncol, 2, &ok) ;
    if (!ok || nrow > Int_max || ncol > Int_max || nzmax > Int_max)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate header */
    /* ---------------------------------------------------------------------- */

    A = CHOLMOD(malloc) (sizeof (cholmod_sparse), 1, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }
    PRINT1 (("cholmod_allocate_sparse %d-by-%d nzmax %d sorted %d packed %d"
		" xtype %d\n", nrow, ncol, nzmax, sorted, packed, xtype)) ;

    nzmax = MAX (1, nzmax) ;

    A->nrow = nrow ;
    A->ncol = ncol ;
    A->nzmax = nzmax ;
    A->packed = packed ;    /* default is packed (A->nz not present) */
    A->stype = stype ;
    A->itype = ITYPE ;
    A->xtype = xtype ;
    A->dtype = DTYPE ;

    A->nz = NULL ;
    A->p = NULL ;
    A->i = NULL ;
    A->x = NULL ;
    A->z = NULL ;

    /* A 1-by-m matrix always has sorted columns */
    A->sorted = (nrow <= 1) ? TRUE : sorted ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix itself */
    /* ---------------------------------------------------------------------- */

    /* allocate O(ncol) space */
    A->p = CHOLMOD(malloc) (((size_t) ncol)+1, sizeof (Int), Common) ;
    if (!packed)
    {
	A->nz = CHOLMOD(malloc) (ncol, sizeof (Int), Common) ;
    }

    /* allocate O(nz) space */
    nzmax0 = 0 ;
    CHOLMOD(realloc_multiple) (nzmax, 1, xtype, &(A->i), NULL, &(A->x), &(A->z),
	    &nzmax0, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free_sparse) (&A, Common) ;
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* initialize A->p and A->nz so that A is an empty matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    for (j = 0 ; j <= (Int) ncol ; j++)
    {
	Ap [j] = 0 ;
    }
    if (!packed)
    {
	Anz = A->nz ;
	for (j = 0 ; j < (Int) ncol ; j++)
	{
	    Anz [j] = 0 ;
	}
    }
    return (A) ;
}


/* ========================================================================== */
/* === cholmod_free_sparse ================================================== */
/* ========================================================================== */

/* free a sparse matrix
 *
 * workspace: none
 */

int CHOLMOD(free_sparse)
(
    /* ---- in/out --- */
    cholmod_sparse **AHandle,	/* matrix to deallocate, NULL on output */
    /* --------------- */
    cholmod_common *Common
)
{
    Int n, nz ;
    cholmod_sparse *A ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (AHandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    A = *AHandle ;
    if (A == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    n = A->ncol ;
    nz = A->nzmax ;
    A->p  = CHOLMOD(free) (n+1, sizeof (Int), A->p,  Common) ;
    A->i  = CHOLMOD(free) (nz,  sizeof (Int), A->i,  Common) ;
    A->nz = CHOLMOD(free) (n,   sizeof (Int), A->nz, Common) ;

    switch (A->xtype)
    {
	case CHOLMOD_REAL:
	    A->x = CHOLMOD(free) (nz, sizeof (double), A->x,  Common) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    A->x = CHOLMOD(free) (nz, 2*sizeof (double), A->x,  Common) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    A->x = CHOLMOD(free) (nz, sizeof (double), A->x,  Common) ;
	    A->z = CHOLMOD(free) (nz, sizeof (double), A->z,  Common) ;
	    break ;
    }

    *AHandle = CHOLMOD(free) (1, sizeof (cholmod_sparse), (*AHandle), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_sparse ============================================ */
/* ========================================================================== */

/* Change the size of A->i, A->x, and A->z, or allocate them if their current
 * size is zero.  A->x and A->z are not modified if A->xtype is CHOLMOD_PATTERN.
 * A->z is not modified unless A->xtype is CHOLMOD_ZOMPLEX.
 * 
 * workspace: none
 */

int CHOLMOD(reallocate_sparse)
(
    /* ---- input ---- */
    size_t nznew,	/* new # of entries in A */
    /* ---- in/out --- */
    cholmod_sparse *A,	/* matrix to reallocate */
    /* --------------- */
    cholmod_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (A, FALSE) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;
    PRINT1 (("realloc matrix %d to %d, xtype: %d\n",
		A->nzmax, nznew, A->xtype)) ;

    /* ---------------------------------------------------------------------- */
    /* resize the matrix */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(realloc_multiple) (MAX (1,nznew), 1, A->xtype, &(A->i), NULL,
	    &(A->x), &(A->z), &(A->nzmax), Common) ;

    return (Common->status == CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_speye ======================================================== */
/* ========================================================================== */

/* Return a sparse identity matrix. */

cholmod_sparse *CHOLMOD(speye)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of A */
    size_t ncol,	/* # of columns of A */
    int xtype,		/* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Az ;
    cholmod_sparse *A ;
    Int *Ap, *Ai ;
    Int j, n ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix */
    /* ---------------------------------------------------------------------- */

    n = MIN (nrow, ncol) ;
    A = CHOLMOD(allocate_sparse) (nrow, ncol, n, TRUE, TRUE, 0, xtype,
	    Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory or inputs invalid */
    }

    /* ---------------------------------------------------------------------- */
    /* create the identity matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;

    for (j = 0 ; j < n ; j++)
    {
	Ap [j] = j ;
    }
    for (j = n ; j <= ((Int) ncol) ; j++)
    {
	Ap [j] = n ;
    }
    for (j = 0 ; j < n ; j++)
    {
	Ai [j] = j ;
    }

    switch (xtype)
    {
	case CHOLMOD_REAL:
	    for (j = 0 ; j < n ; j++)
	    {
		Ax [j] = 1 ;
	    }
	    break ;

	case CHOLMOD_COMPLEX:
	    for (j = 0 ; j < n ; j++)
	    {
		Ax [2*j  ] = 1 ;
		Ax [2*j+1] = 0 ;
	    }
	    break ;

	case CHOLMOD_ZOMPLEX:
	    for (j = 0 ; j < n ; j++)
	    {
		Ax [j] = 1 ;
	    }
	    for (j = 0 ; j < n ; j++)
	    {
		Az [j] = 0 ;
	    }
	    break ;
    }

    return (A) ;
}


/* ========================================================================== */
/* === cholmod_spzeros ====================================================== */
/* ========================================================================== */

/* Return a sparse zero matrix. */

cholmod_sparse *CHOLMOD(spzeros)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of A */
    size_t ncol,	/* # of columns of A */
    size_t nzmax,	/* max # of nonzeros of A */
    int xtype,		/* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix */
    /* ---------------------------------------------------------------------- */

    return (CHOLMOD(allocate_sparse) (nrow, ncol, nzmax, TRUE, TRUE, 0, xtype,
	    Common)) ;
}


/* ========================================================================== */
/* === cholmod_nnz ========================================================== */
/* ========================================================================== */

/* Return the number of entries in a sparse matrix.
 *
 * workspace: none
 * integer overflow cannot occur, since the matrix is already allocated.
 */

SuiteSparse_long CHOLMOD(nnz)
(
    /* ---- input ---- */
    cholmod_sparse *A,
    /* --------------- */
    cholmod_common *Common
)
{
    Int *Ap, *Anz ;
    size_t nz ;
    Int j, ncol ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (EMPTY) ;
    RETURN_IF_NULL (A, EMPTY) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, EMPTY) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* return nnz (A) */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    if (A->packed)
    {
	Ap = A->p ;
	RETURN_IF_NULL (Ap, EMPTY) ;
	nz = Ap [ncol] ;
    }
    else
    {
	Anz = A->nz ;
	RETURN_IF_NULL (Anz, EMPTY) ;
	nz = 0 ;
	for (j = 0 ; j < ncol ; j++)
	{
	    nz += MAX (0, Anz [j]) ;
	}
    }
    return (nz) ;
}


/* ========================================================================== */
/* === cholmod_copy_sparse ================================================== */
/* ========================================================================== */

/* C = A.  Create an exact copy of a sparse matrix, with one exception.
 * Entries in unused space are not copied (they might not be initialized,
 * and copying them would cause program checkers such as purify and
 * valgrind to complain).  The xtype of the resulting matrix C is the same as
 * the xtype of the input matrix A.
 *
 * See also Core/cholmod_copy, which copies a matrix with possible changes
 * in stype, presence of diagonal entries, pattern vs. numerical values,
 * real and/or imaginary parts, and so on.
 */

cholmod_sparse *CHOLMOD(copy_sparse)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to copy */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Cx, *Az, *Cz ;
    Int *Ap, *Ai, *Anz, *Cp, *Ci, *Cnz ;
    cholmod_sparse *C ;
    Int p, pend, j, ncol, packed, nzmax, nz, xtype ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    if (A->stype != 0 && A->nrow != A->ncol)
    {
	ERROR (CHOLMOD_INVALID, "rectangular matrix with stype != 0 invalid") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    ASSERT (CHOLMOD(dump_sparse) (A, "A original", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    ncol = A->ncol ;
    nzmax = A->nzmax ;
    packed = A->packed ;
    Ap = A->p ;
    Ai = A->i ;
    Ax = A->x ;
    Az = A->z ;
    Anz = A->nz ;
    xtype = A->xtype ;

    /* ---------------------------------------------------------------------- */
    /* allocate the copy */
    /* ---------------------------------------------------------------------- */

    C = CHOLMOD(allocate_sparse) (A->nrow, A->ncol, A->nzmax, A->sorted,
	    A->packed, A->stype, A->xtype, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    Cp = C->p ;
    Ci = C->i ;
    Cx = C->x ;
    Cz = C->z ;
    Cnz = C->nz ;

    /* ---------------------------------------------------------------------- */
    /* copy the matrix */
    /* ---------------------------------------------------------------------- */

    for (j = 0 ; j <= ncol ; j++)
    {
	Cp [j] = Ap [j] ;
    }

    if (packed)
    {
	nz = Ap [ncol] ;
	for (p = 0 ; p < nz ; p++)
	{
	    Ci [p] = Ai [p] ;
	}

	switch (xtype)
	{
	    case CHOLMOD_REAL:
		for (p = 0 ; p < nz ; p++)
		{
		    Cx [p] = Ax [p] ;
		}
		break ;

	    case CHOLMOD_COMPLEX:
		for (p = 0 ; p < 2*nz ; p++)
		{
		    Cx [p] = Ax [p] ;
		}
		break ;

	    case CHOLMOD_ZOMPLEX:
		for (p = 0 ; p < nz ; p++)
		{
		    Cx [p] = Ax [p] ;
		    Cz [p] = Az [p] ;
		}
		break ;
	}

    }
    else
    {

	for (j = 0 ; j < ncol ; j++)
	{
	    Cnz [j] = Anz [j] ;
	}

	switch (xtype)
	{
	    case CHOLMOD_PATTERN:
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Ci [p] = Ai [p] ;
		    }
		}
		break ;

	    case CHOLMOD_REAL:
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Ci [p] = Ai [p] ;
			Cx [p] = Ax [p] ;
		    }
		}
		break ;

	    case CHOLMOD_COMPLEX:
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Ci [p] = Ai [p] ;
			Cx [2*p  ] = Ax [2*p  ] ;
			Cx [2*p+1] = Ax [2*p+1] ;
		    }
		}
		break ;

	    case CHOLMOD_ZOMPLEX:
		for (j = 0 ; j < ncol ; j++)
		{
		    p = Ap [j] ;
		    pend = p + Anz [j] ;
		    for ( ; p < pend ; p++)
		    {
			Ci [p] = Ai [p] ;
			Cx [p] = Ax [p] ;
			Cz [p] = Az [p] ;
		    }
		}
		break ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (C, "C copy", Common) >= 0) ;
    return (C) ;
}
