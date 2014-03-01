/* ========================================================================== */
/* === Core/cholmod_triplet ================================================= */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2006,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Core utility routines for the cholmod_triplet object:
 *
 * A sparse matrix held in triplet form is the simplest one for a user to
 * create.  It consists of a list of nz entries in arbitrary order, held in
 * three arrays: i, j, and x, each of length nk.  The kth entry is in row i[k],
 * column j[k], with value x[k].  There may be duplicate values; if A(i,j)
 * appears more than once, its value is the sum of the entries with those row
 * and column indices.
 *
 * Primary routines:
 * -----------------
 * cholmod_allocate_triplet	allocate a triplet matrix
 * cholmod_free_triplet		free a triplet matrix
 *
 * Secondary routines:
 * -------------------
 * cholmod_reallocate_triplet	reallocate a triplet matrix
 * cholmod_sparse_to_triplet	create a triplet matrix copy of a sparse matrix
 * cholmod_triplet_to_sparse	create a sparse matrix copy of a triplet matrix
 * cholmod_copy_triplet		create a copy of a triplet matrix
 *
 * The relationship between an m-by-n cholmod_sparse matrix A and a
 * cholmod_triplet matrix (i, j, and x) is identical to how they are used in
 * the MATLAB "sparse" and "find" functions:
 *
 *	[i j x] = find (A)
 *	[m n] = size (A)
 *	A = sparse (i,j,x,m,n)
 *
 * with the exception that the cholmod_sparse matrix may be "unpacked", may
 * have either sorted or unsorted columns (depending on the option selected),
 * and may be symmetric with just the upper or lower triangular part stored.
 * Likewise, the cholmod_triplet matrix may contain just the entries in the
 * upper or lower triangular part of a symmetric matrix.
 *
 * MATLAB sparse matrices are always "packed", always have sorted columns,
 * and always store both parts of a symmetric matrix.  In some cases, MATLAB
 * behaves like CHOLMOD by ignoring entries in the upper or lower triangular
 * part of a matrix that is otherwise assumed to be symmetric (such as the
 * input to chol).  In CHOLMOD, that option is a characteristic of the object.
 * In MATLAB, that option is based on how a matrix is used as the input to
 * a function.
 *
 * The triplet matrix is provided to give the user a simple way of constructing
 * a sparse matrix.  There are very few operations supported for triplet
 * matrices.  The assumption is that they will be converted to cholmod_sparse
 * matrix form first.
 *
 * Adding two triplet matrices simply involves concatenating the contents of
 * the three arrays (i, j, and x).   To permute a triplet matrix, just replace
 * the row and column indices with their permuted values.  For example, if
 * P is a permutation vector, then P [k] = j means row/column j is the kth
 * row/column in C=P*A*P'.  In MATLAB notation, C=A(p,p).  If Pinv is an array
 * of size n and T is the triplet form of A, then:
 *
 *	Ti = T->i ;
 *	Tj = T->j ;
 *	for (k = 0 ; k < n  ; k++) Pinv [P [k]] = k ;
 *	for (k = 0 ; k < nz ; k++) Ti [k] = Pinv [Ti [k]] ;
 *	for (k = 0 ; k < nz ; k++) Tj [k] = Pinv [Tj [k]] ;
 *
 * overwrites T with the triplet form of C=P*A*P'.  The conversion
 *
 *	C = cholmod_triplet_to_sparse (T, 0, &Common) ;
 *
 * will then return the matrix C = P*A*P'.
 *
 * Note that T->stype > 0 means that entries in the lower triangular part of
 * T are transposed into the upper triangular part when T is converted to
 * sparse matrix (cholmod_sparse) form with cholmod_triplet_to_sparse.  The
 * opposite is true for T->stype < 0.
 *
 * Since the triplet matrix T is so simple to generate, it's quite easy
 * to remove entries that you do not want, prior to converting T to the
 * cholmod_sparse form.  So if you include these entries in T, CHOLMOD
 * assumes that there must be a reason (such as the one above).  Thus,
 * no entry in a triplet matrix is ever ignored.
 *
 * Other operations, such as extacting a submatrix, horizontal and vertical
 * concatenation, multiply a triplet matrix times a dense matrix, are also
 * simple.  Multiplying two triplet matrices is not trivial; the simplest
 * method is to convert them to cholmod_sparse matrices first.
 *
 * Supports all xtypes (pattern, real, complex, and zomplex).
 */
 
#include "cholmod_internal.h"
#include "cholmod_core.h"


/* ========================================================================== */
/* === TEMPLATE ============================================================= */
/* ========================================================================== */

#define PATTERN
#include "t_cholmod_triplet.c"
#define REAL
#include "t_cholmod_triplet.c"
#define COMPLEX
#include "t_cholmod_triplet.c"
#define ZOMPLEX
#include "t_cholmod_triplet.c"


/* ========================================================================== */
/* === cholmod_allocate_triplet ============================================= */
/* ========================================================================== */

/* allocate space for a triplet matrix
 *
 * workspace: none
 */

cholmod_triplet *CHOLMOD(allocate_triplet)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of T */
    size_t ncol,	/* # of columns of T */
    size_t nzmax,	/* max # of nonzeros of T */
    int stype,		/* stype of T */
    int xtype,		/* CHOLMOD_PATTERN, _REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_triplet *T ;
    size_t nzmax0 ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
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

    T = CHOLMOD(malloc) (sizeof (cholmod_triplet), 1, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    PRINT1 (("cholmod_allocate_triplet %d-by-%d nzmax %d xtype %d\n",
		nrow, ncol, nzmax, xtype)) ;

    nzmax = MAX (1, nzmax) ;

    T->nrow = nrow ;
    T->ncol = ncol ;
    T->nzmax = nzmax ;
    T->nnz = 0 ;
    T->stype = stype ;
    T->itype = ITYPE ;
    T->xtype = xtype ;
    T->dtype = DTYPE ;

    T->j = NULL ;
    T->i = NULL ;
    T->x = NULL ;
    T->z = NULL ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix itself */
    /* ---------------------------------------------------------------------- */

    nzmax0 = 0 ;
    CHOLMOD(realloc_multiple) (nzmax, 2, xtype, &(T->i), &(T->j),
		&(T->x), &(T->z), &nzmax0, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free_triplet) (&T, Common) ;
	return (NULL) ;	    /* out of memory */
    }

    return (T) ;
}


/* ========================================================================== */
/* === cholmod_free_triplet ================================================= */
/* ========================================================================== */

/* free a triplet matrix
 *
 * workspace: none
 */

int CHOLMOD(free_triplet)
(
    /* ---- in/out --- */
    cholmod_triplet **THandle,    /* matrix to deallocate, NULL on output */
    /* --------------- */
    cholmod_common *Common
)
{
    Int nz ;
    cholmod_triplet *T ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (THandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    T = *THandle ;
    if (T == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    nz = T->nzmax ;
    T->j = CHOLMOD(free) (nz, sizeof (Int), T->j, Common) ;
    T->i = CHOLMOD(free) (nz, sizeof (Int), T->i, Common) ;
    if (T->xtype == CHOLMOD_REAL)
    {
	T->x = CHOLMOD(free) (nz, sizeof (double), T->x, Common) ;
    }
    else if (T->xtype == CHOLMOD_COMPLEX)
    {
	T->x = CHOLMOD(free) (nz, 2*sizeof (double), T->x, Common) ;
    }
    else if (T->xtype == CHOLMOD_ZOMPLEX)
    {
	T->x = CHOLMOD(free) (nz, sizeof (double), T->x, Common) ;
	T->z = CHOLMOD(free) (nz, sizeof (double), T->z, Common) ;
    }
    *THandle = CHOLMOD(free) (1, sizeof (cholmod_triplet), (*THandle), Common) ;
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_reallocate_triplet =========================================== */
/* ========================================================================== */

/* Change the size of T->i, T->j, and T->x, or allocate them if their current
 * size is zero.  T->x is not modified if T->xtype is CHOLMOD_PATTERN.
 *
 * workspace: none
 */

int CHOLMOD(reallocate_triplet)
(
    /* ---- input ---- */
    size_t nznew,	/* new # of entries in T */
    /* ---- in/out --- */
    cholmod_triplet *T,	/* triplet matrix to modify */
    /* --------------- */
    cholmod_common *Common
)
{

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (T, FALSE) ;
    RETURN_IF_XTYPE_INVALID (T, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, FALSE) ;
    Common->status = CHOLMOD_OK ;
    PRINT1 (("realloc triplet %d to %d, xtype: %d\n",
		T->nzmax, nznew, T->xtype)) ;

    /* ---------------------------------------------------------------------- */
    /* resize the matrix */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(realloc_multiple) (MAX (1,nznew), 2, T->xtype, &(T->i), &(T->j),
	    &(T->x), &(T->z), &(T->nzmax), Common) ;

    return (Common->status == CHOLMOD_OK) ;
}


/* ========================================================================== */
/* === cholmod_triplet_to_sparse ============================================ */
/* ========================================================================== */

/* Convert a set of triplets into a cholmod_sparse matrix.  In MATLAB notation,
 * for unsymmetric matrices:
 *
 *	A = sparse (Ti, Tj, Tx, nrow, ncol, nzmax) ;
 *
 * For the symmetric upper case:
 *
 *	A = sparse (min(Ti,Tj), max(Ti,Tj), Tx, nrow, ncol, nzmax) ;
 *
 * For the symmetric lower case:
 *
 *	A = sparse (max(Ti,Tj), min(Ti,Tj), Tx, nrow, ncol, nzmax) ;
 *
 * If Tx is NULL, then A->x is not allocated, and only the pattern of A is
 * computed.  A is returned in packed form, and can be of any stype
 * (upper/lower/unsymmetric).  It has enough space to hold the values in T,
 * or nzmax, whichever is larger.
 *
 * workspace: Iwork (max (nrow,ncol))
 *	allocates a temporary copy of its output matrix.
 *
 * The resulting sparse matrix has the same xtype as the input triplet matrix.
 */

cholmod_sparse *CHOLMOD(triplet_to_sparse)
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* matrix to copy */
    size_t nzmax,	/* allocate at least this much space in output matrix */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *R, *A = NULL ;
    Int *Wj, *Rp, *Ri, *Rnz, *Ti, *Tj ;
    Int i, j, p, k, stype, nrow, ncol, nz, ok ;
    size_t anz = 0 ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (T, NULL) ;
    Ti = T->i ;
    Tj = T->j ;
    RETURN_IF_NULL (Ti, NULL) ;
    RETURN_IF_NULL (Tj, NULL) ;
    RETURN_IF_XTYPE_INVALID (T, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    stype = SIGN (T->stype) ;
    if (stype && T->nrow != T->ncol)
    {
	/* inputs invalid */
	ERROR (CHOLMOD_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    DEBUG (CHOLMOD(dump_triplet) (T, "T", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    nrow = T->nrow ;
    ncol = T->ncol ;
    nz = T->nnz ;

    /* ---------------------------------------------------------------------- */
    /* allocate workspace */
    /* ---------------------------------------------------------------------- */

    CHOLMOD(allocate_work) (0, MAX (nrow, ncol), 0, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* allocate temporary matrix R */
    /* ---------------------------------------------------------------------- */

    R = CHOLMOD(allocate_sparse) (ncol, nrow, nz, FALSE, FALSE, -stype,
	    T->xtype, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    Rp = R->p ;
    Ri = R->i ;
    Rnz = R->nz ;

    /* ---------------------------------------------------------------------- */
    /* count the entries in each row of A (also counting duplicates) */
    /* ---------------------------------------------------------------------- */

    for (i = 0 ; i < nrow ; i++)
    {
	Rnz [i] = 0 ;	
    }

    if (stype > 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i >= nrow || j < 0 || j >= ncol)
	    {
		ERROR (CHOLMOD_INVALID, "index out of range") ;
		break ;
	    }
	    /* A will be symmetric with just the upper triangular part stored.
	     * Create a matrix R that is lower triangular.  Entries in the
	     * upper part of R are transposed to the lower part. */
	    Rnz [MIN (i,j)]++ ;
	}
    }
    else if (stype < 0)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i >= nrow || j < 0 || j >= ncol)
	    {
		ERROR (CHOLMOD_INVALID, "index out of range") ;
		break ;
	    }
	    /* A will be symmetric with just the lower triangular part stored.
	     * Create a matrix R that is upper triangular.  Entries in the
	     * lower part of R are transposed to the upper part. */
	    Rnz [MAX (i,j)]++ ;
	}
    }
    else
    {
	for (k = 0 ; k < nz ; k++)
	{
	    i = Ti [k] ;
	    j = Tj [k] ;
	    if (i < 0 || i >= nrow || j < 0 || j >= ncol)
	    {
		ERROR (CHOLMOD_INVALID, "index out of range") ;
		break ;
	    }
	    /* constructing an unsymmetric matrix */
	    Rnz [i]++ ;
	}
    }

    if (Common->status < CHOLMOD_OK)
    {
	/* triplet matrix is invalid */
	CHOLMOD(free_sparse) (&R, Common) ;
	return (NULL) ;
    }

    /* ---------------------------------------------------------------------- */
    /* construct the row pointers */
    /* ---------------------------------------------------------------------- */

    p = 0 ;
    for (i = 0 ; i < nrow ; i++)
    {
	Rp [i] = p ;
	p += Rnz [i] ;
    }
    Rp [nrow] = p ;

    /* use Wj (i/l/l) as temporary row pointers */
    Wj = Common->Iwork ;	/* size MAX (nrow,ncol) FUTURE WORK: (i/l/l) */
    for (i = 0 ; i < nrow ; i++)
    {
	Wj [i] = Rp [i] ;
    }

    /* ---------------------------------------------------------------------- */
    /* construct triplet matrix, using template routine */
    /* ---------------------------------------------------------------------- */

    switch (T->xtype)
    {
	case CHOLMOD_PATTERN:
	    anz = p_cholmod_triplet_to_sparse (T, R, Common) ;
	    break ;

	case CHOLMOD_REAL:
	    anz = r_cholmod_triplet_to_sparse (T, R, Common) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    anz = c_cholmod_triplet_to_sparse (T, R, Common) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    anz = z_cholmod_triplet_to_sparse (T, R, Common) ;
	    break ;
    }

    /* ---------------------------------------------------------------------- */
    /* A = R' (array transpose, not complex conjugate transpose) */
    /* ---------------------------------------------------------------------- */

    /* workspace: Iwork (R->nrow), which is A->ncol */

    ASSERT (CHOLMOD(dump_sparse) (R, "R", Common) >= 0) ;

    A = CHOLMOD(allocate_sparse) (nrow, ncol, MAX (anz, nzmax), TRUE, TRUE,
	stype, T->xtype, Common) ;

    if (stype)
    {
	ok = CHOLMOD(transpose_sym) (R, 1, NULL, A, Common) ;
    }
    else
    {
	ok = CHOLMOD(transpose_unsym) (R, 1, NULL, NULL, 0, A, Common) ; 
    }

    CHOLMOD(free_sparse) (&R, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free_sparse) (&A, Common) ;
    }

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_sparse) (A, "A = triplet(T) result", Common) >= 0) ;
    return (A) ;
}


/* ========================================================================== */
/* === cholmod_sparse_to_triplet ============================================ */
/* ========================================================================== */

/* Converts a sparse column-oriented matrix to triplet form.
 * The resulting triplet matrix has the same xtype as the sparse matrix.
 *
 * workspace: none
 */

cholmod_triplet *CHOLMOD(sparse_to_triplet)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to copy */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Ax, *Az, *Tx, *Tz ;
    Int *Ap, *Ai, *Ti, *Tj, *Anz ;
    cholmod_triplet *T ;
    Int i, xtype, p, pend, k, j, nrow, ncol, nz, stype, packed, up, lo,
	both ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    stype = SIGN (A->stype) ;
    nrow = A->nrow ;
    ncol = A->ncol ;
    if (stype && nrow != ncol)
    {
	/* inputs invalid */
	ERROR (CHOLMOD_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Ax = A->x ;
    Az = A->z ;
    xtype = A->xtype ;
    Common->status = CHOLMOD_OK ;

    ASSERT (CHOLMOD(dump_sparse) (A, "A", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* allocate triplet matrix */
    /* ---------------------------------------------------------------------- */

    nz = CHOLMOD(nnz) (A, Common) ;
    T = CHOLMOD(allocate_triplet) (nrow, ncol, nz, A->stype, A->xtype, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* convert to a sparse matrix */
    /* ---------------------------------------------------------------------- */

    Ap = A->p ;
    Ai = A->i ;
    Anz = A->nz ;
    packed = A->packed ;

    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    T->stype = A->stype ;

    both = (A->stype == 0) ;
    up = (A->stype > 0) ;
    lo = (A->stype < 0) ;

    k = 0 ;

    for (j = 0 ; j < ncol ; j++)
    {
	p = Ap [j] ;
	pend = (packed) ? (Ap [j+1]) : (p + Anz [j]) ;
	for ( ; p < pend ; p++)
	{
	    i = Ai [p] ;
	    if (both || (up && i <= j) || (lo && i >= j))
	    {
		Ti [k] = Ai [p] ;
		Tj [k] = j ;

		if (xtype == CHOLMOD_REAL)
		{
		    Tx [k] = Ax [p] ;
		}
		else if (xtype == CHOLMOD_COMPLEX)
		{
		    Tx [2*k  ] = Ax [2*p  ] ;
		    Tx [2*k+1] = Ax [2*p+1] ;
		}
		else if (xtype == CHOLMOD_ZOMPLEX)
		{
		    Tx [k] = Ax [p] ;
		    Tz [k] = Az [p] ;
		}

		k++ ;
		ASSERT (k <= nz) ;
	    }
	}
    }

    T->nnz = k ;

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_triplet) (T, "T", Common)) ;
    return (T) ;
}


/* ========================================================================== */
/* === cholmod_copy_triplet ================================================= */
/* ========================================================================== */

/* Create an exact copy of a triplet matrix, except that entries in unused
 * space are not copied (they might not be initialized, and copying them would
 * cause program checkers such as purify and valgrind to complain).
 * The output triplet matrix has the same xtype as the input triplet matrix.
 */

cholmod_triplet *CHOLMOD(copy_triplet)
(
    /* ---- input ---- */
    cholmod_triplet *T,	/* matrix to copy */
    /* --------------- */
    cholmod_common *Common
)
{
    double *Tx, *Tz, *Cx, *Cz ;
    Int *Ci, *Cj, *Ti, *Tj ;
    cholmod_triplet *C ;
    Int xtype, k, nz ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (T, NULL) ;
    RETURN_IF_XTYPE_INVALID (T, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    nz = T->nnz ;
    Ti = T->i ;
    Tj = T->j ;
    Tx = T->x ;
    Tz = T->z ;
    xtype = T->xtype ;
    RETURN_IF_NULL (Ti, NULL) ;
    RETURN_IF_NULL (Tj, NULL) ;
    Common->status = CHOLMOD_OK ;
    DEBUG (CHOLMOD(dump_triplet) (T, "T input", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* allocate copy */
    /* ---------------------------------------------------------------------- */

    C = CHOLMOD(allocate_triplet) (T->nrow, T->ncol, T->nzmax, T->stype,
	    xtype, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    /* ---------------------------------------------------------------------- */
    /* copy the triplet matrix */
    /* ---------------------------------------------------------------------- */

    Ci = C->i ;
    Cj = C->j ;
    Cx = C->x ;
    Cz = C->z ;
    C->nnz = nz ;

    for (k = 0 ; k < nz ; k++)
    {
	Ci [k] = Ti [k] ;
    }
    for (k = 0 ; k < nz ; k++)
    {
	Cj [k] = Tj [k] ;
    }

    if (xtype == CHOLMOD_REAL)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    Cx [k] = Tx [k] ;
	}
    }
    else if (xtype == CHOLMOD_COMPLEX)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    Cx [2*k  ] = Tx [2*k  ] ;
	    Cx [2*k+1] = Tx [2*k+1] ;
	}
    }
    else if (xtype == CHOLMOD_ZOMPLEX)
    {
	for (k = 0 ; k < nz ; k++)
	{
	    Cx [k] = Tx [k] ;
	    Cz [k] = Tz [k] ;
	}
    }

    /* ---------------------------------------------------------------------- */
    /* return the result */
    /* ---------------------------------------------------------------------- */

    ASSERT (CHOLMOD(dump_triplet) (C, "C triplet copy", Common)) ;
    return (C) ;
}
