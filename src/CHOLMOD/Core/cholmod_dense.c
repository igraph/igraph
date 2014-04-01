/* ========================================================================== */
/* === Core/cholmod_dense =================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Core Module.  Copyright (C) 2005-2013,
 * Univ. of Florida.  Author: Timothy A. Davis
 * The CHOLMOD/Core Module is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* Core utility routines for the cholmod_dense object:
 *
 * The solve routines and some of the MatrixOps and Modify routines use dense
 * matrices as inputs.  These are held in column-major order.  With a leading
 * dimension of d, the entry in row i and column j is held in x [i+j*d].
 *
 * Primary routines:
 * -----------------
 * cholmod_allocate_dense	allocate a dense matrix
 * cholmod_free_dense		free a dense matrix
 *
 * Secondary routines:
 * -------------------
 * cholmod_zeros		allocate a dense matrix of all zeros
 * cholmod_ones			allocate a dense matrix of all ones
 * cholmod_eye			allocate a dense identity matrix 
 * cholmod_sparse_to_dense	create a dense matrix copy of a sparse matrix
 * cholmod_dense_to_sparse	create a sparse matrix copy of a dense matrix
 * cholmod_copy_dense		create a copy of a dense matrix
 * cholmod_copy_dense2		copy a dense matrix (pre-allocated)
 *
 * All routines in this file can handle the real, complex, and zomplex cases.
 * Pattern-only dense matrices are not supported.  cholmod_sparse_to_dense can
 * take a pattern-only input sparse matrix, however, and cholmod_dense_to_sparse
 * can generate a pattern-only output sparse matrix.
 */

#include "cholmod_internal.h"
#include "cholmod_core.h"

/* ========================================================================== */
/* === TEMPLATE ============================================================= */
/* ========================================================================== */

#define PATTERN
#include "t_cholmod_dense.c"
#define REAL
#include "t_cholmod_dense.c"
#define COMPLEX
#include "t_cholmod_dense.c"
#define ZOMPLEX
#include "t_cholmod_dense.c"


/* ========================================================================== */
/* === cholmod_allocate_dense =============================================== */
/* ========================================================================== */

/* Allocate a dense matrix with leading dimension d.  The space is not
 * initialized.
 */

cholmod_dense *CHOLMOD(allocate_dense)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of matrix */
    size_t ncol,	/* # of columns of matrix */
    size_t d,		/* leading dimension */
    int xtype,		/* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *X ;
    size_t nzmax, nzmax0 ;
    int ok = TRUE ;

    /* ---------------------------------------------------------------------- */
    /* get inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    if (d < nrow)
    {
	ERROR (CHOLMOD_INVALID, "leading dimension invalid") ;
	return (NULL) ;
    }
    if (xtype < CHOLMOD_REAL || xtype > CHOLMOD_ZOMPLEX)
    {
	ERROR (CHOLMOD_INVALID, "xtype invalid") ;
	return (NULL) ;
    }

    /* ensure the dimensions do not cause integer overflow */
    (void) CHOLMOD(add_size_t) (ncol, 2, &ok) ;

    /* nzmax = MAX (1, d*ncol) ; */
    nzmax = CHOLMOD(mult_size_t) (d, ncol, &ok) ;
    nzmax = MAX (1, nzmax) ;

    if (!ok || nrow > Int_max || ncol > Int_max || nzmax > Int_max)
    {
	ERROR (CHOLMOD_TOO_LARGE, "problem too large") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate header */
    /* ---------------------------------------------------------------------- */

    X = CHOLMOD(malloc) (sizeof (cholmod_dense), 1, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory */
    }

    PRINT1 (("cholmod_allocate_dense %d-by-%d nzmax %d xtype %d\n",
		nrow, ncol, nzmax, xtype)) ;

    X->nrow = nrow ;
    X->ncol = ncol ;
    X->nzmax = nzmax ;
    X->xtype = xtype ;
    X->dtype = DTYPE ;
    X->x = NULL ;
    X->z = NULL ;
    X->d = d ;

    /* ---------------------------------------------------------------------- */
    /* allocate the matrix itself */
    /* ---------------------------------------------------------------------- */

    nzmax0 = 0 ;
    CHOLMOD(realloc_multiple) (nzmax, 0, xtype, NULL, NULL, &(X->x), &(X->z),
	    &nzmax0, Common) ;

    if (Common->status < CHOLMOD_OK)
    {
	CHOLMOD(free_dense) (&X, Common) ;
	return (NULL) ;	    /* out of memory */
    }

    return (X) ;
}


/* ========================================================================== */
/* === cholmod_zeros ======================================================== */
/* ========================================================================== */

/* Allocate a dense matrix and set it to zero */

cholmod_dense *CHOLMOD(zeros)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of matrix */
    size_t ncol,	/* # of columns of matrix */
    int xtype,		/* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *X ;
    double *Xx, *Xz ;
    Int i, nz ;

    /* ---------------------------------------------------------------------- */
    /* allocate a dense matrix and set it to zero */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    X = CHOLMOD(allocate_dense) (nrow, ncol, nrow, xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* NULL Common, out of memory, or inputs invalid */
    }

    Xx = X->x ;
    Xz = X->z ;
    nz = MAX (1, X->nzmax) ;

    switch (xtype)
    {
	case CHOLMOD_REAL:
	    for (i = 0 ; i < nz ; i++)
	    {
		Xx [i] = 0 ;
	    }
	    break ;

	case CHOLMOD_COMPLEX:
	    for (i = 0 ; i < 2*nz ; i++)
	    {
		Xx [i] = 0 ;
	    }
	    break ;

	case CHOLMOD_ZOMPLEX:
	    for (i = 0 ; i < nz ; i++)
	    {
		Xx [i] = 0 ;
	    }
	    for (i = 0 ; i < nz ; i++)
	    {
		Xz [i] = 0 ;
	    }
	    break ;
    }

    return (X) ;
}


/* ========================================================================== */
/* === cholmod_ones ========================================================= */
/* ========================================================================== */

/* Allocate a dense matrix and set it to zero */

cholmod_dense *CHOLMOD(ones)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of matrix */
    size_t ncol,	/* # of columns of matrix */
    int xtype,		/* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *X ;
    double *Xx, *Xz ;
    Int i, nz ;

    /* ---------------------------------------------------------------------- */
    /* allocate a dense matrix and set it to all ones */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    X = CHOLMOD(allocate_dense) (nrow, ncol, nrow, xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* NULL Common, out of memory, or inputs invalid */
    }

    Xx = X->x ;
    Xz = X->z ;
    nz = MAX (1, X->nzmax) ;

    switch (xtype)
    {
	case CHOLMOD_REAL:
	    for (i = 0 ; i < nz ; i++)
	    {
		Xx [i] = 1 ;
	    }
	    break ;

	case CHOLMOD_COMPLEX:
	    for (i = 0 ; i < nz ; i++)
	    {
		Xx [2*i  ] = 1 ;
		Xx [2*i+1] = 0 ;
	    }
	    break ;

	case CHOLMOD_ZOMPLEX:
	    for (i = 0 ; i < nz ; i++)
	    {
		Xx [i] = 1 ;
	    }
	    for (i = 0 ; i < nz ; i++)
	    {
		Xz [i] = 0 ;
	    }
	    break ;
    }

    return (X) ;
}


/* ========================================================================== */
/* === cholmod_eye ========================================================== */
/* ========================================================================== */

/* Allocate a dense matrix and set it to the identity matrix */

cholmod_dense *CHOLMOD(eye)
(
    /* ---- input ---- */
    size_t nrow,	/* # of rows of matrix */
    size_t ncol,	/* # of columns of matrix */
    int xtype,		/* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *X ;
    double *Xx, *Xz ;
    Int i, n, nz ;

    /* ---------------------------------------------------------------------- */
    /* allocate a dense matrix and set it to the identity matrix */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    X = CHOLMOD(zeros) (nrow, ncol, xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* NULL Common, out of memory, or inputs invalid */
    }

    nz = MAX (1, nrow*ncol) ;
    Xx = X->x ;
    Xz = X->z ;

    n = MIN (nrow, ncol) ;

    switch (xtype)
    {
	case CHOLMOD_REAL:
	case CHOLMOD_ZOMPLEX:
	    for (i = 0 ; i < n ; i++)
	    {
		Xx [i + i*nrow] = 1 ;
	    }
	    break ;

	case CHOLMOD_COMPLEX:
	    for (i = 0 ; i < n ; i++)
	    {
		Xx [2 * (i + i*nrow)] = 1 ;
	    }
	    break ;
    }

    return (X) ;
}

/* ========================================================================== */
/* === cholmod_free_dense =================================================== */
/* ========================================================================== */

/* free a dense matrix
 *
 * workspace: none
 */

int CHOLMOD(free_dense)
(
    /* ---- in/out --- */
    cholmod_dense **XHandle,	/* dense matrix to deallocate, NULL on output */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *X ;

    RETURN_IF_NULL_COMMON (FALSE) ;

    if (XHandle == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }
    X = *XHandle ;
    if (X == NULL)
    {
	/* nothing to do */
	return (TRUE) ;
    }

    switch (X->xtype)
    {
	case CHOLMOD_REAL:
	    X->x = CHOLMOD(free) (X->nzmax, sizeof (double), X->x, Common) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    X->x = CHOLMOD(free) (X->nzmax, 2*sizeof (double), X->x, Common) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    X->x = CHOLMOD(free) (X->nzmax, sizeof (double), X->x, Common) ;
	    X->z = CHOLMOD(free) (X->nzmax, sizeof (double), X->z, Common) ;
	    break ;
    }

    *XHandle = CHOLMOD(free) (1, sizeof (cholmod_dense), (*XHandle), Common) ;
    return (TRUE) ;
}

/* ========================================================================== */
/* === cholmod_ensure_dense ================================================= */
/* ========================================================================== */

/* Ensure that the input matrix has a certain size and type.  If not, free
 * the existing matrix and reallocate one of the right size and type.
 * Returns a pointer to the cholmod_dense matrix, possibly reallocated.
 * Also modifies the input matrix handle, XHandle, if necessary.
 */

cholmod_dense *CHOLMOD(ensure_dense)
(
    /* ---- input/output ---- */
    cholmod_dense **XHandle,    /* matrix handle to check */
    /* ---- input ---- */
    size_t nrow,	/* # of rows of matrix */
    size_t ncol,	/* # of columns of matrix */
    size_t d,		/* leading dimension */
    int xtype,		/* CHOLMOD_REAL, _COMPLEX, or _ZOMPLEX */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *X ;

    RETURN_IF_NULL_COMMON (NULL) ;
    if (XHandle == NULL)
    {
        ERROR (CHOLMOD_INVALID, "matrix invalid") ;
        return (NULL) ;
    }

    X = *XHandle ;
    if (X == NULL || X->nrow != nrow || X->ncol != ncol
        || X->d != d || X->xtype != xtype)
    {
        /* Matrix X is not allocated, or has the wrong size.  Free it and
         * reallocate it in the right size and shape.  If an error occurs
         * (out of memory or inputs nrow, etc invalid), then the error is
         * set in cholmod_allocate_dense and X is returned as NULL. */
#if 0
        if (X == NULL)
        {
            printf ("oops, X was null\n") ;
        }
        else
        {
            printf ("oops, nrow %g %g ncol %g %g d %g %g xtype %g %g\n",
                (double) X->nrow, (double) nrow,
                (double) X->ncol, (double) ncol,
                (double) X->d, (double) d,
                (double) X->xtype, (double) xtype
                ) ;
        }
#endif
        CHOLMOD(free_dense) (XHandle, Common) ;
        X = CHOLMOD(allocate_dense) (nrow, ncol, d, xtype, Common) ;
        *XHandle = X ;
    }
    return (X) ;
}


/* ========================================================================== */
/* === cholmod_sparse_to_dense ============================================== */
/* ========================================================================== */

/* Convert a sparse matrix to a dense matrix.
 * The output dense matrix has the same xtype as the input sparse matrix,
 * except that a pattern-only sparse matrix A is converted into a real dense
 * matrix X, with 1's and 0's.  All xtypes are supported.
 */

cholmod_dense *CHOLMOD(sparse_to_dense)
(
    /* ---- input ---- */
    cholmod_sparse *A,	/* matrix to copy */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *X = NULL ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (A, NULL) ;
    RETURN_IF_XTYPE_INVALID (A, CHOLMOD_PATTERN, CHOLMOD_ZOMPLEX, NULL) ;
    if (A->stype && A->nrow != A->ncol)
    {
	ERROR (CHOLMOD_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;
    ASSERT (CHOLMOD(dump_sparse) (A, "A", Common) >= 0) ;

    /* ---------------------------------------------------------------------- */
    /* convert the matrix, using template routine */
    /* ---------------------------------------------------------------------- */

    switch (A->xtype)
    {
	case CHOLMOD_PATTERN:
	    X = p_cholmod_sparse_to_dense (A, Common) ;
	    break ;

	case CHOLMOD_REAL:
	    X = r_cholmod_sparse_to_dense (A, Common) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    X = c_cholmod_sparse_to_dense (A, Common) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    X = z_cholmod_sparse_to_dense (A, Common) ;
	    break ;
    }
    return (X) ;
}


/* ========================================================================== */
/* === cholmod_dense_to_sparse ============================================== */
/* ========================================================================== */

/* Convert a dense matrix to a sparse matrix, similar to the MATLAB statements:
 *
 * C = sparse (X)			values = TRUE
 * C = spones (sparse (X))		values = FALSE
 *
 * except that X must be double (it can be of many different types in MATLAB)
 *
 * The resulting sparse matrix C has the same numeric xtype as the input dense
 * matrix X, unless "values" is FALSE (in which case C is real, where C(i,j)=1
 * if (i,j) is an entry in X.
 */

cholmod_sparse *CHOLMOD(dense_to_sparse)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* matrix to copy */
    int values,		/* TRUE if values to be copied, FALSE otherwise */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_sparse *C = NULL ;

    DEBUG (CHOLMOD(dump_dense) (X, "X", Common)) ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (X, NULL) ;
    RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, NULL) ;
    if (X->d < X->nrow)
    {
	ERROR (CHOLMOD_INVALID, "matrix invalid") ;
	return (NULL) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* convert the matrix, using template routine */
    /* ---------------------------------------------------------------------- */

    switch (X->xtype)
    {
	case CHOLMOD_REAL:
	    C = r_cholmod_dense_to_sparse (X, values, Common) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    C = c_cholmod_dense_to_sparse (X, values, Common) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    C = z_cholmod_dense_to_sparse (X, values, Common) ;
	    break ;
    }
    return (C) ;
}


/* ========================================================================== */
/* === cholmod_copy_dense2 ================================================== */
/* ========================================================================== */

/* Y = X, where X and Y are both already allocated.  The leading dimensions of
 * X and Y may differ, but both must be >= the # of rows in X and Y.
 * Entries in rows nrow to d-1 are not copied from X, since the space might not
 * be initialized.  Y->nzmax is unchanged.  X->nzmax is typically
 * (X->d)*(X->ncol), but a user might modify that condition outside of any
 * CHOLMOD routine.
 *
 * The two dense matrices X and Y must have the same numeric xtype.
 */

int CHOLMOD(copy_dense2)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* matrix to copy */
    /* ---- output --- */
    cholmod_dense *Y,	/* copy of matrix X */
    /* --------------- */
    cholmod_common *Common
)
{
    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (FALSE) ;
    RETURN_IF_NULL (X, FALSE) ;
    RETURN_IF_NULL (Y, FALSE) ;
    RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    RETURN_IF_XTYPE_INVALID (Y, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, FALSE) ;
    if (X->nrow != Y->nrow || X->ncol != Y->ncol || X->xtype != Y->xtype)
    {
	ERROR (CHOLMOD_INVALID, "X and Y must have same dimensions and xtype") ;
	return (FALSE) ;
    }
    if (X->d < X->nrow || Y->d < Y->nrow
	    || (X->d * X->ncol) > X->nzmax || (Y->d * Y->ncol) > Y->nzmax)
    {
	ERROR (CHOLMOD_INVALID, "X and/or Y invalid") ;
	return (FALSE) ;
    }
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* copy the matrix, using template routine */
    /* ---------------------------------------------------------------------- */

    switch (X->xtype)
    {
	case CHOLMOD_REAL:
	    r_cholmod_copy_dense2 (X, Y) ;
	    break ;

	case CHOLMOD_COMPLEX:
	    c_cholmod_copy_dense2 (X, Y) ;
	    break ;

	case CHOLMOD_ZOMPLEX:
	    z_cholmod_copy_dense2 (X, Y) ;
	    break ;
    }
    return (TRUE) ;
}


/* ========================================================================== */
/* === cholmod_copy_dense =================================================== */
/* ========================================================================== */

/* Y = X, copy a dense matrix */

cholmod_dense *CHOLMOD(copy_dense)
(
    /* ---- input ---- */
    cholmod_dense *X,	/* matrix to copy */
    /* --------------- */
    cholmod_common *Common
)
{
    cholmod_dense *Y ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    RETURN_IF_NULL_COMMON (NULL) ;
    RETURN_IF_NULL (X, NULL) ;
    RETURN_IF_XTYPE_INVALID (X, CHOLMOD_REAL, CHOLMOD_ZOMPLEX, NULL) ;
    Common->status = CHOLMOD_OK ;

    /* ---------------------------------------------------------------------- */
    /* allocate result */
    /* ---------------------------------------------------------------------- */

    Y = CHOLMOD(allocate_dense) (X->nrow, X->ncol, X->d, X->xtype, Common) ;
    if (Common->status < CHOLMOD_OK)
    {
	return (NULL) ;	    /* out of memory or X invalid */
    }

    /* ---------------------------------------------------------------------- */
    /* Y = X */
    /* ---------------------------------------------------------------------- */

    /* This cannot fail (X and Y are allocated, and have the same nrow, ncol
     * d, and xtype) */
    CHOLMOD(copy_dense2) (X, Y, Common) ;

    /* ---------------------------------------------------------------------- */
    /* return result */
    /* ---------------------------------------------------------------------- */

    return (Y) ;
}
