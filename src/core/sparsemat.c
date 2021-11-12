/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <cs.h>

#include "igraph_sparsemat.h"
#include "igraph_error.h"
#include "igraph_interface.h"
#include "igraph_constructors.h"
#include "igraph_memory.h"
#include "igraph_vector_ptr.h"
#include "igraph_attributes.h"

#include <string.h>

/**
 * \section about_sparsemat About sparse matrices
 *
 * <para>
 * The <code>igraph_sparsemat_t</code> data type stores sparse matrices,
 * i.e. matrices in which the majority of the elements are zero.
 * </para>
 *
 * <para>The data type is essentially a wrapper to some of the
 * functions in the CXSparse library, by Tim Davis, see
 * http://faculty.cse.tamu.edu/davis/suitesparse.html
 * </para>
 *
 * <para>
 * Matrices can be stored in two formats: triplet and
 * column-compressed. The triplet format is intended for sparse matrix
 * initialization, as it is easy to add new (non-zero) elements to
 * it. Most of the computations are done on sparse matrices in
 * column-compressed format, after the user has converted the triplet
 * matrix to column-compressed, via \ref igraph_sparsemat_compress().
 * </para>
 *
 * <para>
 * Both formats are dynamic, in the sense that new elements can be
 * added to them, possibly resulting the allocation of more memory.
 * </para>
 *
 * <para>
 * Row and column indices follow the C convention and are zero-based.
 * </para>
 *
 * <para>
 * \example examples/simple/igraph_sparsemat.c
 * \example examples/simple/igraph_sparsemat3.c
 * \example examples/simple/igraph_sparsemat4.c
 * \example examples/simple/igraph_sparsemat6.c
 * \example examples/simple/igraph_sparsemat7.c
 * \example examples/simple/igraph_sparsemat8.c
 * </para>
 */

/**
 * \function igraph_sparsemat_init
 * \brief Initializes a sparse matrix, in triplet format.
 *
 * This is the most common way to create a sparse matrix, together
 * with the \ref igraph_sparsemat_entry() function, which can be used to
 * add the non-zero elements one by one. Once done, the user can call
 * \ref igraph_sparsemat_compress() to convert the matrix to
 * column-compressed, to allow computations with it.
 *
 * </para><para>The user must call \ref igraph_sparsemat_destroy() on
 * the matrix to deallocate the memory, once the matrix is no more
 * needed.
 * \param A Pointer to a not yet initialized sparse matrix.
 * \param rows The number of rows in the matrix.
 * \param cols The number of columns.
 * \param nzmax The maximum number of non-zero elements in the
 *    matrix. It is not compulsory to get this right, but it is
 *    useful for the allocation of the proper amount of memory.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_init(igraph_sparsemat_t *A, int rows, int cols, int nzmax) {

    if (rows < 0) {
        IGRAPH_ERROR("Negative number of rows", IGRAPH_EINVAL);
    }
    if (cols < 0) {
        IGRAPH_ERROR("Negative number of columns", IGRAPH_EINVAL);
    }

    A->cs = cs_spalloc( rows, cols, nzmax, /*values=*/ 1,
                        /*triplet=*/ 1);
    if (!A->cs) {
        IGRAPH_ERROR("Cannot allocate memory for sparse matrix", IGRAPH_ENOMEM);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_copy
 * \brief Copies a sparse matrix.
 *
 * Create a sparse matrix object, by copying another one. The source
 * matrix can be either in triplet or column-compressed format.
 *
 * </para><para>
 * Exactly the same amount of memory will be allocated to the
 * copy matrix, as it is currently for the original one.
 * \param to Pointer to an uninitialized sparse matrix, the copy will
 *    be created here.
 * \param from The sparse matrix to copy.
 * \return Error code.
 *
 * Time complexity: O(n+nzmax), the number of columns plus the maximum
 * number of non-zero elements.
 */

int igraph_sparsemat_copy(igraph_sparsemat_t *to,
                          const igraph_sparsemat_t *from) {

    CS_INT ne = from->cs->nz == -1 ? from->cs->n + 1 : from->cs->nzmax;

    to->cs = cs_spalloc(from->cs->m, from->cs->n, from->cs->nzmax,
                        /*values=*/ 1,
                        /*triplet=*/ igraph_sparsemat_is_triplet(from));

    to->cs->nzmax = from->cs->nzmax;
    to->cs->m     = from->cs->m;
    to->cs->n     = from->cs->n;
    to->cs->nz    = from->cs->nz;

    memcpy(to->cs->p, from->cs->p, sizeof(int) * (size_t) ne);
    memcpy(to->cs->i, from->cs->i, sizeof(int) * (size_t) (from->cs->nzmax));
    memcpy(to->cs->x, from->cs->x, sizeof(double) * (size_t) (from->cs->nzmax));

    return 0;
}

/**
 * \function igraph_sparsemat_destroy
 * \brief Deallocates memory used by a sparse matrix.
 *
 * One destroyed, the sparse matrix must be initialized again, before
 * calling any other operation on it.
 * \param A The sparse matrix to destroy.
 *
 * Time complexity: O(1).
 */

void igraph_sparsemat_destroy(igraph_sparsemat_t *A) {
    cs_spfree(A->cs);
}

/**
 * \function igraph_sparsemat_realloc
 * \brief Allocates more (or less) memory for a sparse matrix.
 *
 * Sparse matrices automatically allocate more memory, as needed. To
 * control memory allocation, the user can call this function, to
 * allocate memory for a given number of non-zero elements.
 *
 * \param A The sparse matrix, it can be in triplet or
 *    column-compressed format.
 * \param nzmax The new maximum number of non-zero elements.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_realloc(igraph_sparsemat_t *A, int nzmax) {
    if (!cs_sprealloc(A->cs, nzmax)) {
        IGRAPH_ERROR("Could not allocate more memory for sparse matrix.", IGRAPH_ENOMEM);
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_sparsemat_nrow
 * \brief Number of rows.
 *
 * \param A The input matrix, in triplet or column-compressed format.
 * \return The number of rows in the \p A matrix.
 *
 * Time complexity: O(1).
 */

long int igraph_sparsemat_nrow(const igraph_sparsemat_t *A) {
    return A->cs->m;
}

/**
 * \function igraph_sparsemat_ncol
 * \brief Number of columns.
 *
 * \param A The input matrix, in triplet or column-compressed format.
 * \return The number of columns in the \p A matrix.
 *
 * Time complexity: O(1).
 */

long int igraph_sparsemat_ncol(const igraph_sparsemat_t *A) {
    return A->cs->n;
}

/**
 * \function igraph_sparsemat_type
 * \brief Type of a sparse matrix (triplet or column-compressed).
 *
 * Gives whether a sparse matrix is stored in the triplet format or in
 * column-compressed format.
 * \param A The input matrix.
 * \return Either \c IGRAPH_SPARSEMAT_CC or \c
 * IGRAPH_SPARSEMAT_TRIPLET.
 *
 * Time complexity: O(1).
 */

igraph_sparsemat_type_t igraph_sparsemat_type(const igraph_sparsemat_t *A) {
    return A->cs->nz < 0 ? IGRAPH_SPARSEMAT_CC : IGRAPH_SPARSEMAT_TRIPLET;
}

/**
 * \function igraph_sparsemat_is_triplet
 * \brief Is this sparse matrix in triplet format?
 *
 * Decides whether a sparse matrix is in triplet format.
 * \param A The input matrix.
 * \return One if the input matrix is in triplet format, zero
 * otherwise.
 *
 * Time complexity: O(1).
 */

igraph_bool_t igraph_sparsemat_is_triplet(const igraph_sparsemat_t *A) {
    return A->cs->nz >= 0;
}

/**
 * \function igraph_sparsemat_is_cc
 * \brief Is this sparse matrix in column-compressed format?
 *
 * Decides whether a sparse matrix is in column-compressed format.
 * \param A The input matrix.
 * \return One if the input matrix is in column-compressed format, zero
 * otherwise.
 *
 * Time complexity: O(1).
 */

igraph_bool_t igraph_sparsemat_is_cc(const igraph_sparsemat_t *A) {
    return A->cs->nz < 0;
}

/**
 * \function igraph_sparsemat_permute
 * \brief Permutes the rows and columns of a sparse matrix.
 *
 * \param A The input matrix, it must be in column-compressed format.
 * \param p Integer vector, giving the permutation of the rows.
 * \param q Integer vector, the permutation of the columns.
 * \param res Pointer to an uninitialized sparse matrix, the result is
 *   stored here.
 * \return Error code.
 *
 * Time complexity: O(m+n+nz), the number of rows plus the number of
 * columns plus the number of non-zero elements in the matrix.
 */

int igraph_sparsemat_permute(const igraph_sparsemat_t *A,
                             const igraph_vector_int_t *p,
                             const igraph_vector_int_t *q,
                             igraph_sparsemat_t *res) {

    CS_INT nrow = A->cs->m, ncol = A->cs->n;
    igraph_vector_int_t pinv;
    CS_INT i;

    if (nrow != igraph_vector_int_size(p)) {
        IGRAPH_ERROR("Invalid row permutation length", IGRAPH_FAILURE);
    }
    if (ncol != igraph_vector_int_size(q)) {
        IGRAPH_ERROR("Invalid column permutation length", IGRAPH_FAILURE);
    }

    /* We invert the permutation by hand */
    IGRAPH_CHECK(igraph_vector_int_init(&pinv, nrow));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &pinv);
    for (i = 0; i < nrow; i++) {
        VECTOR(pinv)[ VECTOR(*p)[i] ] = (int) i;
    }

    /* And call the permutation routine */
    res->cs = cs_permute(A->cs, VECTOR(pinv), VECTOR(*q), /*values=*/ 1);
    if (!res->cs) {
        IGRAPH_ERROR("Cannot index sparse matrix", IGRAPH_FAILURE);
    }

    igraph_vector_int_destroy(&pinv);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_sparsemat_index_rows(const igraph_sparsemat_t *A,
                                         const igraph_vector_int_t *p,
                                         igraph_sparsemat_t *res,
                                         igraph_real_t *constres) {

    igraph_sparsemat_t II, II2;
    CS_INT nrow = A->cs->m;
    long int idx_rows = igraph_vector_int_size(p);
    long int k;

    /* Create index matrix */
    IGRAPH_CHECK(igraph_sparsemat_init(&II2, (int) idx_rows, (int) nrow,
                                       (int) idx_rows));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &II2);
    for (k = 0; k < idx_rows; k++) {
        igraph_sparsemat_entry(&II2, (int) k, VECTOR(*p)[k], 1.0);
    }
    IGRAPH_CHECK(igraph_sparsemat_compress(&II2, &II));
    igraph_sparsemat_destroy(&II2);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &II);

    /* Multiply */
    IGRAPH_CHECK(igraph_sparsemat_multiply(&II, A, res));
    igraph_sparsemat_destroy(&II);
    IGRAPH_FINALLY_CLEAN(1);

    if (constres) {
        if (res->cs->p[1] != 0) {
            *constres = res->cs->x[0];
        } else {
            *constres = 0.0;
        }
    }

    return 0;
}

static int igraph_i_sparsemat_index_cols(const igraph_sparsemat_t *A,
                                         const igraph_vector_int_t *q,
                                         igraph_sparsemat_t *res,
                                         igraph_real_t *constres) {

    igraph_sparsemat_t JJ, JJ2;
    CS_INT ncol = A->cs->n;
    long int idx_cols = igraph_vector_int_size(q);
    long int k;

    /* Create index matrix */
    IGRAPH_CHECK(igraph_sparsemat_init(&JJ2, (int) ncol, (int) idx_cols,
                                       (int) idx_cols));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &JJ2);
    for (k = 0; k < idx_cols; k++) {
        igraph_sparsemat_entry(&JJ2, VECTOR(*q)[k], (int) k, 1.0);
    }
    IGRAPH_CHECK(igraph_sparsemat_compress(&JJ2, &JJ));
    igraph_sparsemat_destroy(&JJ2);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &JJ);

    /* Multiply */
    IGRAPH_CHECK(igraph_sparsemat_multiply(A, &JJ, res));
    igraph_sparsemat_destroy(&JJ);
    IGRAPH_FINALLY_CLEAN(1);

    if (constres) {
        if (res->cs->p [1] != 0) {
            *constres = res->cs->x [0];
        } else {
            *constres = 0.0;
        }
    }

    return 0;
}

/**
 * \function igraph_sparsemat_index
 * \brief Extracts a submatrix or a single element.
 *
 * This function indexes into a sparse matrix.
 * It serves two purposes. First, it can extract
 * submatrices from a sparse matrix. Second, as a special case, it can
 * extract a single element from a sparse matrix.
 *
 * \param A The input matrix, it must be in column-compressed format.
 * \param p An integer vector, or a null pointer. The selected row
 *    index or indices. A null pointer selects all rows.
 * \param q An integer vector, or a null pointer. The selected column
 *    index or indices. A null pointer selects all columns.
 * \param res Pointer to an uninitialized sparse matrix, or a null
 *    pointer. If not a null pointer, then the selected submatrix is
 *    stored here.
 * \param constres Pointer to a real variable or a null pointer. If
 *    not a null pointer, then the first non-zero element in the
 *    selected submatrix is stored here, if there is one. Otherwise
 *    zero is stored here. This behavior is handy if one
 *    wants to select a single entry from the matrix.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_index(const igraph_sparsemat_t *A,
                           const igraph_vector_int_t *p,
                           const igraph_vector_int_t *q,
                           igraph_sparsemat_t *res,
                           igraph_real_t *constres) {

    igraph_sparsemat_t II, JJ, II2, JJ2, tmp;
    CS_INT nrow = A->cs->m;
    CS_INT ncol = A->cs->n;
    long int idx_rows = p ? igraph_vector_int_size(p) : -1;
    long int idx_cols = q ? igraph_vector_int_size(q) : -1;
    long int k;

    igraph_sparsemat_t *myres = res, mres;

    if (!p && !q) {
        IGRAPH_ERROR("No index vectors", IGRAPH_EINVAL);
    }

    if (!res && (idx_rows != 1 || idx_cols != 1)) {
        IGRAPH_ERROR("Sparse matrix indexing: must give `res' if not a "
                     "single element is selected", IGRAPH_EINVAL);
    }

    if (!q) {
        return igraph_i_sparsemat_index_rows(A, p, res, constres);
    }
    if (!p) {
        return igraph_i_sparsemat_index_cols(A, q, res, constres);
    }

    if (!res) {
        myres = &mres;
    }

    /* Create first index matrix */
    IGRAPH_CHECK(igraph_sparsemat_init(&II2, (int) idx_rows, (int) nrow,
                                       (int) idx_rows));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &II2);
    for (k = 0; k < idx_rows; k++) {
        igraph_sparsemat_entry(&II2, (int) k, VECTOR(*p)[k], 1.0);
    }
    IGRAPH_CHECK(igraph_sparsemat_compress(&II2, &II));
    igraph_sparsemat_destroy(&II2);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &II);

    /* Create second index matrix */
    IGRAPH_CHECK(igraph_sparsemat_init(&JJ2, (int) ncol, (int) idx_cols,
                                       (int) idx_cols));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &JJ2);
    for (k = 0; k < idx_cols; k++) {
        igraph_sparsemat_entry(&JJ2, VECTOR(*q)[k], (int) k, 1.0);
    }
    IGRAPH_CHECK(igraph_sparsemat_compress(&JJ2, &JJ));
    igraph_sparsemat_destroy(&JJ2);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &JJ);

    /* Multiply */
    IGRAPH_CHECK(igraph_sparsemat_multiply(&II, A, &tmp));
    igraph_sparsemat_destroy(&II);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
    IGRAPH_CHECK(igraph_sparsemat_multiply(&tmp, &JJ, myres));
    igraph_sparsemat_destroy(&tmp);
    igraph_sparsemat_destroy(&JJ);
    IGRAPH_FINALLY_CLEAN(2);

    if (constres) {
        if (myres->cs->p [1] != 0) {
            *constres = myres->cs->x [0];
        } else {
            *constres = 0.0;
        }
    }

    if (!res) {
        igraph_sparsemat_destroy(myres);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_entry
 * \brief Adds an element to a sparse matrix.
 *
 * This function can be used to add the entries to a sparse matrix,
 * after initializing it with \ref igraph_sparsemat_init(). If you add
 * multiple entries in the same position, they will all be saved, and
 * the resulting value is the sum of all entries in that position.
 *
 * \param A The input matrix, it must be in triplet format.
 * \param row The row index of the entry to add.
 * \param col The column index of the entry to add.
 * \param elem The value of the entry.
 * \return Error code.
 *
 * Time complexity: O(1) on average.
 */

int igraph_sparsemat_entry(igraph_sparsemat_t *A, int row, int col,
                           igraph_real_t elem) {
    if (!igraph_sparsemat_is_triplet(A)) {
        IGRAPH_ERROR("Entries can only be added to sparse matrices that are in triplet format.",
                     IGRAPH_EINVAL);
    }

    if (!cs_entry(A->cs, row, col, elem)) {
        IGRAPH_ERROR("Cannot add entry to sparse matrix.",
                     IGRAPH_FAILURE);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_sparsemat_compress
 * \brief Converts a sparse matrix to column-compressed format.
 *
 * Converts a sparse matrix from triplet format to column-compressed format.
 * Almost all sparse matrix operations require that the matrix is in
 * column-compressed format.
 *
 * \param A The input matrix, it must be in triplet format.
 * \param res Pointer to an uninitialized sparse matrix object, the
 *    compressed version of \p A is stored here.
 * \return Error code.
 *
 * Time complexity: O(nz) where \c nz is the number of non-zero elements.
 */

int igraph_sparsemat_compress(const igraph_sparsemat_t *A,
                              igraph_sparsemat_t *res) {

    if (! igraph_sparsemat_is_triplet(A)) {
        IGRAPH_ERROR("Sparse matrix to compress is not in triplet format.", IGRAPH_EINVAL);
    }
    res->cs = cs_compress(A->cs);
    if (!res->cs) {
        IGRAPH_ERROR("Cannot compress sparse matrix", IGRAPH_FAILURE);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_sparsemat_transpose
 * \brief Transposes a sparse matrix.
 *
 * \param A The input matrix, column-compressed or triple format.
 * \param res Pointer to an uninitialized sparse matrix, the result is
 *    stored here.
 * \param values If this is non-zero, the matrix transpose is
 *    calculated the normal way. If it is zero, then only the pattern
 *    of the input matrix is stored in the result, the values are not.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_transpose(const igraph_sparsemat_t *A,
                               igraph_sparsemat_t *res,
                               int values) {

    if (A->cs->nz < 0) {
        /* column-compressed */
        res->cs = cs_transpose(A->cs, values);
        if (!res->cs) {
            IGRAPH_ERROR("Cannot transpose sparse matrix", IGRAPH_FAILURE);
        }
    } else {
        /* triplets */
        CS_INT *tmp;
        IGRAPH_CHECK(igraph_sparsemat_copy(res, A));
        tmp = res->cs->p;
        res->cs->p = res->cs->i;
        res->cs->i = tmp;
    }
    return 0;
}

static int igraph_i_sparsemat_is_symmetric_cc(const igraph_sparsemat_t *A, igraph_bool_t *result) {
    igraph_sparsemat_t t, tt;
    igraph_bool_t res;
    int nz;

    IGRAPH_CHECK(igraph_sparsemat_transpose(A, &t, /*values=*/ 1));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &t);
    IGRAPH_CHECK(igraph_sparsemat_dupl(&t));
    IGRAPH_CHECK(igraph_sparsemat_transpose(&t, &tt, /*values=*/ 1));
    igraph_sparsemat_destroy(&t);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tt);
    IGRAPH_CHECK(igraph_sparsemat_transpose(&tt, &t, /*values=*/ 1));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &t);

    nz = t.cs->p[t.cs->n];
    res = memcmp(t.cs->i, tt.cs->i, sizeof(int) * (size_t) nz) == 0;
    res = res && memcmp(t.cs->p, tt.cs->p, sizeof(int) *
                        (size_t)(t.cs->n + 1)) == 0;
    res = res && memcmp(t.cs->x, tt.cs->x, sizeof(igraph_real_t) * (size_t)nz) == 0;

    igraph_sparsemat_destroy(&t);
    igraph_sparsemat_destroy(&tt);
    IGRAPH_FINALLY_CLEAN(2);

    *result = res;

    return IGRAPH_SUCCESS;
}

static int igraph_i_sparsemat_is_symmetric_triplet(const igraph_sparsemat_t *A, igraph_bool_t *result) {
    igraph_sparsemat_t tmp;

    IGRAPH_CHECK(igraph_sparsemat_compress(A, &tmp));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
    IGRAPH_CHECK(igraph_i_sparsemat_is_symmetric_cc(&tmp, result));
    igraph_sparsemat_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

igraph_bool_t igraph_sparsemat_is_symmetric(const igraph_sparsemat_t *A) {
    igraph_bool_t res = 0;

    if (A->cs->m != A->cs->n) {
        return 0;
    }

    /* TODO(ntamas): return values from igraph_i_sparsemat_is_symmetric_... are
     * ignored here; this should be fixed. Right now these functions don't
     * change 'res' if they fail so we will report matrices as not being
     * symmetric if an error happens */
    if (A->cs->nz < 0) {
        igraph_i_sparsemat_is_symmetric_cc(A, &res);
    } else {
        igraph_i_sparsemat_is_symmetric_triplet(A, &res);
    }

    return res;
}

/**
 * \function igraph_sparsemat_dupl
 * \brief Removes duplicate elements from a sparse matrix.
 *
 * It is possible that a column-compressed sparse matrix stores a
 * single matrix entry in multiple pieces. The entry is then the sum
 * of all its pieces. (Some functions create matrices like this.) This
 * function eliminates the multiple pieces.
 * \param A The input matrix, in column-compressed format.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_dupl(igraph_sparsemat_t *A) {

    if (!cs_dupl(A->cs)) {
        IGRAPH_ERROR("Cannot remove duplicates from sparse matrix",
                     IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_fkeep
 * \brief Filters the elements of a sparse matrix.
 *
 * This function can be used to filter the (non-zero) elements of a
 * sparse matrix. For all entries, it calls the supplied function and
 * depending on the return values either keeps, or deleted the element
 * from the matrix.
 * \param A The input matrix, in column-compressed format.
 * \param fkeep The filter function. It must take four arguments: the
 *    first is an \c int, the row index of the entry, the second is
 *    another \c int, the column index. The third is \c igraph_real_t,
 *    the value of the entry. The fourth element is a \c void pointer,
 *    the \p other argument is passed here. The function must return
 *    an \c int. If this is zero, then the entry is deleted, otherwise
 *    it is kept.
 * \param other A \c void pointer that is passed to the filtering
 * function.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_fkeep(
    igraph_sparsemat_t *A,
    igraph_integer_t (*fkeep)(igraph_integer_t, igraph_integer_t, igraph_real_t, void*),
    void *other
) {

    IGRAPH_ASSERT(A);
    IGRAPH_ASSERT(fkeep);
    if (!igraph_sparsemat_is_cc(A)) {
        IGRAPH_ERROR("The sparse matrix is not in compressed format.", IGRAPH_EINVAL);
    }
    if (cs_fkeep(A->cs, fkeep, other) < 0) {
        IGRAPH_ERROR("External function cs_keep has returned an unknown error while filtering the matrix.", IGRAPH_FAILURE);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_sparsemat_dropzeros
 * \brief Drops the zero elements from a sparse matrix.
 *
 * As a result of matrix operations, some of the entries in a sparse
 * matrix might be zero. This function removes these entries.
 * \param A The input matrix, it must be in column-compressed format.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_dropzeros(igraph_sparsemat_t *A) {

    if (!cs_dropzeros(A->cs)) {
        IGRAPH_ERROR("Cannot drop zeros from sparse matrix", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_droptol
 * \brief Drops the almost zero elements from a sparse matrix.
 *
 * This function is similar to \ref igraph_sparsemat_dropzeros(), but it
 * also drops entries that are closer to zero than the given tolerance
 * threshold.
 * \param A The input matrix, it must be in column-compressed format.
 * \param tol Real number, giving the tolerance threshold.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_droptol(igraph_sparsemat_t *A, igraph_real_t tol) {

    IGRAPH_ASSERT(A);
    if (!igraph_sparsemat_is_cc(A)) {
        IGRAPH_ERROR("The sparse matrix is not in compressed format.", IGRAPH_EINVAL);
    }
    if (cs_droptol(A->cs, tol) < 0) {
        IGRAPH_ERROR("External function cs_droptol has returned an unknown error.", IGRAPH_FAILURE);
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_sparsemat_multiply
 * \brief Matrix multiplication.
 *
 * Multiplies two sparse matrices.
 * \param A The first input matrix (left hand side), in
 *   column-compressed format.
 * \param B The second input matrix (right hand side), in
 *   column-compressed format.
 * \param res Pointer to an uninitialized sparse matrix, the result is
 *   stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_multiply(const igraph_sparsemat_t *A,
                              const igraph_sparsemat_t *B,
                              igraph_sparsemat_t *res) {

    res->cs = cs_multiply(A->cs, B->cs);
    if (!res->cs) {
        IGRAPH_ERROR("Cannot multiply matrices", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_add
 * \brief Sum of two sparse matrices.
 *
 * \param A The first input matrix, in column-compressed format.
 * \param B The second input matrix, in column-compressed format.
 * \param alpha Real scalar, \p A is multiplied by \p alpha before the
 *    addition.
 * \param beta Real scalar, \p B is multiplied by \p beta before the
 *    addition.
 * \param res Pointer to an uninitialized sparse matrix, the result
 *    is stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_add(const igraph_sparsemat_t *A,
                         const igraph_sparsemat_t *B,
                         igraph_real_t alpha,
                         igraph_real_t beta,
                         igraph_sparsemat_t *res) {

    res->cs = cs_add(A->cs, B->cs, alpha, beta);
    if (!res->cs) {
        IGRAPH_ERROR("Cannot add matrices", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_gaxpy
 * \brief Matrix-vector product, added to another vector.
 *
 * \param A The input matrix, in column-compressed format.
 * \param x The input vector, its size must match the number of
 *    columns in \p A.
 * \param res This vector is added to the matrix-vector product
 *    and it is overwritten by the result.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_gaxpy(const igraph_sparsemat_t *A,
                           const igraph_vector_t *x,
                           igraph_vector_t *res) {

    if (A->cs->n != igraph_vector_size(x) ||
        A->cs->m != igraph_vector_size(res)) {
        IGRAPH_ERROR("Invalid matrix/vector size for multiplication",
                     IGRAPH_EINVAL);
    }

    if (! (cs_gaxpy(A->cs, VECTOR(*x), VECTOR(*res)))) {
        IGRAPH_ERROR("Cannot perform sparse matrix vector multiplication",
                     IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_lsolve
 * \brief Solves a lower-triangular linear system.
 *
 * Solve the Lx=b linear equation system, where the L coefficient
 * matrix is square and lower-triangular, with a zero-free diagonal.
 * \param L The input matrix, in column-compressed format.
 * \param b The right hand side of the linear system.
 * \param res An initialized vector, the result is stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_lsolve(const igraph_sparsemat_t *L,
                            const igraph_vector_t *b,
                            igraph_vector_t *res) {

    if (L->cs->m != L->cs->n) {
        IGRAPH_ERROR("Cannot perform lower triangular solve", IGRAPH_NONSQUARE);
    }

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    if (! cs_lsolve(L->cs, VECTOR(*res))) {
        IGRAPH_ERROR("Cannot perform lower triangular solve", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_ltsolve
 * \brief Solves an upper-triangular linear system.
 *
 * Solve the L'x=b linear equation system, where the L
 * matrix is square and lower-triangular, with a zero-free diagonal.
 * \param L The input matrix, in column-compressed format.
 * \param b The right hand side of the linear system.
 * \param res An initialized vector, the result is stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_ltsolve(const igraph_sparsemat_t *L,
                             const igraph_vector_t *b,
                             igraph_vector_t *res) {

    if (L->cs->m != L->cs->n) {
        IGRAPH_ERROR("Cannot perform transposed lower triangular solve",
                     IGRAPH_NONSQUARE);
    }

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    if (!cs_ltsolve(L->cs, VECTOR(*res))) {
        IGRAPH_ERROR("Cannot perform lower triangular solve", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_usolve
 * \brief Solves an upper-triangular linear system.
 *
 * Solves the Ux=b upper triangular system.
 * \param U The input matrix, in column-compressed format.
 * \param b The right hand side of the linear system.
 * \param res An initialized vector, the result is stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_usolve(const igraph_sparsemat_t *U,
                            const igraph_vector_t *b,
                            igraph_vector_t *res) {

    if (U->cs->m != U->cs->n) {
        IGRAPH_ERROR("Cannot perform upper triangular solve", IGRAPH_NONSQUARE);
    }

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    if (! cs_usolve(U->cs, VECTOR(*res))) {
        IGRAPH_ERROR("Cannot perform upper triangular solve", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_utsolve
 * \brief Solves a lower-triangular linear system.
 *
 * This is the same as \ref igraph_sparsemat_usolve(), but U'x=b is
 * solved, where the apostrophe denotes the transpose.
 * \param U The input matrix, in column-compressed format.
 * \param b The right hand side of the linear system.
 * \param res An initialized vector, the result is stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_utsolve(const igraph_sparsemat_t *U,
                             const igraph_vector_t *b,
                             igraph_vector_t *res) {

    if (U->cs->m != U->cs->n) {
        IGRAPH_ERROR("Cannot perform transposed upper triangular solve",
                     IGRAPH_NONSQUARE);
    }

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    if (!cs_utsolve(U->cs, VECTOR(*res))) {
        IGRAPH_ERROR("Cannot perform transposed upper triangular solve",
                     IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_cholsol
 * \brief Solves a symmetric linear system via Cholesky decomposition.
 *
 * Solve Ax=b, where A is a symmetric positive definite matrix.
 * \param A The input matrix, in column-compressed format.
 * \param v The right hand side.
 * \param res An initialized vector, the result is stored here.
 * \param order An integer giving the ordering method to use for the
 *    factorization. Zero is the natural ordering; if it is one, then
 *    the fill-reducing minimum-degree ordering of A+A' is used.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_cholsol(const igraph_sparsemat_t *A,
                             const igraph_vector_t *b,
                             igraph_vector_t *res,
                             int order) {

    if (A->cs->m != A->cs->n) {
        IGRAPH_ERROR("Cannot perform sparse symmetric solve",
                     IGRAPH_NONSQUARE);
    }

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    if (! cs_cholsol(order, A->cs, VECTOR(*res))) {
        IGRAPH_ERROR("Cannot perform sparse symmetric solve", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_lusol
 * \brief Solves a linear system via LU decomposition.
 *
 * Solve Ax=b, via LU factorization of A.
 * \param A The input matrix, in column-compressed format.
 * \param b The right hand side of the equation.
 * \param res An initialized vector, the result is stored here.
 * \param order The ordering method to use, zero means the natural
 *    ordering, one means the fill-reducing minimum-degree ordering of
 *    A+A', two means the ordering of A'*A, after removing the dense
 *    rows from A. Three means the ordering of A'*A.
 * \param tol Real number, the tolerance limit to use for the numeric
 *    LU factorization.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_lusol(const igraph_sparsemat_t *A,
                           const igraph_vector_t *b,
                           igraph_vector_t *res,
                           int order,
                           igraph_real_t tol) {

    if (A->cs->m != A->cs->n) {
        IGRAPH_ERROR("Cannot perform LU solve",
                     IGRAPH_NONSQUARE);
    }

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    if (! cs_lusol(order, A->cs, VECTOR(*res), tol)) {
        IGRAPH_ERROR("Cannot perform LU solve", IGRAPH_FAILURE);
    }

    return 0;
}

static int igraph_i_sparsemat_cc(igraph_t *graph, const igraph_sparsemat_t *A,
                                 igraph_bool_t directed) {

    igraph_vector_t edges;
    CS_INT no_of_nodes = A->cs->m;
    CS_INT no_of_edges = A->cs->p[A->cs->n];
    CS_INT *p = A->cs->p;
    CS_INT *i = A->cs->i;
    long int from = 0;
    long int to = 0;
    long int e = 0;

    if (no_of_nodes != A->cs->n) {
        IGRAPH_ERROR("Cannot create graph object", IGRAPH_NONSQUARE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    while (*p < no_of_edges) {
        while (to < * (p + 1)) {
            if (directed || from >= *i) {
                VECTOR(edges)[e++] = from;
                VECTOR(edges)[e++] = (*i);
            }
            to++;
            i++;
        }
        from++;
        p++;
    }
    igraph_vector_resize(&edges, e);

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_sparsemat_triplet(igraph_t *graph, const igraph_sparsemat_t *A,
                                      igraph_bool_t directed) {

    igraph_vector_t edges;
    CS_INT no_of_nodes = A->cs->m;
    CS_INT no_of_edges = A->cs->nz;
    CS_INT *i = A->cs->p;
    CS_INT *j = A->cs->i;
    long int e;

    if (no_of_nodes != A->cs->n) {
        IGRAPH_ERROR("Cannot create graph object", IGRAPH_NONSQUARE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);

    for (e = 0; e < 2 * no_of_edges; i++, j++) {
        if (directed || *i >= *j) {
            VECTOR(edges)[e++] = (*i);
            VECTOR(edges)[e++] = (*j);
        }
    }
    igraph_vector_resize(&edges, e);

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_sparsemat
 * \brief Creates an igraph graph from a sparse matrix.
 *
 * One edge is created for each non-zero entry in the matrix. If you
 * have a symmetric matrix, and want to create an undirected graph,
 * then delete the entries in the upper diagonal first, or call \ref
 * igraph_simplify() on the result graph to eliminate the multiple
 * edges.
 * \param graph Pointer to an uninitialized igraph_t object, the
 *    graphs is stored here.
 * \param A The input matrix, in triplet or column-compressed format.
 * \param directed Boolean scalar, whether to create a directed
 *    graph.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat(igraph_t *graph, const igraph_sparsemat_t *A,
                     igraph_bool_t directed) {

    if (A->cs->nz < 0) {
        return (igraph_i_sparsemat_cc(graph, A, directed));
    } else {
        return (igraph_i_sparsemat_triplet(graph, A, directed));
    }
}

static int igraph_i_weighted_sparsemat_cc(const igraph_sparsemat_t *A,
                                          igraph_bool_t directed, const char *attr,
                                          igraph_bool_t loops,
                                          igraph_vector_t *edges,
                                          igraph_vector_t *weights) {

    CS_INT no_of_edges = A->cs->p[A->cs->n];
    CS_INT *p = A->cs->p;
    CS_INT *i = A->cs->i;
    CS_ENTRY *x = A->cs->x;
    long int from = 0;
    long int to = 0;
    long int e = 0, w = 0;

    IGRAPH_UNUSED(attr);

    igraph_vector_resize(edges, no_of_edges * 2);
    igraph_vector_resize(weights, no_of_edges);

    while (*p < no_of_edges) {
        while (to < * (p + 1)) {
            if ( (loops || from != *i) && (directed || from >= *i) && *x != 0) {
                VECTOR(*edges)[e++] = (*i);
                VECTOR(*edges)[e++] = from;
                VECTOR(*weights)[w++] = (*x);
            }
            to++;
            i++;
            x++;
        }
        from++;
        p++;
    }

    igraph_vector_resize(edges, e);
    igraph_vector_resize(weights, w);

    return 0;
}

static int igraph_i_weighted_sparsemat_triplet(const igraph_sparsemat_t *A,
                                               igraph_bool_t directed,
                                               const char *attr,
                                               igraph_bool_t loops,
                                               igraph_vector_t *edges,
                                               igraph_vector_t *weights) {

    IGRAPH_UNUSED(A); IGRAPH_UNUSED(directed); IGRAPH_UNUSED(attr);
    IGRAPH_UNUSED(loops); IGRAPH_UNUSED(edges); IGRAPH_UNUSED(weights);

    /* TODO */
    IGRAPH_ERROR("Triplet matrices are not implemented",
                 IGRAPH_UNIMPLEMENTED);
}

int igraph_weighted_sparsemat(igraph_t *graph, const igraph_sparsemat_t *A,
                              igraph_bool_t directed, const char *attr,
                              igraph_bool_t loops) {

    igraph_vector_t edges, weights;
    CS_INT pot_edges = A->cs->nz < 0 ? A->cs->p[A->cs->n] : A->cs->nz;
    const char* default_attr = "weight";
    igraph_vector_ptr_t attr_vec;
    igraph_attribute_record_t attr_rec;
    CS_INT no_of_nodes = A->cs->m;

    if (no_of_nodes != A->cs->n) {
        IGRAPH_ERROR("Cannot create graph object", IGRAPH_NONSQUARE);
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, pot_edges * 2);
    IGRAPH_VECTOR_INIT_FINALLY(&weights, pot_edges);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&attr_vec, 1);

    if (A->cs->nz < 0) {
        IGRAPH_CHECK(igraph_i_weighted_sparsemat_cc(A, directed, attr, loops,
                     &edges, &weights));
    } else {
        IGRAPH_CHECK(igraph_i_weighted_sparsemat_triplet(A, directed, attr,
                     loops, &edges,
                     &weights));
    }

    /* Prepare attribute record */
    attr_rec.name = attr ? attr : default_attr;
    attr_rec.type = IGRAPH_ATTRIBUTE_NUMERIC;
    attr_rec.value = &weights;
    VECTOR(attr_vec)[0] = &attr_rec;

    /* Create graph */
    IGRAPH_CHECK(igraph_empty(graph, (igraph_integer_t) no_of_nodes, directed));
    IGRAPH_FINALLY(igraph_destroy, graph);
    if (igraph_vector_size(&edges) > 0) {
        IGRAPH_CHECK(igraph_add_edges(graph, &edges, &attr_vec));
    }
    IGRAPH_FINALLY_CLEAN(1);

    /* Cleanup */
    igraph_vector_destroy(&edges);
    igraph_vector_destroy(&weights);
    igraph_vector_ptr_destroy(&attr_vec);
    IGRAPH_FINALLY_CLEAN(3);

    return 0;
}

/**
 * \function igraph_get_sparsemat
 * \brief Converts an igraph graph to a sparse matrix.
 *
 * If the graph is undirected, then a symmetric matrix is created.
 * \param graph The input graph.
 * \param res Pointer to an uninitialized sparse matrix. The result
 *    will be stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_get_sparsemat(const igraph_t *graph, igraph_sparsemat_t *res) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    igraph_bool_t directed = igraph_is_directed(graph);
    long int nzmax = directed ? no_of_edges : no_of_edges * 2;
    long int i;

    IGRAPH_CHECK(igraph_sparsemat_init(res, (igraph_integer_t) no_of_nodes,
                                       (igraph_integer_t) no_of_nodes,
                                       (igraph_integer_t) nzmax));

    for (i = 0; i < no_of_edges; i++) {
        long int from = IGRAPH_FROM(graph, i);
        long int to = IGRAPH_TO(graph, i);
        IGRAPH_CHECK(igraph_sparsemat_entry(res, (int) from, (int) to, 1.0));
        if (!directed && from != to) {
            IGRAPH_CHECK(igraph_sparsemat_entry(res, (int) to, (int) from, 1.0));
        }
    }

    return 0;
}

#define CHECK(x) if ((x)<0) { IGRAPH_ERROR("Cannot write to file", IGRAPH_EFILE); }

/**
 * \function igraph_sparsemat_print
 * \brief Prints a sparse matrix to a file.
 *
 * Only the non-zero entries are printed. This function serves more as
 * a debugging utility, as currently there is no function that could
 * read back the printed matrix from the file.
 * \param A The input matrix, triplet or column-compressed format.
 * \param outstream The stream to print it to.
 * \return Error code.
 *
 * Time complexity: O(nz) for triplet matrices, O(n+nz) for
 * column-compressed matrices. nz is the number of non-zero elements,
 * n is the number columns in the matrix.
 */

int igraph_sparsemat_print(const igraph_sparsemat_t *A,
                           FILE *outstream) {

    if (A->cs->nz < 0) {
        /* CC */
        CS_INT j, p;
        for (j = 0; j < A->cs->n; j++) {
            CHECK(fprintf(outstream, "col " CS_ID ": locations " CS_ID " to " CS_ID "\n",
                          j, A->cs->p[j], A->cs->p[j + 1] - 1));
            for (p = A->cs->p[j]; p < A->cs->p[j + 1]; p++) {
                CHECK(fprintf(outstream, CS_ID " : %g\n", A->cs->i[p], A->cs->x[p]));
            }
        }
    } else {
        /* Triplet */
        CS_INT p;
        for (p = 0; p < A->cs->nz; p++) {
            CHECK(fprintf(outstream, CS_ID " " CS_ID " : %g\n",
                          A->cs->i[p], A->cs->p[p], A->cs->x[p]));
        }
    }

    return 0;
}

#undef CHECK

static int igraph_i_sparsemat_eye_triplet(igraph_sparsemat_t *A, int n, int nzmax,
                                          igraph_real_t value) {
    long int i;

    IGRAPH_CHECK(igraph_sparsemat_init(A, n, n, nzmax));

    for (i = 0; i < n; i++) {
        igraph_sparsemat_entry(A, (int) i, (int) i, value);
    }

    return 0;
}

static int igraph_i_sparsemat_eye_cc(igraph_sparsemat_t *A, int n,
                                     igraph_real_t value) {
    CS_INT i;

    A->cs = cs_spalloc(n, n, n, /*values=*/ 1, /*triplet=*/ 0);
    if (!A->cs) {
        IGRAPH_ERROR("Cannot create eye sparse matrix", IGRAPH_FAILURE);
    }

    for (i = 0; i < n; i++) {
        A->cs->p [i] = i;
        A->cs->i [i] = i;
        A->cs->x [i] = value;
    }
    A->cs->p [n] = n;

    return 0;
}

/**
 * \function igraph_sparsemat_eye
 * \brief Creates a sparse identity matrix.
 *
 * \param A An uninitialized sparse matrix, the result is stored
 *   here.
 * \param n The number of rows and number of columns in the matrix.
 * \param nzmax The maximum number of non-zero elements, this
 *   essentially gives the amount of memory that will be allocated for
 *   matrix elements.
 * \param value The value to store in the diagonal.
 * \param compress Whether to create a column-compressed matrix. If
 *   false, then a triplet matrix is created.
 * \return Error code.
 *
 * Time complexity: O(n).
 */

int igraph_sparsemat_eye(igraph_sparsemat_t *A, int n, int nzmax,
                         igraph_real_t value,
                         igraph_bool_t compress) {
    if (compress) {
        return (igraph_i_sparsemat_eye_cc(A, n, value));
    } else {
        return (igraph_i_sparsemat_eye_triplet(A, n, nzmax, value));
    }
}

static int igraph_i_sparsemat_diag_triplet(igraph_sparsemat_t *A, int nzmax,
                                           const igraph_vector_t *values) {

    int i, n = (int) igraph_vector_size(values);

    IGRAPH_CHECK(igraph_sparsemat_init(A, n, n, nzmax));

    for (i = 0; i < n; i++) {
        igraph_sparsemat_entry(A, i, i, VECTOR(*values)[i]);
    }

    return 0;

}

static int igraph_i_sparsemat_diag_cc(igraph_sparsemat_t *A,
                                      const igraph_vector_t *values) {

    CS_INT i, n = igraph_vector_size(values);

    A->cs = cs_spalloc(n, n, n, /*values=*/ 1, /*triplet=*/ 0);
    if (!A->cs) {
        IGRAPH_ERROR("Cannot create eye sparse matrix", IGRAPH_FAILURE);
    }

    for (i = 0; i < n; i++) {
        A->cs->p [i] = i;
        A->cs->i [i] = i;
        A->cs->x [i] = VECTOR(*values)[i];
    }
    A->cs->p [n] = n;

    return 0;

}

/**
 * \function igraph_sparsemat_diag
 * \brief Creates a sparse diagonal matrix.
 *
 * \param A An uninitialized sparse matrix, the result is stored
 *    here.
 * \param nzmax The maximum number of non-zero elements, this
 *   essentially gives the amount of memory that will be allocated for
 *   matrix elements.
 * \param values The values to store in the diagonal, the size of the
 *    matrix defined by the length of this vector.
 * \param compress Whether to create a column-compressed matrix. If
 *   false, then a triplet matrix is created.
 * \return Error code.
 *
 * Time complexity: O(n), the length of the diagonal vector.
 */

int igraph_sparsemat_diag(igraph_sparsemat_t *A, int nzmax,
                          const igraph_vector_t *values,
                          igraph_bool_t compress) {

    if (compress) {
        return (igraph_i_sparsemat_diag_cc(A, values));
    } else {
        return (igraph_i_sparsemat_diag_triplet(A, nzmax, values));
    }
}

static int igraph_i_sparsemat_arpack_multiply(igraph_real_t *to,
                                              const igraph_real_t *from,
                                              int n,
                                              void *extra) {
    igraph_sparsemat_t *A = extra;
    igraph_vector_t vto, vfrom;
    igraph_vector_view(&vto, to, n);
    igraph_vector_view(&vfrom, from, n);
    igraph_vector_null(&vto);
    IGRAPH_CHECK(igraph_sparsemat_gaxpy(A, &vfrom, &vto));
    return 0;
}

typedef struct igraph_i_sparsemat_arpack_rssolve_data_t {
    igraph_sparsemat_symbolic_t *dis;
    igraph_sparsemat_numeric_t *din;
    igraph_real_t tol;
    igraph_sparsemat_solve_t method;
} igraph_i_sparsemat_arpack_rssolve_data_t;

static int igraph_i_sparsemat_arpack_solve(igraph_real_t *to,
                                           const igraph_real_t *from,
                                           int n,
                                           void *extra) {

    igraph_i_sparsemat_arpack_rssolve_data_t *data = extra;
    igraph_vector_t vfrom, vto;

    igraph_vector_view(&vfrom, from, n);
    igraph_vector_view(&vto, to, n);

    if (data->method == IGRAPH_SPARSEMAT_SOLVE_LU) {
        IGRAPH_CHECK(igraph_sparsemat_luresol(data->dis, data->din, &vfrom,
                                              &vto));
    } else if (data->method == IGRAPH_SPARSEMAT_SOLVE_QR) {
        IGRAPH_CHECK(igraph_sparsemat_qrresol(data->dis, data->din, &vfrom,
                                              &vto));

    }

    return 0;
}

/**
 * \function igraph_sparsemat_arpack_rssolve
 * \brief Eigenvalues and eigenvectors of a symmetric sparse matrix via ARPACK.
 *
 * \param The input matrix, must be column-compressed.
 * \param options It is passed to \ref igraph_arpack_rssolve(). See
 *    \ref igraph_arpack_options_t for the details. If \c mode is 1,
 *    then ARPACK uses regular mode, if \c mode is 3, then shift and
 *    invert mode is used and the \c sigma structure member defines
 *    the shift.
 * \param storage Storage for ARPACK. See \ref
 *    igraph_arpack_rssolve() and \ref igraph_arpack_storage_t for
 *    details.
 * \param values An initialized vector or a null pointer, the
 *    eigenvalues are stored here.
 * \param vectors An initialised matrix, or a null pointer, the
 *    eigenvectors are stored here, in the columns.
 * \param solvemethod The method to solve the linear system, if \c
 *    mode is 3, i.e. the shift and invert mode is used.
 *    Possible values:
 *    \clist
 *      \cli IGRAPH_SPARSEMAT_SOLVE_LU
 *           The linear system is solved using LU decomposition.
 *      \cli IGRAPH_SPARSEMAT_SOLVE_QR
 *           The linear system is solved using QR decomposition.
 *    \endclist
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_arpack_rssolve(const igraph_sparsemat_t *A,
                                    igraph_arpack_options_t *options,
                                    igraph_arpack_storage_t *storage,
                                    igraph_vector_t *values,
                                    igraph_matrix_t *vectors,
                                    igraph_sparsemat_solve_t solvemethod) {

    int n = (int) igraph_sparsemat_nrow(A);

    if (n != igraph_sparsemat_ncol(A)) {
        IGRAPH_ERROR("Non-square matrix for ARPACK", IGRAPH_NONSQUARE);
    }

    options->n = n;

    if (options->mode == 1) {
        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_sparsemat_arpack_multiply,
                                           (void*) A, options, storage,
                                           values, vectors));
    } else if (options->mode == 3) {
        igraph_real_t sigma = options->sigma;
        igraph_sparsemat_t OP, eye;
        igraph_sparsemat_symbolic_t symb;
        igraph_sparsemat_numeric_t num;
        igraph_i_sparsemat_arpack_rssolve_data_t data;
        /*-----------------------------------*/
        /* We need to factor the (A-sigma*I) */
        /*-----------------------------------*/

        /* Create (A-sigma*I) */
        IGRAPH_CHECK(igraph_sparsemat_eye(&eye, /*n=*/ n, /*nzmax=*/ n,
                                          /*value=*/ -sigma, /*compress=*/ 1));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, &eye);
        IGRAPH_CHECK(igraph_sparsemat_add(/*A=*/ A, /*B=*/ &eye, /*alpha=*/ 1.0,
                     /*beta=*/ 1.0, /*res=*/ &OP));
        igraph_sparsemat_destroy(&eye);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_FINALLY(igraph_sparsemat_destroy, &OP);

        if (solvemethod == IGRAPH_SPARSEMAT_SOLVE_LU) {
            /* Symbolic analysis */
            IGRAPH_CHECK(igraph_sparsemat_symblu(/*order=*/ 0, &OP, &symb));
            IGRAPH_FINALLY(igraph_sparsemat_symbolic_destroy, &symb);
            /* Numeric LU factorization */
            IGRAPH_CHECK(igraph_sparsemat_lu(&OP, &symb, &num, /*tol=*/ 0));
            IGRAPH_FINALLY(igraph_sparsemat_numeric_destroy, &num);
        } else if (solvemethod == IGRAPH_SPARSEMAT_SOLVE_QR) {
            /* Symbolic analysis */
            IGRAPH_CHECK(igraph_sparsemat_symbqr(/*order=*/ 0, &OP, &symb));
            IGRAPH_FINALLY(igraph_sparsemat_symbolic_destroy, &symb);
            /* Numeric QR factorization */
            IGRAPH_CHECK(igraph_sparsemat_qr(&OP, &symb, &num));
            IGRAPH_FINALLY(igraph_sparsemat_numeric_destroy, &num);
        }

        data.dis = &symb;
        data.din = &num;
        data.tol = options->tol;
        data.method = solvemethod;
        IGRAPH_CHECK(igraph_arpack_rssolve(igraph_i_sparsemat_arpack_solve,
                                           (void*) &data, options, storage,
                                           values, vectors));

        igraph_sparsemat_numeric_destroy(&num);
        igraph_sparsemat_symbolic_destroy(&symb);
        igraph_sparsemat_destroy(&OP);
        IGRAPH_FINALLY_CLEAN(3);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_arpack_rnsolve
 * \brief Eigenvalues and eigenvectors of a nonsymmetric sparse matrix via ARPACK.
 *
 * Eigenvalues and/or eigenvectors of a nonsymmetric sparse matrix.
 * \param A The input matrix, in column-compressed mode.
 * \param options ARPACK options, it is passed to \ref
 *    igraph_arpack_rnsolve(). See also \ref igraph_arpack_options_t
 *    for details.
 * \param storage Storage for ARPACK, this is passed to \ref
 *    igraph_arpack_rnsolve(). See \ref igraph_arpack_storage_t for
 *    details.
 * \param values An initialized matrix, or a null pointer. If not a
 *    null pointer, then the eigenvalues are stored here, the first
 *    column is the real part, the second column is the imaginary
 *    part.
 * \param vectors An initialized matrix, or a null pointer. If not a
 *    null pointer, then the eigenvectors are stored here, please see
 *    \ref igraph_arpack_rnsolve() for the format.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_arpack_rnsolve(const igraph_sparsemat_t *A,
                                    igraph_arpack_options_t *options,
                                    igraph_arpack_storage_t *storage,
                                    igraph_matrix_t *values,
                                    igraph_matrix_t *vectors) {

    int n = (int) igraph_sparsemat_nrow(A);

    if (n != igraph_sparsemat_ncol(A)) {
        IGRAPH_ERROR("Non-square matrix for ARPACK", IGRAPH_NONSQUARE);
    }

    options->n = n;

    return igraph_arpack_rnsolve(igraph_i_sparsemat_arpack_multiply,
                                 (void*) A, options, storage,
                                 values, vectors);
}

/**
 * \function igraph_sparsemat_symbqr
 * \brief Symbolic QR decomposition.
 *
 * QR decomposition of sparse matrices involves two steps, the first
 * is calling this function, and then \ref
 * igraph_sparsemat_qr().
 * \param order The ordering to use: 0 means natural ordering, 1 means
 *   minimum degree ordering of A+A', 2 is minimum degree ordering of
 *   A'A after removing the dense rows from A, and 3 is the minimum
 *   degree ordering of A'A.
 * \param A The input matrix, in column-compressed format.
 * \param dis The result of the symbolic analysis is stored here. Once
 *    not needed anymore, it must be destroyed by calling \ref
 *    igraph_sparsemat_symbolic_destroy().
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_symbqr(long int order, const igraph_sparsemat_t *A,
                            igraph_sparsemat_symbolic_t *dis) {

    dis->symbolic = cs_sqr((int) order, A->cs, /*qr=*/ 1);
    if (!dis->symbolic) {
        IGRAPH_ERROR("Cannot do symbolic QR decomposition", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_symblu
 * \brief Symbolic LU decomposition.
 *
 * LU decomposition of sparse matrices involves two steps, the first
 * is calling this function, and then \ref igraph_sparsemat_lu().
 * \param order The ordering to use: 0 means natural ordering, 1 means
 *   minimum degree ordering of A+A', 2 is minimum degree ordering of
 *   A'A after removing the dense rows from A, and 3 is the minimum
 *   degree ordering of A'A.
 * \param A The input matrix, in column-compressed format.
 * \param dis The result of the symbolic analysis is stored here. Once
 *    not needed anymore, it must be destroyed by calling \ref
 *    igraph_sparsemat_symbolic_destroy().
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_symblu(long int order, const igraph_sparsemat_t *A,
                            igraph_sparsemat_symbolic_t *dis) {

    dis->symbolic = cs_sqr((int) order, A->cs, /*qr=*/ 0);
    if (!dis->symbolic) {
        IGRAPH_ERROR("Cannot do symbolic LU decomposition", IGRAPH_FAILURE);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_lu
 * \brief LU decomposition of a sparse matrix.
 *
 * Performs numeric sparse LU decomposition of a matrix.
 * \param A The input matrix, in column-compressed format.
 * \param dis The symbolic analysis for LU decomposition, coming from
 *    a call to the \ref igraph_sparsemat_symblu() function.
 * \param din The numeric decomposition, the result is stored here. It
 *    can be used to solve linear systems with changing right hand
 *    side vectors, by calling \ref igraph_sparsemat_luresol(). Once
 *    not needed any more, it must be destroyed by calling \ref
 *    igraph_sparsemat_symbolic_destroy() on it.
 * \param tol The tolerance for the numeric LU decomposition.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_lu(const igraph_sparsemat_t *A,
                        const igraph_sparsemat_symbolic_t *dis,
                        igraph_sparsemat_numeric_t *din, double tol) {
    din->numeric = cs_lu(A->cs, dis->symbolic, tol);
    if (!din->numeric) {
        IGRAPH_ERROR("Cannot do LU decomposition", IGRAPH_FAILURE);
    }
    return 0;
}

/**
 * \function igraph_sparsemat_qr
 * \brief QR decomposition of a sparse matrix.
 *
 * Numeric QR decomposition of a sparse matrix.
 * \param A The input matrix, in column-compressed format.
 * \param dis The result of the symbolic QR analysis, from the
 *    function \ref igraph_sparsemat_symbqr().
 * \param din The result of the decomposition is stored here, it can
 *    be used to solve many linear systems with the same coefficient
 *    matrix and changing right hand sides, using the \ref
 *    igraph_sparsemat_qrresol() function. Once not needed any more,
 *    one should call \ref igraph_sparsemat_numeric_destroy() on it to
 *    free the allocated memory.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_qr(const igraph_sparsemat_t *A,
                        const igraph_sparsemat_symbolic_t *dis,
                        igraph_sparsemat_numeric_t *din) {
    din->numeric = cs_qr(A->cs, dis->symbolic);
    if (!din->numeric) {
        IGRAPH_ERROR("Cannot do QR decomposition", IGRAPH_FAILURE);
    }
    return 0;
}

/**
 * \function igraph_sparsemat_luresol
 * \brief Solves a linear system using a precomputed LU decomposition.
 *
 * Uses the LU decomposition of a matrix to solve linear systems.
 * \param dis The symbolic analysis of the coefficient matrix, the
 *    result of \ref igraph_sparsemat_symblu().
 * \param din The LU decomposition, the result of a call to \ref
 *    igraph_sparsemat_lu().
 * \param b A vector that defines the right hand side of the linear
 *    equation system.
 * \param res An initialized vector, the solution of the linear system
 *    is stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_luresol(const igraph_sparsemat_symbolic_t *dis,
                             const igraph_sparsemat_numeric_t *din,
                             const igraph_vector_t *b,
                             igraph_vector_t *res) {
    int n = din->numeric->L->n;
    igraph_real_t *workspace;

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    workspace = IGRAPH_CALLOC(n, igraph_real_t);
    if (!workspace) {
        IGRAPH_ERROR("Cannot LU (re)solve sparse matrix", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, workspace);

    if (!cs_ipvec(din->numeric->pinv, VECTOR(*res), workspace, n)) {
        IGRAPH_ERROR("Cannot LU (re)solve sparse matrix", IGRAPH_FAILURE);
    }
    if (!cs_lsolve(din->numeric->L, workspace)) {
        IGRAPH_ERROR("Cannot LU (re)solve sparse matrix", IGRAPH_FAILURE);
    }
    if (!cs_usolve(din->numeric->U, workspace)) {
        IGRAPH_ERROR("Cannot LU (re)solve sparse matrix", IGRAPH_FAILURE);
    }
    if (!cs_ipvec(dis->symbolic->q, workspace, VECTOR(*res), n)) {
        IGRAPH_ERROR("Cannot LU (re)solve sparse matrix", IGRAPH_FAILURE);
    }

    IGRAPH_FREE(workspace);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_sparsemat_qrresol
 * \brief Solves a linear system using a precomputed QR decomposition.
 *
 * Solves a linear system using a QR decomposition of its coefficient
 * matrix.
 * \param dis Symbolic analysis of the coefficient matrix, the result
 *    of \ref igraph_sparsemat_symbqr().
 * \param din The QR decomposition of the coefficient matrix, the
 *    result of \ref igraph_sparsemat_qr().
 * \param b Vector, giving the right hand side of the linear equation
 *    system.
 * \param res An initialized vector, the solution is stored here. It
 *    is resized as needed.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_qrresol(const igraph_sparsemat_symbolic_t *dis,
                             const igraph_sparsemat_numeric_t *din,
                             const igraph_vector_t *b,
                             igraph_vector_t *res) {
    int n = din->numeric->L->n;
    igraph_real_t *workspace;
    int k;

    if (res != b) {
        IGRAPH_CHECK(igraph_vector_update(res, b));
    }

    workspace = IGRAPH_CALLOC(dis->symbolic ? dis->symbolic->m2 : 1,
                              igraph_real_t);
    if (!workspace) {
        IGRAPH_ERROR("Cannot QR (re)solve sparse matrix", IGRAPH_FAILURE);
    }
    IGRAPH_FINALLY(igraph_free, workspace);

    if (!cs_ipvec(dis->symbolic->pinv, VECTOR(*res), workspace, n)) {
        IGRAPH_ERROR("Cannot QR (re)solve sparse matrix", IGRAPH_FAILURE);
    }
    for (k = 0; k < n; k++) {
        if (!cs_happly(din->numeric->L, k, din->numeric->B[k], workspace)) {
            IGRAPH_ERROR("Cannot QR (re)solve sparse matrix", IGRAPH_FAILURE);
        }
    }
    if (!cs_usolve(din->numeric->U, workspace)) {
        IGRAPH_ERROR("Cannot QR (re)solve sparse matrix", IGRAPH_FAILURE);
    }
    if (!cs_ipvec(dis->symbolic->q, workspace, VECTOR(*res), n)) {
        IGRAPH_ERROR("Cannot QR (re)solve sparse matrix", IGRAPH_FAILURE);
    }

    IGRAPH_FREE(workspace);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_sparsemat_symbolic_destroy
 * \brief Deallocates memory after a symbolic decomposition.
 *
 * Frees the memory allocated by \ref igraph_sparsemat_symbqr() or
 * \ref igraph_sparsemat_symblu().
 * \param dis The symbolic analysis.
 *
 * Time complexity: O(1).
 */

void igraph_sparsemat_symbolic_destroy(igraph_sparsemat_symbolic_t *dis) {
    cs_sfree(dis->symbolic);
    dis->symbolic = 0;
}

/**
 * \function igraph_sparsemat_numeric_destroy
 * \brief Deallocates memory after a numeric decomposition.
 *
 * Frees the memoty allocated by \ref igraph_sparsemat_qr() or \ref
 * igraph_sparsemat_lu().
 * \param din The LU or QR decomposition.
 *
 * Time complexity: O(1).
 */

void igraph_sparsemat_numeric_destroy(igraph_sparsemat_numeric_t *din) {
    cs_nfree(din->numeric);
    din->numeric = 0;
}

/**
 * \function igraph_matrix_as_sparsemat
 * \brief Converts a dense matrix to a sparse matrix.
 *
 * \param res An uninitialized sparse matrix, the result is stored
 *    here.
 * \param mat The dense input matrix.
 * \param tol Real scalar, the tolerance. Values closer than \p tol to
 *    zero are considered as zero, and will not be included in the
 *    sparse matrix.
 * \return Error code.
 *
 * Time complexity: O(mn), the number of elements in the dense
 * matrix.
 */

int igraph_matrix_as_sparsemat(igraph_sparsemat_t *res,
                               const igraph_matrix_t *mat,
                               igraph_real_t tol) {
    int nrow = (int) igraph_matrix_nrow(mat);
    int ncol = (int) igraph_matrix_ncol(mat);
    int i, j, nzmax = 0;

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            if (fabs(MATRIX(*mat, i, j)) > tol) {
                nzmax++;
            }
        }
    }

    IGRAPH_CHECK(igraph_sparsemat_init(res, nrow, ncol, nzmax));

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            if (fabs(MATRIX(*mat, i, j)) > tol) {
                IGRAPH_CHECK(igraph_sparsemat_entry(res, i, j, MATRIX(*mat, i, j)));
            }
        }
    }

    return 0;
}

static int igraph_i_sparsemat_as_matrix_cc(igraph_matrix_t *res,
                                           const igraph_sparsemat_t *spmat) {

    long int nrow = igraph_sparsemat_nrow(spmat);
    long int ncol = igraph_sparsemat_ncol(spmat);
    CS_INT from = 0, to = 0;
    CS_INT *p = spmat->cs->p;
    CS_INT *i = spmat->cs->i;
    CS_ENTRY *x = spmat->cs->x;
    CS_INT nzmax = spmat->cs->nzmax;

    IGRAPH_CHECK(igraph_matrix_resize(res, nrow, ncol));
    igraph_matrix_null(res);

    while (*p < nzmax) {
        while (to < * (p + 1)) {
            MATRIX(*res, *i, from) += *x;
            to++;
            i++;
            x++;
        }
        from++;
        p++;
    }

    return 0;
}

static int igraph_i_sparsemat_as_matrix_triplet(igraph_matrix_t *res,
                                                const igraph_sparsemat_t *spmat) {
    long int nrow = igraph_sparsemat_nrow(spmat);
    long int ncol = igraph_sparsemat_ncol(spmat);
    CS_INT *i = spmat->cs->p;
    CS_INT *j = spmat->cs->i;
    CS_ENTRY *x = spmat->cs->x;
    CS_INT nz = spmat->cs->nz;
    CS_INT e;

    IGRAPH_CHECK(igraph_matrix_resize(res, nrow, ncol));
    igraph_matrix_null(res);

    for (e = 0; e < nz; e++, i++, j++, x++) {
        MATRIX(*res, *j, *i) += *x;
    }

    return 0;
}

/**
 * \function igraph_sparsemat_as_matrix
 * \brief Converts a sparse matrix to a dense matrix.
 *
 * \param res Pointer to an initialized matrix, the result is stored
 *    here. It will be resized to the required size.
 * \param spmat The input sparse matrix, in triplet or
 *    column-compressed format.
 * \return Error code.
 *
 * Time complexity: O(mn), the number of elements in the dense
 * matrix.
 */

int igraph_sparsemat_as_matrix(igraph_matrix_t *res,
                               const igraph_sparsemat_t *spmat) {
    if (spmat->cs->nz < 0) {
        return (igraph_i_sparsemat_as_matrix_cc(res, spmat));
    } else {
        return (igraph_i_sparsemat_as_matrix_triplet(res, spmat));
    }
}

/**
 * \function igraph_sparsemat_max
 * \brief Maximum of a sparse matrix.
 *
 * \param A The input matrix, column-compressed.
 * \return The maximum in the input matrix, or \c IGRAPH_NEGINFINITY
 *    if the matrix has zero elements.
 *
 * Time complexity: TODO.
 */

igraph_real_t igraph_sparsemat_max(igraph_sparsemat_t *A) {
    CS_INT i, n;
    CS_ENTRY *ptr;
    igraph_real_t res;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    ptr = A->cs->x;
    n = A->cs->nz == -1 ? A->cs->p[A->cs->n] : A->cs->nz;
    if (n == 0) {
        return IGRAPH_NEGINFINITY;
    }
    res = *ptr;
    for (i = 1; i < n; i++, ptr++) {
        if (*ptr > res) {
            res = *ptr;
        }
    }
    return res;
}

/* TODO: CC matrix don't actually need _dupl,
   because the elements are right beside each other.
   Same for max and minmax. */

/**
 * \function igraph_sparsemat_min
 * \brief Minimum of a sparse matrix.
 *
 * \param A The input matrix, column-compressed.
 * \return The minimum in the input matrix, or \c IGRAPH_POSINFINITY
 *    if the matrix has zero elements.
 *
 * Time complexity: TODO.
 */

igraph_real_t igraph_sparsemat_min(igraph_sparsemat_t *A) {
    CS_INT i, n;
    CS_ENTRY *ptr;
    igraph_real_t res;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    ptr = A->cs->x;
    n = A->cs->nz == -1 ? A->cs->p[A->cs->n] : A->cs->nz;
    if (n == 0) {
        return IGRAPH_POSINFINITY;
    }
    res = *ptr;
    for (i = 1; i < n; i++, ptr++) {
        if (*ptr < res) {
            res = *ptr;
        }
    }
    return res;
}

/**
 * \function igraph_sparsemat_minmax
 * \brief Minimum and maximum of a sparse matrix.
 *
 * \param A The input matrix, column-compressed.
 * \param min The minimum in the input matrix is stored here, or \c
 *    IGRAPH_POSINFINITY if the matrix has zero elements.
 * \param max The maximum in the input matrix is stored here, or \c
 *    IGRAPH_NEGINFINITY if the matrix has zero elements.
 * \return Error code.
 *
 * Time complexity: TODO.
 */


int igraph_sparsemat_minmax(igraph_sparsemat_t *A,
                            igraph_real_t *min, igraph_real_t *max) {
    CS_INT i, n;
    CS_ENTRY *ptr;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    ptr = A->cs->x;
    n = A->cs->nz == -1 ? A->cs->p[A->cs->n] : A->cs->nz;
    if (n == 0) {
        *min = IGRAPH_POSINFINITY;
        *max = IGRAPH_NEGINFINITY;
        return 0;
    }
    *min = *max = *ptr;
    for (i = 1; i < n; i++, ptr++) {
        if (*ptr > *max) {
            *max = *ptr;
        } else if (*ptr < *min) {
            *min = *ptr;
        }
    }
    return 0;
}

/**
 * \function igraph_sparsemat_count_nonzero
 * \brief Counts nonzero elements of a sparse matrix.
 *
 * \param A The input matrix, column-compressed.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

long int igraph_sparsemat_count_nonzero(igraph_sparsemat_t *A) {
    CS_INT i, n;
    CS_ENTRY *ptr;
    int res = 0;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    ptr = A->cs->x;
    n = A->cs->nz == -1 ? A->cs->p[A->cs->n] : A->cs->nz;
    if (n == 0) {
        return 0;
    }
    for (i = 0; i < n; i++, ptr++) {
        if (*ptr) {
            res++;
        }
    }
    return res;
}

/**
 * \function igraph_sparsemat_count_nonzerotol
 * \brief Counts nonzero elements of a sparse matrix, ignoring elements close to zero.
 *
 * Count the number of matrix entries that are closer to zero than \p
 * tol.
 * \param The input matrix, column-compressed.
 * \param Real scalar, the tolerance.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

long int igraph_sparsemat_count_nonzerotol(igraph_sparsemat_t *A,
        igraph_real_t tol) {
    CS_INT i, n;
    CS_ENTRY *ptr;
    int res = 0;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    ptr = A->cs->x;
    n = A->cs->nz == -1 ? A->cs->p[A->cs->n] : A->cs->nz;
    if (n == 0) {
        return 0;
    }
    for (i = 0; i < n; i++, ptr++) {
        if (*ptr < - tol || *ptr > tol) {
            res++;
        }
    }
    return res;
}

static int igraph_i_sparsemat_rowsums_triplet(const igraph_sparsemat_t *A,
                                              igraph_vector_t *res) {
    CS_INT i;
    CS_INT *pi = A->cs->i;
    CS_ENTRY *px = A->cs->x;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    igraph_vector_null(res);

    for (i = 0; i < A->cs->nz; i++, pi++, px++) {
        VECTOR(*res)[ *pi ] += *px;
    }

    return 0;
}

static int igraph_i_sparsemat_rowsums_cc(const igraph_sparsemat_t *A,
                                         igraph_vector_t *res) {
    CS_INT ne = A->cs->p[A->cs->n];
    CS_ENTRY *px = A->cs->x;
    CS_INT *pi = A->cs->i;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    igraph_vector_null(res);

    for (; pi < A->cs->i + ne; pi++, px++) {
        VECTOR(*res)[ *pi ] += *px;
    }

    return 0;
}

/**
 * \function igraph_sparsemat_rowsums
 * \brief Row-wise sums.
 *
 * \param A The input matrix, in triplet or column-compressed format.
 * \param res An initialized vector, the result is stored here. It
 *    will be resized as needed.
 * \return Error code.
 *
 * Time complexity: O(nz), the number of non-zero elements.
 */

int igraph_sparsemat_rowsums(const igraph_sparsemat_t *A,
                             igraph_vector_t *res) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_rowsums_triplet(A, res);
    } else {
        return igraph_i_sparsemat_rowsums_cc(A, res);
    }
}

static int igraph_i_sparsemat_rowmins_triplet(const igraph_sparsemat_t *A,
                                              igraph_vector_t *res) {
    CS_INT i;
    CS_INT *pi = A->cs->i;
    CS_ENTRY *px = A->cs->x;
    double inf = IGRAPH_INFINITY;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    igraph_vector_fill(res, inf);

    for (i = 0; i < A->cs->nz; i++, pi++, px++) {
        if (*px < VECTOR(*res)[ *pi ]) {
            VECTOR(*res)[ *pi ] = *px;
        }
    }

    return 0;
}

static int igraph_i_sparsemat_rowmins_cc(igraph_sparsemat_t *A,
                                         igraph_vector_t *res) {
    CS_INT ne;
    CS_ENTRY *px;
    CS_INT *pi;
    double inf = IGRAPH_INFINITY;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    ne = A->cs->p[A->cs->n];
    px = A->cs->x;
    pi = A->cs->i;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    igraph_vector_fill(res, inf);

    for (; pi < A->cs->i + ne; pi++, px++) {
        if (*px < VECTOR(*res)[ *pi ]) {
            VECTOR(*res)[ *pi ] = *px;
        }
    }

    return 0;
}

int igraph_sparsemat_rowmins(igraph_sparsemat_t *A,
                             igraph_vector_t *res) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_rowmins_triplet(A, res);
    } else {
        return igraph_i_sparsemat_rowmins_cc(A, res);
    }
}


static int igraph_i_sparsemat_rowmaxs_triplet(const igraph_sparsemat_t *A,
                                              igraph_vector_t *res) {
    CS_INT i;
    CS_INT *pi = A->cs->i;
    CS_ENTRY *px = A->cs->x;
    double inf = IGRAPH_NEGINFINITY;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    igraph_vector_fill(res, inf);

    for (i = 0; i < A->cs->nz; i++, pi++, px++) {
        if (*px > VECTOR(*res)[ *pi ]) {
            VECTOR(*res)[ *pi ] = *px;
        }
    }

    return 0;
}

static int igraph_i_sparsemat_rowmaxs_cc(igraph_sparsemat_t *A,
                                         igraph_vector_t *res) {
    CS_INT ne;
    CS_ENTRY *px;
    CS_INT *pi;
    double inf = IGRAPH_NEGINFINITY;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    ne = A->cs->p[A->cs->n];
    px = A->cs->x;
    pi = A->cs->i;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    igraph_vector_fill(res, inf);

    for (; pi < A->cs->i + ne; pi++, px++) {
        if (*px > VECTOR(*res)[ *pi ]) {
            VECTOR(*res)[ *pi ] = *px;
        }
    }

    return 0;
}

int igraph_sparsemat_rowmaxs(igraph_sparsemat_t *A,
                             igraph_vector_t *res) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_rowmaxs_triplet(A, res);
    } else {
        return igraph_i_sparsemat_rowmaxs_cc(A, res);
    }
}

static int igraph_i_sparsemat_colmins_triplet(const igraph_sparsemat_t *A,
                                              igraph_vector_t *res) {
    CS_INT i;
    CS_INT *pp = A->cs->p;
    CS_ENTRY *px = A->cs->x;
    double inf = IGRAPH_INFINITY;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->n));
    igraph_vector_fill(res, inf);

    for (i = 0; i < A->cs->nz; i++, pp++, px++) {
        if (*px < VECTOR(*res)[ *pp ]) {
            VECTOR(*res)[ *pp ] = *px;
        }
    }

    return 0;
}

static int igraph_i_sparsemat_colmins_cc(igraph_sparsemat_t *A,
                                         igraph_vector_t *res) {
    CS_INT n;
    CS_ENTRY *px;
    CS_INT *pp;
    CS_INT *pi;
    double *pr;
    double inf = IGRAPH_INFINITY;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    n = A->cs->n;
    px = A->cs->x;
    pp = A->cs->p;
    pi = A->cs->i;

    IGRAPH_CHECK(igraph_vector_resize(res, n));
    igraph_vector_fill(res, inf);
    pr = VECTOR(*res);

    for (; pp < A->cs->p + n; pp++, pr++) {
        for (; pi < A->cs->i + * (pp + 1); pi++, px++) {
            if (*px < *pr) {
                *pr = *px;
            }
        }
    }
    return 0;
}

int igraph_sparsemat_colmins(igraph_sparsemat_t *A,
                             igraph_vector_t *res) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_colmins_triplet(A, res);
    } else {
        return igraph_i_sparsemat_colmins_cc(A, res);
    }
}

static int igraph_i_sparsemat_colmaxs_triplet(const igraph_sparsemat_t *A,
                                              igraph_vector_t *res) {
    CS_INT i;
    CS_INT *pp = A->cs->p;
    CS_ENTRY *px = A->cs->x;
    double inf = IGRAPH_NEGINFINITY;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->n));
    igraph_vector_fill(res, inf);

    for (i = 0; i < A->cs->nz; i++, pp++, px++) {
        if (*px > VECTOR(*res)[ *pp ]) {
            VECTOR(*res)[ *pp ] = *px;
        }
    }

    return 0;
}

static int igraph_i_sparsemat_colmaxs_cc(igraph_sparsemat_t *A,
                                         igraph_vector_t *res) {
    CS_INT n;
    CS_ENTRY *px;
    CS_INT *pp;
    CS_INT *pi;
    double *pr;
    double inf = IGRAPH_NEGINFINITY;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    n = A->cs->n;
    px = A->cs->x;
    pp = A->cs->p;
    pi = A->cs->i;

    IGRAPH_CHECK(igraph_vector_resize(res, n));
    igraph_vector_fill(res, inf);
    pr = VECTOR(*res);

    for (; pp < A->cs->p + n; pp++, pr++) {
        for (; pi < A->cs->i + * (pp + 1); pi++, px++) {
            if (*px > *pr) {
                *pr = *px;
            }
        }
    }
    return 0;
}

int igraph_sparsemat_colmaxs(igraph_sparsemat_t *A,
                             igraph_vector_t *res) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_colmaxs_triplet(A, res);
    } else {
        return igraph_i_sparsemat_colmaxs_cc(A, res);
    }
}

static int igraph_i_sparsemat_which_min_rows_triplet(igraph_sparsemat_t *A,
                                                     igraph_vector_t *res,
                                                     igraph_vector_int_t *pos) {
    CS_INT i;
    CS_INT *pi = A->cs->i;
    CS_INT *pp = A->cs->p;
    CS_ENTRY *px = A->cs->x;
    double inf = IGRAPH_INFINITY;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    IGRAPH_CHECK(igraph_vector_int_resize(pos, A->cs->m));
    igraph_vector_fill(res, inf);
    igraph_vector_int_null(pos);

    for (i = 0; i < A->cs->nz; i++, pi++, px++, pp++) {
        if (*px < VECTOR(*res)[ *pi ]) {
            VECTOR(*res)[ *pi ] = *px;
            VECTOR(*pos)[ *pi ] = *pp;
        }
    }

    return 0;
}

static int igraph_i_sparsemat_which_min_rows_cc(igraph_sparsemat_t *A,
                                                igraph_vector_t *res,
                                                igraph_vector_int_t *pos) {
    CS_INT n;
    CS_ENTRY *px;
    CS_INT *pp;
    CS_INT *pi;
    double inf = IGRAPH_INFINITY;
    int j;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    n = A->cs->n;
    px = A->cs->x;
    pp = A->cs->p;
    pi = A->cs->i;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->m));
    IGRAPH_CHECK(igraph_vector_int_resize(pos, A->cs->m));
    igraph_vector_fill(res, inf);
    igraph_vector_int_null(pos);

    for (j = 0; pp < A->cs->p + n; pp++, j++) {
        for (; pi < A->cs->i + * (pp + 1); pi++, px++) {
            if (*px < VECTOR(*res)[ *pi ]) {
                VECTOR(*res)[ *pi ] = *px;
                VECTOR(*pos)[ *pi ] = j;
            }
        }
    }

    return 0;
}

int igraph_sparsemat_which_min_rows(igraph_sparsemat_t *A,
                                    igraph_vector_t *res,
                                    igraph_vector_int_t *pos) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_which_min_rows_triplet(A, res, pos);
    } else {
        return igraph_i_sparsemat_which_min_rows_cc(A, res, pos);
    }
}

static int igraph_i_sparsemat_which_min_cols_triplet(igraph_sparsemat_t *A,
                                                     igraph_vector_t *res,
                                                     igraph_vector_int_t *pos) {

    CS_INT i;
    CS_INT *pi = A->cs->i;
    CS_INT *pp = A->cs->p;
    CS_ENTRY *px = A->cs->x;
    double inf = IGRAPH_INFINITY;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->n));
    IGRAPH_CHECK(igraph_vector_int_resize(pos, A->cs->n));
    igraph_vector_fill(res, inf);
    igraph_vector_int_null(pos);

    for (i = 0; i < A->cs->nz; i++, pi++, pp++, px++) {
        if (*px < VECTOR(*res)[ *pp ]) {
            VECTOR(*res)[ *pp ] = *px;
            VECTOR(*pos)[ *pp ] = *pi;
        }
    }

    return 0;
}

static int igraph_i_sparsemat_which_min_cols_cc(igraph_sparsemat_t *A,
                                                igraph_vector_t *res,
                                                igraph_vector_int_t *pos) {
    CS_INT n, j, p;
    CS_ENTRY *px;
    double *pr;
    igraph_integer_t *ppos;
    double inf = IGRAPH_INFINITY;

    IGRAPH_CHECK(igraph_sparsemat_dupl(A));

    n = A->cs->n;
    px = A->cs->x;

    IGRAPH_CHECK(igraph_vector_resize(res, n));
    igraph_vector_fill(res, inf);
    pr = VECTOR(*res);
    IGRAPH_CHECK(igraph_vector_int_resize(pos, n));
    igraph_vector_int_null(pos);
    ppos = VECTOR(*pos);

    for (j = 0; j < A->cs->n; j++, pr++, ppos++) {
        for (p = A->cs->p[j]; p < A->cs->p[j + 1]; p++, px++) {
            if (*px < *pr) {
                *pr = *px;
                *ppos = A->cs->i[p];
            }
        }
    }
    return 0;
}

int igraph_sparsemat_which_min_cols(igraph_sparsemat_t *A,
                                    igraph_vector_t *res,
                                    igraph_vector_int_t *pos) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_which_min_cols_triplet(A, res, pos);
    } else {
        return igraph_i_sparsemat_which_min_cols_cc(A, res, pos);
    }
}

static int igraph_i_sparsemat_colsums_triplet(const igraph_sparsemat_t *A,
                                              igraph_vector_t *res) {
    CS_INT i;
    CS_INT *pp = A->cs->p;
    CS_ENTRY *px = A->cs->x;

    IGRAPH_CHECK(igraph_vector_resize(res, A->cs->n));
    igraph_vector_null(res);

    for (i = 0; i < A->cs->nz; i++, pp++, px++) {
        VECTOR(*res)[ *pp ] += *px;
    }

    return 0;
}

static int igraph_i_sparsemat_colsums_cc(const igraph_sparsemat_t *A,
                                         igraph_vector_t *res) {
    CS_INT n = A->cs->n;
    CS_ENTRY *px = A->cs->x;
    CS_INT *pp = A->cs->p;
    CS_INT *pi = A->cs->i;
    double *pr;

    IGRAPH_CHECK(igraph_vector_resize(res, n));
    igraph_vector_null(res);
    pr = VECTOR(*res);

    for (; pp < A->cs->p + n; pp++, pr++) {
        for (; pi < A->cs->i + * (pp + 1); pi++, px++) {
            *pr += *px;
        }
    }
    return 0;
}

/**
 * \function igraph_sparsemat_colsums
 * \brief Column-wise sums.
 *
 * \param A The input matrix, in triplet or column-compressed format.
 * \param res An initialized vector, the result is stored here. It
 *    will be resized as needed.
 * \return Error code.
 *
 * Time complexity: O(nz) for triplet matrices, O(nz+n) for
 * column-compressed ones, nz is the number of non-zero elements, n is
 * the number of columns.
 */

int igraph_sparsemat_colsums(const igraph_sparsemat_t *A,
                             igraph_vector_t *res) {
    if (igraph_sparsemat_is_triplet(A)) {
        return igraph_i_sparsemat_colsums_triplet(A, res);
    } else {
        return igraph_i_sparsemat_colsums_cc(A, res);
    }
}

/**
 * \function igraph_sparsemat_scale
 * \brief Scales a sparse matrix.
 *
 * Multiplies all elements of a sparse matrix, by the given scalar.
 * \param A The input matrix.
 * \param by The scaling factor.
 * \return Error code.
 *
 * Time complexity: O(nz), the number of non-zero elements in the
 * matrix.
 */

int igraph_sparsemat_scale(igraph_sparsemat_t *A, igraph_real_t by) {

    CS_ENTRY *px = A->cs->x;
    CS_INT n = A->cs->nz == -1 ? A->cs->p[A->cs->n] : A->cs->nz;
    CS_ENTRY *stop = px + n;

    for (; px < stop; px++) {
        *px *= by;
    }

    return 0;
}

/**
 * \function igraph_sparsemat_add_rows
 * \brief Adds rows to a sparse matrix.
 *
 * The current matrix elements are retained and all elements in the
 * new rows are zero.
 * \param A The input matrix, in triplet or column-compressed format.
 * \param n The number of rows to add.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_sparsemat_add_rows(igraph_sparsemat_t *A, long int n) {
    A->cs->m += n;
    return 0;
}

/**
 * \function igraph_sparsemat_add_cols
 * \brief Adds columns to a sparse matrix.
 *
 * The current matrix elements are retained, and all elements in the
 * new columns are zero.
 * \param A The input matrix, in triplet or column-compressed format.
 * \param n The number of columns to add.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

int igraph_sparsemat_add_cols(igraph_sparsemat_t *A, long int n) {
    if (igraph_sparsemat_is_triplet(A)) {
        A->cs->n += n;
    } else {
        CS_INT realloc_ok = 0, i;
        CS_INT *newp = cs_realloc(A->cs->p, (A->cs->n + n + 1), sizeof(int), &realloc_ok);
        if (!realloc_ok) {
            IGRAPH_ERROR("Cannot add columns to sparse matrix", IGRAPH_ENOMEM);
        }
        if (newp != A->cs->p) {
            A->cs->p = newp;
        }
        for (i = A->cs->n + 1; i < A->cs->n + n + 1; i++) {
            A->cs->p[i] = A->cs->p[i - 1];
        }
        A->cs->n += n;
    }
    return 0;
}

/**
 * \function igraph_sparsemat_resize
 * \brief Resizes a sparse matrix.
 *
 * This function resizes a sparse matrix. The resized sparse matrix
 * will be empty.
 *
 * \param A The initialized sparse matrix to resize.
 * \param nrow The new number of rows.
 * \param ncol The new number of columns.
 * \param nzmax The new maximum number of elements.
 * \return Error code.
 *
 * Time complexity: O(nzmax), the maximum number of non-zero elements.
 */

int igraph_sparsemat_resize(igraph_sparsemat_t *A, long int nrow,
                            long int ncol, int nzmax) {

    if (A->cs->nz < 0) {
        igraph_sparsemat_t tmp;
        IGRAPH_CHECK(igraph_sparsemat_init(&tmp, (int) nrow, (int) ncol, nzmax));
        igraph_sparsemat_destroy(A);
        *A = tmp;
    } else {
        IGRAPH_CHECK(igraph_sparsemat_realloc(A, nzmax));
        A->cs->m = (int) nrow;
        A->cs->n = (int) ncol;
        A->cs->nz = 0;
    }
    return 0;
}

/**
 * \function igraph_sparsemat_nonzero_storage
 * \brief Returns number of stored entries of a sparse matrix.
 *
 * This function will return the number of stored entries of a sparse
 * matrix. These entries can be zero, and multiple entries can be
 * at the same position. Use \ref igraph_sparsemat_dupl() to sum
 * duplicate entries, and \ref igraph_sparsemat_dropzeros() to remove
 * zeros.
 *
 * \param A A sparse matrix in either triplet or compressed form.
 * \return Number of stored entries.
 *
 * Time complexity: O(1).
 */

int igraph_sparsemat_nonzero_storage(const igraph_sparsemat_t *A) {
    if (A->cs->nz < 0) {
        return A->cs->p[A->cs->n];
    } else {
        return A->cs->nz;
    }
}

int igraph_sparsemat_getelements(const igraph_sparsemat_t *A,
                                 igraph_vector_int_t *i,
                                 igraph_vector_int_t *j,
                                 igraph_vector_t *x) {
    CS_INT nz = A->cs->nz;
    if (nz < 0) {
        nz = A->cs->p[A->cs->n];
        IGRAPH_CHECK(igraph_vector_int_resize(i, nz));
        IGRAPH_CHECK(igraph_vector_int_resize(j, A->cs->n + 1));
        IGRAPH_CHECK(igraph_vector_resize(x, nz));
        memcpy(VECTOR(*i), A->cs->i, (size_t) nz * sizeof(int));
        memcpy(VECTOR(*j), A->cs->p, (size_t) (A->cs->n + 1) * sizeof(int));
        memcpy(VECTOR(*x), A->cs->x, (size_t) nz * sizeof(igraph_real_t));
    } else {
        IGRAPH_CHECK(igraph_vector_int_resize(i, nz));
        IGRAPH_CHECK(igraph_vector_int_resize(j, nz));
        IGRAPH_CHECK(igraph_vector_resize(x, nz));
        memcpy(VECTOR(*i), A->cs->i, (size_t) nz * sizeof(int));
        memcpy(VECTOR(*j), A->cs->p, (size_t) nz * sizeof(int));
        memcpy(VECTOR(*x), A->cs->x, (size_t) nz * sizeof(igraph_real_t));
    }
    return 0;
}

int igraph_sparsemat_scale_rows(igraph_sparsemat_t *A,
                                const igraph_vector_t *fact) {
    CS_INT *i = A->cs->i;
    CS_ENTRY *x = A->cs->x;
    CS_INT no_of_edges = A->cs->nz < 0 ? A->cs->p[A->cs->n] : A->cs->nz;
    CS_INT e;

    for (e = 0; e < no_of_edges; e++, x++, i++) {
        igraph_real_t f = VECTOR(*fact)[*i];
        (*x) *= f;
    }

    return 0;
}

static int igraph_i_sparsemat_scale_cols_cc(igraph_sparsemat_t *A,
                                            const igraph_vector_t *fact) {
    CS_INT *i = A->cs->i;
    CS_ENTRY *x = A->cs->x;
    CS_INT no_of_edges = A->cs->p[A->cs->n];
    CS_INT e;
    CS_INT c = 0;        /* actual column */

    for (e = 0; e < no_of_edges; e++, x++, i++) {
        igraph_real_t f;
        while (c < A->cs->n && A->cs->p[c + 1] == e) {
            c++;
        }
        f = VECTOR(*fact)[c];
        (*x) *= f;
    }

    return 0;
}

static int igraph_i_sparsemat_scale_cols_triplet(igraph_sparsemat_t *A,
                                                 const igraph_vector_t *fact) {
    CS_INT *j = A->cs->p;
    CS_ENTRY *x = A->cs->x;
    CS_INT no_of_edges = A->cs->nz;
    CS_INT e;

    for (e = 0; e < no_of_edges; e++, x++, j++) {
        igraph_real_t f = VECTOR(*fact)[*j];
        (*x) *= f;
    }

    return 0;
}

int igraph_sparsemat_scale_cols(igraph_sparsemat_t *A,
                                const igraph_vector_t *fact) {
    if (A->cs->nz < 0) {
        return igraph_i_sparsemat_scale_cols_cc(A, fact);
    } else {
        return igraph_i_sparsemat_scale_cols_triplet(A, fact);
    }
}

int igraph_sparsemat_multiply_by_dense(const igraph_sparsemat_t *A,
                                       const igraph_matrix_t *B,
                                       igraph_matrix_t *res) {

    int m = (int) igraph_sparsemat_nrow(A);
    int n = (int) igraph_sparsemat_ncol(A);
    int p = (int) igraph_matrix_ncol(B);
    int i;

    if (igraph_matrix_nrow(B) != n) {
        IGRAPH_ERROR("Invalid dimensions in sparse-dense matrix product",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, m, p));
    igraph_matrix_null(res);

    for (i = 0; i < p; i++) {
        if (!(cs_gaxpy(A->cs, &MATRIX(*B, 0, i), &MATRIX(*res, 0, i)))) {
            IGRAPH_ERROR("Cannot perform sparse-dense matrix multiplication",
                         IGRAPH_FAILURE);
        }
    }

    return 0;
}

int igraph_sparsemat_dense_multiply(const igraph_matrix_t *A,
                                    const igraph_sparsemat_t *B,
                                    igraph_matrix_t *res) {
    int m = (int) igraph_matrix_nrow(A);
    int n = (int) igraph_matrix_ncol(A);
    int p = (int) igraph_sparsemat_ncol(B);
    int r, c;
    CS_INT *Bp = B->cs->p;

    if (igraph_sparsemat_nrow(B) != n) {
        IGRAPH_ERROR("Invalid dimensions in dense-sparse matrix product",
                     IGRAPH_EINVAL);
    }

    if (!igraph_sparsemat_is_cc(B)) {
        IGRAPH_ERROR("Dense-sparse product is only implemented for "
                     "column-compressed sparse matrices", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_matrix_resize(res, m, p));
    igraph_matrix_null(res);

    for (c = 0; c < p; c++) {
        for (r = 0; r < m; r++) {
            int idx = *Bp;
            while (idx < * (Bp + 1)) {
                MATRIX(*res, r, c) += MATRIX(*A, r, B->cs->i[idx]) * B->cs->x[idx];
                idx++;
            }
        }
        Bp++;
    }

    return 0;
}

/**
 * \function igraph_sparsemat_view
 * \brief Initialize a sparse matrix and set all parameters.
 *
 * This function can be used to temporarily handle existing sparse matrix data,
 * usually created by another software library, as an \c igraph_sparsemat_t object,
 * and thus avoid unnecessary copying. It supports data stored in either the
 * compressed sparse column format, or the <code>(i, j, x)</code> triplet format
 * where \c i and \c j are the matrix indices of a non-zero element, and \c x
 * is its value.
 *
 * </para><para>
 * The compressed sparse column (or row) format is commonly used to represent
 * sparse matrix data. It consists of three vectors, the \p p column pointers, the
 * \p i row indices, and the \p x values. <code>p[k]</code> is the number
 * of non-zero entires in matrix columns <code>k-1</code> and lower.
 * <code>p[0]</code> is always zero and <code>p[n]</code> is always the total
 * number of non-zero entires in the matrix. <code>i[l]</code> is the row index
 * of the \c l-th stored element, while <code>x[l]</code> is its value.
 * If a matrix element with indices <code>(j, k)</code> is explicitly stored,
 * it must be located between positions <code>p[k]</code> and <code>p[k+1] - 1</code>
 * (inclusive) in the \p i and \p x vectors.
 *
 * </para><para>
 * Do not call \ref igraph_sparsemat_destroy() on a sparse matrix created with
 * this function. Instead, \ref igraph_free() must be called on the \c cs
 * field of \p A to free the storage allocated by this function.
 *
 * </para><para>
 * Warning: Matrices created with this function must not be used with functions
 * that may reallocate the underlying storage, such as \ref igraph_sparsemat_entry().
 *
 * \param A The non-initialized sparse matrix.
 * \param nzmax The maximum number of entries, typically the actual number of entries.
 * \param m The number of matrix rows.
 * \param n The number of matrix columns.
 * \param p For a compressed matrix, this is the column pointer vector, and
 *          must be of size <code>n+1</code>. For a triplet format matrix, it
 *          is a vector of column indices and must be of size \p nzmax.
 * \param i The row vector. This should contain the row indices of the
 *          elements in \p x. It must be of size \p nzmax.
 * \param x The values of the non-zero elements of the sparse matrix.
 *          It must be of size \p nzmax.
 * \param nz For a compressed matrix, is must be -1. For a triplet format
 *           matrix, is must contain the number of entries.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_sparsemat_view(igraph_sparsemat_t *A, int nzmax, int m, int n,
                          int *p, int *i, double *x, int nz) {

    A->cs = IGRAPH_CALLOC(1, cs_di);
    A->cs->nzmax = nzmax;
    A->cs->m = m;
    A->cs->n = n;
    A->cs->p = p;
    A->cs->i = i;
    A->cs->x = x;
    A->cs->nz = nz;

    return IGRAPH_SUCCESS;
}

int igraph_i_sparsemat_view(igraph_sparsemat_t *A, int nzmax, int m, int n,
                            int *p, int *i, double *x, int nz) {
    IGRAPH_WARNING("igraph_i_sparsemat_view() is deprecated, use igraph_sparsemat_view()");
    return igraph_sparsemat_view(A, nzmax, m, n, p, i, x, nz);
}

int igraph_sparsemat_sort(const igraph_sparsemat_t *A,
                          igraph_sparsemat_t *sorted) {

    igraph_sparsemat_t tmp;

    IGRAPH_CHECK(igraph_sparsemat_transpose(A, &tmp, /*values=*/ 1));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
    IGRAPH_CHECK(igraph_sparsemat_transpose(&tmp, sorted, /*values=*/ 1));
    igraph_sparsemat_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

/**
 * \function igraph_sparsemat_getelements_sorted
 * \brief Returns the sorted elements of a sparse matrix.
 *
 * This function will sort a sparse matrix and return the elements in
 * 3 vectors. Two vectors will indicate where the elements are located,
 * and one will give the elements.
 *
 * \param A A sparse matrix in either triplet or compressed form.
 * \param i An initialized int vector. This will store the rows of the
 *          returned elements.
 * \param j An initialized int vector. For a triplet matrix this will
 *          store the columns of the returned elements. For a compressed
 *          matrix, if the column index is \c k, then <code>j[k]</code>
 *          is the index in \p x of the start of the \c k-th column, and
 *          the last element of \c j is the total number of elements.
 *          The total number of elements in the \c k-th column is
 *          therefore <code>j[k+1] - j[k]</code>. For example, if there
 *          is one element in the first column, and five in the second,
 *          \c j will be set to <code>{0, 1, 6}</code>.
 * \param x An initialized vector. The elements will be placed here.
 * \return Error code.
 *
 * Time complexity: O(n), the number of stored elements in the sparse matrix.
 */

int igraph_sparsemat_getelements_sorted(const igraph_sparsemat_t *A,
                                        igraph_vector_int_t *i,
                                        igraph_vector_int_t *j,
                                        igraph_vector_t *x) {
    if (A->cs->nz < 0) {
        igraph_sparsemat_t tmp;
        IGRAPH_CHECK(igraph_sparsemat_sort(A, &tmp));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
        IGRAPH_CHECK(igraph_sparsemat_getelements(&tmp, i, j, x));
        igraph_sparsemat_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        IGRAPH_CHECK(igraph_sparsemat_getelements(A, i, j, x));
    }

    return IGRAPH_SUCCESS;
}

int igraph_sparsemat_nzmax(const igraph_sparsemat_t *A) {
    return A->cs->nzmax;
}

int igraph_sparsemat_neg(igraph_sparsemat_t *A) {
    CS_INT i, nz = A->cs->nz == -1 ? A->cs->p[A->cs->n] : A->cs->nz;
    CS_ENTRY *px = A->cs->x;

    for (i = 0; i < nz; i++, px++) {
        *px = - (*px);
    }

    return 0;
}

/**
 * \function igraph_sparsemat_iterator_init
 * \brief Initialize a sparse matrix iterator.
 *
 * \param it A pointer to an uninitialized sparse matrix iterator.
 * \param sparsemat Pointer to the sparse matrix.
 * \return Error code. This will always return \c IGRAPH_SUCCESS
 *
 * Time complexity: O(n), the number of columns of the sparse matrix.
 */

int igraph_sparsemat_iterator_init(igraph_sparsemat_iterator_t *it,
                                   igraph_sparsemat_t *sparsemat) {

    it->mat = sparsemat;
    igraph_sparsemat_iterator_reset(it);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_sparsemat_iterator_reset
 * \brief Reset a sparse matrix iterator to the first element.
 *
 * \param it A pointer to the sparse matrix iterator.
 * \return Error code. This will always return \c IGRAPH_SUCCESS
 *
 * Time complexity: O(n), the number of columns of the sparse matrix.
 */

int igraph_sparsemat_iterator_reset(igraph_sparsemat_iterator_t *it) {
    it->pos = 0;
    it->col = 0;
    if (!igraph_sparsemat_is_triplet(it->mat)) {
        while (it->col < it->mat->cs->n &&
               it->mat->cs->p[it->col + 1] == it->pos) {
            it->col ++;
        }
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_sparsemat_iterator_end
 * \brief Query if the iterator is past the last element.
 *
 * \param it A pointer to the sparse matrix iterator.
 * \return true if the iterator is past the last element, false if it
 *         points to an element in a sparse matrix.
 *
 * Time complexity: O(1).
 */

igraph_bool_t
igraph_sparsemat_iterator_end(const igraph_sparsemat_iterator_t *it) {
    CS_INT nz = it->mat->cs->nz == -1 ? it->mat->cs->p[it->mat->cs->n] :
                it->mat->cs->nz;
    return it->pos >= nz;
}

/**
 * \function igraph_sparsemat_iterator_row
 * \brief Return the row of the iterator.
 *
 * \param it A pointer to the sparse matrix iterator.
 * \return The row of the element at the current iterator position.
 *
 * Time complexity: O(1).
 */

int igraph_sparsemat_iterator_row(const igraph_sparsemat_iterator_t *it) {
    return it->mat->cs->i[it->pos];
}

/**
 * \function igraph_sparsemat_iterator_col
 * \brief Return the column of the iterator.
 *
 * \param it A pointer to the sparse matrix iterator.
 * \return The column of the element at the current iterator position.
 *
 * Time complexity: O(1).
 */

int igraph_sparsemat_iterator_col(const igraph_sparsemat_iterator_t *it) {
    if (igraph_sparsemat_is_triplet(it->mat)) {
        return it->mat->cs->p[it->pos];
    } else {
        return it->col;
    }
}

/**
 * \function igraph_sparsemat_iterator_get
 * \brief Return the element at the current iterator position.
 *
 * \param it A pointer to the sparse matrix iterator.
 * \return The value of the element at the current iterator position.
 *
 * Time complexity: O(1).
 */

igraph_real_t
igraph_sparsemat_iterator_get(const igraph_sparsemat_iterator_t *it) {
    return it->mat->cs->x[it->pos];
}

/**
 * \function igraph_sparsemat_iterator_next
 * \brief Let a sparse matrix iterator go to the next element.
 *
 * \param it A pointer to the sparse matrix iterator.
 * \return The position of the iterator in the element vector.
 *
 * Time complexity: O(n), the number of columns of the sparse matrix.
 */

int igraph_sparsemat_iterator_next(igraph_sparsemat_iterator_t *it) {
    it->pos += 1;
    while (it->col < it->mat->cs->n &&
           it->mat->cs->p[it->col + 1] == it->pos) {
        it->col++;
    }
    return it->pos;
}

/**
 * \function igraph_sparsemat_iterator_idx
 * \brief Returns the element vector index of a sparse matrix iterator.
 *
 * \param it A pointer to the sparse matrix iterator.
 * \return The position of the iterator in the element vector.
 *
 * Time complexity: O(1).
 */

int igraph_sparsemat_iterator_idx(const igraph_sparsemat_iterator_t *it) {
    return it->pos;
}
