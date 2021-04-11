/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_types.h"
#include "igraph_spmatrix.h"
#include "igraph_error.h"

#include <string.h>     /* memcpy & co. */

/**
 * \section igraph_spmatrix_constructor_and_destructor Sparse matrix constructors
 * and destructors.
 */

/**
 * \ingroup matrix
 * \function igraph_spmatrix_init
 * \brief Initializes a sparse matrix.
 *
 * </para><para>
 * Every sparse matrix needs to be initialized before using it, this is done
 * by calling this function. A matrix has to be destroyed if it is not
 * needed any more, see \ref igraph_spmatrix_destroy().
 * \param m Pointer to a not yet initialized sparse matrix object to be
 *        initialized.
 * \param nrow The number of rows in the matrix.
 * \param ncol The number of columns in the matrix.
 * \return Error code.
 *
 * Time complexity: operating system dependent.
 */

int igraph_spmatrix_init(igraph_spmatrix_t *m, long int nrow, long int ncol) {
    IGRAPH_ASSERT(m != NULL);
    IGRAPH_VECTOR_INIT_FINALLY(&m->ridx, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&m->cidx, ncol + 1);
    IGRAPH_VECTOR_INIT_FINALLY(&m->data, 0);
    IGRAPH_FINALLY_CLEAN(3);
    m->nrow = nrow;
    m->ncol = ncol;
    return 0;
}

/**
 * \ingroup matrix
 * \function igraph_spmatrix_destroy
 * \brief Destroys a sparse matrix object.
 *
 * </para><para>
 * This function frees all the memory allocated for a sparse matrix
 * object. The destroyed object needs to be reinitialized before using
 * it again.
 * \param m The matrix to destroy.
 *
 * Time complexity: operating system dependent.
 */

void igraph_spmatrix_destroy(igraph_spmatrix_t *m) {
    IGRAPH_ASSERT(m != NULL);
    igraph_vector_destroy(&m->ridx);
    igraph_vector_destroy(&m->cidx);
    igraph_vector_destroy(&m->data);
}

/**
 * \ingroup matrix
 * \function igraph_spmatrix_copy
 * \brief Copies a sparse matrix.
 *
 * </para><para>
 * Creates a sparse matrix object by copying another one.
 * \param to Pointer to an uninitialized sparse matrix object.
 * \param from The initialized sparse matrix object to copy.
 * \return Error code, \c IGRAPH_ENOMEM if there
 *   isn't enough memory to allocate the new sparse matrix.
 *
 * Time complexity: O(n), the number
 * of elements in the matrix.
 */

int igraph_spmatrix_copy(igraph_spmatrix_t *to, const igraph_spmatrix_t *from) {
    IGRAPH_ASSERT(from != NULL);
    IGRAPH_ASSERT(to != NULL);
    to->nrow = from->nrow;
    to->ncol = from->ncol;
    IGRAPH_CHECK(igraph_vector_copy(&to->ridx, &from->ridx));
    IGRAPH_CHECK(igraph_vector_copy(&to->cidx, &from->cidx));
    IGRAPH_CHECK(igraph_vector_copy(&to->data, &from->data));
    return 0;
}

/**
 * \section igraph_spmatrix_accessing_elements Accessing elements of a sparse matrix
 */

/**
 * \ingroup matrix
 * \function igraph_spmatrix_e
 * \brief Accessing an element of a sparse matrix.
 *
 * Note that there are no range checks right now.
 * \param m The matrix object.
 * \param row The index of the row, starting with zero.
 * \param col The index of the column, starting with zero.
 *
 * Time complexity: O(log n), where n is the number of nonzero elements in
 * the requested column.
 */
igraph_real_t igraph_spmatrix_e(const igraph_spmatrix_t *m,
                                long int row, long int col) {
    long int start, end;

    IGRAPH_ASSERT(m != NULL);
    start = (long) VECTOR(m->cidx)[col];
    end = (long) VECTOR(m->cidx)[col + 1] - 1;

    if (end < start) {
        return 0;
    }
    /* Elements residing in column col are between m->data[start] and
     * m->data[end], inclusive, ordered by row index */
    while (start < end - 1) {
        long int mid = (start + end) / 2;
        if (VECTOR(m->ridx)[mid] > row) {
            end = mid;
        } else if (VECTOR(m->ridx)[mid] < row) {
            start = mid;
        } else {
            start = mid;
            break;
        }
    }

    if (VECTOR(m->ridx)[start] == row) {
        return VECTOR(m->data)[start];
    }
    if (VECTOR(m->ridx)[start] != row && VECTOR(m->ridx)[end] == row) {
        return VECTOR(m->data)[end];
    }
    return 0;
}


/**
 * \ingroup matrix
 * \function igraph_spmatrix_set
 * \brief Setting an element of a sparse matrix.
 *
 * Note that there are no range checks right now.
 * \param m The matrix object.
 * \param row The index of the row, starting with zero.
 * \param col The index of the column, starting with zero.
 * \param value The new value.
 *
 * Time complexity: O(log n), where n is the number of nonzero elements in
 * the requested column.
 */
int igraph_spmatrix_set(igraph_spmatrix_t *m, long int row, long int col,
                        igraph_real_t value) {
    long int start, end;

    IGRAPH_ASSERT(m != NULL);
    start = (long) VECTOR(m->cidx)[col];
    end = (long) VECTOR(m->cidx)[col + 1] - 1;

    if (end < start) {
        /* First element in the column */
        if (value == 0.0) {
            return 0;
        }
        IGRAPH_CHECK(igraph_vector_insert(&m->ridx, start, row));
        IGRAPH_CHECK(igraph_vector_insert(&m->data, start, value));
        for (start = col + 1; start < m->ncol + 1; start++) {
            VECTOR(m->cidx)[start]++;
        }
        return 0;
    }

    /* Elements residing in column col are between m->data[start] and
     * m->data[end], inclusive, ordered by row index */
    while (start < end - 1) {
        long int mid = (start + end) / 2;
        if (VECTOR(m->ridx)[mid] > row) {
            end = mid;
        } else if (VECTOR(m->ridx)[mid] < row) {
            start = mid;
        } else {
            start = mid;
            break;
        }
    }

    if (VECTOR(m->ridx)[start] == row) {
        /* Overwriting a value - or deleting it if it has been overwritten by zero */
        if (value == 0) {
            igraph_vector_remove(&m->ridx, start);
            igraph_vector_remove(&m->data, start);
            for (start = col + 1; start < m->ncol + 1; start++) {
                VECTOR(m->cidx)[start]--;
            }
        } else {
            VECTOR(m->data)[start] = value;
        }
        return 0;
    } else if (VECTOR(m->ridx)[end] == row) {
        /* Overwriting a value - or deleting it if it has been overwritten by zero */
        if (value == 0) {
            igraph_vector_remove(&m->ridx, end);
            igraph_vector_remove(&m->data, end);
            for (start = col + 1; start < m->ncol + 1; start++) {
                VECTOR(m->cidx)[start]--;
            }
        } else {
            VECTOR(m->data)[end] = value;
        }
        return 0;
    }

    /* New element has to be inserted, but only if not a zero is
     * being written into the matrix */
    if (value != 0.0) {
        if (VECTOR(m->ridx)[end] < row) {
            IGRAPH_CHECK(igraph_vector_insert(&m->ridx, end + 1, row));
            IGRAPH_CHECK(igraph_vector_insert(&m->data, end + 1, value));
        } else if (VECTOR(m->ridx)[start] < row) {
            IGRAPH_CHECK(igraph_vector_insert(&m->ridx, start + 1, row));
            IGRAPH_CHECK(igraph_vector_insert(&m->data, start + 1, value));
        } else {
            IGRAPH_CHECK(igraph_vector_insert(&m->ridx, start, row));
            IGRAPH_CHECK(igraph_vector_insert(&m->data, start, value));
        }
        for (start = col + 1; start < m->ncol + 1; start++) {
            VECTOR(m->cidx)[start]++;
        }
    }
    return 0;
}


/**
 * \ingroup matrix
 * \function igraph_spmatrix_add_e
 * \brief Adding a real value to an element of a sparse matrix.
 *
 * Note that there are no range checks right now. This is implemented to avoid
 * double lookup of a given element in the matrix by using \ref igraph_spmatrix_e()
 * and \ref igraph_spmatrix_set() consecutively.
 *
 * \param m The matrix object.
 * \param row The index of the row, starting with zero.
 * \param col The index of the column, starting with zero.
 * \param value The value to add.
 *
 * Time complexity: O(log n), where n is the number of nonzero elements in
 * the requested column.
 */
int igraph_spmatrix_add_e(igraph_spmatrix_t *m, long int row, long int col,
                          igraph_real_t value) {
    long int start, end;

    IGRAPH_ASSERT(m != NULL);
    start = (long) VECTOR(m->cidx)[col];
    end = (long) VECTOR(m->cidx)[col + 1] - 1;

    if (end < start) {
        /* First element in the column */
        if (value == 0.0) {
            return 0;
        }
        IGRAPH_CHECK(igraph_vector_insert(&m->ridx, start, row));
        IGRAPH_CHECK(igraph_vector_insert(&m->data, start, value));
        for (start = col + 1; start < m->ncol + 1; start++) {
            VECTOR(m->cidx)[start]++;
        }
        return 0;
    }

    /* Elements residing in column col are between m->data[start] and
     * m->data[end], inclusive, ordered by row index */
    while (start < end - 1) {
        long int mid = (start + end) / 2;
        if (VECTOR(m->ridx)[mid] > row) {
            end = mid;
        } else if (VECTOR(m->ridx)[mid] < row) {
            start = mid;
        } else {
            start = mid;
            break;
        }
    }

    if (VECTOR(m->ridx)[start] == row) {
        /* Overwriting a value */
        if (VECTOR(m->data)[start] == -1) {
            igraph_vector_remove(&m->ridx, start);
            igraph_vector_remove(&m->data, start);
            for (start = col + 1; start < m->ncol + 1; start++) {
                VECTOR(m->cidx)[start]--;
            }
        } else {
            VECTOR(m->data)[start] += value;
        }
        return 0;
    } else if (VECTOR(m->ridx)[end] == row) {
        /* Overwriting a value */
        if (VECTOR(m->data)[end] == -1) {
            igraph_vector_remove(&m->ridx, end);
            igraph_vector_remove(&m->data, end);
            for (start = col + 1; start < m->ncol + 1; start++) {
                VECTOR(m->cidx)[start]--;
            }
        } else {
            VECTOR(m->data)[end] += value;
        }
        return 0;
    }

    /* New element has to be inserted, but only if not a zero is
     * being added to a zero element of the matrix */
    if (value != 0.0) {
        if (VECTOR(m->ridx)[end] < row) {
            IGRAPH_CHECK(igraph_vector_insert(&m->ridx, end + 1, row));
            IGRAPH_CHECK(igraph_vector_insert(&m->data, end + 1, value));
        } else if (VECTOR(m->ridx)[start] < row) {
            IGRAPH_CHECK(igraph_vector_insert(&m->ridx, start + 1, row));
            IGRAPH_CHECK(igraph_vector_insert(&m->data, start + 1, value));
        } else {
            IGRAPH_CHECK(igraph_vector_insert(&m->ridx, start, row));
            IGRAPH_CHECK(igraph_vector_insert(&m->data, start, value));
        }
        for (start = col + 1; start < m->ncol + 1; start++) {
            VECTOR(m->cidx)[start]++;
        }
    }
    return 0;
}

/**
 * \function igraph_spmatrix_add_col_values
 * \brief Adds the values of a column to another column.
 *
 * \param to The index of the column to be added to.
 * \param from The index of the column to be added.
 * \return Error code.
 */
int igraph_spmatrix_add_col_values(igraph_spmatrix_t *m, long int to, long int from) {
    long int i;
    if (to < 0 || to >= m->ncol) {
       IGRAPH_ERROR("The 'to' column does not exist.", IGRAPH_EINVAL);
    }
    if (from < 0 || from >= m->ncol) {
       IGRAPH_ERROR("The 'from' column does not exist.", IGRAPH_EINVAL);
    }
    /* TODO: I think this implementation could be speeded up if I don't use
     * igraph_spmatrix_add_e directly -- but maybe it's not worth the fuss */
    for (i = (long int) VECTOR(m->cidx)[from]; i < VECTOR(m->cidx)[from + 1]; i++) {
        IGRAPH_CHECK(igraph_spmatrix_add_e(m, (long int) VECTOR(m->ridx)[i],
                                           to, VECTOR(m->data)[i]));
    }

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup matrix
 * \function igraph_spmatrix_resize
 * \brief Resizes a sparse matrix.
 *
 * </para><para>
 * This function resizes a sparse matrix by adding more elements to it.
 * The matrix retains its data even after resizing it, except for the data
 * which lies outside the new boundaries (if the new size is smaller).
 * \param m Pointer to an already initialized sparse matrix object.
 * \param nrow The number of rows in the resized matrix.
 * \param ncol The number of columns in the resized matrix.
 * \return Error code.
 *
 * Time complexity: O(n).
 * n is the number of elements in the old matrix.
 */

int igraph_spmatrix_resize(igraph_spmatrix_t *m, long int nrow, long int ncol) {
    long int i, j, ci, ei, mincol;
    IGRAPH_ASSERT(m != NULL);
    /* Iterating through the matrix data and deleting unnecessary data. */
    /* At the same time, we create the new indices as well */
    if (nrow < m->nrow) {
        ei = j = 0;
        mincol = (m->ncol < ncol) ? m->ncol : ncol;
        for (ci = 0; ci < mincol; ci++) {
            for (; ei < VECTOR(m->cidx)[ci + 1]; ei++) {
                if (VECTOR(m->ridx)[ei] < nrow) {
                    VECTOR(m->ridx)[j] = VECTOR(m->ridx)[ei];
                    VECTOR(m->data)[j] = VECTOR(m->data)[ei];
                    j++;
                }
            }
            VECTOR(m->cidx)[ci] = j;
        }
        /* Contract the row index and the data vector */
        IGRAPH_CHECK(igraph_vector_resize(&m->ridx, j));
        IGRAPH_CHECK(igraph_vector_resize(&m->cidx, j));
    }
    /* Updating cidx */
    IGRAPH_CHECK(igraph_vector_resize(&m->cidx, ncol + 1));
    for (i = m->ncol + 1; i < ncol + 1; i++) {
        VECTOR(m->cidx)[i] = VECTOR(m->cidx)[m->ncol];
    }
    m->nrow = nrow;
    m->ncol = ncol;
    return 0;
}

/**
 * \ingroup matrix
 * \function igraph_spmatrix_count_nonzero
 * \brief The number of non-zero elements in a sparse matrix.
 *
 * \param m Pointer to an initialized sparse matrix object.
 * \return The size of the matrix.
 *
 * Time complexity: O(1).
 */

long int igraph_spmatrix_count_nonzero(const igraph_spmatrix_t *m) {
    IGRAPH_ASSERT(m != NULL);
    return igraph_vector_size(&m->data);
}


/**
 * \ingroup matrix
 * \function igraph_spmatrix_size
 * \brief The number of elements in a sparse matrix.
 *
 * \param m Pointer to an initialized sparse matrix object.
 * \return The size of the matrix.
 *
 * Time complexity: O(1).
 */

long int igraph_spmatrix_size(const igraph_spmatrix_t *m) {
    IGRAPH_ASSERT(m != NULL);
    return (m->nrow) * (m->ncol);
}

/**
 * \ingroup matrix
 * \function igraph_spmatrix_nrow
 * \brief The number of rows in a sparse matrix.
 *
 * \param m Pointer to an initialized sparse matrix object.
 * \return The number of rows in the matrix.
 *
 * Time complexity: O(1).
 */

long int igraph_spmatrix_nrow(const igraph_spmatrix_t *m) {
    IGRAPH_ASSERT(m != NULL);
    return m->nrow;
}

/**
 * \ingroup matrix
 * \function igraph_spmatrix_ncol
 * \brief The number of columns in a sparse matrix.
 *
 * \param m Pointer to an initialized sparse matrix object.
 * \return The number of columns in the sparse matrix.
 *
 * Time complexity: O(1).
 */

long int igraph_spmatrix_ncol(const igraph_spmatrix_t *m) {
    IGRAPH_ASSERT(m != NULL);
    return m->ncol;
}

/**
 * \ingroup matrix
 * \brief Copies a sparse matrix to a regular C array.
 *
 * </para><para>
 * The matrix is copied columnwise, as this is the format most
 * programs and languages use.
 * The C array should be of sufficient size, there are (of course) no
 * range checks done.
 * \param m Pointer to an initialized sparse matrix object.
 * \param to Pointer to a C array, the place to copy the data to.
 * \return Error code.
 *
 * Time complexity: O(n),
 * n is the number of
 * elements in the matrix.
 */

int igraph_spmatrix_copy_to(const igraph_spmatrix_t *m, igraph_real_t *to) {
    long int c, dest_idx, idx;

    memset(to, 0, sizeof(igraph_real_t) * (size_t) igraph_spmatrix_size(m));
    for (c = 0, dest_idx = 0; c < m->ncol; c++, dest_idx += m->nrow) {
        for (idx = (long int) VECTOR(m->cidx)[c]; idx < VECTOR(m->cidx)[c + 1]; idx++) {
            to[dest_idx + (long)VECTOR(m->ridx)[idx]] = VECTOR(m->data)[idx];
        }
    }
    return 0;
}

/**
 * \ingroup matrix
 * \brief Sets all element in a sparse matrix to zero.
 *
 * \param m Pointer to an initialized matrix object.
 * \return Error code, always returns with success.
 *
 * Time complexity: O(n),
 * n is the number of columns in the matrix
 */

int igraph_spmatrix_null(igraph_spmatrix_t *m) {
    IGRAPH_ASSERT(m != NULL);
    igraph_vector_clear(&m->data);
    igraph_vector_clear(&m->ridx);
    igraph_vector_null(&m->cidx);
    return 0;
}

/**
 * \ingroup matrix
 * \function igraph_spmatrix_add_cols
 * \brief Adds columns to a sparse matrix.
 * \param m The sparse matrix object.
 * \param n The number of columns to add.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_spmatrix_add_cols(igraph_spmatrix_t *m, long int n) {
    igraph_spmatrix_resize(m, m->nrow, m->ncol + n);
    return 0;
}

/**
 * \ingroup matrix
 * \function igraph_spmatrix_add_rows
 * \brief Adds rows to a sparse matrix.
 * \param m The sparse matrix object.
 * \param n The number of rows to add.
 * \return Error code.
 *
 * Time complexity: O(1).
 */

int igraph_spmatrix_add_rows(igraph_spmatrix_t *m, long int n) {
    igraph_spmatrix_resize(m, m->nrow + n, m->ncol);
    return 0;
}

/**
 * \function igraph_spmatrix_clear_row
 * \brief Clears a row in the matrix (sets all of its elements to zero).
 * \param m The matrix.
 * \param row The index of the row to be cleared.
 * \return Error code.
 *
 * Time complexity: O(n), the number of nonzero elements in the matrix.
 */

int igraph_spmatrix_clear_row(igraph_spmatrix_t *m, long int row) {
    if (row < 0 || row >= m->nrow) {
       IGRAPH_ERROR("The row does not exist.", IGRAPH_EINVAL);
    }
    long int ci, ei, i, j, nremove = 0, nremove_old = 0;
    igraph_vector_t permvec;

    IGRAPH_ASSERT(m != NULL);
    IGRAPH_VECTOR_INIT_FINALLY(&permvec, igraph_vector_size(&m->data));
    for (ci = 0, i = 0, j = 1; ci < m->ncol; ci++) {
        for (ei = (long int) VECTOR(m->cidx)[ci]; ei < VECTOR(m->cidx)[ci + 1]; ei++) {
            if (VECTOR(m->ridx)[ei] == row) {
                /* this element will be deleted, so all elements in cidx from the
                 * column index of this element will have to be decreased by one */
                nremove++;
            } else {
                /* this element will be kept */
                VECTOR(permvec)[i] = j;
                j++;
            }
            i++;
        }
        if (ci > 0) {
            VECTOR(m->cidx)[ci] -= nremove_old;
        }
        nremove_old = nremove;
    }
    VECTOR(m->cidx)[m->ncol] -= nremove;
    igraph_vector_permdelete(&m->ridx, &permvec, nremove);
    igraph_vector_permdelete(&m->data, &permvec, nremove);
    igraph_vector_destroy(&permvec);
    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

/* Unused local functions---temporarily disabled */
#if 0
static int igraph_i_spmatrix_clear_row_fast(igraph_spmatrix_t *m, long int row) {
    long int ei, n;

    IGRAPH_ASSERT(m != NULL);
    n = igraph_vector_size(&m->data);
    for (ei = 0; ei < n; ei++) {
        if (VECTOR(m->ridx)[ei] == row) {
            VECTOR(m->data)[ei] = 0.0;
        }
    }
    return 0;
}

static int igraph_i_spmatrix_cleanup(igraph_spmatrix_t *m) {
    long int ci, ei, i, j, nremove = 0, nremove_old = 0;
    igraph_vector_t permvec;

    IGRAPH_ASSERT(m != NULL);
    IGRAPH_VECTOR_INIT_FINALLY(&permvec, igraph_vector_size(&m->data));
    for (ci = 0, i = 0, j = 1; ci < m->ncol; ci++) {
        for (ei = (long int) VECTOR(m->cidx)[ci]; ei < VECTOR(m->cidx)[ci + 1]; ei++) {
            if (VECTOR(m->data)[ei] == 0.0) {
                /* this element will be deleted, so all elements in cidx from the
                 * column index of this element will have to be decreased by one */
                nremove++;
            } else {
                /* this element will be kept */
                VECTOR(permvec)[i] = j;
                j++;
            }
            i++;
        }
        if (ci > 0) {
            VECTOR(m->cidx)[ci] -= nremove_old;
        }
        nremove_old = nremove;
    }
    VECTOR(m->cidx)[m->ncol] -= nremove;
    igraph_vector_permdelete(&m->ridx, &permvec, nremove);
    igraph_vector_permdelete(&m->data, &permvec, nremove);
    igraph_vector_destroy(&permvec);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}
#endif

/**
 * \function igraph_spmatrix_clear_col
 * \brief Clears a column in the matrix (sets all of its elements to zero).
 * \param m The matrix.
 * \param col The index of the column to be cleared.
 * \return Error code.
 *
 * Time complexity: TODO
 */

int igraph_spmatrix_clear_col(igraph_spmatrix_t *m, long int col) {
    if (col < 0 || col >= m->ncol) {
       IGRAPH_ERROR("The column does not exist.", IGRAPH_EINVAL);
    }
    long int i, n;
    IGRAPH_ASSERT(m != NULL);
    n = (long)VECTOR(m->cidx)[col + 1] - (long)VECTOR(m->cidx)[col];
    if (n == 0) {
        return 0;
    }
    igraph_vector_remove_section(&m->ridx, (long int) VECTOR(m->cidx)[col],
                                 (long int) VECTOR(m->cidx)[col + 1]);
    igraph_vector_remove_section(&m->data, (long int) VECTOR(m->cidx)[col],
                                 (long int) VECTOR(m->cidx)[col + 1]);
    for (i = col + 1; i <= m->ncol; i++) {
        VECTOR(m->cidx)[i] -= n;
    }
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_spmatrix_scale
 * \brief Multiplies each element of the sparse matrix by a constant.
 * \param m The matrix.
 * \param by The constant.
 *
 * Time complexity: O(n), the number of elements in the matrix.
 */

void igraph_spmatrix_scale(igraph_spmatrix_t *m, igraph_real_t by) {
    IGRAPH_ASSERT(m != NULL);
    igraph_vector_scale(&m->data, by);
}

/**
 * \function igraph_spmatrix_colsums
 * \brief Calculates the column sums of the matrix.
 * \param m The matrix.
 * \param res An initialized \c igraph_vector_t, the result will be stored here.
 *   The vector will be resized as needed.
 *
 * Time complexity: O(n), the number of nonzero elements in the matrix.
 */

int igraph_spmatrix_colsums(const igraph_spmatrix_t *m, igraph_vector_t *res) {
    long int i, c;
    IGRAPH_ASSERT(m != NULL);
    IGRAPH_CHECK(igraph_vector_resize(res, m->ncol));
    igraph_vector_null(res);
    for (c = 0; c < m->ncol; c++) {
        for (i = (long int) VECTOR(m->cidx)[c]; i < VECTOR(m->cidx)[c + 1]; i++) {
            VECTOR(*res)[c] += VECTOR(m->data)[i];
        }
    }
    return 0;
}

/**
 * \function igraph_spmatrix_rowsums
 * \brief Calculates the row sums of the matrix.
 * \param m The matrix.
 * \param res An initialized \c igraph_vector_t, the result will be stored here.
 *   The vector will be resized as needed.
 *
 * Time complexity: O(n), the number of nonzero elements in the matrix.
 */

int igraph_spmatrix_rowsums(const igraph_spmatrix_t *m, igraph_vector_t *res) {
    long int i, n;
    IGRAPH_ASSERT(m != NULL);

    IGRAPH_CHECK(igraph_vector_resize(res, m->nrow));
    n = igraph_vector_size(&m->data);
    igraph_vector_null(res);
    for (i = 0; i < n; i++) {
        VECTOR(*res)[(long int)VECTOR(m->ridx)[i]] += VECTOR(m->data)[i];
    }
    return 0;
}

/**
 * \function igraph_spmatrix_max_nonzero
 * \brief Returns the maximum nonzero element of a matrix.
 * If the matrix is empty, zero is returned.
 *
 * \param m the matrix object.
 * \param ridx the row index of the maximum element if not \c NULL.
 * \param cidx the column index of the maximum element if not \c NULL.
 *
 * Time complexity: O(n), the number of nonzero elements in the matrix.
 */
igraph_real_t igraph_spmatrix_max_nonzero(const igraph_spmatrix_t *m,
        igraph_real_t *ridx, igraph_real_t *cidx) {
    igraph_real_t res;
    long int i, n, maxidx;

    IGRAPH_ASSERT(m != NULL);
    n = igraph_vector_size(&m->data);
    if (n == 0) {
        return 0.0;
    }

    maxidx = -1;
    for (i = 0; i < n; i++)
        if (VECTOR(m->data)[i] != 0.0 &&
            (maxidx == -1 || VECTOR(m->data)[i] >= VECTOR(m->data)[maxidx])) {
            maxidx = i;
        }

    if (maxidx == -1) {
        return 0.0;
    }

    res = VECTOR(m->data)[maxidx];
    if (ridx != 0) {
        *ridx = VECTOR(m->ridx)[maxidx];
    }
    if (cidx != 0) {
        igraph_vector_binsearch(&m->cidx, maxidx, &i);
        while (VECTOR(m->cidx)[i + 1] == VECTOR(m->cidx)[i]) {
            i++;
        }
        *cidx = (igraph_real_t)i;
    }
    return res;
}

/**
 * \function igraph_spmatrix_max
 * \brief Returns the maximum element of a matrix.
 * If the matrix is empty, zero is returned.
 *
 * \param m the matrix object.
 * \param ridx the row index of the maximum element if not \c NULL.
 * \param cidx the column index of the maximum element if not \c NULL.
 *
 * Time complexity: O(n), the number of nonzero elements in the matrix.
 */
igraph_real_t igraph_spmatrix_max(const igraph_spmatrix_t *m,
                                  igraph_real_t *ridx, igraph_real_t *cidx) {
    igraph_real_t res;
    long int i, j, k, maxidx;

    IGRAPH_ASSERT(m != NULL);
    i = igraph_vector_size(&m->data);
    if (i == 0) {
        return 0.0;
    }

    maxidx = (long)igraph_vector_which_max(&m->data);
    res = VECTOR(m->data)[maxidx];
    if (res >= 0.0 || i == m->nrow * m->ncol) {
        if (ridx != 0) {
            *ridx = VECTOR(m->ridx)[maxidx];
        }
        if (cidx != 0) {
            igraph_vector_binsearch(&m->cidx, maxidx, &i);
            i--;
            while (i < m->ncol - 1 && VECTOR(m->cidx)[i + 1] == VECTOR(m->cidx)[i]) {
                i++;
            }
            *cidx = (igraph_real_t)i;
        }
        return res;
    }
    /* the maximal nonzero element is negative and there is at least a
     * single zero
     */
    res = 0.0;
    if (cidx != 0 || ridx != 0) {
        for (i = 0; i < m->ncol; i++) {
            if (VECTOR(m->cidx)[i + 1] - VECTOR(m->cidx)[i] < m->nrow) {
                if (cidx != 0) {
                    *cidx = i;
                }
                if (ridx != 0) {
                    for (j = (long int) VECTOR(m->cidx)[i], k = 0;
                         j < VECTOR(m->cidx)[i + 1]; j++, k++) {
                        if (VECTOR(m->ridx)[j] != k) {
                            *ridx = k;
                            break;
                        }
                    }
                }
                break;
            }
        }
    }

    return res;
}


/* Unused function, temporarily disabled */
/*
static int igraph_i_spmatrix_get_col_nonzero_indices(const igraph_spmatrix_t *m,
                                                     igraph_vector_t *res, long int col) {
    long int i, n;
    IGRAPH_ASSERT(m != NULL);
    n = (long int) (VECTOR(m->cidx)[col + 1] - VECTOR(m->cidx)[col]);
    IGRAPH_CHECK(igraph_vector_resize(res, n));
    for (i = (long int) VECTOR(m->cidx)[col], n = 0;
         i < VECTOR(m->cidx)[col + 1]; i++, n++)
        if (VECTOR(m->data)[i] != 0.0) {
            VECTOR(*res)[n] = VECTOR(m->ridx)[i];
        }
    return 0;
}
*/


/**
 * \section igraph_spmatrix_iterating Iterating over the non-zero elements of a sparse matrix
 *
 * <para>The \type igraph_spmatrix_iter_t type represents an iterator that can
 * be used to step over the non-zero elements of a sparse matrix in columnwise
 * order efficiently. In general, you shouldn't modify the elements of the matrix
 * while iterating over it; doing so will probably invalidate the iterator, but
 * there are no checks to prevent you from doing this.</para>
 *
 * <para>To access the row index of the current element of the iterator, use its
 * \c ri field. Similarly, the \c ci field stores the column index of the current
 * element and the \c value field stores the value of the element.</para>
 */

/**
 * \function igraph_spmatrix_iter_create
 * \brief Creates a sparse matrix iterator corresponding to the given matrix.
 *
 * \param  mit  pointer to the matrix iterator being initialized
 * \param  m    pointer to the matrix we will be iterating over
 * \return  Error code. The current implementation is always successful.
 *
 * Time complexity: O(1).
 */
int igraph_spmatrix_iter_create(igraph_spmatrix_iter_t *mit, const igraph_spmatrix_t *m) {
    mit->m = m;
    IGRAPH_CHECK(igraph_spmatrix_iter_reset(mit));
    return 0;
}

/**
 * \function igraph_spmatrix_iter_reset
 * \brief Resets a sparse matrix iterator.
 *
 * </para><para>
 * After resetting, the iterator will point to the first nonzero element (if any).
 *
 * \param  mit  pointer to the matrix iterator being reset
 * \return  Error code. The current implementation is always successful.
 *
 * Time complexity: O(1).
 */
int igraph_spmatrix_iter_reset(igraph_spmatrix_iter_t *mit) {
    IGRAPH_ASSERT(mit->m);

    if (igraph_spmatrix_count_nonzero(mit->m) == 0) {
        mit->pos = mit->ri = mit->ci = -1L;
        mit->value = -1;
        return 0;
    }

    mit->ci = 0;
    mit->pos = -1;

    IGRAPH_CHECK(igraph_spmatrix_iter_next(mit));

    return 0;
}

/**
 * \function igraph_spmatrix_iter_next
 * \brief Moves a sparse matrix iterator to the next nonzero element.
 *
 * </para><para>
 * You should call this function only if \ref igraph_spmatrix_iter_end()
 * returns FALSE (0).
 *
 * \param  mit  pointer to the matrix iterator being moved
 * \return  Error code. The current implementation is always successful.
 *
 * Time complexity: O(1).
 */
int igraph_spmatrix_iter_next(igraph_spmatrix_iter_t *mit) {
    mit->pos++;

    if (igraph_spmatrix_iter_end(mit)) {
        return 0;
    }

    mit->ri = (long int)VECTOR(mit->m->ridx)[mit->pos];
    mit->value = VECTOR(mit->m->data)[mit->pos];

    while (VECTOR(mit->m->cidx)[mit->ci + 1] <= mit->pos) {
        mit->ci++;
    }

    return 0;
}

/**
 * \function igraph_spmatrix_iter_end
 * \brief Checks whether there are more elements in the iterator.
 *
 * </para><para>
 * You should call this function before calling \ref igraph_spmatrix_iter_next()
 * to make sure you have more elements in the iterator.
 *
 * \param  mit  pointer to the matrix iterator being checked
 * \return   TRUE (1) if there are more elements in the iterator,
 *           FALSE (0) otherwise.
 *
 * Time complexity: O(1).
 */
igraph_bool_t igraph_spmatrix_iter_end(igraph_spmatrix_iter_t *mit) {
    return mit->pos >= igraph_spmatrix_count_nonzero(mit->m);
}

/**
 * \function igraph_spmatrix_iter_destroy
 * \brief Frees the memory used by the iterator.
 *
 * </para><para>
 * The current implementation does not allocate any memory upon
 * creation, so this function does nothing. However, since there is
 * no guarantee that future implementations will not allocate any
 * memory in \ref igraph_spmatrix_iter_create(), you are still
 * required to call this function whenever you are done with the
 * iterator.
 *
 * \param  mit  pointer to the matrix iterator being destroyed
 *
 * Time complexity: O(1).
 */
void igraph_spmatrix_iter_destroy(igraph_spmatrix_iter_t *mit) {
    IGRAPH_UNUSED(mit);
    /* Nothing to do at the moment */
}

#ifndef USING_R
/**
 * \function igraph_spmatrix_print
 * \brief Prints a sparse matrix.
 *
 * Prints a sparse matrix to the standard output. Only the non-zero entries
 * are printed.
 *
 * \return Error code.
 *
 * Time complexity: O(n), the number of non-zero elements.
 */
int igraph_spmatrix_print(const igraph_spmatrix_t* matrix) {
    return igraph_spmatrix_fprint(matrix, stdout);
}
#endif

/**
 * \function igraph_spmatrix_fprint
 * \brief Prints a sparse matrix to the given file.
 *
 * Prints a sparse matrix to the given file. Only the non-zero entries
 * are printed.
 *
 * \return Error code.
 *
 * Time complexity: O(n), the number of non-zero elements.
 */
int igraph_spmatrix_fprint(const igraph_spmatrix_t* matrix, FILE *file) {
    igraph_spmatrix_iter_t mit;

    IGRAPH_CHECK(igraph_spmatrix_iter_create(&mit, matrix));
    IGRAPH_FINALLY(igraph_spmatrix_iter_destroy, &mit);
    while (!igraph_spmatrix_iter_end(&mit)) {
        fprintf(file, "[%ld, %ld] = %.4f\n", (long int)mit.ri,
                (long int)mit.ci, mit.value);
        igraph_spmatrix_iter_next(&mit);
    }
    igraph_spmatrix_iter_destroy(&mit);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}


