/*
   igraph library.
   Copyright (C) 2007-2025  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

typedef struct TYPE(igraph_matrix) {
    TYPE(igraph_vector) data;
    igraph_int_t nrow, ncol;
} TYPE(igraph_matrix);

/*---------------*/
/* Allocation    */
/*---------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, init)(
    TYPE(igraph_matrix) *m, igraph_int_t nrow, igraph_int_t ncol);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, init_array)(
    TYPE(igraph_matrix)* m, const BASE* data, igraph_int_t nrow, igraph_int_t ncol, igraph_matrix_storage_t storage);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, init_copy)(
    TYPE(igraph_matrix) *to, const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, destroy)(TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_matrix, capacity)(const TYPE(igraph_matrix) *m);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

/* MATRIX */
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_matrix, get)(
    const TYPE(igraph_matrix) *m, igraph_int_t row, igraph_int_t col);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE* FUNCTION(igraph_matrix, get_ptr)(
    const TYPE(igraph_matrix) *m, igraph_int_t row, igraph_int_t col);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, set)(
    TYPE(igraph_matrix)* m, igraph_int_t row, igraph_int_t col, BASE value);

/*------------------------------*/
/* Initializing matrix elements */
/*------------------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_matrix, null)(TYPE(igraph_matrix) *m);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, fill)(TYPE(igraph_matrix) *m, BASE e);

/*-----------------------*/
/* Matrix views          */
/*-----------------------*/

IGRAPH_EXPORT TYPE(igraph_matrix) FUNCTION(igraph_matrix, view)(
    const BASE *data,
    igraph_int_t nrow, igraph_int_t ncol);
IGRAPH_EXPORT TYPE(igraph_matrix) FUNCTION(igraph_matrix, view_from_vector)(
    const TYPE(igraph_vector) *v,
    igraph_int_t ncol
);

/*------------------*/
/* Copying matrices */
/*------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_matrix, copy_to)(const TYPE(igraph_matrix) *m, BASE *to, igraph_matrix_storage_t storage);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, update)(TYPE(igraph_matrix) *to,
                                                  const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, rbind)(TYPE(igraph_matrix) *to,
                                                 const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, cbind)(TYPE(igraph_matrix) *to,
                                                 const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, swap)(TYPE(igraph_matrix) *m1, TYPE(igraph_matrix) *m2);

/*--------------------------*/
/* Copying rows and columns */
/*--------------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, get_row)(
    const TYPE(igraph_matrix) *m, TYPE(igraph_vector) *res, igraph_int_t index);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, get_col)(
    const TYPE(igraph_matrix) *m, TYPE(igraph_vector) *res, igraph_int_t index);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, set_row)(
    TYPE(igraph_matrix) *m, const TYPE(igraph_vector) *v, igraph_int_t index);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, set_col)(
    TYPE(igraph_matrix) *m, const TYPE(igraph_vector) *v, igraph_int_t index);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, select_rows)(
    const TYPE(igraph_matrix) *m, TYPE(igraph_matrix) *res, const igraph_vector_int_t *rows);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, select_cols)(
    const TYPE(igraph_matrix) *m, TYPE(igraph_matrix) *res, const igraph_vector_int_t *cols);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, select_rows_cols)(
    const TYPE(igraph_matrix) *m, TYPE(igraph_matrix) *res,
    const igraph_vector_int_t *rows, const igraph_vector_int_t *cols);

/*-----------------------------*/
/* Exchanging rows and columns */
/*-----------------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, swap_rows)(
    TYPE(igraph_matrix) *m, igraph_int_t i, igraph_int_t j);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, swap_cols)(
    TYPE(igraph_matrix) *m, igraph_int_t i, igraph_int_t j);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, swap_rowcol)(
    TYPE(igraph_matrix) *m, igraph_int_t i, igraph_int_t j);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, transpose)(TYPE(igraph_matrix) *m);

/*-----------------------------*/
/* Matrix operations           */
/*-----------------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, add)(TYPE(igraph_matrix) *m1,
                                               const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, sub)(TYPE(igraph_matrix) *m1,
                                               const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, mul_elements)(TYPE(igraph_matrix) *m1,
                                                        const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, div_elements)(TYPE(igraph_matrix) *m1,
                                                        const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, scale)(TYPE(igraph_matrix) *m, BASE by);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, add_constant)(TYPE(igraph_matrix) *m, BASE plus);

/*-----------------------------*/
/* Finding minimum and maximum */
/*-----------------------------*/

#ifndef NOTORDERED
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t FUNCTION(igraph_matrix, min)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t FUNCTION(igraph_matrix, max)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, which_min)(
    const TYPE(igraph_matrix) *m, igraph_int_t *i, igraph_int_t *j);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, which_max)(
    const TYPE(igraph_matrix) *m, igraph_int_t *i, igraph_int_t *j);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, minmax)(
    const TYPE(igraph_matrix) *m, BASE *min, BASE *max);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, which_minmax)(
    const TYPE(igraph_matrix) *m, igraph_int_t *imin, igraph_int_t *jmin,
    igraph_int_t *imax, igraph_int_t *jmax);
#endif

/*------------------------------*/
/* Comparison                   */
/*------------------------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, all_e)(const TYPE(igraph_matrix) *lhs,
                                                           const TYPE(igraph_matrix) *rhs);
#ifndef NOTORDERED
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, all_l)(const TYPE(igraph_matrix) *lhs,
                                                           const TYPE(igraph_matrix) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, all_g)(const TYPE(igraph_matrix) *lhs,
                                                           const TYPE(igraph_matrix) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, all_le)(const TYPE(igraph_matrix) *lhs,
                                                            const TYPE(igraph_matrix) *rhs);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, all_ge)(const TYPE(igraph_matrix) *lhs,
                                                            const TYPE(igraph_matrix) *rhs);
#endif

/*-------------------*/
/* Matrix properties */
/*-------------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, isnull)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, empty)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_matrix, size)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_matrix, nrow)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_int_t FUNCTION(igraph_matrix, ncol)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, is_symmetric)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_matrix, sum)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE BASE FUNCTION(igraph_matrix, prod)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, rowsum)(const TYPE(igraph_matrix) *m,
                                                  TYPE(igraph_vector) *res);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, colsum)(const TYPE(igraph_matrix) *m,
                                                  TYPE(igraph_vector) *res);
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, is_equal)(const TYPE(igraph_matrix) *m1,
                                                              const TYPE(igraph_matrix) *m2);
#ifndef NOTORDERED
IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_real_t FUNCTION(igraph_matrix, maxdifference)(const TYPE(igraph_matrix) *m1,
                                                                   const TYPE(igraph_matrix) *m2);
#endif

/*------------------------*/
/* Searching for elements */
/*------------------------*/

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t FUNCTION(igraph_matrix, contains)(
    const TYPE(igraph_matrix) *m, BASE e);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, search)(
    const TYPE(igraph_matrix) *m, igraph_int_t from, BASE what,
    igraph_int_t *pos, igraph_int_t *row, igraph_int_t *col);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, resize)(
    TYPE(igraph_matrix) *m, igraph_int_t nrow, igraph_int_t ncol);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, resize_min)(
    TYPE(igraph_matrix) *m);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, add_cols)(
    TYPE(igraph_matrix) *m, igraph_int_t n);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, add_rows)(
    TYPE(igraph_matrix) *m, igraph_int_t n);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, remove_col)(
    TYPE(igraph_matrix) *m, igraph_int_t col);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, remove_row)(
    TYPE(igraph_matrix) *m, igraph_int_t row);

/*------------------------*/
/* Print as text          */
/*------------------------*/

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, print)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, fprint)(const TYPE(igraph_matrix) *m, FILE *file);

#ifdef OUT_FORMAT
IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, printf)(const TYPE(igraph_matrix) *m, const char *format);
#endif /* OUT_FORMAT */

/*-----------------------------------------*/
/* Operations specific to complex matrices */
/*-----------------------------------------*/

#ifdef BASE_COMPLEX

IGRAPH_EXPORT igraph_error_t igraph_matrix_complex_real(const igraph_matrix_complex_t *v,
                                                        igraph_matrix_t *real);
IGRAPH_EXPORT igraph_error_t igraph_matrix_complex_imag(const igraph_matrix_complex_t *v,
                                                        igraph_matrix_t *imag);
IGRAPH_EXPORT igraph_error_t igraph_matrix_complex_realimag(const igraph_matrix_complex_t *v,
                                                            igraph_matrix_t *real,
                                                            igraph_matrix_t *imag);
IGRAPH_EXPORT igraph_error_t igraph_matrix_complex_create(igraph_matrix_complex_t *v,
                                                          const igraph_matrix_t *real,
                                                          const igraph_matrix_t *imag);
IGRAPH_EXPORT igraph_error_t igraph_matrix_complex_create_polar(igraph_matrix_complex_t *v,
                                                                const igraph_matrix_t *r,
                                                                const igraph_matrix_t *theta);

IGRAPH_EXPORT IGRAPH_FUNCATTR_PURE igraph_bool_t igraph_matrix_complex_all_almost_e(igraph_matrix_complex_t *lhs,
                                                                igraph_matrix_complex_t *rhs,
                                                                igraph_real_t eps);

#endif /* BASE_COMPLEX */

IGRAPH_EXPORT igraph_error_t FUNCTION(igraph_matrix, permdelete_rows)(
    TYPE(igraph_matrix) *m, igraph_int_t *index, igraph_int_t nremove);
