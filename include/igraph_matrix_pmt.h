/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2007-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

typedef struct TYPE(igraph_matrix) {
    TYPE(igraph_vector) data;
    igraph_long_t nrow, ncol;
} TYPE(igraph_matrix);

/*---------------*/
/* Allocation    */
/*---------------*/

DECLDIR igraph_error_t FUNCTION(igraph_matrix, init)(TYPE(igraph_matrix) *m,
        igraph_long_t nrow, igraph_long_t ncol);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, copy)(TYPE(igraph_matrix) *to,
        const TYPE(igraph_matrix) *from);
DECLDIR void FUNCTION(igraph_matrix, destroy)(TYPE(igraph_matrix) *m);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, capacity)(const TYPE(igraph_matrix) *m);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

/* MATRIX */
DECLDIR BASE FUNCTION(igraph_matrix, e)(const TYPE(igraph_matrix) *m,
                                        igraph_long_t row, igraph_long_t col);
DECLDIR BASE* FUNCTION(igraph_matrix, e_ptr)(const TYPE(igraph_matrix) *m,
                                     igraph_long_t row, igraph_long_t col);
DECLDIR void FUNCTION(igraph_matrix, set)(TYPE(igraph_matrix)* m, igraph_long_t row, igraph_long_t col,
        BASE value);

/*------------------------------*/
/* Initializing matrix elements */
/*------------------------------*/

DECLDIR void FUNCTION(igraph_matrix, null)(TYPE(igraph_matrix) *m);
DECLDIR void FUNCTION(igraph_matrix, fill)(TYPE(igraph_matrix) *m, BASE e);

/*-----------------------*/
/* Matrix views          */
/*-----------------------*/

DECLDIR const TYPE(igraph_matrix) *FUNCTION(igraph_matrix, view)(const TYPE(igraph_matrix) *m,
        const BASE *data,
        igraph_long_t nrow,
        igraph_long_t ncol);

/*------------------*/
/* Copying matrices */
/*------------------*/

DECLDIR void FUNCTION(igraph_matrix, copy_to)(const TYPE(igraph_matrix) *m, BASE *to);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, update)(TYPE(igraph_matrix) *to,
        const TYPE(igraph_matrix) *from);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, rbind)(TYPE(igraph_matrix) *to,
        const TYPE(igraph_matrix) *from);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, cbind)(TYPE(igraph_matrix) *to,
        const TYPE(igraph_matrix) *from);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, swap)(TYPE(igraph_matrix) *m1, TYPE(igraph_matrix) *m2);

/*--------------------------*/
/* Copying rows and columns */
/*--------------------------*/

DECLDIR igraph_error_t FUNCTION(igraph_matrix, get_row)(const TYPE(igraph_matrix) *m,
        TYPE(igraph_vector) *res, igraph_long_t index);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, get_col)(const TYPE(igraph_matrix) *m,
        TYPE(igraph_vector) *res, igraph_long_t index);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, set_row)(TYPE(igraph_matrix) *m,
        const TYPE(igraph_vector) *v, igraph_long_t index);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, set_col)(TYPE(igraph_matrix) *m,
        const TYPE(igraph_vector) *v, igraph_long_t index);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, select_rows)(const TYPE(igraph_matrix) *m,
        TYPE(igraph_matrix) *res,
        const igraph_vector_t *rows);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, select_cols)(const TYPE(igraph_matrix) *m,
        TYPE(igraph_matrix) *res,
        const igraph_vector_t *cols);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, select_rows_cols)(const TYPE(igraph_matrix) *m,
        TYPE(igraph_matrix) *res,
        const igraph_vector_t *rows,
        const igraph_vector_t *cols);

/*-----------------------------*/
/* Exchanging rows and columns */
/*-----------------------------*/

DECLDIR igraph_error_t FUNCTION(igraph_matrix, swap_rows)(TYPE(igraph_matrix) *m,
        igraph_long_t i, igraph_long_t j);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, swap_cols)(TYPE(igraph_matrix) *m,
        igraph_long_t i, igraph_long_t j);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, swap_rowcol)(TYPE(igraph_matrix) *m,
        igraph_long_t i, igraph_long_t j);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, transpose)(TYPE(igraph_matrix) *m);

/*-----------------------------*/
/* Matrix operations           */
/*-----------------------------*/

DECLDIR igraph_error_t FUNCTION(igraph_matrix, add)(TYPE(igraph_matrix) *m1,
        const TYPE(igraph_matrix) *m2);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, sub)(TYPE(igraph_matrix) *m1,
        const TYPE(igraph_matrix) *m2);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, mul_elements)(TYPE(igraph_matrix) *m1,
        const TYPE(igraph_matrix) *m2);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, div_elements)(TYPE(igraph_matrix) *m1,
        const TYPE(igraph_matrix) *m2);
DECLDIR void FUNCTION(igraph_matrix, scale)(TYPE(igraph_matrix) *m, BASE by);
DECLDIR void FUNCTION(igraph_matrix, add_constant)(TYPE(igraph_matrix) *m, BASE plus);

/*-----------------------------*/
/* Finding minimum and maximum */
/*-----------------------------*/

DECLDIR igraph_real_t FUNCTION(igraph_matrix, min)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_real_t FUNCTION(igraph_matrix, max)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, which_min)(const TYPE(igraph_matrix) *m,
        igraph_long_t *i, igraph_long_t *j);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, which_max)(const TYPE(igraph_matrix) *m,
        igraph_long_t *i, igraph_long_t *j);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, minmax)(const TYPE(igraph_matrix) *m,
        BASE *min, BASE *max);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, which_minmax)(const TYPE(igraph_matrix) *m,
        igraph_long_t *imin, igraph_long_t *jmin,
        igraph_long_t *imax, igraph_long_t *jmax);

/*------------------------------*/
/* Comparison                   */
/*------------------------------*/

DECLDIR igraph_bool_t FUNCTION(igraph_matrix, all_e)(const TYPE(igraph_matrix) *lhs,
        const TYPE(igraph_matrix) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, all_l)(const TYPE(igraph_matrix) *lhs,
        const TYPE(igraph_matrix) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, all_g)(const TYPE(igraph_matrix) *lhs,
        const TYPE(igraph_matrix) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, all_le)(const TYPE(igraph_matrix) *lhs,
        const TYPE(igraph_matrix) *rhs);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, all_ge)(const TYPE(igraph_matrix) *lhs,
        const TYPE(igraph_matrix) *rhs);

/*-------------------*/
/* Matrix properties */
/*-------------------*/

DECLDIR igraph_bool_t FUNCTION(igraph_matrix, isnull)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, empty)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_long_t FUNCTION(igraph_matrix, size)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, nrow)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, ncol)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, is_symmetric)(const TYPE(igraph_matrix) *m);
DECLDIR BASE FUNCTION(igraph_matrix, sum)(const TYPE(igraph_matrix) *m);
DECLDIR BASE FUNCTION(igraph_matrix, prod)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, rowsum)(const TYPE(igraph_matrix) *m,
        TYPE(igraph_vector) *res);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, colsum)(const TYPE(igraph_matrix) *m,
        TYPE(igraph_vector) *res);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, is_equal)(const TYPE(igraph_matrix) *m1,
        const TYPE(igraph_matrix) *m2);
DECLDIR igraph_real_t FUNCTION(igraph_matrix, maxdifference)(const TYPE(igraph_matrix) *m1,
        const TYPE(igraph_matrix) *m2);

/*------------------------*/
/* Searching for elements */
/*------------------------*/

DECLDIR igraph_bool_t FUNCTION(igraph_matrix, contains)(const TYPE(igraph_matrix) *m,
        BASE e);
DECLDIR igraph_bool_t FUNCTION(igraph_matrix, search)(const TYPE(igraph_matrix) *m,
        igraph_long_t from, BASE what,
        igraph_long_t *pos,
        igraph_long_t *row, igraph_long_t *col);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

DECLDIR igraph_error_t FUNCTION(igraph_matrix, resize)(TYPE(igraph_matrix) *m,
        igraph_long_t nrow, igraph_long_t ncol);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, resize_min)(TYPE(igraph_matrix) *m);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, add_cols)(TYPE(igraph_matrix) *m, igraph_long_t n);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, add_rows)(TYPE(igraph_matrix) *m, igraph_long_t n);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, remove_col)(TYPE(igraph_matrix) *m, igraph_long_t col);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, remove_row)(TYPE(igraph_matrix) *m, igraph_long_t row);

/*------------------------*/
/* Print as text          */
/*------------------------*/

DECLDIR igraph_error_t FUNCTION(igraph_matrix, print)(const TYPE(igraph_matrix) *m);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, printf)(const TYPE(igraph_matrix) *m,
                                    const char *format);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, fprint)(const TYPE(igraph_matrix) *m,
                                    FILE *file);

#ifdef BASE_COMPLEX

DECLDIR igraph_error_t igraph_matrix_complex_real(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *real);
DECLDIR igraph_error_t igraph_matrix_complex_imag(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *imag);
DECLDIR igraph_error_t igraph_matrix_complex_realimag(const igraph_matrix_complex_t *v,
                                   igraph_matrix_t *real,
                                   igraph_matrix_t *imag);
DECLDIR igraph_error_t igraph_matrix_complex_create(igraph_matrix_complex_t *v,
                                 const igraph_matrix_t *real,
                                 const igraph_matrix_t *imag);
DECLDIR igraph_error_t igraph_matrix_complex_create_polar(igraph_matrix_complex_t *v,
                                       const igraph_matrix_t *r,
                                       const igraph_matrix_t *theta);

#endif

DECLDIR igraph_error_t FUNCTION(igraph_matrix, permdelete_rows)(TYPE(igraph_matrix) *m,
        igraph_long_t *index, igraph_long_t nremove);
DECLDIR igraph_error_t FUNCTION(igraph_matrix, delete_rows_neg)(TYPE(igraph_matrix) *m,
        const igraph_vector_t *neg,
        igraph_long_t nremove);

