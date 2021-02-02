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
    long int nrow, ncol;
} TYPE(igraph_matrix);

/*---------------*/
/* Allocation    */
/*---------------*/

IGRAPH_EXPORT int FUNCTION(igraph_matrix, init)(TYPE(igraph_matrix) *m,
                                                long int nrow, long int ncol);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, copy)(TYPE(igraph_matrix) *to,
                                                const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, destroy)(TYPE(igraph_matrix) *m);
IGRAPH_EXPORT long int FUNCTION(igraph_matrix, capacity)(const TYPE(igraph_matrix) *m);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

/* MATRIX */
IGRAPH_EXPORT BASE FUNCTION(igraph_matrix, e)(const TYPE(igraph_matrix) *m,
                                              long int row, long int col);
IGRAPH_EXPORT BASE* FUNCTION(igraph_matrix, e_ptr)(const TYPE(igraph_matrix) *m,
                                                   long int row, long int col);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, set)(TYPE(igraph_matrix)* m, long int row, long int col,
                                                BASE value);

/*------------------------------*/
/* Initializing matrix elements */
/*------------------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_matrix, null)(TYPE(igraph_matrix) *m);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, fill)(TYPE(igraph_matrix) *m, BASE e);

/*-----------------------*/
/* Matrix views          */
/*-----------------------*/

IGRAPH_EXPORT const TYPE(igraph_matrix) *FUNCTION(igraph_matrix, view)(const TYPE(igraph_matrix) *m,
                                                                       const BASE *data,
                                                                       long int nrow,
                                                                       long int ncol);

/*------------------*/
/* Copying matrices */
/*------------------*/

IGRAPH_EXPORT void FUNCTION(igraph_matrix, copy_to)(const TYPE(igraph_matrix) *m, BASE *to);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, update)(TYPE(igraph_matrix) *to,
                                                  const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, rbind)(TYPE(igraph_matrix) *to,
                                                 const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, cbind)(TYPE(igraph_matrix) *to,
                                                 const TYPE(igraph_matrix) *from);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, swap)(TYPE(igraph_matrix) *m1, TYPE(igraph_matrix) *m2);

/*--------------------------*/
/* Copying rows and columns */
/*--------------------------*/

IGRAPH_EXPORT int FUNCTION(igraph_matrix, get_row)(const TYPE(igraph_matrix) *m,
                                                   TYPE(igraph_vector) *res, long int index);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, get_col)(const TYPE(igraph_matrix) *m,
                                                   TYPE(igraph_vector) *res, long int index);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, set_row)(TYPE(igraph_matrix) *m,
                                                   const TYPE(igraph_vector) *v, long int index);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, set_col)(TYPE(igraph_matrix) *m,
                                                   const TYPE(igraph_vector) *v, long int index);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, select_rows)(const TYPE(igraph_matrix) *m,
                                                       TYPE(igraph_matrix) *res,
                                                       const igraph_vector_t *rows);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, select_cols)(const TYPE(igraph_matrix) *m,
                                                       TYPE(igraph_matrix) *res,
                                                       const igraph_vector_t *cols);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, select_rows_cols)(const TYPE(igraph_matrix) *m,
                                                            TYPE(igraph_matrix) *res,
                                                            const igraph_vector_t *rows,
                                                            const igraph_vector_t *cols);

/*-----------------------------*/
/* Exchanging rows and columns */
/*-----------------------------*/

IGRAPH_EXPORT int FUNCTION(igraph_matrix, swap_rows)(TYPE(igraph_matrix) *m,
                                                     long int i, long int j);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, swap_cols)(TYPE(igraph_matrix) *m,
                                                     long int i, long int j);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, swap_rowcol)(TYPE(igraph_matrix) *m,
                                                       long int i, long int j);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, transpose)(TYPE(igraph_matrix) *m);

/*-----------------------------*/
/* Matrix operations           */
/*-----------------------------*/

IGRAPH_EXPORT int FUNCTION(igraph_matrix, add)(TYPE(igraph_matrix) *m1,
                                               const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, sub)(TYPE(igraph_matrix) *m1,
                                               const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, mul_elements)(TYPE(igraph_matrix) *m1,
                                                        const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, div_elements)(TYPE(igraph_matrix) *m1,
                                                        const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, scale)(TYPE(igraph_matrix) *m, BASE by);
IGRAPH_EXPORT void FUNCTION(igraph_matrix, add_constant)(TYPE(igraph_matrix) *m, BASE plus);

/*-----------------------------*/
/* Finding minimum and maximum */
/*-----------------------------*/

IGRAPH_EXPORT igraph_real_t FUNCTION(igraph_matrix, min)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT igraph_real_t FUNCTION(igraph_matrix, max)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, which_min)(const TYPE(igraph_matrix) *m,
                                                     long int *i, long int *j);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, which_max)(const TYPE(igraph_matrix) *m,
                                                     long int *i, long int *j);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, minmax)(const TYPE(igraph_matrix) *m,
                                                  BASE *min, BASE *max);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, which_minmax)(const TYPE(igraph_matrix) *m,
                                                        long int *imin, long int *jmin,
                                                        long int *imax, long int *jmax);

/*------------------------------*/
/* Comparison                   */
/*------------------------------*/

IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, all_e)(const TYPE(igraph_matrix) *lhs,
                                                           const TYPE(igraph_matrix) *rhs);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, all_l)(const TYPE(igraph_matrix) *lhs,
                                                           const TYPE(igraph_matrix) *rhs);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, all_g)(const TYPE(igraph_matrix) *lhs,
                                                           const TYPE(igraph_matrix) *rhs);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, all_le)(const TYPE(igraph_matrix) *lhs,
                                                            const TYPE(igraph_matrix) *rhs);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, all_ge)(const TYPE(igraph_matrix) *lhs,
                                                            const TYPE(igraph_matrix) *rhs);

/*-------------------*/
/* Matrix properties */
/*-------------------*/

IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, isnull)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, empty)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT long int FUNCTION(igraph_matrix, size)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT long int FUNCTION(igraph_matrix, nrow)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT long int FUNCTION(igraph_matrix, ncol)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, is_symmetric)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT BASE FUNCTION(igraph_matrix, sum)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT BASE FUNCTION(igraph_matrix, prod)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, rowsum)(const TYPE(igraph_matrix) *m,
                                                  TYPE(igraph_vector) *res);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, colsum)(const TYPE(igraph_matrix) *m,
                                                  TYPE(igraph_vector) *res);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, is_equal)(const TYPE(igraph_matrix) *m1,
                                                              const TYPE(igraph_matrix) *m2);
IGRAPH_EXPORT igraph_real_t FUNCTION(igraph_matrix, maxdifference)(const TYPE(igraph_matrix) *m1,
                                                                   const TYPE(igraph_matrix) *m2);

/*------------------------*/
/* Searching for elements */
/*------------------------*/

IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, contains)(const TYPE(igraph_matrix) *m,
                                                              BASE e);
IGRAPH_EXPORT igraph_bool_t FUNCTION(igraph_matrix, search)(const TYPE(igraph_matrix) *m,
                                                            long int from, BASE what,
                                                            long int *pos,
                                                            long int *row, long int *col);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

IGRAPH_EXPORT int FUNCTION(igraph_matrix, resize)(TYPE(igraph_matrix) *m,
                                                  long int nrow, long int ncol);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, resize_min)(TYPE(igraph_matrix) *m);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, add_cols)(TYPE(igraph_matrix) *m, long int n);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, add_rows)(TYPE(igraph_matrix) *m, long int n);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, remove_col)(TYPE(igraph_matrix) *m, long int col);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, remove_row)(TYPE(igraph_matrix) *m, long int row);

/*------------------------*/
/* Print as text          */
/*------------------------*/

IGRAPH_EXPORT int FUNCTION(igraph_matrix, print)(const TYPE(igraph_matrix) *m);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, printf)(const TYPE(igraph_matrix) *m,
                                                  const char *format);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, fprint)(const TYPE(igraph_matrix) *m,
                                                  FILE *file);

#ifdef BASE_COMPLEX

IGRAPH_EXPORT int igraph_matrix_complex_real(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *real);
IGRAPH_EXPORT int igraph_matrix_complex_imag(const igraph_matrix_complex_t *v,
                               igraph_matrix_t *imag);
IGRAPH_EXPORT int igraph_matrix_complex_realimag(const igraph_matrix_complex_t *v,
                                   igraph_matrix_t *real,
                                   igraph_matrix_t *imag);
IGRAPH_EXPORT int igraph_matrix_complex_create(igraph_matrix_complex_t *v,
                                 const igraph_matrix_t *real,
                                 const igraph_matrix_t *imag);
IGRAPH_EXPORT int igraph_matrix_complex_create_polar(igraph_matrix_complex_t *v,
                                       const igraph_matrix_t *r,
                                       const igraph_matrix_t *theta);

#endif

IGRAPH_EXPORT int FUNCTION(igraph_matrix, permdelete_rows)(TYPE(igraph_matrix) *m,
                                                           long int *index, long int nremove);
IGRAPH_EXPORT int FUNCTION(igraph_matrix, delete_rows_neg)(TYPE(igraph_matrix) *m,
                                                           const igraph_vector_t *neg,
                                                           long int nremove);

