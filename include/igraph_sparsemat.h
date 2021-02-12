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

#ifndef IGRAPH_SPARSEMAT_H
#define IGRAPH_SPARSEMAT_H

#include "igraph_decls.h"
#include "igraph_types.h"
#include "igraph_vector.h"
#include "igraph_datatype.h"
#include "igraph_arpack.h"

#include <stdio.h>

__BEGIN_DECLS

struct cs_di_sparse;
struct cs_di_symbolic;
struct cs_di_numeric;

typedef struct {
    struct cs_di_sparse *cs;
} igraph_sparsemat_t;

typedef struct {
    struct cs_di_symbolic *symbolic;
} igraph_sparsemat_symbolic_t;

typedef struct {
    struct cs_di_numeric *numeric;
} igraph_sparsemat_numeric_t;

typedef enum { IGRAPH_SPARSEMAT_TRIPLET,
               IGRAPH_SPARSEMAT_CC
             } igraph_sparsemat_type_t;

typedef struct {
    igraph_sparsemat_t *mat;
    int pos;
    int col;
} igraph_sparsemat_iterator_t;

IGRAPH_EXPORT int igraph_sparsemat_init(igraph_sparsemat_t *A, int rows, int cols, int nzmax);
IGRAPH_EXPORT int igraph_sparsemat_copy(igraph_sparsemat_t *to,
                                        const igraph_sparsemat_t *from);
IGRAPH_EXPORT void igraph_sparsemat_destroy(igraph_sparsemat_t *A);
IGRAPH_EXPORT int igraph_sparsemat_realloc(igraph_sparsemat_t *A, int nzmax);

IGRAPH_EXPORT long int igraph_sparsemat_nrow(const igraph_sparsemat_t *A);
IGRAPH_EXPORT long int igraph_sparsemat_ncol(const igraph_sparsemat_t *B);
IGRAPH_EXPORT igraph_sparsemat_type_t igraph_sparsemat_type(const igraph_sparsemat_t *A);
IGRAPH_EXPORT igraph_bool_t igraph_sparsemat_is_triplet(const igraph_sparsemat_t *A);
IGRAPH_EXPORT igraph_bool_t igraph_sparsemat_is_cc(const igraph_sparsemat_t *A);

IGRAPH_EXPORT int igraph_sparsemat_permute(const igraph_sparsemat_t *A,
                                           const igraph_vector_int_t *p,
                                           const igraph_vector_int_t *q,
                                           igraph_sparsemat_t *res);

IGRAPH_EXPORT int igraph_sparsemat_index(const igraph_sparsemat_t *A,
                                         const igraph_vector_int_t *p,
                                         const igraph_vector_int_t *q,
                                         igraph_sparsemat_t *res,
                                         igraph_real_t *constres);

IGRAPH_EXPORT int igraph_sparsemat_entry(igraph_sparsemat_t *A, int row, int col,
                                         igraph_real_t elem);
IGRAPH_EXPORT int igraph_sparsemat_compress(const igraph_sparsemat_t *A,
                                            igraph_sparsemat_t *res);
IGRAPH_EXPORT int igraph_sparsemat_transpose(const igraph_sparsemat_t *A,
                                             igraph_sparsemat_t *res, int values);
IGRAPH_EXPORT igraph_bool_t igraph_sparsemat_is_symmetric(const igraph_sparsemat_t *A);
IGRAPH_EXPORT int igraph_sparsemat_dupl(igraph_sparsemat_t *A);
IGRAPH_EXPORT int igraph_sparsemat_fkeep(igraph_sparsemat_t *A,
                                         int (*fkeep)(int, int, igraph_real_t, void*),
                                         void *other);
IGRAPH_EXPORT int igraph_sparsemat_dropzeros(igraph_sparsemat_t *A);
IGRAPH_EXPORT int igraph_sparsemat_droptol(igraph_sparsemat_t *A, igraph_real_t tol);
IGRAPH_EXPORT int igraph_sparsemat_multiply(const igraph_sparsemat_t *A,
                                            const igraph_sparsemat_t *B,
                                            igraph_sparsemat_t *res);
IGRAPH_EXPORT int igraph_sparsemat_add(const igraph_sparsemat_t *A,
                                       const igraph_sparsemat_t *B,
                                       igraph_real_t alpha,
                                       igraph_real_t beta,
                                       igraph_sparsemat_t *res);
IGRAPH_EXPORT int igraph_sparsemat_gaxpy(const igraph_sparsemat_t *A,
                                         const igraph_vector_t *x,
                                         igraph_vector_t *res);

IGRAPH_EXPORT int igraph_sparsemat_lsolve(const igraph_sparsemat_t *A,
                                          const igraph_vector_t *b,
                                          igraph_vector_t *res);
IGRAPH_EXPORT int igraph_sparsemat_ltsolve(const igraph_sparsemat_t *A,
                                           const igraph_vector_t *b,
                                           igraph_vector_t *res);
IGRAPH_EXPORT int igraph_sparsemat_usolve(const igraph_sparsemat_t *A,
                                          const igraph_vector_t *b,
                                          igraph_vector_t *res);
IGRAPH_EXPORT int igraph_sparsemat_utsolve(const igraph_sparsemat_t *A,
                                           const igraph_vector_t *b,
                                           igraph_vector_t *res);

IGRAPH_EXPORT int igraph_sparsemat_cholsol(const igraph_sparsemat_t *A,
                                           const igraph_vector_t *b,
                                           igraph_vector_t *res,
                                           int order);

IGRAPH_EXPORT int igraph_sparsemat_lusol(const igraph_sparsemat_t *A,
                                         const igraph_vector_t *b,
                                         igraph_vector_t *res,
                                         int order,
                                         igraph_real_t tol);

IGRAPH_EXPORT int igraph_sparsemat_print(const igraph_sparsemat_t *A,
                                         FILE *outstream);

IGRAPH_EXPORT int igraph_sparsemat_eye(igraph_sparsemat_t *A, int n, int nzmax,
                                       igraph_real_t value,
                                       igraph_bool_t compress);

IGRAPH_EXPORT int igraph_sparsemat_diag(igraph_sparsemat_t *A, int nzmax,
                                        const igraph_vector_t *values,
                                        igraph_bool_t compress);

IGRAPH_EXPORT int igraph_sparsemat(igraph_t *graph, const igraph_sparsemat_t *A,
                                   igraph_bool_t directed);

IGRAPH_EXPORT int igraph_weighted_sparsemat(igraph_t *graph, const igraph_sparsemat_t *A,
                                            igraph_bool_t directed, const char *attr,
                                            igraph_bool_t loops);

IGRAPH_EXPORT int igraph_get_sparsemat(const igraph_t *graph, igraph_sparsemat_t *res);

IGRAPH_EXPORT int igraph_matrix_as_sparsemat(igraph_sparsemat_t *res,
                                             const igraph_matrix_t *mat,
                                             igraph_real_t tol);

IGRAPH_EXPORT int igraph_sparsemat_as_matrix(igraph_matrix_t *res,
                                             const igraph_sparsemat_t *spmat);

typedef enum { IGRAPH_SPARSEMAT_SOLVE_LU,
               IGRAPH_SPARSEMAT_SOLVE_QR
             } igraph_sparsemat_solve_t;

IGRAPH_EXPORT int igraph_sparsemat_arpack_rssolve(const igraph_sparsemat_t *A,
                                                  igraph_arpack_options_t *options,
                                                  igraph_arpack_storage_t *storage,
                                                  igraph_vector_t *values,
                                                  igraph_matrix_t *vectors,
                                                  igraph_sparsemat_solve_t solvemethod);

IGRAPH_EXPORT int igraph_sparsemat_arpack_rnsolve(const igraph_sparsemat_t *A,
                                                  igraph_arpack_options_t *options,
                                                  igraph_arpack_storage_t *storage,
                                                  igraph_matrix_t *values,
                                                  igraph_matrix_t *vectors);

IGRAPH_EXPORT int igraph_sparsemat_lu(const igraph_sparsemat_t *A,
                                      const igraph_sparsemat_symbolic_t *dis,
                                      igraph_sparsemat_numeric_t *din, double tol);

IGRAPH_EXPORT int igraph_sparsemat_qr(const igraph_sparsemat_t *A,
                                      const igraph_sparsemat_symbolic_t *dis,
                                      igraph_sparsemat_numeric_t *din);

IGRAPH_EXPORT int igraph_sparsemat_luresol(const igraph_sparsemat_symbolic_t *dis,
                                           const igraph_sparsemat_numeric_t *din,
                                           const igraph_vector_t *b,
                                           igraph_vector_t *res);

IGRAPH_EXPORT int igraph_sparsemat_qrresol(const igraph_sparsemat_symbolic_t *dis,
                                           const igraph_sparsemat_numeric_t *din,
                                           const igraph_vector_t *b,
                                           igraph_vector_t *res);

IGRAPH_EXPORT int igraph_sparsemat_symbqr(long int order, const igraph_sparsemat_t *A,
                                          igraph_sparsemat_symbolic_t *dis);

IGRAPH_EXPORT int igraph_sparsemat_symblu(long int order, const igraph_sparsemat_t *A,
                                          igraph_sparsemat_symbolic_t *dis);


IGRAPH_EXPORT void igraph_sparsemat_symbolic_destroy(igraph_sparsemat_symbolic_t *dis);
IGRAPH_EXPORT void igraph_sparsemat_numeric_destroy(igraph_sparsemat_numeric_t *din);

IGRAPH_EXPORT igraph_real_t igraph_sparsemat_max(igraph_sparsemat_t *A);
IGRAPH_EXPORT igraph_real_t igraph_sparsemat_min(igraph_sparsemat_t *A);
IGRAPH_EXPORT int igraph_sparsemat_minmax(igraph_sparsemat_t *A,
                                          igraph_real_t *min, igraph_real_t *max);

IGRAPH_EXPORT long int igraph_sparsemat_count_nonzero(igraph_sparsemat_t *A);
IGRAPH_EXPORT long int igraph_sparsemat_count_nonzerotol(igraph_sparsemat_t *A,
                                                         igraph_real_t tol);
IGRAPH_EXPORT int igraph_sparsemat_rowsums(const igraph_sparsemat_t *A,
                                           igraph_vector_t *res);
IGRAPH_EXPORT int igraph_sparsemat_colsums(const igraph_sparsemat_t *A,
                                           igraph_vector_t *res);

IGRAPH_EXPORT int igraph_sparsemat_rowmins(igraph_sparsemat_t *A,
                                           igraph_vector_t *res);
IGRAPH_EXPORT int igraph_sparsemat_colmins(igraph_sparsemat_t *A,
                                           igraph_vector_t *res);

IGRAPH_EXPORT int igraph_sparsemat_rowmaxs(igraph_sparsemat_t *A,
                                           igraph_vector_t *res);
IGRAPH_EXPORT int igraph_sparsemat_colmaxs(igraph_sparsemat_t *A,
                                           igraph_vector_t *res);

IGRAPH_EXPORT int igraph_sparsemat_which_min_rows(igraph_sparsemat_t *A,
                                                  igraph_vector_t *res,
                                                  igraph_vector_int_t *pos);
IGRAPH_EXPORT int igraph_sparsemat_which_min_cols(igraph_sparsemat_t *A,
                                                  igraph_vector_t *res,
                                                  igraph_vector_int_t *pos);

IGRAPH_EXPORT int igraph_sparsemat_scale(igraph_sparsemat_t *A, igraph_real_t by);


IGRAPH_EXPORT int igraph_sparsemat_add_rows(igraph_sparsemat_t *A, long int n);
IGRAPH_EXPORT int igraph_sparsemat_add_cols(igraph_sparsemat_t *A, long int n);
IGRAPH_EXPORT int igraph_sparsemat_resize(igraph_sparsemat_t *A, long int nrow,
                                          long int ncol, int nzmax);
IGRAPH_EXPORT int igraph_sparsemat_nonzero_storage(const igraph_sparsemat_t *A);
IGRAPH_EXPORT int igraph_sparsemat_getelements(const igraph_sparsemat_t *A,
                                               igraph_vector_int_t *i,
                                               igraph_vector_int_t *j,
                                               igraph_vector_t *x);
IGRAPH_EXPORT int igraph_sparsemat_getelements_sorted(const igraph_sparsemat_t *A,
                                                      igraph_vector_int_t *i,
                                                      igraph_vector_int_t *j,
                                                      igraph_vector_t *x);
IGRAPH_EXPORT int igraph_sparsemat_scale_rows(igraph_sparsemat_t *A,
                                              const igraph_vector_t *fact);
IGRAPH_EXPORT int igraph_sparsemat_scale_cols(igraph_sparsemat_t *A,
                                              const igraph_vector_t *fact);
IGRAPH_EXPORT int igraph_sparsemat_multiply_by_dense(const igraph_sparsemat_t *A,
                                                     const igraph_matrix_t *B,
                                                     igraph_matrix_t *res);
IGRAPH_EXPORT int igraph_sparsemat_dense_multiply(const igraph_matrix_t *A,
                                                  const igraph_sparsemat_t *B,
                                                  igraph_matrix_t *res);

IGRAPH_EXPORT int igraph_sparsemat_view(igraph_sparsemat_t *A, int nzmax, int m, int n,
                                        int *p, int *i, double *x, int nz);
IGRAPH_EXPORT IGRAPH_DEPRECATED int igraph_i_sparsemat_view(igraph_sparsemat_t *A, int nzmax, int m, int n,
                                                            int *p, int *i, double *x, int nz);

IGRAPH_EXPORT int igraph_sparsemat_sort(const igraph_sparsemat_t *A,
                                        igraph_sparsemat_t *sorted);

IGRAPH_EXPORT int igraph_sparsemat_nzmax(const igraph_sparsemat_t *A);

IGRAPH_EXPORT int igraph_sparsemat_neg(igraph_sparsemat_t *A);

IGRAPH_EXPORT int igraph_sparsemat_iterator_init(igraph_sparsemat_iterator_t *it,
                                                 igraph_sparsemat_t *sparsemat);
IGRAPH_EXPORT int igraph_sparsemat_iterator_reset(igraph_sparsemat_iterator_t *it);
IGRAPH_EXPORT igraph_bool_t igraph_sparsemat_iterator_end(const igraph_sparsemat_iterator_t *it);
IGRAPH_EXPORT int igraph_sparsemat_iterator_row(const igraph_sparsemat_iterator_t *it);
IGRAPH_EXPORT int igraph_sparsemat_iterator_col(const igraph_sparsemat_iterator_t *it);
IGRAPH_EXPORT int igraph_sparsemat_iterator_idx(const igraph_sparsemat_iterator_t *it);
IGRAPH_EXPORT igraph_real_t igraph_sparsemat_iterator_get(const igraph_sparsemat_iterator_t *it);
IGRAPH_EXPORT int igraph_sparsemat_iterator_next(igraph_sparsemat_iterator_t *it);

__END_DECLS

#endif
