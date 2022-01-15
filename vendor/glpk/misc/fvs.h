/* fvs.h (sparse vector in FVS format) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2016 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#ifndef FVS_H
#define FVS_H

typedef struct FVS FVS;

struct FVS
{     /* sparse vector in FVS (Full Vector Storage) format */
      int n;
      /* vector dimension (total number of elements) */
      int nnz;
      /* number of non-zero elements, 0 <= nnz <= n */
      int *ind; /* int ind[1+n]; */
      /* ind[0] is not used;
       * ind[k] = j, 1 <= k <= nnz, means that vec[j] != 0
       * non-zero indices in the array ind are stored in arbitrary
       * order; if vec[j] = 0, its index j SHOULD NOT be presented in
       * the array ind */
      double *vec; /* double vec[1+n]; */
      /* vec[0] is not used;
       * vec[j], 1 <= j <= n, is a numeric value of j-th element */
};

#define fvs_alloc_vec _glp_fvs_alloc_vec
void fvs_alloc_vec(FVS *x, int n);
/* allocate sparse vector */

#define fvs_check_vec _glp_fvs_check_vec
void fvs_check_vec(const FVS *x);
/* check sparse vector */

#define fvs_gather_vec _glp_fvs_gather_vec
void fvs_gather_vec(FVS *x, double eps);
/* gather sparse vector */

#define fvs_clear_vec _glp_fvs_clear_vec
void fvs_clear_vec(FVS *x);
/* clear sparse vector */

#define fvs_copy_vec _glp_fvs_copy_vec
void fvs_copy_vec(FVS *x, const FVS *y);
/* copy sparse vector */

#define fvs_adjust_vec _glp_fvs_adjust_vec
void fvs_adjust_vec(FVS *x, double eps);
/* replace tiny vector elements by exact zeros */

#define fvs_free_vec _glp_fvs_free_vec
void fvs_free_vec(FVS *x);
/* deallocate sparse vector */

#endif

/* eof */
