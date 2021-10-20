/* spv.h (operations on sparse vectors) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2007-2017 Free Software Foundation, Inc.
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

#ifndef SPV_H
#define SPV_H

typedef struct SPV SPV;

struct SPV
{     /* sparse vector v = (v[j]) */
      int n;
      /* dimension, n >= 0 */
      int nnz;
      /* number of non-zero components, 0 <= nnz <= n */
      int *pos; /* int pos[1+n]; */
      /* pos[j] = k, 1 <= j <= n, is position of (non-zero) v[j] in the
       * arrays ind and val, where 1 <= k <= nnz; pos[j] = 0 means that
       * v[j] is structural zero */
      int *ind; /* int ind[1+n]; */
      /* ind[k] = j, 1 <= k <= nnz, is index of v[j] */
      double *val; /* double val[1+n]; */
      /* val[k], 1 <= k <= nnz, is a numeric value of v[j] */
};

#define spv_create_vec _glp_spv_create_vec
SPV *spv_create_vec(int n);
/* create sparse vector */

#define spv_check_vec _glp_spv_check_vec
void spv_check_vec(SPV *v);
/* check that sparse vector has correct representation */

#define spv_get_vj _glp_spv_get_vj
double spv_get_vj(SPV *v, int j);
/* retrieve component of sparse vector */

#define spv_set_vj _glp_spv_set_vj
void spv_set_vj(SPV *v, int j, double val);
/* set/change component of sparse vector */

#define spv_clear_vec _glp_spv_clear_vec
void spv_clear_vec(SPV *v);
/* set all components of sparse vector to zero */

#define spv_clean_vec _glp_spv_clean_vec
void spv_clean_vec(SPV *v, double eps);
/* remove zero or small components from sparse vector */

#define spv_copy_vec _glp_spv_copy_vec
void spv_copy_vec(SPV *x, SPV *y);
/* copy sparse vector (x := y) */

#define spv_linear_comb _glp_spv_linear_comb
void spv_linear_comb(SPV *x, double a, SPV *y);
/* compute linear combination (x := x + a * y) */

#define spv_delete_vec _glp_spv_delete_vec
void spv_delete_vec(SPV *v);
/* delete sparse vector */

#endif

/* eof */
