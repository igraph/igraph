/* glplux.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
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

#include "glplux.h"
#define xfault xerror
#define dmp_create_poolx(size) dmp_create_pool()

/*----------------------------------------------------------------------
// lux_create - create LU-factorization.
//
// SYNOPSIS
//
// #include "glplux.h"
// LUX *lux_create(int n);
//
// DESCRIPTION
//
// The routine lux_create creates LU-factorization data structure for
// a matrix of the order n. Initially the factorization corresponds to
// the unity matrix (F = V = P = Q = I, so A = I).
//
// RETURNS
//
// The routine returns a pointer to the created LU-factorization data
// structure, which represents the unity matrix of the order n. */

LUX *lux_create(int n)
{     LUX *lux;
      int k;
      if (n < 1)
         xfault("lux_create: n = %d; invalid parameter\n", n);
      lux = xmalloc(sizeof(LUX));
      lux->n = n;
      lux->pool = dmp_create_poolx(sizeof(LUXELM));
      lux->F_row = xcalloc(1+n, sizeof(LUXELM *));
      lux->F_col = xcalloc(1+n, sizeof(LUXELM *));
      lux->V_piv = xcalloc(1+n, sizeof(mpq_t));
      lux->V_row = xcalloc(1+n, sizeof(LUXELM *));
      lux->V_col = xcalloc(1+n, sizeof(LUXELM *));
      lux->P_row = xcalloc(1+n, sizeof(int));
      lux->P_col = xcalloc(1+n, sizeof(int));
      lux->Q_row = xcalloc(1+n, sizeof(int));
      lux->Q_col = xcalloc(1+n, sizeof(int));
      for (k = 1; k <= n; k++)
      {  lux->F_row[k] = lux->F_col[k] = NULL;
         mpq_init(lux->V_piv[k]);
         mpq_set_si(lux->V_piv[k], 1, 1);
         lux->V_row[k] = lux->V_col[k] = NULL;
         lux->P_row[k] = lux->P_col[k] = k;
         lux->Q_row[k] = lux->Q_col[k] = k;
      }
      lux->rank = n;
      return lux;
}

/*----------------------------------------------------------------------
// initialize - initialize LU-factorization data structures.
//
// This routine initializes data structures for subsequent computing
// the LU-factorization of a given matrix A, which is specified by the
// formal routine col. On exit V = A and F = P = Q = I, where I is the
// unity matrix. */

static void initialize(LUX *lux, int (*col)(void *info, int j,
      int ind[], mpq_t val[]), void *info, LUXWKA *wka)
{     int n = lux->n;
      DMP *pool = lux->pool;
      LUXELM **F_row = lux->F_row;
      LUXELM **F_col = lux->F_col;
      mpq_t *V_piv = lux->V_piv;
      LUXELM **V_row = lux->V_row;
      LUXELM **V_col = lux->V_col;
      int *P_row = lux->P_row;
      int *P_col = lux->P_col;
      int *Q_row = lux->Q_row;
      int *Q_col = lux->Q_col;
      int *R_len = wka->R_len;
      int *R_head = wka->R_head;
      int *R_prev = wka->R_prev;
      int *R_next = wka->R_next;
      int *C_len = wka->C_len;
      int *C_head = wka->C_head;
      int *C_prev = wka->C_prev;
      int *C_next = wka->C_next;
      LUXELM *fij, *vij;
      int i, j, k, len, *ind;
      mpq_t *val;
      /* F := I */
      for (i = 1; i <= n; i++)
      {  while (F_row[i] != NULL)
         {  fij = F_row[i], F_row[i] = fij->r_next;
            mpq_clear(fij->val);
            dmp_free_atom(pool, fij, sizeof(LUXELM));
         }
      }
      for (j = 1; j <= n; j++) F_col[j] = NULL;
      /* V := 0 */
      for (k = 1; k <= n; k++) mpq_set_si(V_piv[k], 0, 1);
      for (i = 1; i <= n; i++)
      {  while (V_row[i] != NULL)
         {  vij = V_row[i], V_row[i] = vij->r_next;
            mpq_clear(vij->val);
            dmp_free_atom(pool, vij, sizeof(LUXELM));
         }
      }
      for (j = 1; j <= n; j++) V_col[j] = NULL;
      /* V := A */
      ind = xcalloc(1+n, sizeof(int));
      val = xcalloc(1+n, sizeof(mpq_t));
      for (k = 1; k <= n; k++) mpq_init(val[k]);
      for (j = 1; j <= n; j++)
      {  /* obtain j-th column of matrix A */
         len = col(info, j, ind, val);
         if (!(0 <= len && len <= n))
            xfault("lux_decomp: j = %d: len = %d; invalid column length"
               "\n", j, len);
         /* copy elements of j-th column to matrix V */
         for (k = 1; k <= len; k++)
         {  /* get row index of a[i,j] */
            i = ind[k];
            if (!(1 <= i && i <= n))
               xfault("lux_decomp: j = %d: i = %d; row index out of ran"
                  "ge\n", j, i);
            /* check for duplicate indices */
            if (V_row[i] != NULL && V_row[i]->j == j)
               xfault("lux_decomp: j = %d: i = %d; duplicate row indice"
                  "s not allowed\n", j, i);
            /* check for zero value */
            if (mpq_sgn(val[k]) == 0)
               xfault("lux_decomp: j = %d: i = %d; zero elements not al"
                  "lowed\n", j, i);
            /* add new element v[i,j] = a[i,j] to V */
            vij = dmp_get_atom(pool, sizeof(LUXELM));
            vij->i = i, vij->j = j;
            mpq_init(vij->val);
            mpq_set(vij->val, val[k]);
            vij->r_prev = NULL;
            vij->r_next = V_row[i];
            vij->c_prev = NULL;
            vij->c_next = V_col[j];
            if (vij->r_next != NULL) vij->r_next->r_prev = vij;
            if (vij->c_next != NULL) vij->c_next->c_prev = vij;
            V_row[i] = V_col[j] = vij;
         }
      }
      xfree(ind);
      for (k = 1; k <= n; k++) mpq_clear(val[k]);
      xfree(val);
      /* P := Q := I */
      for (k = 1; k <= n; k++)
         P_row[k] = P_col[k] = Q_row[k] = Q_col[k] = k;
      /* the rank of A and V is not determined yet */
      lux->rank = -1;
      /* initially the entire matrix V is active */
      /* determine its row lengths */
      for (i = 1; i <= n; i++)
      {  len = 0;
         for (vij = V_row[i]; vij != NULL; vij = vij->r_next) len++;
         R_len[i] = len;
      }
      /* build linked lists of active rows */
      for (len = 0; len <= n; len++) R_head[len] = 0;
      for (i = 1; i <= n; i++)
      {  len = R_len[i];
         R_prev[i] = 0;
         R_next[i] = R_head[len];
         if (R_next[i] != 0) R_prev[R_next[i]] = i;
         R_head[len] = i;
      }
      /* determine its column lengths */
      for (j = 1; j <= n; j++)
      {  len = 0;
         for (vij = V_col[j]; vij != NULL; vij = vij->c_next) len++;
         C_len[j] = len;
      }
      /* build linked lists of active columns */
      for (len = 0; len <= n; len++) C_head[len] = 0;
      for (j = 1; j <= n; j++)
      {  len = C_len[j];
         C_prev[j] = 0;
         C_next[j] = C_head[len];
         if (C_next[j] != 0) C_prev[C_next[j]] = j;
         C_head[len] = j;
      }
      return;
}

/*----------------------------------------------------------------------
// find_pivot - choose a pivot element.
//
// This routine chooses a pivot element v[p,q] in the active submatrix
// of matrix U = P*V*Q.
//
// It is assumed that on entry the matrix U has the following partially
// triangularized form:
//
//       1       k         n
//    1  x x x x x x x x x x
//       . x x x x x x x x x
//       . . x x x x x x x x
//       . . . x x x x x x x
//    k  . . . . * * * * * *
//       . . . . * * * * * *
//       . . . . * * * * * *
//       . . . . * * * * * *
//       . . . . * * * * * *
//    n  . . . . * * * * * *
//
// where rows and columns k, k+1, ..., n belong to the active submatrix
// (elements of the active submatrix are marked by '*').
//
// Since the matrix U = P*V*Q is not stored, the routine works with the
// matrix V. It is assumed that the row-wise representation corresponds
// to the matrix V, but the column-wise representation corresponds to
// the active submatrix of the matrix V, i.e. elements of the matrix V,
// which does not belong to the active submatrix, are missing from the
// column linked lists. It is also assumed that each active row of the
// matrix V is in the set R[len], where len is number of non-zeros in
// the row, and each active column of the matrix V is in the set C[len],
// where len is number of non-zeros in the column (in the latter case
// only elements of the active submatrix are counted; such elements are
// marked by '*' on the figure above).
//
// Due to exact arithmetic any non-zero element of the active submatrix
// can be chosen as a pivot. However, to keep sparsity of the matrix V
// the routine uses Markowitz strategy, trying to choose such element
// v[p,q], which has smallest Markowitz cost (nr[p]-1) * (nc[q]-1),
// where nr[p] and nc[q] are the number of non-zero elements, resp., in
// p-th row and in q-th column of the active submatrix.
//
// In order to reduce the search, i.e. not to walk through all elements
// of the active submatrix, the routine exploits a technique proposed by
// I.Duff. This technique is based on using the sets R[len] and C[len]
// of active rows and columns.
//
// On exit the routine returns a pointer to a pivot v[p,q] chosen, or
// NULL, if the active submatrix is empty. */

static LUXELM *find_pivot(LUX *lux, LUXWKA *wka)
{     int n = lux->n;
      LUXELM **V_row = lux->V_row;
      LUXELM **V_col = lux->V_col;
      int *R_len = wka->R_len;
      int *R_head = wka->R_head;
      int *R_next = wka->R_next;
      int *C_len = wka->C_len;
      int *C_head = wka->C_head;
      int *C_next = wka->C_next;
      LUXELM *piv, *some, *vij;
      int i, j, len, min_len, ncand, piv_lim = 5;
      double best, cost;
      /* nothing is chosen so far */
      piv = NULL, best = DBL_MAX, ncand = 0;
      /* if in the active submatrix there is a column that has the only
         non-zero (column singleton), choose it as a pivot */
      j = C_head[1];
      if (j != 0)
      {  xassert(C_len[j] == 1);
         piv = V_col[j];
         xassert(piv != NULL && piv->c_next == NULL);
         goto done;
      }
      /* if in the active submatrix there is a row that has the only
         non-zero (row singleton), choose it as a pivot */
      i = R_head[1];
      if (i != 0)
      {  xassert(R_len[i] == 1);
         piv = V_row[i];
         xassert(piv != NULL && piv->r_next == NULL);
         goto done;
      }
      /* there are no singletons in the active submatrix; walk through
         other non-empty rows and columns */
      for (len = 2; len <= n; len++)
      {  /* consider active columns having len non-zeros */
         for (j = C_head[len]; j != 0; j = C_next[j])
         {  /* j-th column has len non-zeros */
            /* find an element in the row of minimal length */
            some = NULL, min_len = INT_MAX;
            for (vij = V_col[j]; vij != NULL; vij = vij->c_next)
            {  if (min_len > R_len[vij->i])
                  some = vij, min_len = R_len[vij->i];
               /* if Markowitz cost of this element is not greater than
                  (len-1)**2, it can be chosen right now; this heuristic
                  reduces the search and works well in many cases */
               if (min_len <= len)
               {  piv = some;
                  goto done;
               }
            }
            /* j-th column has been scanned */
            /* the minimal element found is a next pivot candidate */
            xassert(some != NULL);
            ncand++;
            /* compute its Markowitz cost */
            cost = (double)(min_len - 1) * (double)(len - 1);
            /* choose between the current candidate and this element */
            if (cost < best) piv = some, best = cost;
            /* if piv_lim candidates have been considered, there is a
               doubt that a much better candidate exists; therefore it
               is the time to terminate the search */
            if (ncand == piv_lim) goto done;
         }
         /* now consider active rows having len non-zeros */
         for (i = R_head[len]; i != 0; i = R_next[i])
         {  /* i-th row has len non-zeros */
            /* find an element in the column of minimal length */
            some = NULL, min_len = INT_MAX;
            for (vij = V_row[i]; vij != NULL; vij = vij->r_next)
            {  if (min_len > C_len[vij->j])
                  some = vij, min_len = C_len[vij->j];
               /* if Markowitz cost of this element is not greater than
                  (len-1)**2, it can be chosen right now; this heuristic
                  reduces the search and works well in many cases */
               if (min_len <= len)
               {  piv = some;
                  goto done;
               }
            }
            /* i-th row has been scanned */
            /* the minimal element found is a next pivot candidate */
            xassert(some != NULL);
            ncand++;
            /* compute its Markowitz cost */
            cost = (double)(len - 1) * (double)(min_len - 1);
            /* choose between the current candidate and this element */
            if (cost < best) piv = some, best = cost;
            /* if piv_lim candidates have been considered, there is a
               doubt that a much better candidate exists; therefore it
               is the time to terminate the search */
            if (ncand == piv_lim) goto done;
         }
      }
done: /* bring the pivot v[p,q] to the factorizing routine */
      return piv;
}

/*----------------------------------------------------------------------
// eliminate - perform gaussian elimination.
//
// This routine performs elementary gaussian transformations in order
// to eliminate subdiagonal elements in the k-th column of the matrix
// U = P*V*Q using the pivot element u[k,k], where k is the number of
// the current elimination step.
//
// The parameter piv specifies the pivot element v[p,q] = u[k,k].
//
// Each time when the routine applies the elementary transformation to
// a non-pivot row of the matrix V, it stores the corresponding element
// to the matrix F in order to keep the main equality A = F*V.
//
// The routine assumes that on entry the matrices L = P*F*inv(P) and
// U = P*V*Q are the following:
//
//       1       k                  1       k         n
//    1  1 . . . . . . . . .     1  x x x x x x x x x x
//       x 1 . . . . . . . .        . x x x x x x x x x
//       x x 1 . . . . . . .        . . x x x x x x x x
//       x x x 1 . . . . . .        . . . x x x x x x x
//    k  x x x x 1 . . . . .     k  . . . . * * * * * *
//       x x x x _ 1 . . . .        . . . . # * * * * *
//       x x x x _ . 1 . . .        . . . . # * * * * *
//       x x x x _ . . 1 . .        . . . . # * * * * *
//       x x x x _ . . . 1 .        . . . . # * * * * *
//    n  x x x x _ . . . . 1     n  . . . . # * * * * *
//
//            matrix L                   matrix U
//
// where rows and columns of the matrix U with numbers k, k+1, ..., n
// form the active submatrix (eliminated elements are marked by '#' and
// other elements of the active submatrix are marked by '*'). Note that
// each eliminated non-zero element u[i,k] of the matrix U gives the
// corresponding element l[i,k] of the matrix L (marked by '_').
//
// Actually all operations are performed on the matrix V. Should note
// that the row-wise representation corresponds to the matrix V, but the
// column-wise representation corresponds to the active submatrix of the
// matrix V, i.e. elements of the matrix V, which doesn't belong to the
// active submatrix, are missing from the column linked lists.
//
// Let u[k,k] = v[p,q] be the pivot. In order to eliminate subdiagonal
// elements u[i',k] = v[i,q], i' = k+1, k+2, ..., n, the routine applies
// the following elementary gaussian transformations:
//
//    (i-th row of V) := (i-th row of V) - f[i,p] * (p-th row of V),
//
// where f[i,p] = v[i,q] / v[p,q] is a gaussian multiplier.
//
// Additionally, in order to keep the main equality A = F*V, each time
// when the routine applies the transformation to i-th row of the matrix
// V, it also adds f[i,p] as a new element to the matrix F.
//
// IMPORTANT: On entry the working arrays flag and work should contain
// zeros. This status is provided by the routine on exit. */

static void eliminate(LUX *lux, LUXWKA *wka, LUXELM *piv, int flag[],
      mpq_t work[])
{     DMP *pool = lux->pool;
      LUXELM **F_row = lux->F_row;
      LUXELM **F_col = lux->F_col;
      mpq_t *V_piv = lux->V_piv;
      LUXELM **V_row = lux->V_row;
      LUXELM **V_col = lux->V_col;
      int *R_len = wka->R_len;
      int *R_head = wka->R_head;
      int *R_prev = wka->R_prev;
      int *R_next = wka->R_next;
      int *C_len = wka->C_len;
      int *C_head = wka->C_head;
      int *C_prev = wka->C_prev;
      int *C_next = wka->C_next;
      LUXELM *fip, *vij, *vpj, *viq, *next;
      mpq_t temp;
      int i, j, p, q;
      mpq_init(temp);
      /* determine row and column indices of the pivot v[p,q] */
      xassert(piv != NULL);
      p = piv->i, q = piv->j;
      /* remove p-th (pivot) row from the active set; it will never
         return there */
      if (R_prev[p] == 0)
         R_head[R_len[p]] = R_next[p];
      else
         R_next[R_prev[p]] = R_next[p];
      if (R_next[p] == 0)
         ;
      else
         R_prev[R_next[p]] = R_prev[p];
      /* remove q-th (pivot) column from the active set; it will never
         return there */
      if (C_prev[q] == 0)
         C_head[C_len[q]] = C_next[q];
      else
         C_next[C_prev[q]] = C_next[q];
      if (C_next[q] == 0)
         ;
      else
         C_prev[C_next[q]] = C_prev[q];
      /* store the pivot value in a separate array */
      mpq_set(V_piv[p], piv->val);
      /* remove the pivot from p-th row */
      if (piv->r_prev == NULL)
         V_row[p] = piv->r_next;
      else
         piv->r_prev->r_next = piv->r_next;
      if (piv->r_next == NULL)
         ;
      else
         piv->r_next->r_prev = piv->r_prev;
      R_len[p]--;
      /* remove the pivot from q-th column */
      if (piv->c_prev == NULL)
         V_col[q] = piv->c_next;
      else
         piv->c_prev->c_next = piv->c_next;
      if (piv->c_next == NULL)
         ;
      else
         piv->c_next->c_prev = piv->c_prev;
      C_len[q]--;
      /* free the space occupied by the pivot */
      mpq_clear(piv->val);
      dmp_free_atom(pool, piv, sizeof(LUXELM));
      /* walk through p-th (pivot) row, which already does not contain
         the pivot v[p,q], and do the following... */
      for (vpj = V_row[p]; vpj != NULL; vpj = vpj->r_next)
      {  /* get column index of v[p,j] */
         j = vpj->j;
         /* store v[p,j] in the working array */
         flag[j] = 1;
         mpq_set(work[j], vpj->val);
         /* remove j-th column from the active set; it will return there
            later with a new length */
         if (C_prev[j] == 0)
            C_head[C_len[j]] = C_next[j];
         else
            C_next[C_prev[j]] = C_next[j];
         if (C_next[j] == 0)
            ;
         else
            C_prev[C_next[j]] = C_prev[j];
         /* v[p,j] leaves the active submatrix, so remove it from j-th
            column; however, v[p,j] is kept in p-th row */
         if (vpj->c_prev == NULL)
            V_col[j] = vpj->c_next;
         else
            vpj->c_prev->c_next = vpj->c_next;
         if (vpj->c_next == NULL)
            ;
         else
            vpj->c_next->c_prev = vpj->c_prev;
         C_len[j]--;
      }
      /* now walk through q-th (pivot) column, which already does not
         contain the pivot v[p,q], and perform gaussian elimination */
      while (V_col[q] != NULL)
      {  /* element v[i,q] has to be eliminated */
         viq = V_col[q];
         /* get row index of v[i,q] */
         i = viq->i;
         /* remove i-th row from the active set; later it will return
            there with a new length */
         if (R_prev[i] == 0)
            R_head[R_len[i]] = R_next[i];
         else
            R_next[R_prev[i]] = R_next[i];
         if (R_next[i] == 0)
            ;
         else
            R_prev[R_next[i]] = R_prev[i];
         /* compute gaussian multiplier f[i,p] = v[i,q] / v[p,q] and
            store it in the matrix F */
         fip = dmp_get_atom(pool, sizeof(LUXELM));
         fip->i = i, fip->j = p;
         mpq_init(fip->val);
         mpq_div(fip->val, viq->val, V_piv[p]);
         fip->r_prev = NULL;
         fip->r_next = F_row[i];
         fip->c_prev = NULL;
         fip->c_next = F_col[p];
         if (fip->r_next != NULL) fip->r_next->r_prev = fip;
         if (fip->c_next != NULL) fip->c_next->c_prev = fip;
         F_row[i] = F_col[p] = fip;
         /* v[i,q] has to be eliminated, so remove it from i-th row */
         if (viq->r_prev == NULL)
            V_row[i] = viq->r_next;
         else
            viq->r_prev->r_next = viq->r_next;
         if (viq->r_next == NULL)
            ;
         else
            viq->r_next->r_prev = viq->r_prev;
         R_len[i]--;
         /* and also from q-th column */
         V_col[q] = viq->c_next;
         C_len[q]--;
         /* free the space occupied by v[i,q] */
         mpq_clear(viq->val);
         dmp_free_atom(pool, viq, sizeof(LUXELM));
         /* perform gaussian transformation:
            (i-th row) := (i-th row) - f[i,p] * (p-th row)
            note that now p-th row, which is in the working array,
            does not contain the pivot v[p,q], and i-th row does not
            contain the element v[i,q] to be eliminated */
         /* walk through i-th row and transform existing non-zero
            elements */
         for (vij = V_row[i]; vij != NULL; vij = next)
         {  next = vij->r_next;
            /* get column index of v[i,j] */
            j = vij->j;
            /* v[i,j] := v[i,j] - f[i,p] * v[p,j] */
            if (flag[j])
            {  /* v[p,j] != 0 */
               flag[j] = 0;
               mpq_mul(temp, fip->val, work[j]);
               mpq_sub(vij->val, vij->val, temp);
               if (mpq_sgn(vij->val) == 0)
               {  /* new v[i,j] is zero, so remove it from the active
                     submatrix */
                  /* remove v[i,j] from i-th row */
                  if (vij->r_prev == NULL)
                     V_row[i] = vij->r_next;
                  else
                     vij->r_prev->r_next = vij->r_next;
                  if (vij->r_next == NULL)
                     ;
                  else
                     vij->r_next->r_prev = vij->r_prev;
                  R_len[i]--;
                  /* remove v[i,j] from j-th column */
                  if (vij->c_prev == NULL)
                     V_col[j] = vij->c_next;
                  else
                     vij->c_prev->c_next = vij->c_next;
                  if (vij->c_next == NULL)
                     ;
                  else
                     vij->c_next->c_prev = vij->c_prev;
                  C_len[j]--;
                  /* free the space occupied by v[i,j] */
                  mpq_clear(vij->val);
                  dmp_free_atom(pool, vij, sizeof(LUXELM));
               }
            }
         }
         /* now flag is the pattern of the set v[p,*] \ v[i,*] */
         /* walk through p-th (pivot) row and create new elements in
            i-th row, which appear due to fill-in */
         for (vpj = V_row[p]; vpj != NULL; vpj = vpj->r_next)
         {  j = vpj->j;
            if (flag[j])
            {  /* create new non-zero v[i,j] = 0 - f[i,p] * v[p,j] and
                  add it to i-th row and j-th column */
               vij = dmp_get_atom(pool, sizeof(LUXELM));
               vij->i = i, vij->j = j;
               mpq_init(vij->val);
               mpq_mul(vij->val, fip->val, work[j]);
               mpq_neg(vij->val, vij->val);
               vij->r_prev = NULL;
               vij->r_next = V_row[i];
               vij->c_prev = NULL;
               vij->c_next = V_col[j];
               if (vij->r_next != NULL) vij->r_next->r_prev = vij;
               if (vij->c_next != NULL) vij->c_next->c_prev = vij;
               V_row[i] = V_col[j] = vij;
               R_len[i]++, C_len[j]++;
            }
            else
            {  /* there is no fill-in, because v[i,j] already exists in
                  i-th row; restore the flag, which was reset before */
               flag[j] = 1;
            }
         }
         /* now i-th row has been completely transformed and can return
            to the active set with a new length */
         R_prev[i] = 0;
         R_next[i] = R_head[R_len[i]];
         if (R_next[i] != 0) R_prev[R_next[i]] = i;
         R_head[R_len[i]] = i;
      }
      /* at this point q-th (pivot) column must be empty */
      xassert(C_len[q] == 0);
      /* walk through p-th (pivot) row again and do the following... */
      for (vpj = V_row[p]; vpj != NULL; vpj = vpj->r_next)
      {  /* get column index of v[p,j] */
         j = vpj->j;
         /* erase v[p,j] from the working array */
         flag[j] = 0;
         mpq_set_si(work[j], 0, 1);
         /* now j-th column has been completely transformed, so it can
            return to the active list with a new length */
         C_prev[j] = 0;
         C_next[j] = C_head[C_len[j]];
         if (C_next[j] != 0) C_prev[C_next[j]] = j;
         C_head[C_len[j]] = j;
      }
      mpq_clear(temp);
      /* return to the factorizing routine */
      return;
}

/*----------------------------------------------------------------------
// lux_decomp - compute LU-factorization.
//
// SYNOPSIS
//
// #include "glplux.h"
// int lux_decomp(LUX *lux, int (*col)(void *info, int j, int ind[],
//    mpq_t val[]), void *info);
//
// DESCRIPTION
//
// The routine lux_decomp computes LU-factorization of a given square
// matrix A.
//
// The parameter lux specifies LU-factorization data structure built by
// means of the routine lux_create.
//
// The formal routine col specifies the original matrix A. In order to
// obtain j-th column of the matrix A the routine lux_decomp calls the
// routine col with the parameter j (1 <= j <= n, where n is the order
// of A). In response the routine col should store row indices and
// numerical values of non-zero elements of j-th column of A to the
// locations ind[1], ..., ind[len] and val[1], ..., val[len], resp.,
// where len is the number of non-zeros in j-th column, which should be
// returned on exit. Neiter zero nor duplicate elements are allowed.
//
// The parameter info is a transit pointer passed to the formal routine
// col; it can be used for various purposes.
//
// RETURNS
//
// The routine lux_decomp returns the singularity flag. Zero flag means
// that the original matrix A is non-singular while non-zero flag means
// that A is (exactly!) singular.
//
// Note that LU-factorization is valid in both cases, however, in case
// of singularity some rows of the matrix V (including pivot elements)
// will be empty.
//
// REPAIRING SINGULAR MATRIX
//
// If the routine lux_decomp returns non-zero flag, it provides all
// necessary information that can be used for "repairing" the matrix A,
// where "repairing" means replacing linearly dependent columns of the
// matrix A by appropriate columns of the unity matrix. This feature is
// needed when the routine lux_decomp is used for reinverting the basis
// matrix within the simplex method procedure.
//
// On exit linearly dependent columns of the matrix U have the numbers
// rank+1, rank+2, ..., n, where rank is the exact rank of the matrix A
// stored by the routine to the member lux->rank. The correspondence
// between columns of A and U is the same as between columns of V and U.
// Thus, linearly dependent columns of the matrix A have the numbers
// Q_col[rank+1], Q_col[rank+2], ..., Q_col[n], where Q_col is an array
// representing the permutation matrix Q in column-like format. It is
// understood that each j-th linearly dependent column of the matrix U
// should be replaced by the unity vector, where all elements are zero
// except the unity diagonal element u[j,j]. On the other hand j-th row
// of the matrix U corresponds to the row of the matrix V (and therefore
// of the matrix A) with the number P_row[j], where P_row is an array
// representing the permutation matrix P in row-like format. Thus, each
// j-th linearly dependent column of the matrix U should be replaced by
// a column of the unity matrix with the number P_row[j].
//
// The code that repairs the matrix A may look like follows:
//
//    for (j = rank+1; j <= n; j++)
//    {  replace column Q_col[j] of the matrix A by column P_row[j] of
//       the unity matrix;
//    }
//
// where rank, P_row, and Q_col are members of the structure LUX. */

int lux_decomp(LUX *lux, int (*col)(void *info, int j, int ind[],
      mpq_t val[]), void *info)
{     int n = lux->n;
      LUXELM **V_row = lux->V_row;
      LUXELM **V_col = lux->V_col;
      int *P_row = lux->P_row;
      int *P_col = lux->P_col;
      int *Q_row = lux->Q_row;
      int *Q_col = lux->Q_col;
      LUXELM *piv, *vij;
      LUXWKA *wka;
      int i, j, k, p, q, t, *flag;
      mpq_t *work;
      /* allocate working area */
      wka = xmalloc(sizeof(LUXWKA));
      wka->R_len = xcalloc(1+n, sizeof(int));
      wka->R_head = xcalloc(1+n, sizeof(int));
      wka->R_prev = xcalloc(1+n, sizeof(int));
      wka->R_next = xcalloc(1+n, sizeof(int));
      wka->C_len = xcalloc(1+n, sizeof(int));
      wka->C_head = xcalloc(1+n, sizeof(int));
      wka->C_prev = xcalloc(1+n, sizeof(int));
      wka->C_next = xcalloc(1+n, sizeof(int));
      /* initialize LU-factorization data structures */
      initialize(lux, col, info, wka);
      /* allocate working arrays */
      flag = xcalloc(1+n, sizeof(int));
      work = xcalloc(1+n, sizeof(mpq_t));
      for (k = 1; k <= n; k++)
      {  flag[k] = 0;
         mpq_init(work[k]);
      }
      /* main elimination loop */
      for (k = 1; k <= n; k++)
      {  /* choose a pivot element v[p,q] */
         piv = find_pivot(lux, wka);
         if (piv == NULL)
         {  /* no pivot can be chosen, because the active submatrix is
               empty */
            break;
         }
         /* determine row and column indices of the pivot element */
         p = piv->i, q = piv->j;
         /* let v[p,q] correspond to u[i',j']; permute k-th and i'-th
            rows and k-th and j'-th columns of the matrix U = P*V*Q to
            move the element u[i',j'] to the position u[k,k] */
         i = P_col[p], j = Q_row[q];
         xassert(k <= i && i <= n && k <= j && j <= n);
         /* permute k-th and i-th rows of the matrix U */
         t = P_row[k];
         P_row[i] = t, P_col[t] = i;
         P_row[k] = p, P_col[p] = k;
         /* permute k-th and j-th columns of the matrix U */
         t = Q_col[k];
         Q_col[j] = t, Q_row[t] = j;
         Q_col[k] = q, Q_row[q] = k;
         /* eliminate subdiagonal elements of k-th column of the matrix
            U = P*V*Q using the pivot element u[k,k] = v[p,q] */
         eliminate(lux, wka, piv, flag, work);
      }
      /* determine the rank of A (and V) */
      lux->rank = k - 1;
      /* free working arrays */
      xfree(flag);
      for (k = 1; k <= n; k++) mpq_clear(work[k]);
      xfree(work);
      /* build column lists of the matrix V using its row lists */
      for (j = 1; j <= n; j++)
         xassert(V_col[j] == NULL);
      for (i = 1; i <= n; i++)
      {  for (vij = V_row[i]; vij != NULL; vij = vij->r_next)
         {  j = vij->j;
            vij->c_prev = NULL;
            vij->c_next = V_col[j];
            if (vij->c_next != NULL) vij->c_next->c_prev = vij;
            V_col[j] = vij;
         }
      }
      /* free working area */
      xfree(wka->R_len);
      xfree(wka->R_head);
      xfree(wka->R_prev);
      xfree(wka->R_next);
      xfree(wka->C_len);
      xfree(wka->C_head);
      xfree(wka->C_prev);
      xfree(wka->C_next);
      xfree(wka);
      /* return to the calling program */
      return (lux->rank < n);
}

/*----------------------------------------------------------------------
// lux_f_solve - solve system F*x = b or F'*x = b.
//
// SYNOPSIS
//
// #include "glplux.h"
// void lux_f_solve(LUX *lux, int tr, mpq_t x[]);
//
// DESCRIPTION
//
// The routine lux_f_solve solves either the system F*x = b (if the
// flag tr is zero) or the system F'*x = b (if the flag tr is non-zero),
// where the matrix F is a component of LU-factorization specified by
// the parameter lux, F' is a matrix transposed to F.
//
// On entry the array x should contain elements of the right-hand side
// vector b in locations x[1], ..., x[n], where n is the order of the
// matrix F. On exit this array will contain elements of the solution
// vector x in the same locations. */

void lux_f_solve(LUX *lux, int tr, mpq_t x[])
{     int n = lux->n;
      LUXELM **F_row = lux->F_row;
      LUXELM **F_col = lux->F_col;
      int *P_row = lux->P_row;
      LUXELM *fik, *fkj;
      int i, j, k;
      mpq_t temp;
      mpq_init(temp);
      if (!tr)
      {  /* solve the system F*x = b */
         for (j = 1; j <= n; j++)
         {  k = P_row[j];
            if (mpq_sgn(x[k]) != 0)
            {  for (fik = F_col[k]; fik != NULL; fik = fik->c_next)
               {  mpq_mul(temp, fik->val, x[k]);
                  mpq_sub(x[fik->i], x[fik->i], temp);
               }
            }
         }
      }
      else
      {  /* solve the system F'*x = b */
         for (i = n; i >= 1; i--)
         {  k = P_row[i];
            if (mpq_sgn(x[k]) != 0)
            {  for (fkj = F_row[k]; fkj != NULL; fkj = fkj->r_next)
               {  mpq_mul(temp, fkj->val, x[k]);
                  mpq_sub(x[fkj->j], x[fkj->j], temp);
               }
            }
         }
      }
      mpq_clear(temp);
      return;
}

/*----------------------------------------------------------------------
// lux_v_solve - solve system V*x = b or V'*x = b.
//
// SYNOPSIS
//
// #include "glplux.h"
// void lux_v_solve(LUX *lux, int tr, double x[]);
//
// DESCRIPTION
//
// The routine lux_v_solve solves either the system V*x = b (if the
// flag tr is zero) or the system V'*x = b (if the flag tr is non-zero),
// where the matrix V is a component of LU-factorization specified by
// the parameter lux, V' is a matrix transposed to V.
//
// On entry the array x should contain elements of the right-hand side
// vector b in locations x[1], ..., x[n], where n is the order of the
// matrix V. On exit this array will contain elements of the solution
// vector x in the same locations. */

void lux_v_solve(LUX *lux, int tr, mpq_t x[])
{     int n = lux->n;
      mpq_t *V_piv = lux->V_piv;
      LUXELM **V_row = lux->V_row;
      LUXELM **V_col = lux->V_col;
      int *P_row = lux->P_row;
      int *Q_col = lux->Q_col;
      LUXELM *vij;
      int i, j, k;
      mpq_t *b, temp;
      b = xcalloc(1+n, sizeof(mpq_t));
      for (k = 1; k <= n; k++)
         mpq_init(b[k]), mpq_set(b[k], x[k]), mpq_set_si(x[k], 0, 1);
      mpq_init(temp);
      if (!tr)
      {  /* solve the system V*x = b */
         for (k = n; k >= 1; k--)
         {  i = P_row[k], j = Q_col[k];
            if (mpq_sgn(b[i]) != 0)
            {  mpq_set(x[j], b[i]);
               mpq_div(x[j], x[j], V_piv[i]);
               for (vij = V_col[j]; vij != NULL; vij = vij->c_next)
               {  mpq_mul(temp, vij->val, x[j]);
                  mpq_sub(b[vij->i], b[vij->i], temp);
               }
            }
         }
      }
      else
      {  /* solve the system V'*x = b */
         for (k = 1; k <= n; k++)
         {  i = P_row[k], j = Q_col[k];
            if (mpq_sgn(b[j]) != 0)
            {  mpq_set(x[i], b[j]);
               mpq_div(x[i], x[i], V_piv[i]);
               for (vij = V_row[i]; vij != NULL; vij = vij->r_next)
               {  mpq_mul(temp, vij->val, x[i]);
                  mpq_sub(b[vij->j], b[vij->j], temp);
               }
            }
         }
      }
      for (k = 1; k <= n; k++) mpq_clear(b[k]);
      mpq_clear(temp);
      xfree(b);
      return;
}

/*----------------------------------------------------------------------
// lux_solve - solve system A*x = b or A'*x = b.
//
// SYNOPSIS
//
// #include "glplux.h"
// void lux_solve(LUX *lux, int tr, mpq_t x[]);
//
// DESCRIPTION
//
// The routine lux_solve solves either the system A*x = b (if the flag
// tr is zero) or the system A'*x = b (if the flag tr is non-zero),
// where the parameter lux specifies LU-factorization of the matrix A,
// A' is a matrix transposed to A.
//
// On entry the array x should contain elements of the right-hand side
// vector b in locations x[1], ..., x[n], where n is the order of the
// matrix A. On exit this array will contain elements of the solution
// vector x in the same locations. */

void lux_solve(LUX *lux, int tr, mpq_t x[])
{     if (lux->rank < lux->n)
         xfault("lux_solve: LU-factorization has incomplete rank\n");
      if (!tr)
      {  /* A = F*V, therefore inv(A) = inv(V)*inv(F) */
         lux_f_solve(lux, 0, x);
         lux_v_solve(lux, 0, x);
      }
      else
      {  /* A' = V'*F', therefore inv(A') = inv(F')*inv(V') */
         lux_v_solve(lux, 1, x);
         lux_f_solve(lux, 1, x);
      }
      return;
}

/*----------------------------------------------------------------------
// lux_delete - delete LU-factorization.
//
// SYNOPSIS
//
// #include "glplux.h"
// void lux_delete(LUX *lux);
//
// DESCRIPTION
//
// The routine lux_delete deletes LU-factorization data structure,
// which the parameter lux points to, freeing all the memory allocated
// to this object. */

void lux_delete(LUX *lux)
{     int n = lux->n;
      LUXELM *fij, *vij;
      int i;
      for (i = 1; i <= n; i++)
      {  for (fij = lux->F_row[i]; fij != NULL; fij = fij->r_next)
            mpq_clear(fij->val);
         mpq_clear(lux->V_piv[i]);
         for (vij = lux->V_row[i]; vij != NULL; vij = vij->r_next)
            mpq_clear(vij->val);
      }
      dmp_delete_pool(lux->pool);
      xfree(lux->F_row);
      xfree(lux->F_col);
      xfree(lux->V_piv);
      xfree(lux->V_row);
      xfree(lux->V_col);
      xfree(lux->P_row);
      xfree(lux->P_col);
      xfree(lux->Q_row);
      xfree(lux->Q_col);
      xfree(lux);
      return;
}

/* eof */
