/* npp6.c (translate feasibility problem to CNF-SAT) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2011-2017 Free Software Foundation, Inc.
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

#include "env.h"
#include "npp.h"

/***********************************************************************
*  npp_sat_free_row - process free (unbounded) row
*
*  This routine processes row p, which is free (i.e. has no finite
*  bounds):
*
*     -inf < sum a[p,j] x[j] < +inf.                                 (1)
*
*  The constraint (1) cannot be active and therefore it is redundant,
*  so the routine simply removes it from the original problem. */

void npp_sat_free_row(NPP *npp, NPPROW *p)
{     /* the row should be free */
      xassert(p->lb == -DBL_MAX && p->ub == +DBL_MAX);
      /* remove the row from the problem */
      npp_del_row(npp, p);
      return;
}

/***********************************************************************
*  npp_sat_fixed_col - process fixed column
*
*  This routine processes column q, which is fixed:
*
*     x[q] = s[q],                                                   (1)
*
*  where s[q] is a fixed column value.
*
*  The routine substitutes fixed value s[q] into constraint rows and
*  then removes column x[q] from the original problem.
*
*  Substitution of x[q] = s[q] into row i gives:
*
*     L[i] <= sum a[i,j] x[j] <= U[i]   ==>
*              j
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] x[q] <= U[i]   ==>
*            j!=q
*
*     L[i] <= sum a[i,j] x[j] + a[i,q] s[q] <= U[i]   ==>
*            j!=q
*
*     L~[i] <= sum a[i,j] x[j] <= U~[i],
*             j!=q
*
*  where
*
*     L~[i] = L[i] - a[i,q] s[q],                                    (2)
*
*     U~[i] = U[i] - a[i,q] s[q]                                     (3)
*
*  are, respectively, lower and upper bound of row i in the transformed
*  problem.
*
*  On recovering solution x[q] is assigned the value of s[q]. */

struct sat_fixed_col
{     /* fixed column */
      int q;
      /* column reference number for variable x[q] */
      int s;
      /* value, at which x[q] is fixed */
};

static int rcv_sat_fixed_col(NPP *, void *);

int npp_sat_fixed_col(NPP *npp, NPPCOL *q)
{     struct sat_fixed_col *info;
      NPPROW *i;
      NPPAIJ *aij;
      int temp;
      /* the column should be fixed */
      xassert(q->lb == q->ub);
      /* create transformation stack entry */
      info = npp_push_tse(npp,
         rcv_sat_fixed_col, sizeof(struct sat_fixed_col));
      info->q = q->j;
      info->s = (int)q->lb;
      xassert((double)info->s == q->lb);
      /* substitute x[q] = s[q] into constraint rows */
      if (info->s == 0)
         goto skip;
      for (aij = q->ptr; aij != NULL; aij = aij->c_next)
      {  i = aij->row;
         if (i->lb != -DBL_MAX)
         {  i->lb -= aij->val * (double)info->s;
            temp = (int)i->lb;
            if ((double)temp != i->lb)
               return 1; /* integer arithmetic error */
         }
         if (i->ub != +DBL_MAX)
         {  i->ub -= aij->val * (double)info->s;
            temp = (int)i->ub;
            if ((double)temp != i->ub)
               return 2; /* integer arithmetic error */
         }
      }
skip: /* remove the column from the problem */
      npp_del_col(npp, q);
      return 0;
}

static int rcv_sat_fixed_col(NPP *npp, void *info_)
{     struct sat_fixed_col *info = info_;
      npp->c_value[info->q] = (double)info->s;
      return 0;
}

/***********************************************************************
*  npp_sat_is_bin_comb - test if row is binary combination
*
*  This routine tests if the specified row is a binary combination,
*  i.e. all its constraint coefficients are +1 and -1 and all variables
*  are binary. If the test was passed, the routine returns non-zero,
*  otherwise zero. */

int npp_sat_is_bin_comb(NPP *npp, NPPROW *row)
{     NPPCOL *col;
      NPPAIJ *aij;
      xassert(npp == npp);
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  if (!(aij->val == +1.0 || aij->val == -1.0))
            return 0; /* non-unity coefficient */
         col = aij->col;
         if (!(col->is_int && col->lb == 0.0 && col->ub == 1.0))
            return 0; /* non-binary column */
      }
      return 1; /* test was passed */
}

/***********************************************************************
*  npp_sat_num_pos_coef - determine number of positive coefficients
*
*  This routine returns the number of positive coefficients in the
*  specified row. */

int npp_sat_num_pos_coef(NPP *npp, NPPROW *row)
{     NPPAIJ *aij;
      int num = 0;
      xassert(npp == npp);
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  if (aij->val > 0.0)
            num++;
      }
      return num;
}

/***********************************************************************
*  npp_sat_num_neg_coef - determine number of negative coefficients
*
*  This routine returns the number of negative coefficients in the
*  specified row. */

int npp_sat_num_neg_coef(NPP *npp, NPPROW *row)
{     NPPAIJ *aij;
      int num = 0;
      xassert(npp == npp);
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  if (aij->val < 0.0)
            num++;
      }
      return num;
}

/***********************************************************************
*  npp_sat_is_cover_ineq - test if row is covering inequality
*
*  The canonical form of a covering inequality is the following:
*
*     sum x[j] >= 1,                                                 (1)
*   j in J
*
*  where all x[j] are binary variables.
*
*  In general case a covering inequality may have one of the following
*  two forms:
*
*     sum  x[j] -  sum  x[j] >= 1 - |J-|,                            (2)
*   j in J+      j in J-
*
*
*     sum  x[j] -  sum  x[j] <= |J+| - 1.                            (3)
*   j in J+      j in J-
*
*  Obviously, the inequality (2) can be transformed to the form (1) by
*  substitution x[j] = 1 - x'[j] for all j in J-, where x'[j] is the
*  negation of variable x[j]. And the inequality (3) can be transformed
*  to (2) by multiplying both left- and right-hand sides by -1.
*
*  This routine returns one of the following codes:
*
*  0, if the specified row is not a covering inequality;
*
*  1, if the specified row has the form (2);
*
*  2, if the specified row has the form (3). */

int npp_sat_is_cover_ineq(NPP *npp, NPPROW *row)
{     xassert(npp == npp);
      if (row->lb != -DBL_MAX && row->ub == +DBL_MAX)
      {  /* row is inequality of '>=' type */
         if (npp_sat_is_bin_comb(npp, row))
         {  /* row is a binary combination */
            if (row->lb == 1.0 - npp_sat_num_neg_coef(npp, row))
            {  /* row has the form (2) */
               return 1;
            }
         }
      }
      else if (row->lb == -DBL_MAX && row->ub != +DBL_MAX)
      {  /* row is inequality of '<=' type */
         if (npp_sat_is_bin_comb(npp, row))
         {  /* row is a binary combination */
            if (row->ub == npp_sat_num_pos_coef(npp, row) - 1.0)
            {  /* row has the form (3) */
               return 2;
            }
         }
      }
      /* row is not a covering inequality */
      return 0;
}

/***********************************************************************
*  npp_sat_is_pack_ineq - test if row is packing inequality
*
*  The canonical form of a packing inequality is the following:
*
*     sum x[j] <= 1,                                                 (1)
*   j in J
*
*  where all x[j] are binary variables.
*
*  In general case a packing inequality may have one of the following
*  two forms:
*
*     sum  x[j] -  sum  x[j] <= 1 - |J-|,                            (2)
*   j in J+      j in J-
*
*
*     sum  x[j] -  sum  x[j] >= |J+| - 1.                            (3)
*   j in J+      j in J-
*
*  Obviously, the inequality (2) can be transformed to the form (1) by
*  substitution x[j] = 1 - x'[j] for all j in J-, where x'[j] is the
*  negation of variable x[j]. And the inequality (3) can be transformed
*  to (2) by multiplying both left- and right-hand sides by -1.
*
*  This routine returns one of the following codes:
*
*  0, if the specified row is not a packing inequality;
*
*  1, if the specified row has the form (2);
*
*  2, if the specified row has the form (3). */

int npp_sat_is_pack_ineq(NPP *npp, NPPROW *row)
{     xassert(npp == npp);
      if (row->lb == -DBL_MAX && row->ub != +DBL_MAX)
      {  /* row is inequality of '<=' type */
         if (npp_sat_is_bin_comb(npp, row))
         {  /* row is a binary combination */
            if (row->ub == 1.0 - npp_sat_num_neg_coef(npp, row))
            {  /* row has the form (2) */
               return 1;
            }
         }
      }
      else if (row->lb != -DBL_MAX && row->ub == +DBL_MAX)
      {  /* row is inequality of '>=' type */
         if (npp_sat_is_bin_comb(npp, row))
         {  /* row is a binary combination */
            if (row->lb == npp_sat_num_pos_coef(npp, row) - 1.0)
            {  /* row has the form (3) */
               return 2;
            }
         }
      }
      /* row is not a packing inequality */
      return 0;
}

/***********************************************************************
*  npp_sat_is_partn_eq - test if row is partitioning equality
*
*  The canonical form of a partitioning equality is the following:
*
*     sum x[j] = 1,                                                  (1)
*   j in J
*
*  where all x[j] are binary variables.
*
*  In general case a partitioning equality may have one of the following
*  two forms:
*
*     sum  x[j] -  sum  x[j] = 1 - |J-|,                             (2)
*   j in J+      j in J-
*
*
*     sum  x[j] -  sum  x[j] = |J+| - 1.                             (3)
*   j in J+      j in J-
*
*  Obviously, the equality (2) can be transformed to the form (1) by
*  substitution x[j] = 1 - x'[j] for all j in J-, where x'[j] is the
*  negation of variable x[j]. And the equality (3) can be transformed
*  to (2) by multiplying both left- and right-hand sides by -1.
*
*  This routine returns one of the following codes:
*
*  0, if the specified row is not a partitioning equality;
*
*  1, if the specified row has the form (2);
*
*  2, if the specified row has the form (3). */

int npp_sat_is_partn_eq(NPP *npp, NPPROW *row)
{     xassert(npp == npp);
      if (row->lb == row->ub)
      {  /* row is equality constraint */
         if (npp_sat_is_bin_comb(npp, row))
         {  /* row is a binary combination */
            if (row->lb == 1.0 - npp_sat_num_neg_coef(npp, row))
            {  /* row has the form (2) */
               return 1;
            }
            if (row->ub == npp_sat_num_pos_coef(npp, row) - 1.0)
            {  /* row has the form (3) */
               return 2;
            }
         }
      }
      /* row is not a partitioning equality */
      return 0;
}

/***********************************************************************
*  npp_sat_reverse_row - multiply both sides of row by -1
*
*  This routines multiplies by -1 both left- and right-hand sides of
*  the specified row:
*
*     L <= sum x[j] <= U,
*
*  that results in the following row:
*
*     -U <= sum (-x[j]) <= -L.
*
*  If no integer overflow occured, the routine returns zero, otherwise
*  non-zero. */

int npp_sat_reverse_row(NPP *npp, NPPROW *row)
{     NPPAIJ *aij;
      int temp, ret = 0;
      double old_lb, old_ub;
      xassert(npp == npp);
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  aij->val = -aij->val;
         temp = (int)aij->val;
         if ((double)temp != aij->val)
            ret = 1;
      }
      old_lb = row->lb, old_ub = row->ub;
      if (old_ub == +DBL_MAX)
         row->lb = -DBL_MAX;
      else
      {  row->lb = -old_ub;
         temp = (int)row->lb;
         if ((double)temp != row->lb)
            ret = 2;
      }
      if (old_lb == -DBL_MAX)
         row->ub = +DBL_MAX;
      else
      {  row->ub = -old_lb;
         temp = (int)row->ub;
         if ((double)temp != row->ub)
            ret = 3;
      }
      return ret;
}

/***********************************************************************
*  npp_sat_split_pack - split packing inequality
*
*  Let there be given a packing inequality in canonical form:
*
*     sum  t[j] <= 1,                                                (1)
*   j in J
*
*  where t[j] = x[j] or t[j] = 1 - x[j], x[j] is a binary variable.
*  And let J = J1 U J2 is a partition of the set of literals. Then the
*  inequality (1) is obviously equivalent to the following two packing
*  inequalities:
*
*     sum  t[j] <= y       <-->   sum  t[j] + (1 - y) <= 1,          (2)
*   j in J1                     j in J1
*
*     sum  t[j] <= 1 - y   <-->   sum  t[j] + y       <= 1,          (3)
*   j in J2                     j in J2
*
*  where y is a new binary variable added to the transformed problem.
*
*  Assuming that the specified row is a packing inequality (1), this
*  routine constructs the set J1 by including there first nlit literals
*  (terms) from the specified row, and the set J2 = J \ J1. Then the
*  routine creates a new row, which corresponds to inequality (2), and
*  replaces the specified row with inequality (3). */

NPPROW *npp_sat_split_pack(NPP *npp, NPPROW *row, int nlit)
{     NPPROW *rrr;
      NPPCOL *col;
      NPPAIJ *aij;
      int k;
      /* original row should be packing inequality (1) */
      xassert(npp_sat_is_pack_ineq(npp, row) == 1);
      /* and nlit should be less than the number of literals (terms)
         in the original row */
      xassert(0 < nlit && nlit < npp_row_nnz(npp, row));
      /* create new row corresponding to inequality (2) */
      rrr = npp_add_row(npp);
      rrr->lb = -DBL_MAX, rrr->ub = 1.0;
      /* move first nlit literals (terms) from the original row to the
         new row; the original row becomes inequality (3) */
      for (k = 1; k <= nlit; k++)
      {  aij = row->ptr;
         xassert(aij != NULL);
         /* add literal to the new row */
         npp_add_aij(npp, rrr, aij->col, aij->val);
         /* correct rhs */
         if (aij->val < 0.0)
            rrr->ub -= 1.0, row->ub += 1.0;
         /* remove literal from the original row */
         npp_del_aij(npp, aij);
      }
      /* create new binary variable y */
      col = npp_add_col(npp);
      col->is_int = 1, col->lb = 0.0, col->ub = 1.0;
      /* include literal (1 - y) in the new row */
      npp_add_aij(npp, rrr, col, -1.0);
      rrr->ub -= 1.0;
      /* include literal y in the original row */
      npp_add_aij(npp, row, col, +1.0);
      return rrr;
}

/***********************************************************************
*  npp_sat_encode_pack - encode packing inequality
*
*  Given a packing inequality in canonical form:
*
*     sum  t[j] <= 1,                                                (1)
*   j in J
*
*  where t[j] = x[j] or t[j] = 1 - x[j], x[j] is a binary variable,
*  this routine translates it to CNF by replacing it with the following
*  equivalent set of edge packing inequalities:
*
*     t[j] + t[k] <= 1   for all j, k in J, j != k.                  (2)
*
*  Then the routine transforms each edge packing inequality (2) to
*  corresponding covering inequality (that encodes two-literal clause)
*  by multiplying both its part by -1:
*
*     - t[j] - t[k] >= -1   <-->   (1 - t[j]) + (1 - t[k]) >= 1.     (3)
*
*  On exit the routine removes the original row from the problem. */

void npp_sat_encode_pack(NPP *npp, NPPROW *row)
{     NPPROW *rrr;
      NPPAIJ *aij, *aik;
      /* original row should be packing inequality (1) */
      xassert(npp_sat_is_pack_ineq(npp, row) == 1);
      /* create equivalent system of covering inequalities (3) */
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  /* due to symmetry only one of inequalities t[j] + t[k] <= 1
            and t[k] <= t[j] <= 1 can be considered */
         for (aik = aij->r_next; aik != NULL; aik = aik->r_next)
         {  /* create edge packing inequality (2) */
            rrr = npp_add_row(npp);
            rrr->lb = -DBL_MAX, rrr->ub = 1.0;
            npp_add_aij(npp, rrr, aij->col, aij->val);
            if (aij->val < 0.0)
               rrr->ub -= 1.0;
            npp_add_aij(npp, rrr, aik->col, aik->val);
            if (aik->val < 0.0)
               rrr->ub -= 1.0;
            /* and transform it to covering inequality (3) */
            npp_sat_reverse_row(npp, rrr);
            xassert(npp_sat_is_cover_ineq(npp, rrr) == 1);
         }
      }
      /* remove the original row from the problem */
      npp_del_row(npp, row);
      return;
}

/***********************************************************************
*  npp_sat_encode_sum2 - encode 2-bit summation
*
*  Given a set containing two literals x and y this routine encodes
*  the equality
*
*     x + y = s + 2 * c,                                             (1)
*
*  where
*
*     s = (x + y) % 2                                                (2)
*
*  is a binary variable modeling the low sum bit, and
*
*     c = (x + y) / 2                                                (3)
*
*  is a binary variable modeling the high (carry) sum bit. */

void npp_sat_encode_sum2(NPP *npp, NPPLSE *set, NPPSED *sed)
{     NPPROW *row;
      int x, y, s, c;
      /* the set should contain exactly two literals */
      xassert(set != NULL);
      xassert(set->next != NULL);
      xassert(set->next->next == NULL);
      sed->x = set->lit;
      xassert(sed->x.neg == 0 || sed->x.neg == 1);
      sed->y = set->next->lit;
      xassert(sed->y.neg == 0 || sed->y.neg == 1);
      sed->z.col = NULL, sed->z.neg = 0;
      /* perform encoding s = (x + y) % 2 */
      sed->s = npp_add_col(npp);
      sed->s->is_int = 1, sed->s->lb = 0.0, sed->s->ub = 1.0;
      for (x = 0; x <= 1; x++)
      {  for (y = 0; y <= 1; y++)
         {  for (s = 0; s <= 1; s++)
            {  if ((x + y) % 2 != s)
               {  /* generate CNF clause to disable infeasible
                     combination */
                  row = npp_add_row(npp);
                  row->lb = 1.0, row->ub = +DBL_MAX;
                  if (x == sed->x.neg)
                     npp_add_aij(npp, row, sed->x.col, +1.0);
                  else
                  {  npp_add_aij(npp, row, sed->x.col, -1.0);
                     row->lb -= 1.0;
                  }
                  if (y == sed->y.neg)
                     npp_add_aij(npp, row, sed->y.col, +1.0);
                  else
                  {  npp_add_aij(npp, row, sed->y.col, -1.0);
                     row->lb -= 1.0;
                  }
                  if (s == 0)
                     npp_add_aij(npp, row, sed->s, +1.0);
                  else
                  {  npp_add_aij(npp, row, sed->s, -1.0);
                     row->lb -= 1.0;
                  }
               }
            }
         }
      }
      /* perform encoding c = (x + y) / 2 */
      sed->c = npp_add_col(npp);
      sed->c->is_int = 1, sed->c->lb = 0.0, sed->c->ub = 1.0;
      for (x = 0; x <= 1; x++)
      {  for (y = 0; y <= 1; y++)
         {  for (c = 0; c <= 1; c++)
            {  if ((x + y) / 2 != c)
               {  /* generate CNF clause to disable infeasible
                     combination */
                  row = npp_add_row(npp);
                  row->lb = 1.0, row->ub = +DBL_MAX;
                  if (x == sed->x.neg)
                     npp_add_aij(npp, row, sed->x.col, +1.0);
                  else
                  {  npp_add_aij(npp, row, sed->x.col, -1.0);
                     row->lb -= 1.0;
                  }
                  if (y == sed->y.neg)
                     npp_add_aij(npp, row, sed->y.col, +1.0);
                  else
                  {  npp_add_aij(npp, row, sed->y.col, -1.0);
                     row->lb -= 1.0;
                  }
                  if (c == 0)
                     npp_add_aij(npp, row, sed->c, +1.0);
                  else
                  {  npp_add_aij(npp, row, sed->c, -1.0);
                     row->lb -= 1.0;
                  }
               }
            }
         }
      }
      return;
}

/***********************************************************************
*  npp_sat_encode_sum3 - encode 3-bit summation
*
*  Given a set containing at least three literals this routine chooses
*  some literals x, y, z from that set and encodes the equality
*
*     x + y + z = s + 2 * c,                                         (1)
*
*  where
*
*     s = (x + y + z) % 2                                            (2)
*
*  is a binary variable modeling the low sum bit, and
*
*     c = (x + y + z) / 2                                            (3)
*
*  is a binary variable modeling the high (carry) sum bit. */

void npp_sat_encode_sum3(NPP *npp, NPPLSE *set, NPPSED *sed)
{     NPPROW *row;
      int x, y, z, s, c;
      /* the set should contain at least three literals */
      xassert(set != NULL);
      xassert(set->next != NULL);
      xassert(set->next->next != NULL);
      sed->x = set->lit;
      xassert(sed->x.neg == 0 || sed->x.neg == 1);
      sed->y = set->next->lit;
      xassert(sed->y.neg == 0 || sed->y.neg == 1);
      sed->z = set->next->next->lit;
      xassert(sed->z.neg == 0 || sed->z.neg == 1);
      /* perform encoding s = (x + y + z) % 2 */
      sed->s = npp_add_col(npp);
      sed->s->is_int = 1, sed->s->lb = 0.0, sed->s->ub = 1.0;
      for (x = 0; x <= 1; x++)
      {  for (y = 0; y <= 1; y++)
         {  for (z = 0; z <= 1; z++)
            {  for (s = 0; s <= 1; s++)
               {  if ((x + y + z) % 2 != s)
                  {  /* generate CNF clause to disable infeasible
                        combination */
                     row = npp_add_row(npp);
                     row->lb = 1.0, row->ub = +DBL_MAX;
                     if (x == sed->x.neg)
                        npp_add_aij(npp, row, sed->x.col, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->x.col, -1.0);
                        row->lb -= 1.0;
                     }
                     if (y == sed->y.neg)
                        npp_add_aij(npp, row, sed->y.col, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->y.col, -1.0);
                        row->lb -= 1.0;
                     }
                     if (z == sed->z.neg)
                        npp_add_aij(npp, row, sed->z.col, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->z.col, -1.0);
                        row->lb -= 1.0;
                     }
                     if (s == 0)
                        npp_add_aij(npp, row, sed->s, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->s, -1.0);
                        row->lb -= 1.0;
                     }
                  }
               }
            }
         }
      }
      /* perform encoding c = (x + y + z) / 2 */
      sed->c = npp_add_col(npp);
      sed->c->is_int = 1, sed->c->lb = 0.0, sed->c->ub = 1.0;
      for (x = 0; x <= 1; x++)
      {  for (y = 0; y <= 1; y++)
         {  for (z = 0; z <= 1; z++)
            {  for (c = 0; c <= 1; c++)
               {  if ((x + y + z) / 2 != c)
                  {  /* generate CNF clause to disable infeasible
                        combination */
                     row = npp_add_row(npp);
                     row->lb = 1.0, row->ub = +DBL_MAX;
                     if (x == sed->x.neg)
                        npp_add_aij(npp, row, sed->x.col, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->x.col, -1.0);
                        row->lb -= 1.0;
                     }
                     if (y == sed->y.neg)
                        npp_add_aij(npp, row, sed->y.col, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->y.col, -1.0);
                        row->lb -= 1.0;
                     }
                     if (z == sed->z.neg)
                        npp_add_aij(npp, row, sed->z.col, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->z.col, -1.0);
                        row->lb -= 1.0;
                     }
                     if (c == 0)
                        npp_add_aij(npp, row, sed->c, +1.0);
                     else
                     {  npp_add_aij(npp, row, sed->c, -1.0);
                        row->lb -= 1.0;
                     }
                  }
               }
            }
         }
      }
      return;
}

/***********************************************************************
*  npp_sat_encode_sum_ax - encode linear combination of 0-1 variables
*
*  PURPOSE
*
*  Given a linear combination of binary variables:
*
*     sum a[j] x[j],                                                 (1)
*      j
*
*  which is the linear form of the specified row, this routine encodes
*  (i.e. translates to CNF) the following equality:
*
*                        n
*     sum |a[j]| t[j] = sum 2**(k-1) * y[k],                         (2)
*      j                k=1
*
*  where t[j] = x[j] (if a[j] > 0) or t[j] = 1 - x[j] (if a[j] < 0),
*  and y[k] is either t[j] or a new literal created by the routine or
*  a constant zero. Note that the sum in the right-hand side of (2) can
*  be thought as a n-bit representation of the sum in the left-hand
*  side, which is a non-negative integer number.
*
*  ALGORITHM
*
*  First, the number of bits, n, sufficient to represent any value in
*  the left-hand side of (2) is determined. Obviously, n is the number
*  of bits sufficient to represent the sum (sum |a[j]|).
*
*  Let
*
*               n
*     |a[j]| = sum 2**(k-1) b[j,k],                                  (3)
*              k=1
*
*  where b[j,k] is k-th bit in a n-bit representation of |a[j]|. Then
*
*                          m            n
*     sum |a[j]| * t[j] = sum 2**(k-1) sum b[j,k] * t[j].            (4)
*      j                  k=1          j=1
*
*  Introducing the set
*
*     J[k] = { j : b[j,k] = 1 }                                      (5)
*
*  allows rewriting (4) as follows:
*
*                          n
*     sum |a[j]| * t[j] = sum 2**(k-1)  sum    t[j].                 (6)
*      j                  k=1         j in J[k]
*
*  Thus, our goal is to provide |J[k]| <= 1 for all k, in which case
*  we will have the representation (1).
*
*  Let |J[k]| = 2, i.e. J[k] has exactly two literals u and v. In this
*  case we can apply the following transformation:
*
*     u + v = s + 2 * c,                                             (7)
*
*  where s and c are, respectively, low (sum) and high (carry) bits of
*  the sum of two bits. This allows to replace two literals u and v in
*  J[k] by one literal s, and carry out literal c to J[k+1].
*
*  If |J[k]| >= 3, i.e. J[k] has at least three literals u, v, and w,
*  we can apply the following transformation:
*
*     u + v + w = s + 2 * c.                                         (8)
*
*  Again, literal s replaces literals u, v, and w in J[k], and literal
*  c goes into J[k+1].
*
*  On exit the routine stores each literal from J[k] in element y[k],
*  1 <= k <= n. If J[k] is empty, y[k] is set to constant false.
*
*  RETURNS
*
*  The routine returns n, the number of literals in the right-hand side
*  of (2), 0 <= n <= NBIT_MAX. If the sum (sum |a[j]|) is too large, so
*  more than NBIT_MAX (= 31) literals are needed to encode the original
*  linear combination, the routine returns a negative value. */

#define NBIT_MAX 31
/* maximal number of literals in the right hand-side of (2) */

static NPPLSE *remove_lse(NPP *npp, NPPLSE *set, NPPCOL *col)
{     /* remove specified literal from specified literal set */
      NPPLSE *lse, *prev = NULL;
      for (lse = set; lse != NULL; prev = lse, lse = lse->next)
         if (lse->lit.col == col) break;
      xassert(lse != NULL);
      if (prev == NULL)
         set = lse->next;
      else
         prev->next = lse->next;
      dmp_free_atom(npp->pool, lse, sizeof(NPPLSE));
      return set;
}

int npp_sat_encode_sum_ax(NPP *npp, NPPROW *row, NPPLIT y[])
{     NPPAIJ *aij;
      NPPLSE *set[1+NBIT_MAX], *lse;
      NPPSED sed;
      int k, n, temp;
      double sum;
      /* compute the sum (sum |a[j]|) */
      sum = 0.0;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
         sum += fabs(aij->val);
      /* determine n, the number of bits in the sum */
      temp = (int)sum;
      if ((double)temp != sum)
         return -1; /* integer arithmetic error */
      for (n = 0; temp > 0; n++, temp >>= 1);
      xassert(0 <= n && n <= NBIT_MAX);
      /* build initial sets J[k], 1 <= k <= n; see (5) */
      /* set[k] is a pointer to the list of literals in J[k] */
      for (k = 1; k <= n; k++)
         set[k] = NULL;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  temp = (int)fabs(aij->val);
         xassert((int)temp == fabs(aij->val));
         for (k = 1; temp > 0; k++, temp >>= 1)
         {  if (temp & 1)
            {  xassert(k <= n);
               lse = dmp_get_atom(npp->pool, sizeof(NPPLSE));
               lse->lit.col = aij->col;
               lse->lit.neg = (aij->val > 0.0 ? 0 : 1);
               lse->next = set[k];
               set[k] = lse;
            }
         }
      }
      /* main transformation loop */
      for (k = 1; k <= n; k++)
      {  /* reduce J[k] and set y[k] */
         for (;;)
         {  if (set[k] == NULL)
            {  /* J[k] is empty */
               /* set y[k] to constant false */
               y[k].col = NULL, y[k].neg = 0;
               break;
            }
            if (set[k]->next == NULL)
            {  /* J[k] contains one literal */
               /* set y[k] to that literal */
               y[k] = set[k]->lit;
               dmp_free_atom(npp->pool, set[k], sizeof(NPPLSE));
               break;
            }
            if (set[k]->next->next == NULL)
            {  /* J[k] contains two literals */
               /* apply transformation (7) */
               npp_sat_encode_sum2(npp, set[k], &sed);
            }
            else
            {  /* J[k] contains at least three literals */
               /* apply transformation (8) */
               npp_sat_encode_sum3(npp, set[k], &sed);
               /* remove third literal from set[k] */
               set[k] = remove_lse(npp, set[k], sed.z.col);
            }
            /* remove second literal from set[k] */
            set[k] = remove_lse(npp, set[k], sed.y.col);
            /* remove first literal from set[k] */
            set[k] = remove_lse(npp, set[k], sed.x.col);
            /* include new literal s to set[k] */
            lse = dmp_get_atom(npp->pool, sizeof(NPPLSE));
            lse->lit.col = sed.s, lse->lit.neg = 0;
            lse->next = set[k];
            set[k] = lse;
            /* include new literal c to set[k+1] */
            xassert(k < n); /* FIXME: can "overflow" happen? */
            lse = dmp_get_atom(npp->pool, sizeof(NPPLSE));
            lse->lit.col = sed.c, lse->lit.neg = 0;
            lse->next = set[k+1];
            set[k+1] = lse;
         }
      }
      return n;
}

/***********************************************************************
*  npp_sat_normalize_clause - normalize clause
*
*  This routine normalizes the specified clause, which is a disjunction
*  of literals, by replacing multiple literals, which refer to the same
*  binary variable, with a single literal.
*
*  On exit the routine returns the number of literals in the resulting
*  clause. However, if the specified clause includes both a literal and
*  its negation, the routine returns a negative value meaning that the
*  clause is equivalent to the value true. */

int npp_sat_normalize_clause(NPP *npp, int size, NPPLIT lit[])
{     int j, k, new_size;
      xassert(npp == npp);
      xassert(size >= 0);
      new_size = 0;
      for (k = 1; k <= size; k++)
      {  for (j = 1; j <= new_size; j++)
         {  if (lit[k].col == lit[j].col)
            {  /* lit[k] refers to the same variable as lit[j], which
                  is already included in the resulting clause */
               if (lit[k].neg == lit[j].neg)
               {  /* ignore lit[k] due to the idempotent law */
                  goto skip;
               }
               else
               {  /* lit[k] is NOT lit[j]; the clause is equivalent to
                     the value true */
                  return -1;
               }
            }
         }
         /* include lit[k] in the resulting clause */
         lit[++new_size] = lit[k];
skip:    ;
      }
      return new_size;
}

/***********************************************************************
*  npp_sat_encode_clause - translate clause to cover inequality
*
*  Given a clause
*
*     OR  t[j],                                                      (1)
*   j in J
*
*  where t[j] is a literal, i.e. t[j] = x[j] or t[j] = NOT x[j], this
*  routine translates it to the following equivalent cover inequality,
*  which is added to the transformed problem:
*
*     sum t[j] >= 1,                                                 (2)
*   j in J
*
*  where t[j] = x[j] or t[j] = 1 - x[j].
*
*  If necessary, the clause should be normalized before a call to this
*  routine. */

NPPROW *npp_sat_encode_clause(NPP *npp, int size, NPPLIT lit[])
{     NPPROW *row;
      int k;
      xassert(size >= 1);
      row = npp_add_row(npp);
      row->lb = 1.0, row->ub = +DBL_MAX;
      for (k = 1; k <= size; k++)
      {  xassert(lit[k].col != NULL);
         if (lit[k].neg == 0)
            npp_add_aij(npp, row, lit[k].col, +1.0);
         else if (lit[k].neg == 1)
         {  npp_add_aij(npp, row, lit[k].col, -1.0);
            row->lb -= 1.0;
         }
         else
            xassert(lit != lit);
      }
      return row;
}

/***********************************************************************
*  npp_sat_encode_geq - encode "not less than" constraint
*
*  PURPOSE
*
*  This routine translates to CNF the following constraint:
*
*      n
*     sum 2**(k-1) * y[k] >= b,                                      (1)
*     k=1
*
*  where y[k] is either a literal (i.e. y[k] = x[k] or y[k] = 1 - x[k])
*  or constant false (zero), b is a given lower bound.
*
*  ALGORITHM
*
*  If b < 0, the constraint is redundant, so assume that b >= 0. Let
*
*          n
*     b = sum 2**(k-1) b[k],                                         (2)
*         k=1
*
*  where b[k] is k-th binary digit of b. (Note that if b >= 2**n and
*  therefore cannot be represented in the form (2), the constraint (1)
*  is infeasible.) In this case the condition (1) is equivalent to the
*  following condition:
*
*     y[n] y[n-1] ... y[2] y[1] >= b[n] b[n-1] ... b[2] b[1],        (3)
*
*  where ">=" is understood lexicographically.
*
*  Algorithmically the condition (3) can be tested as follows:
*
*     for (k = n; k >= 1; k--)
*     {  if (y[k] < b[k])
*           y is less than b;
*        if (y[k] > b[k])
*           y is greater than b;
*     }
*     y is equal to b;
*
*  Thus, y is less than b iff there exists k, 1 <= k <= n, for which
*  the following condition is satisfied:
*
*     y[n] = b[n] AND ... AND y[k+1] = b[k+1] AND y[k] < b[k].       (4)
*
*  Negating the condition (4) we have that y is not less than b iff for
*  all k, 1 <= k <= n, the following condition is satisfied:
*
*     y[n] != b[n] OR ... OR y[k+1] != b[k+1] OR y[k] >= b[k].       (5)
*
*  Note that if b[k] = 0, the literal y[k] >= b[k] is always true, in
*  which case the entire clause (5) is true and can be omitted.
*
*  RETURNS
*
*  Normally the routine returns zero. However, if the constraint (1) is
*  infeasible, the routine returns non-zero. */

int npp_sat_encode_geq(NPP *npp, int n, NPPLIT y[], int rhs)
{     NPPLIT lit[1+NBIT_MAX];
      int j, k, size, temp, b[1+NBIT_MAX];
      xassert(0 <= n && n <= NBIT_MAX);
      /* if the constraint (1) is redundant, do nothing */
      if (rhs < 0)
         return 0;
      /* determine binary digits of b according to (2) */
      for (k = 1, temp = rhs; k <= n; k++, temp >>= 1)
         b[k] = temp & 1;
      if (temp != 0)
      {  /* b >= 2**n; the constraint (1) is infeasible */
         return 1;
      }
      /* main transformation loop */
      for (k = 1; k <= n; k++)
      {  /* build the clause (5) for current k */
         size = 0; /* clause size = number of literals */
         /* add literal y[k] >= b[k] */
         if (b[k] == 0)
         {  /* b[k] = 0 -> the literal is true */
            goto skip;
         }
         else if (y[k].col == NULL)
         {  /* y[k] = 0, b[k] = 1 -> the literal is false */
            xassert(y[k].neg == 0);
         }
         else
         {  /* add literal y[k] = 1 */
            lit[++size] = y[k];
         }
         for (j = k+1; j <= n; j++)
         {  /* add literal y[j] != b[j] */
            if (y[j].col == NULL)
            {  xassert(y[j].neg == 0);
               if (b[j] == 0)
               {  /* y[j] = 0, b[j] = 0 -> the literal is false */
                  continue;
               }
               else
               {  /* y[j] = 0, b[j] = 1 -> the literal is true */
                  goto skip;
               }
            }
            else
            {  lit[++size] = y[j];
               if (b[j] != 0)
                  lit[size].neg = 1 - lit[size].neg;
            }
         }
         /* normalize the clause */
         size = npp_sat_normalize_clause(npp, size, lit);
         if (size < 0)
         {  /* the clause is equivalent to the value true */
            goto skip;
         }
         if (size == 0)
         {  /* the clause is equivalent to the value false; this means
               that the constraint (1) is infeasible */
            return 2;
         }
         /* translate the clause to corresponding cover inequality */
         npp_sat_encode_clause(npp, size, lit);
skip:    ;
      }
      return 0;
}

/***********************************************************************
*  npp_sat_encode_leq - encode "not greater than" constraint
*
*  PURPOSE
*
*  This routine translates to CNF the following constraint:
*
*      n
*     sum 2**(k-1) * y[k] <= b,                                      (1)
*     k=1
*
*  where y[k] is either a literal (i.e. y[k] = x[k] or y[k] = 1 - x[k])
*  or constant false (zero), b is a given upper bound.
*
*  ALGORITHM
*
*  If b < 0, the constraint is infeasible, so assume that b >= 0. Let
*
*          n
*     b = sum 2**(k-1) b[k],                                         (2)
*         k=1
*
*  where b[k] is k-th binary digit of b. (Note that if b >= 2**n and
*  therefore cannot be represented in the form (2), the constraint (1)
*  is redundant.) In this case the condition (1) is equivalent to the
*  following condition:
*
*     y[n] y[n-1] ... y[2] y[1] <= b[n] b[n-1] ... b[2] b[1],        (3)
*
*  where "<=" is understood lexicographically.
*
*  Algorithmically the condition (3) can be tested as follows:
*
*     for (k = n; k >= 1; k--)
*     {  if (y[k] < b[k])
*           y is less than b;
*        if (y[k] > b[k])
*           y is greater than b;
*     }
*     y is equal to b;
*
*  Thus, y is greater than b iff there exists k, 1 <= k <= n, for which
*  the following condition is satisfied:
*
*     y[n] = b[n] AND ... AND y[k+1] = b[k+1] AND y[k] > b[k].       (4)
*
*  Negating the condition (4) we have that y is not greater than b iff
*  for all k, 1 <= k <= n, the following condition is satisfied:
*
*     y[n] != b[n] OR ... OR y[k+1] != b[k+1] OR y[k] <= b[k].       (5)
*
*  Note that if b[k] = 1, the literal y[k] <= b[k] is always true, in
*  which case the entire clause (5) is true and can be omitted.
*
*  RETURNS
*
*  Normally the routine returns zero. However, if the constraint (1) is
*  infeasible, the routine returns non-zero. */

int npp_sat_encode_leq(NPP *npp, int n, NPPLIT y[], int rhs)
{     NPPLIT lit[1+NBIT_MAX];
      int j, k, size, temp, b[1+NBIT_MAX];
      xassert(0 <= n && n <= NBIT_MAX);
      /* check if the constraint (1) is infeasible */
      if (rhs < 0)
         return 1;
      /* determine binary digits of b according to (2) */
      for (k = 1, temp = rhs; k <= n; k++, temp >>= 1)
         b[k] = temp & 1;
      if (temp != 0)
      {  /* b >= 2**n; the constraint (1) is redundant */
         return 0;
      }
      /* main transformation loop */
      for (k = 1; k <= n; k++)
      {  /* build the clause (5) for current k */
         size = 0; /* clause size = number of literals */
         /* add literal y[k] <= b[k] */
         if (b[k] == 1)
         {  /* b[k] = 1 -> the literal is true */
            goto skip;
         }
         else if (y[k].col == NULL)
         {  /* y[k] = 0, b[k] = 0 -> the literal is true */
            xassert(y[k].neg == 0);
            goto skip;
         }
         else
         {  /* add literal y[k] = 0 */
            lit[++size] = y[k];
            lit[size].neg = 1 - lit[size].neg;
         }
         for (j = k+1; j <= n; j++)
         {  /* add literal y[j] != b[j] */
            if (y[j].col == NULL)
            {  xassert(y[j].neg == 0);
               if (b[j] == 0)
               {  /* y[j] = 0, b[j] = 0 -> the literal is false */
                  continue;
               }
               else
               {  /* y[j] = 0, b[j] = 1 -> the literal is true */
                  goto skip;
               }
            }
            else
            {  lit[++size] = y[j];
               if (b[j] != 0)
                  lit[size].neg = 1 - lit[size].neg;
            }
         }
         /* normalize the clause */
         size = npp_sat_normalize_clause(npp, size, lit);
         if (size < 0)
         {  /* the clause is equivalent to the value true */
            goto skip;
         }
         if (size == 0)
         {  /* the clause is equivalent to the value false; this means
               that the constraint (1) is infeasible */
            return 2;
         }
         /* translate the clause to corresponding cover inequality */
         npp_sat_encode_clause(npp, size, lit);
skip:    ;
      }
      return 0;
}

/***********************************************************************
*  npp_sat_encode_row - encode constraint (row) of general type
*
*  PURPOSE
*
*  This routine translates to CNF the following constraint (row):
*
*     L <= sum a[j] x[j] <= U,                                       (1)
*           j
*
*  where all x[j] are binary variables.
*
*  ALGORITHM
*
*  First, the routine performs substitution x[j] = t[j] for j in J+
*  and x[j] = 1 - t[j] for j in J-, where J+ = { j : a[j] > 0 } and
*  J- = { j : a[j] < 0 }. This gives:
*
*     L <=  sum  a[j] t[j] +   sum  a[j] (1 - t[j]) <= U  ==>
*         j in J+            j in J-
*
*     L' <= sum |a[j]| t[j] <= U',                                   (2)
*            j
*
*  where
*
*     L' = L -   sum  a[j],   U' = U -   sum  a[j].                  (3)
*              j in J-                 j in J-
*
*  (Actually only new bounds L' and U' are computed.)
*
*  Then the routine translates to CNF the following equality:
*
*                        n
*     sum |a[j]| t[j] = sum 2**(k-1) * y[k],                         (4)
*      j                k=1
*
*  where y[k] is either some t[j] or a new literal or a constant zero
*  (see the routine npp_sat_encode_sum_ax).
*
*  Finally, the routine translates to CNF the following conditions:
*
*      n
*     sum 2**(k-1) * y[k] >= L'                                      (5)
*     k=1
*
*  and
*
*      n
*     sum 2**(k-1) * y[k] <= U'                                      (6)
*     k=1
*
*  (see the routines npp_sat_encode_geq and npp_sat_encode_leq).
*
*  All resulting clauses are encoded as cover inequalities and included
*  into the transformed problem.
*
*  Note that on exit the routine removes the specified constraint (row)
*  from the original problem.
*
*  RETURNS
*
*  The routine returns one of the following codes:
*
*  0 - translation was successful;
*  1 - constraint (1) was found infeasible;
*  2 - integer arithmetic error occured. */

int npp_sat_encode_row(NPP *npp, NPPROW *row)
{     NPPAIJ *aij;
      NPPLIT y[1+NBIT_MAX];
      int n, rhs;
      double lb, ub;
      /* the row should not be free */
      xassert(!(row->lb == -DBL_MAX && row->ub == +DBL_MAX));
      /* compute new bounds L' and U' (3) */
      lb = row->lb;
      ub = row->ub;
      for (aij = row->ptr; aij != NULL; aij = aij->r_next)
      {  if (aij->val < 0.0)
         {  if (lb != -DBL_MAX)
               lb -= aij->val;
            if (ub != -DBL_MAX)
               ub -= aij->val;
         }
      }
      /* encode the equality (4) */
      n = npp_sat_encode_sum_ax(npp, row, y);
      if (n < 0)
         return 2; /* integer arithmetic error */
      /* encode the condition (5) */
      if (lb != -DBL_MAX)
      {  rhs = (int)lb;
         if ((double)rhs != lb)
            return 2; /* integer arithmetic error */
         if (npp_sat_encode_geq(npp, n, y, rhs) != 0)
            return 1; /* original constraint is infeasible */
      }
      /* encode the condition (6) */
      if (ub != +DBL_MAX)
      {  rhs = (int)ub;
         if ((double)rhs != ub)
            return 2; /* integer arithmetic error */
         if (npp_sat_encode_leq(npp, n, y, rhs) != 0)
            return 1; /* original constraint is infeasible */
      }
      /* remove the specified row from the problem */
      npp_del_row(npp, row);
      return 0;
}

/***********************************************************************
*  npp_sat_encode_prob - encode 0-1 feasibility problem
*
*  This routine translates the specified 0-1 feasibility problem to an
*  equivalent SAT-CNF problem.
*
*  N.B. Currently this is a very crude implementation.
*
*  RETURNS
*
*  0           success;
*
*  GLP_ENOPFS  primal/integer infeasibility detected;
*
*  GLP_ERANGE  integer overflow occured. */

int npp_sat_encode_prob(NPP *npp)
{     NPPROW *row, *next_row, *prev_row;
      NPPCOL *col, *next_col;
      int cover = 0, pack = 0, partn = 0, ret;
      /* process and remove free rows */
      for (row = npp->r_head; row != NULL; row = next_row)
      {  next_row = row->next;
         if (row->lb == -DBL_MAX && row->ub == +DBL_MAX)
            npp_sat_free_row(npp, row);
      }
      /* process and remove fixed columns */
      for (col = npp->c_head; col != NULL; col = next_col)
      {  next_col = col->next;
         if (col->lb == col->ub)
            xassert(npp_sat_fixed_col(npp, col) == 0);
      }
      /* only binary variables should remain */
      for (col = npp->c_head; col != NULL; col = col->next)
         xassert(col->is_int && col->lb == 0.0 && col->ub == 1.0);
      /* new rows may be added to the end of the row list, so we walk
         from the end to beginning of the list */
      for (row = npp->r_tail; row != NULL; row = prev_row)
      {  prev_row = row->prev;
         /* process special cases */
         ret = npp_sat_is_cover_ineq(npp, row);
         if (ret != 0)
         {  /* row is covering inequality */
            cover++;
            /* since it already encodes a clause, just transform it to
               canonical form */
            if (ret == 2)
            {  xassert(npp_sat_reverse_row(npp, row) == 0);
               ret = npp_sat_is_cover_ineq(npp, row);
            }
            xassert(ret == 1);
            continue;
         }
         ret = npp_sat_is_partn_eq(npp, row);
         if (ret != 0)
         {  /* row is partitioning equality */
            NPPROW *cov;
            NPPAIJ *aij;
            partn++;
            /* transform it to canonical form */
            if (ret == 2)
            {  xassert(npp_sat_reverse_row(npp, row) == 0);
               ret = npp_sat_is_partn_eq(npp, row);
            }
            xassert(ret == 1);
            /* and split it into covering and packing inequalities,
               both in canonical forms */
            cov = npp_add_row(npp);
            cov->lb = row->lb, cov->ub = +DBL_MAX;
            for (aij = row->ptr; aij != NULL; aij = aij->r_next)
               npp_add_aij(npp, cov, aij->col, aij->val);
            xassert(npp_sat_is_cover_ineq(npp, cov) == 1);
            /* the cover inequality already encodes a clause and do
               not need any further processing */
            row->lb = -DBL_MAX;
            xassert(npp_sat_is_pack_ineq(npp, row) == 1);
            /* the packing inequality will be processed below */
            pack--;
         }
         ret = npp_sat_is_pack_ineq(npp, row);
         if (ret != 0)
         {  /* row is packing inequality */
            NPPROW *rrr;
            int nlit, desired_nlit = 4;
            pack++;
            /* transform it to canonical form */
            if (ret == 2)
            {  xassert(npp_sat_reverse_row(npp, row) == 0);
               ret = npp_sat_is_pack_ineq(npp, row);
            }
            xassert(ret == 1);
            /* process the packing inequality */
            for (;;)
            {  /* determine the number of literals in the remaining
                  inequality */
               nlit = npp_row_nnz(npp, row);
               if (nlit <= desired_nlit)
                  break;
               /* split the current inequality into one having not more
                  than desired_nlit literals and remaining one */
               rrr = npp_sat_split_pack(npp, row, desired_nlit-1);
               /* translate the former inequality to CNF and remove it
                  from the original problem */
               npp_sat_encode_pack(npp, rrr);
            }
            /* translate the remaining inequality to CNF and remove it
               from the original problem */
            npp_sat_encode_pack(npp, row);
            continue;
         }
         /* translate row of general type to CNF and remove it from the
            original problem */
         ret = npp_sat_encode_row(npp, row);
         if (ret == 0)
            ;
         else if (ret == 1)
            ret = GLP_ENOPFS;
         else if (ret == 2)
            ret = GLP_ERANGE;
         else
            xassert(ret != ret);
         if (ret != 0)
            goto done;
      }
      ret = 0;
      if (cover != 0)
         xprintf("%d covering inequalities\n", cover);
      if (pack != 0)
         xprintf("%d packing inequalities\n", pack);
      if (partn != 0)
         xprintf("%d partitioning equalities\n", partn);
done: return ret;
}

/* eof */
