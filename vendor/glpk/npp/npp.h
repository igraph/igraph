/* npp.h (LP/MIP preprocessor) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2009-2017 Free Software Foundation, Inc.
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

#ifndef NPP_H
#define NPP_H

#include "prob.h"

#if 0 /* 20/XI-2017 */
typedef struct NPP NPP;
#else
typedef struct glp_prep NPP;
#endif
typedef struct NPPROW NPPROW;
typedef struct NPPCOL NPPCOL;
typedef struct NPPAIJ NPPAIJ;
typedef struct NPPTSE NPPTSE;
typedef struct NPPLFE NPPLFE;

#if 0 /* 20/XI-2017 */
struct NPP
#else
struct glp_prep
#endif
{     /* LP/MIP preprocessor workspace */
      /*--------------------------------------------------------------*/
      /* original problem segment */
      int orig_dir;
      /* optimization direction flag:
         GLP_MIN - minimization
         GLP_MAX - maximization */
      int orig_m;
      /* number of rows */
      int orig_n;
      /* number of columns */
      int orig_nnz;
      /* number of non-zero constraint coefficients */
      /*--------------------------------------------------------------*/
      /* transformed problem segment (always minimization) */
      DMP *pool;
      /* memory pool to store problem components */
      char *name;
      /* problem name (1 to 255 chars); NULL means no name is assigned
         to the problem */
      char *obj;
      /* objective function name (1 to 255 chars); NULL means no name
         is assigned to the objective function */
      double c0;
      /* constant term of the objective function */
      int nrows;
      /* number of rows introduced into the problem; this count
         increases by one every time a new row is added and never
         decreases; thus, actual number of rows may be less than nrows
         due to row deletions */
      int ncols;
      /* number of columns introduced into the problem; this count
         increases by one every time a new column is added and never
         decreases; thus, actual number of column may be less than
         ncols due to column deletions */
      NPPROW *r_head;
      /* pointer to the beginning of the row list */
      NPPROW *r_tail;
      /* pointer to the end of the row list */
      NPPCOL *c_head;
      /* pointer to the beginning of the column list */
      NPPCOL *c_tail;
      /* pointer to the end of the column list */
      /*--------------------------------------------------------------*/
      /* transformation history */
      DMP *stack;
      /* memory pool to store transformation entries */
      NPPTSE *top;
      /* pointer to most recent transformation entry */
#if 0 /* 16/XII-2009 */
      int count[1+25];
      /* transformation statistics */
#endif
      /*--------------------------------------------------------------*/
      /* resultant (preprocessed) problem segment */
      int m;
      /* number of rows */
      int n;
      /* number of columns */
      int nnz;
      /* number of non-zero constraint coefficients */
      int *row_ref; /* int row_ref[1+m]; */
      /* row_ref[i], 1 <= i <= m, is the reference number assigned to
         a row, which is i-th row of the resultant problem */
      int *col_ref; /* int col_ref[1+n]; */
      /* col_ref[j], 1 <= j <= n, is the reference number assigned to
         a column, which is j-th column of the resultant problem */
      /*--------------------------------------------------------------*/
      /* recovered solution segment */
      int sol;
      /* solution indicator:
         GLP_SOL - basic solution
         GLP_IPT - interior-point solution
         GLP_MIP - mixed integer solution */
      int scaling;
      /* scaling option:
         GLP_OFF - scaling is disabled
         GLP_ON  - scaling is enabled */
      int p_stat;
      /* status of primal basic solution:
         GLP_UNDEF  - primal solution is undefined
         GLP_FEAS   - primal solution is feasible
         GLP_INFEAS - primal solution is infeasible
         GLP_NOFEAS - no primal feasible solution exists */
      int d_stat;
      /* status of dual basic solution:
         GLP_UNDEF  - dual solution is undefined
         GLP_FEAS   - dual solution is feasible
         GLP_INFEAS - dual solution is infeasible
         GLP_NOFEAS - no dual feasible solution exists */
      int t_stat;
      /* status of interior-point solution:
         GLP_UNDEF  - interior solution is undefined
         GLP_OPT    - interior solution is optimal */
      int i_stat;
      /* status of mixed integer solution:
         GLP_UNDEF  - integer solution is undefined
         GLP_OPT    - integer solution is optimal
         GLP_FEAS   - integer solution is feasible
         GLP_NOFEAS - no integer solution exists */
      char *r_stat; /* char r_stat[1+nrows]; */
      /* r_stat[i], 1 <= i <= nrows, is status of i-th row:
         GLP_BS - inactive constraint
         GLP_NL - active constraint on lower bound
         GLP_NU - active constraint on upper bound
         GLP_NF - active free row
         GLP_NS - active equality constraint */
      char *c_stat; /* char c_stat[1+nrows]; */
      /* c_stat[j], 1 <= j <= nrows, is status of j-th column:
         GLP_BS - basic variable
         GLP_NL - non-basic variable on lower bound
         GLP_NU - non-basic variable on upper bound
         GLP_NF - non-basic free variable
         GLP_NS - non-basic fixed variable */
      double *r_pi; /* double r_pi[1+nrows]; */
      /* r_pi[i], 1 <= i <= nrows, is Lagrange multiplier (dual value)
         for i-th row (constraint) */
      double *c_value; /* double c_value[1+ncols]; */
      /* c_value[j], 1 <= j <= ncols, is primal value of j-th column
         (structural variable) */
};

struct NPPROW
{     /* row (constraint) */
      int i;
      /* reference number assigned to the row, 1 <= i <= nrows */
      char *name;
      /* row name (1 to 255 chars); NULL means no name is assigned to
         the row */
      double lb;
      /* lower bound; -DBL_MAX means the row has no lower bound */
      double ub;
      /* upper bound; +DBL_MAX means the row has no upper bound */
      NPPAIJ *ptr;
      /* pointer to the linked list of constraint coefficients */
      int temp;
      /* working field used by preprocessor routines */
      NPPROW *prev;
      /* pointer to previous row in the row list */
      NPPROW *next;
      /* pointer to next row in the row list */
};

struct NPPCOL
{     /* column (variable) */
      int j;
      /* reference number assigned to the column, 1 <= j <= ncols */
      char *name;
      /* column name (1 to 255 chars); NULL means no name is assigned
         to the column */
      char is_int;
      /* 0 means continuous variable; 1 means integer variable */
      double lb;
      /* lower bound; -DBL_MAX means the column has no lower bound */
      double ub;
      /* upper bound; +DBL_MAX means the column has no upper bound */
      double coef;
      /* objective coefficient */
      NPPAIJ *ptr;
      /* pointer to the linked list of constraint coefficients */
      int temp;
      /* working field used by preprocessor routines */
#if 1 /* 28/XII-2009 */
      union
      {  double ll;
         /* implied column lower bound */
         int pos;
         /* vertex ordinal number corresponding to this binary column
            in the conflict graph (0, if the vertex does not exist) */
      }  ll;
      union
      {  double uu;
         /* implied column upper bound */
         int neg;
         /* vertex ordinal number corresponding to complement of this
            binary column in the conflict graph (0, if the vertex does
            not exist) */
      }  uu;
#endif
      NPPCOL *prev;
      /* pointer to previous column in the column list */
      NPPCOL *next;
      /* pointer to next column in the column list */
};

struct NPPAIJ
{     /* constraint coefficient */
      NPPROW *row;
      /* pointer to corresponding row */
      NPPCOL *col;
      /* pointer to corresponding column */
      double val;
      /* (non-zero) coefficient value */
      NPPAIJ *r_prev;
      /* pointer to previous coefficient in the same row */
      NPPAIJ *r_next;
      /* pointer to next coefficient in the same row */
      NPPAIJ *c_prev;
      /* pointer to previous coefficient in the same column */
      NPPAIJ *c_next;
      /* pointer to next coefficient in the same column */
};

struct NPPTSE
{     /* transformation stack entry */
      int (*func)(NPP *npp, void *info);
      /* pointer to routine performing back transformation */
      void *info;
      /* pointer to specific info (depends on the transformation) */
      NPPTSE *link;
      /* pointer to another entry created *before* this entry */
};

struct NPPLFE
{     /* linear form element */
      int ref;
      /* row/column reference number */
      double val;
      /* (non-zero) coefficient value */
      NPPLFE *next;
      /* pointer to another element */
};

#define npp_create_wksp _glp_npp_create_wksp
NPP *npp_create_wksp(void);
/* create LP/MIP preprocessor workspace */

#define npp_insert_row _glp_npp_insert_row
void npp_insert_row(NPP *npp, NPPROW *row, int where);
/* insert row to the row list */

#define npp_remove_row _glp_npp_remove_row
void npp_remove_row(NPP *npp, NPPROW *row);
/* remove row from the row list */

#define npp_activate_row _glp_npp_activate_row
void npp_activate_row(NPP *npp, NPPROW *row);
/* make row active */

#define npp_deactivate_row _glp_npp_deactivate_row
void npp_deactivate_row(NPP *npp, NPPROW *row);
/* make row inactive */

#define npp_insert_col _glp_npp_insert_col
void npp_insert_col(NPP *npp, NPPCOL *col, int where);
/* insert column to the column list */

#define npp_remove_col _glp_npp_remove_col
void npp_remove_col(NPP *npp, NPPCOL *col);
/* remove column from the column list */

#define npp_activate_col _glp_npp_activate_col
void npp_activate_col(NPP *npp, NPPCOL *col);
/* make column active */

#define npp_deactivate_col _glp_npp_deactivate_col
void npp_deactivate_col(NPP *npp, NPPCOL *col);
/* make column inactive */

#define npp_add_row _glp_npp_add_row
NPPROW *npp_add_row(NPP *npp);
/* add new row to the current problem */

#define npp_add_col _glp_npp_add_col
NPPCOL *npp_add_col(NPP *npp);
/* add new column to the current problem */

#define npp_add_aij _glp_npp_add_aij
NPPAIJ *npp_add_aij(NPP *npp, NPPROW *row, NPPCOL *col, double val);
/* add new element to the constraint matrix */

#define npp_row_nnz _glp_npp_row_nnz
int npp_row_nnz(NPP *npp, NPPROW *row);
/* count number of non-zero coefficients in row */

#define npp_col_nnz _glp_npp_col_nnz
int npp_col_nnz(NPP *npp, NPPCOL *col);
/* count number of non-zero coefficients in column */

#define npp_push_tse _glp_npp_push_tse
void *npp_push_tse(NPP *npp, int (*func)(NPP *npp, void *info),
      int size);
/* push new entry to the transformation stack */

#define npp_erase_row _glp_npp_erase_row
void npp_erase_row(NPP *npp, NPPROW *row);
/* erase row content to make it empty */

#define npp_del_row _glp_npp_del_row
void npp_del_row(NPP *npp, NPPROW *row);
/* remove row from the current problem */

#define npp_del_col _glp_npp_del_col
void npp_del_col(NPP *npp, NPPCOL *col);
/* remove column from the current problem */

#define npp_del_aij _glp_npp_del_aij
void npp_del_aij(NPP *npp, NPPAIJ *aij);
/* remove element from the constraint matrix */

#define npp_load_prob _glp_npp_load_prob
void npp_load_prob(NPP *npp, glp_prob *orig, int names, int sol,
      int scaling);
/* load original problem into the preprocessor workspace */

#define npp_build_prob _glp_npp_build_prob
void npp_build_prob(NPP *npp, glp_prob *prob);
/* build resultant (preprocessed) problem */

#define npp_postprocess _glp_npp_postprocess
void npp_postprocess(NPP *npp, glp_prob *prob);
/* postprocess solution from the resultant problem */

#define npp_unload_sol _glp_npp_unload_sol
void npp_unload_sol(NPP *npp, glp_prob *orig);
/* store solution to the original problem */

#define npp_delete_wksp _glp_npp_delete_wksp
void npp_delete_wksp(NPP *npp);
/* delete LP/MIP preprocessor workspace */

#define npp_error()

#define npp_free_row _glp_npp_free_row
void npp_free_row(NPP *npp, NPPROW *p);
/* process free (unbounded) row */

#define npp_geq_row _glp_npp_geq_row
void npp_geq_row(NPP *npp, NPPROW *p);
/* process row of 'not less than' type */

#define npp_leq_row _glp_npp_leq_row
void npp_leq_row(NPP *npp, NPPROW *p);
/* process row of 'not greater than' type */

#define npp_free_col _glp_npp_free_col
void npp_free_col(NPP *npp, NPPCOL *q);
/* process free (unbounded) column */

#define npp_lbnd_col _glp_npp_lbnd_col
void npp_lbnd_col(NPP *npp, NPPCOL *q);
/* process column with (non-zero) lower bound */

#define npp_ubnd_col _glp_npp_ubnd_col
void npp_ubnd_col(NPP *npp, NPPCOL *q);
/* process column with upper bound */

#define npp_dbnd_col _glp_npp_dbnd_col
void npp_dbnd_col(NPP *npp, NPPCOL *q);
/* process non-negative column with upper bound */

#define npp_fixed_col _glp_npp_fixed_col
void npp_fixed_col(NPP *npp, NPPCOL *q);
/* process fixed column */

#define npp_make_equality _glp_npp_make_equality
int npp_make_equality(NPP *npp, NPPROW *p);
/* process row with almost identical bounds */

#define npp_make_fixed _glp_npp_make_fixed
int npp_make_fixed(NPP *npp, NPPCOL *q);
/* process column with almost identical bounds */

#define npp_empty_row _glp_npp_empty_row
int npp_empty_row(NPP *npp, NPPROW *p);
/* process empty row */

#define npp_empty_col _glp_npp_empty_col
int npp_empty_col(NPP *npp, NPPCOL *q);
/* process empty column */

#define npp_implied_value _glp_npp_implied_value
int npp_implied_value(NPP *npp, NPPCOL *q, double s);
/* process implied column value */

#define npp_eq_singlet _glp_npp_eq_singlet
int npp_eq_singlet(NPP *npp, NPPROW *p);
/* process row singleton (equality constraint) */

#define npp_implied_lower _glp_npp_implied_lower
int npp_implied_lower(NPP *npp, NPPCOL *q, double l);
/* process implied column lower bound */

#define npp_implied_upper _glp_npp_implied_upper
int npp_implied_upper(NPP *npp, NPPCOL *q, double u);
/* process implied upper bound of column */

#define npp_ineq_singlet _glp_npp_ineq_singlet
int npp_ineq_singlet(NPP *npp, NPPROW *p);
/* process row singleton (inequality constraint) */

#define npp_implied_slack _glp_npp_implied_slack
void npp_implied_slack(NPP *npp, NPPCOL *q);
/* process column singleton (implied slack variable) */

#define npp_implied_free _glp_npp_implied_free
int npp_implied_free(NPP *npp, NPPCOL *q);
/* process column singleton (implied free variable) */

#define npp_eq_doublet _glp_npp_eq_doublet
NPPCOL *npp_eq_doublet(NPP *npp, NPPROW *p);
/* process row doubleton (equality constraint) */

#define npp_forcing_row _glp_npp_forcing_row
int npp_forcing_row(NPP *npp, NPPROW *p, int at);
/* process forcing row */

#define npp_analyze_row _glp_npp_analyze_row
int npp_analyze_row(NPP *npp, NPPROW *p);
/* perform general row analysis */

#define npp_inactive_bound _glp_npp_inactive_bound
void npp_inactive_bound(NPP *npp, NPPROW *p, int which);
/* remove row lower/upper inactive bound */

#define npp_implied_bounds _glp_npp_implied_bounds
void npp_implied_bounds(NPP *npp, NPPROW *p);
/* determine implied column bounds */

#define npp_binarize_prob _glp_npp_binarize_prob
int npp_binarize_prob(NPP *npp);
/* binarize MIP problem */

#define npp_is_packing _glp_npp_is_packing
int npp_is_packing(NPP *npp, NPPROW *row);
/* test if constraint is packing inequality */

#define npp_hidden_packing _glp_npp_hidden_packing
int npp_hidden_packing(NPP *npp, NPPROW *row);
/* identify hidden packing inequality */

#define npp_implied_packing _glp_npp_implied_packing
int npp_implied_packing(NPP *npp, NPPROW *row, int which,
      NPPCOL *var[], char set[]);
/* identify implied packing inequality */

#define npp_is_covering _glp_npp_is_covering
int npp_is_covering(NPP *npp, NPPROW *row);
/* test if constraint is covering inequality */

#define npp_hidden_covering _glp_npp_hidden_covering
int npp_hidden_covering(NPP *npp, NPPROW *row);
/* identify hidden covering inequality */

#define npp_is_partitioning _glp_npp_is_partitioning
int npp_is_partitioning(NPP *npp, NPPROW *row);
/* test if constraint is partitioning equality */

#define npp_reduce_ineq_coef _glp_npp_reduce_ineq_coef
int npp_reduce_ineq_coef(NPP *npp, NPPROW *row);
/* reduce inequality constraint coefficients */

#define npp_clean_prob _glp_npp_clean_prob
void npp_clean_prob(NPP *npp);
/* perform initial LP/MIP processing */

#define npp_process_row _glp_npp_process_row
int npp_process_row(NPP *npp, NPPROW *row, int hard);
/* perform basic row processing */

#define npp_improve_bounds _glp_npp_improve_bounds
int npp_improve_bounds(NPP *npp, NPPROW *row, int flag);
/* improve current column bounds */

#define npp_process_col _glp_npp_process_col
int npp_process_col(NPP *npp, NPPCOL *col);
/* perform basic column processing */

#define npp_process_prob _glp_npp_process_prob
int npp_process_prob(NPP *npp, int hard);
/* perform basic LP/MIP processing */

#define npp_simplex _glp_npp_simplex
int npp_simplex(NPP *npp, const glp_smcp *parm);
/* process LP prior to applying primal/dual simplex method */

#define npp_integer _glp_npp_integer
int npp_integer(NPP *npp, const glp_iocp *parm);
/* process MIP prior to applying branch-and-bound method */

/**********************************************************************/

#define npp_sat_free_row _glp_npp_sat_free_row
void npp_sat_free_row(NPP *npp, NPPROW *p);
/* process free (unbounded) row */

#define npp_sat_fixed_col _glp_npp_sat_fixed_col
int npp_sat_fixed_col(NPP *npp, NPPCOL *q);
/* process fixed column */

#define npp_sat_is_bin_comb _glp_npp_sat_is_bin_comb
int npp_sat_is_bin_comb(NPP *npp, NPPROW *row);
/* test if row is binary combination */

#define npp_sat_num_pos_coef _glp_npp_sat_num_pos_coef
int npp_sat_num_pos_coef(NPP *npp, NPPROW *row);
/* determine number of positive coefficients */

#define npp_sat_num_neg_coef _glp_npp_sat_num_neg_coef
int npp_sat_num_neg_coef(NPP *npp, NPPROW *row);
/* determine number of negative coefficients */

#define npp_sat_is_cover_ineq _glp_npp_sat_is_cover_ineq
int npp_sat_is_cover_ineq(NPP *npp, NPPROW *row);
/* test if row is covering inequality */

#define npp_sat_is_pack_ineq _glp_npp_sat_is_pack_ineq
int npp_sat_is_pack_ineq(NPP *npp, NPPROW *row);
/* test if row is packing inequality */

#define npp_sat_is_partn_eq _glp_npp_sat_is_partn_eq
int npp_sat_is_partn_eq(NPP *npp, NPPROW *row);
/* test if row is partitioning equality */

#define npp_sat_reverse_row _glp_npp_sat_reverse_row
int npp_sat_reverse_row(NPP *npp, NPPROW *row);
/* multiply both sides of row by -1 */

#define npp_sat_split_pack _glp_npp_sat_split_pack
NPPROW *npp_sat_split_pack(NPP *npp, NPPROW *row, int nnn);
/* split packing inequality */

#define npp_sat_encode_pack _glp_npp_sat_encode_pack
void npp_sat_encode_pack(NPP *npp, NPPROW *row);
/* encode packing inequality */

typedef struct NPPLIT NPPLIT;
typedef struct NPPLSE NPPLSE;
typedef struct NPPSED NPPSED;

struct NPPLIT
{     /* literal (binary variable or its negation) */
      NPPCOL *col;
      /* pointer to binary variable; NULL means constant false */
      int neg;
      /* negation flag:
         0 - literal is variable (or constant false)
         1 - literal is negation of variable (or constant true) */
};

struct NPPLSE
{     /* literal set element */
      NPPLIT lit;
      /* literal */
      NPPLSE *next;
      /* pointer to another element */
};

struct NPPSED
{     /* summation encoding descriptor */
      /* this struct describes the equality
            x + y + z = s + 2 * c,
         which was encoded as CNF and included into the transformed
         problem; here x and y are literals, z is either a literal or
         constant zero, s and c are binary variables modeling, resp.,
         the low and high (carry) sum bits */
      NPPLIT x, y, z;
      /* literals; if z.col = NULL, z is constant zero */
      NPPCOL *s, *c;
      /* binary variables modeling the sum bits */
};

#define npp_sat_encode_sum2 _glp_npp_sat_encode_sum2
void npp_sat_encode_sum2(NPP *npp, NPPLSE *set, NPPSED *sed);
/* encode 2-bit summation */

#define npp_sat_encode_sum3 _glp_npp_sat_encode_sum3
void npp_sat_encode_sum3(NPP *npp, NPPLSE *set, NPPSED *sed);
/* encode 3-bit summation */

#define npp_sat_encode_sum_ax _glp_npp_sat_encode_sum_ax
int npp_sat_encode_sum_ax(NPP *npp, NPPROW *row, NPPLIT y[]);
/* encode linear combination of 0-1 variables */

#define npp_sat_normalize_clause _glp_npp_sat_normalize_clause
int npp_sat_normalize_clause(NPP *npp, int size, NPPLIT lit[]);
/* normalize clause */

#define npp_sat_encode_clause _glp_npp_sat_encode_clause
NPPROW *npp_sat_encode_clause(NPP *npp, int size, NPPLIT lit[]);
/* translate clause to cover inequality */

#define npp_sat_encode_geq _glp_npp_sat_encode_geq
int npp_sat_encode_geq(NPP *npp, int n, NPPLIT y[], int rhs);
/* encode "not less than" constraint */

#define npp_sat_encode_leq _glp_npp_sat_encode_leq
int npp_sat_encode_leq(NPP *npp, int n, NPPLIT y[], int rhs);
/* encode "not greater than" constraint */

#define npp_sat_encode_row _glp_npp_sat_encode_row
int npp_sat_encode_row(NPP *npp, NPPROW *row);
/* encode constraint (row) of general type */

#define npp_sat_encode_prob _glp_npp_sat_encode_prob
int npp_sat_encode_prob(NPP *npp);
/* encode 0-1 feasibility problem */

#endif

/* eof */
