/* glpk.h */

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

#ifndef GLPK_H
#define GLPK_H

#include <stdarg.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* library version numbers: */
#define GLP_MAJOR_VERSION  4
#define GLP_MINOR_VERSION  45

#ifndef GLP_PROB_DEFINED
#define GLP_PROB_DEFINED
typedef struct glp_prob glp_prob;
/* LP/MIP problem object */
#endif

/* optimization direction flag: */
#define GLP_MIN            1  /* minimization */
#define GLP_MAX            2  /* maximization */

/* kind of structural variable: */
#define GLP_CV             1  /* continuous variable */
#define GLP_IV             2  /* integer variable */
#define GLP_BV             3  /* binary variable */

/* type of auxiliary/structural variable: */
#define GLP_FR             1  /* free variable */
#define GLP_LO             2  /* variable with lower bound */
#define GLP_UP             3  /* variable with upper bound */
#define GLP_DB             4  /* double-bounded variable */
#define GLP_FX             5  /* fixed variable */

/* status of auxiliary/structural variable: */
#define GLP_BS             1  /* basic variable */
#define GLP_NL             2  /* non-basic variable on lower bound */
#define GLP_NU             3  /* non-basic variable on upper bound */
#define GLP_NF             4  /* non-basic free variable */
#define GLP_NS             5  /* non-basic fixed variable */

/* scaling options: */
#define GLP_SF_GM       0x01  /* perform geometric mean scaling */
#define GLP_SF_EQ       0x10  /* perform equilibration scaling */
#define GLP_SF_2N       0x20  /* round scale factors to power of two */
#define GLP_SF_SKIP     0x40  /* skip if problem is well scaled */
#define GLP_SF_AUTO     0x80  /* choose scaling options automatically */

/* solution indicator: */
#define GLP_SOL            1  /* basic solution */
#define GLP_IPT            2  /* interior-point solution */
#define GLP_MIP            3  /* mixed integer solution */

/* solution status: */
#define GLP_UNDEF          1  /* solution is undefined */
#define GLP_FEAS           2  /* solution is feasible */
#define GLP_INFEAS         3  /* solution is infeasible */
#define GLP_NOFEAS         4  /* no feasible solution exists */
#define GLP_OPT            5  /* solution is optimal */
#define GLP_UNBND          6  /* solution is unbounded */

typedef struct
{     /* basis factorization control parameters */
      int msg_lev;            /* (reserved) */
      int type;               /* factorization type: */
#define GLP_BF_FT          1  /* LUF + Forrest-Tomlin */
#define GLP_BF_BG          2  /* LUF + Schur compl. + Bartels-Golub */
#define GLP_BF_GR          3  /* LUF + Schur compl. + Givens rotation */
      int lu_size;            /* luf.sv_size */
      double piv_tol;         /* luf.piv_tol */
      int piv_lim;            /* luf.piv_lim */
      int suhl;               /* luf.suhl */
      double eps_tol;         /* luf.eps_tol */
      double max_gro;         /* luf.max_gro */
      int nfs_max;            /* fhv.hh_max */
      double upd_tol;         /* fhv.upd_tol */
      int nrs_max;            /* lpf.n_max */
      int rs_size;            /* lpf.v_size */
      double foo_bar[38];     /* (reserved) */
} glp_bfcp;

typedef struct
{     /* simplex method control parameters */
      int msg_lev;            /* message level: */
#define GLP_MSG_OFF        0  /* no output */
#define GLP_MSG_ERR        1  /* warning and error messages only */
#define GLP_MSG_ON         2  /* normal output */
#define GLP_MSG_ALL        3  /* full output */
#define GLP_MSG_DBG        4  /* debug output */
      int meth;               /* simplex method option: */
#define GLP_PRIMAL         1  /* use primal simplex */
#define GLP_DUALP          2  /* use dual; if it fails, use primal */
#define GLP_DUAL           3  /* use dual simplex */
      int pricing;            /* pricing technique: */
#define GLP_PT_STD      0x11  /* standard (Dantzig rule) */
#define GLP_PT_PSE      0x22  /* projected steepest edge */
      int r_test;             /* ratio test technique: */
#define GLP_RT_STD      0x11  /* standard (textbook) */
#define GLP_RT_HAR      0x22  /* two-pass Harris' ratio test */
      double tol_bnd;         /* spx.tol_bnd */
      double tol_dj;          /* spx.tol_dj */
      double tol_piv;         /* spx.tol_piv */
      double obj_ll;          /* spx.obj_ll */
      double obj_ul;          /* spx.obj_ul */
      int it_lim;             /* spx.it_lim */
      int tm_lim;             /* spx.tm_lim (milliseconds) */
      int out_frq;            /* spx.out_frq */
      int out_dly;            /* spx.out_dly (milliseconds) */
      int presolve;           /* enable/disable using LP presolver */
      double foo_bar[36];     /* (reserved) */
} glp_smcp;

typedef struct
{     /* interior-point solver control parameters */
      int msg_lev;            /* message level (see glp_smcp) */
      int ord_alg;            /* ordering algorithm: */
#define GLP_ORD_NONE       0  /* natural (original) ordering */
#define GLP_ORD_QMD        1  /* quotient minimum degree (QMD) */
#define GLP_ORD_AMD        2  /* approx. minimum degree (AMD) */
#define GLP_ORD_SYMAMD     3  /* approx. minimum degree (SYMAMD) */
      double foo_bar[48];     /* (reserved) */
} glp_iptcp;

#ifndef GLP_TREE_DEFINED
#define GLP_TREE_DEFINED
typedef struct glp_tree glp_tree;
/* branch-and-bound tree */
#endif

typedef struct
{     /* integer optimizer control parameters */
      int msg_lev;            /* message level (see glp_smcp) */
      int br_tech;            /* branching technique: */
#define GLP_BR_FFV         1  /* first fractional variable */
#define GLP_BR_LFV         2  /* last fractional variable */
#define GLP_BR_MFV         3  /* most fractional variable */
#define GLP_BR_DTH         4  /* heuristic by Driebeck and Tomlin */
#define GLP_BR_PCH         5  /* hybrid pseudocost heuristic */
      int bt_tech;            /* backtracking technique: */
#define GLP_BT_DFS         1  /* depth first search */
#define GLP_BT_BFS         2  /* breadth first search */
#define GLP_BT_BLB         3  /* best local bound */
#define GLP_BT_BPH         4  /* best projection heuristic */
      double tol_int;         /* mip.tol_int */
      double tol_obj;         /* mip.tol_obj */
      int tm_lim;             /* mip.tm_lim (milliseconds) */
      int out_frq;            /* mip.out_frq (milliseconds) */
      int out_dly;            /* mip.out_dly (milliseconds) */
      void (*cb_func)(glp_tree *T, void *info);
                              /* mip.cb_func */
      void *cb_info;          /* mip.cb_info */
      int cb_size;            /* mip.cb_size */
      int pp_tech;            /* preprocessing technique: */
#define GLP_PP_NONE        0  /* disable preprocessing */
#define GLP_PP_ROOT        1  /* preprocessing only on root level */
#define GLP_PP_ALL         2  /* preprocessing on all levels */
      double mip_gap;         /* relative MIP gap tolerance */
      int mir_cuts;           /* MIR cuts       (GLP_ON/GLP_OFF) */
      int gmi_cuts;           /* Gomory's cuts  (GLP_ON/GLP_OFF) */
      int cov_cuts;           /* cover cuts     (GLP_ON/GLP_OFF) */
      int clq_cuts;           /* clique cuts    (GLP_ON/GLP_OFF) */
      int presolve;           /* enable/disable using MIP presolver */
      int binarize;           /* try to binarize integer variables */
      int fp_heur;            /* feasibility pump heuristic */
#if 1 /* 28/V-2010 */
      int alien;              /* use alien solver */
#endif
      double foo_bar[29];     /* (reserved) */
} glp_iocp;

typedef struct
{     /* additional row attributes */
      int level;
      /* subproblem level at which the row was added */
      int origin;
      /* row origin flag: */
#define GLP_RF_REG         0  /* regular constraint */
#define GLP_RF_LAZY        1  /* "lazy" constraint */
#define GLP_RF_CUT         2  /* cutting plane constraint */
      int klass;
      /* row class descriptor: */
#define GLP_RF_GMI         1  /* Gomory's mixed integer cut */
#define GLP_RF_MIR         2  /* mixed integer rounding cut */
#define GLP_RF_COV         3  /* mixed cover cut */
#define GLP_RF_CLQ         4  /* clique cut */
      double foo_bar[7];
      /* (reserved) */
} glp_attr;

/* enable/disable flag: */
#define GLP_ON             1  /* enable something */
#define GLP_OFF            0  /* disable something */

/* reason codes: */
#define GLP_IROWGEN     0x01  /* request for row generation */
#define GLP_IBINGO      0x02  /* better integer solution found */
#define GLP_IHEUR       0x03  /* request for heuristic solution */
#define GLP_ICUTGEN     0x04  /* request for cut generation */
#define GLP_IBRANCH     0x05  /* request for branching */
#define GLP_ISELECT     0x06  /* request for subproblem selection */
#define GLP_IPREPRO     0x07  /* request for preprocessing */

/* branch selection indicator: */
#define GLP_NO_BRNCH       0  /* select no branch */
#define GLP_DN_BRNCH       1  /* select down-branch */
#define GLP_UP_BRNCH       2  /* select up-branch */

/* return codes: */
#define GLP_EBADB       0x01  /* invalid basis */
#define GLP_ESING       0x02  /* singular matrix */
#define GLP_ECOND       0x03  /* ill-conditioned matrix */
#define GLP_EBOUND      0x04  /* invalid bounds */
#define GLP_EFAIL       0x05  /* solver failed */
#define GLP_EOBJLL      0x06  /* objective lower limit reached */
#define GLP_EOBJUL      0x07  /* objective upper limit reached */
#define GLP_EITLIM      0x08  /* iteration limit exceeded */
#define GLP_ETMLIM      0x09  /* time limit exceeded */
#define GLP_ENOPFS      0x0A  /* no primal feasible solution */
#define GLP_ENODFS      0x0B  /* no dual feasible solution */
#define GLP_EROOT       0x0C  /* root LP optimum not provided */
#define GLP_ESTOP       0x0D  /* search terminated by application */
#define GLP_EMIPGAP     0x0E  /* relative mip gap tolerance reached */
#define GLP_ENOFEAS     0x0F  /* no primal/dual feasible solution */
#define GLP_ENOCVG      0x10  /* no convergence */
#define GLP_EINSTAB     0x11  /* numerical instability */
#define GLP_EDATA       0x12  /* invalid data */
#define GLP_ERANGE      0x13  /* result out of range */

/* condition indicator: */
#define GLP_KKT_PE         1  /* primal equalities */
#define GLP_KKT_PB         2  /* primal bounds */
#define GLP_KKT_DE         3  /* dual equalities */
#define GLP_KKT_DB         4  /* dual bounds */
#define GLP_KKT_CS         5  /* complementary slackness */

/* MPS file format: */
#define GLP_MPS_DECK       1  /* fixed (ancient) */
#define GLP_MPS_FILE       2  /* free (modern) */

typedef struct
{     /* MPS format control parameters */
      int blank;
      /* character code to replace blanks in symbolic names */
      char *obj_name;
      /* objective row name */
      double tol_mps;
      /* zero tolerance for MPS data */
      double foo_bar[17];
      /* (reserved for use in the future) */
} glp_mpscp;

typedef struct
{     /* CPLEX LP format control parameters */
      double foo_bar[20];
      /* (reserved for use in the future) */
} glp_cpxcp;

#ifndef GLP_TRAN_DEFINED
#define GLP_TRAN_DEFINED
typedef struct glp_tran glp_tran;
/* MathProg translator workspace */
#endif

glp_prob *glp_create_prob(void);
/* create problem object */

void glp_set_prob_name(glp_prob *P, const char *name);
/* assign (change) problem name */

void glp_set_obj_name(glp_prob *P, const char *name);
/* assign (change) objective function name */

void glp_set_obj_dir(glp_prob *P, int dir);
/* set (change) optimization direction flag */

int glp_add_rows(glp_prob *P, int nrs);
/* add new rows to problem object */

int glp_add_cols(glp_prob *P, int ncs);
/* add new columns to problem object */

void glp_set_row_name(glp_prob *P, int i, const char *name);
/* assign (change) row name */

void glp_set_col_name(glp_prob *P, int j, const char *name);
/* assign (change) column name */

void glp_set_row_bnds(glp_prob *P, int i, int type, double lb,
      double ub);
/* set (change) row bounds */

void glp_set_col_bnds(glp_prob *P, int j, int type, double lb,
      double ub);
/* set (change) column bounds */

void glp_set_obj_coef(glp_prob *P, int j, double coef);
/* set (change) obj. coefficient or constant term */

void glp_set_mat_row(glp_prob *P, int i, int len, const int ind[],
      const double val[]);
/* set (replace) row of the constraint matrix */

void glp_set_mat_col(glp_prob *P, int j, int len, const int ind[],
      const double val[]);
/* set (replace) column of the constraint matrix */

void glp_load_matrix(glp_prob *P, int ne, const int ia[],
      const int ja[], const double ar[]);
/* load (replace) the whole constraint matrix */

int glp_check_dup(int m, int n, int ne, const int ia[], const int ja[]);
/* check for duplicate elements in sparse matrix */

void glp_sort_matrix(glp_prob *P);
/* sort elements of the constraint matrix */

void glp_del_rows(glp_prob *P, int nrs, const int num[]);
/* delete specified rows from problem object */

void glp_del_cols(glp_prob *P, int ncs, const int num[]);
/* delete specified columns from problem object */

void glp_copy_prob(glp_prob *dest, glp_prob *prob, int names);
/* copy problem object content */

void glp_erase_prob(glp_prob *P);
/* erase problem object content */

void glp_delete_prob(glp_prob *P);
/* delete problem object */

const char *glp_get_prob_name(glp_prob *P);
/* retrieve problem name */

const char *glp_get_obj_name(glp_prob *P);
/* retrieve objective function name */

int glp_get_obj_dir(glp_prob *P);
/* retrieve optimization direction flag */

int glp_get_num_rows(glp_prob *P);
/* retrieve number of rows */

int glp_get_num_cols(glp_prob *P);
/* retrieve number of columns */

const char *glp_get_row_name(glp_prob *P, int i);
/* retrieve row name */

const char *glp_get_col_name(glp_prob *P, int j);
/* retrieve column name */

int glp_get_row_type(glp_prob *P, int i);
/* retrieve row type */

double glp_get_row_lb(glp_prob *P, int i);
/* retrieve row lower bound */

double glp_get_row_ub(glp_prob *P, int i);
/* retrieve row upper bound */

int glp_get_col_type(glp_prob *P, int j);
/* retrieve column type */

double glp_get_col_lb(glp_prob *P, int j);
/* retrieve column lower bound */

double glp_get_col_ub(glp_prob *P, int j);
/* retrieve column upper bound */

double glp_get_obj_coef(glp_prob *P, int j);
/* retrieve obj. coefficient or constant term */

int glp_get_num_nz(glp_prob *P);
/* retrieve number of constraint coefficients */

int glp_get_mat_row(glp_prob *P, int i, int ind[], double val[]);
/* retrieve row of the constraint matrix */

int glp_get_mat_col(glp_prob *P, int j, int ind[], double val[]);
/* retrieve column of the constraint matrix */

void glp_create_index(glp_prob *P);
/* create the name index */

int glp_find_row(glp_prob *P, const char *name);
/* find row by its name */

int glp_find_col(glp_prob *P, const char *name);
/* find column by its name */

void glp_delete_index(glp_prob *P);
/* delete the name index */

void glp_set_rii(glp_prob *P, int i, double rii);
/* set (change) row scale factor */

void glp_set_sjj(glp_prob *P, int j, double sjj);
/* set (change) column scale factor */

double glp_get_rii(glp_prob *P, int i);
/* retrieve row scale factor */

double glp_get_sjj(glp_prob *P, int j);
/* retrieve column scale factor */

void glp_scale_prob(glp_prob *P, int flags);
/* scale problem data */

void glp_unscale_prob(glp_prob *P);
/* unscale problem data */

void glp_set_row_stat(glp_prob *P, int i, int stat);
/* set (change) row status */

void glp_set_col_stat(glp_prob *P, int j, int stat);
/* set (change) column status */

void glp_std_basis(glp_prob *P);
/* construct standard initial LP basis */

void glp_adv_basis(glp_prob *P, int flags);
/* construct advanced initial LP basis */

void glp_cpx_basis(glp_prob *P);
/* construct Bixby's initial LP basis */

int glp_simplex(glp_prob *P, const glp_smcp *parm);
/* solve LP problem with the simplex method */

int glp_exact(glp_prob *P, const glp_smcp *parm);
/* solve LP problem in exact arithmetic */

void glp_init_smcp(glp_smcp *parm);
/* initialize simplex method control parameters */

int glp_get_status(glp_prob *P);
/* retrieve generic status of basic solution */

int glp_get_prim_stat(glp_prob *P);
/* retrieve status of primal basic solution */

int glp_get_dual_stat(glp_prob *P);
/* retrieve status of dual basic solution */

double glp_get_obj_val(glp_prob *P);
/* retrieve objective value (basic solution) */

int glp_get_row_stat(glp_prob *P, int i);
/* retrieve row status */

double glp_get_row_prim(glp_prob *P, int i);
/* retrieve row primal value (basic solution) */

double glp_get_row_dual(glp_prob *P, int i);
/* retrieve row dual value (basic solution) */

int glp_get_col_stat(glp_prob *P, int j);
/* retrieve column status */

double glp_get_col_prim(glp_prob *P, int j);
/* retrieve column primal value (basic solution) */

double glp_get_col_dual(glp_prob *P, int j);
/* retrieve column dual value (basic solution) */

int glp_get_unbnd_ray(glp_prob *P);
/* determine variable causing unboundedness */

int glp_interior(glp_prob *P, const glp_iptcp *parm);
/* solve LP problem with the interior-point method */

void glp_init_iptcp(glp_iptcp *parm);
/* initialize interior-point solver control parameters */

int glp_ipt_status(glp_prob *P);
/* retrieve status of interior-point solution */

double glp_ipt_obj_val(glp_prob *P);
/* retrieve objective value (interior point) */

double glp_ipt_row_prim(glp_prob *P, int i);
/* retrieve row primal value (interior point) */

double glp_ipt_row_dual(glp_prob *P, int i);
/* retrieve row dual value (interior point) */

double glp_ipt_col_prim(glp_prob *P, int j);
/* retrieve column primal value (interior point) */

double glp_ipt_col_dual(glp_prob *P, int j);
/* retrieve column dual value (interior point) */

void glp_set_col_kind(glp_prob *P, int j, int kind);
/* set (change) column kind */

int glp_get_col_kind(glp_prob *P, int j);
/* retrieve column kind */

int glp_get_num_int(glp_prob *P);
/* retrieve number of integer columns */

int glp_get_num_bin(glp_prob *P);
/* retrieve number of binary columns */

int glp_intopt(glp_prob *P, const glp_iocp *parm);
/* solve MIP problem with the branch-and-bound method */

void glp_init_iocp(glp_iocp *parm);
/* initialize integer optimizer control parameters */

int glp_mip_status(glp_prob *P);
/* retrieve status of MIP solution */

double glp_mip_obj_val(glp_prob *P);
/* retrieve objective value (MIP solution) */

double glp_mip_row_val(glp_prob *P, int i);
/* retrieve row value (MIP solution) */

double glp_mip_col_val(glp_prob *P, int j);
/* retrieve column value (MIP solution) */

int glp_print_sol(glp_prob *P, const char *fname);
/* write basic solution in printable format */

int glp_read_sol(glp_prob *P, const char *fname);
/* read basic solution from text file */

int glp_write_sol(glp_prob *P, const char *fname);
/* write basic solution to text file */

int glp_print_ranges(glp_prob *P, int len, const int list[],
      int flags, const char *fname);
/* print sensitivity analysis report */

int glp_print_ipt(glp_prob *P, const char *fname);
/* write interior-point solution in printable format */

int glp_read_ipt(glp_prob *P, const char *fname);
/* read interior-point solution from text file */

int glp_write_ipt(glp_prob *P, const char *fname);
/* write interior-point solution to text file */

int glp_print_mip(glp_prob *P, const char *fname);
/* write MIP solution in printable format */

int glp_read_mip(glp_prob *P, const char *fname);
/* read MIP solution from text file */

int glp_write_mip(glp_prob *P, const char *fname);
/* write MIP solution to text file */

int glp_bf_exists(glp_prob *P);
/* check if the basis factorization exists */

int glp_factorize(glp_prob *P);
/* compute the basis factorization */

int glp_bf_updated(glp_prob *P);
/* check if the basis factorization has been updated */

void glp_get_bfcp(glp_prob *P, glp_bfcp *parm);
/* retrieve basis factorization control parameters */

void glp_set_bfcp(glp_prob *P, const glp_bfcp *parm);
/* change basis factorization control parameters */

int glp_get_bhead(glp_prob *P, int k);
/* retrieve the basis header information */

int glp_get_row_bind(glp_prob *P, int i);
/* retrieve row index in the basis header */

int glp_get_col_bind(glp_prob *P, int j);
/* retrieve column index in the basis header */

void glp_ftran(glp_prob *P, double x[]);
/* perform forward transformation (solve system B*x = b) */

void glp_btran(glp_prob *P, double x[]);
/* perform backward transformation (solve system B'*x = b) */

int glp_warm_up(glp_prob *P);
/* "warm up" LP basis */

int glp_eval_tab_row(glp_prob *P, int k, int ind[], double val[]);
/* compute row of the simplex tableau */

int glp_eval_tab_col(glp_prob *P, int k, int ind[], double val[]);
/* compute column of the simplex tableau */

int glp_transform_row(glp_prob *P, int len, int ind[], double val[]);
/* transform explicitly specified row */

int glp_transform_col(glp_prob *P, int len, int ind[], double val[]);
/* transform explicitly specified column */

int glp_prim_rtest(glp_prob *P, int len, const int ind[],
      const double val[], int dir, double eps);
/* perform primal ratio test */

int glp_dual_rtest(glp_prob *P, int len, const int ind[],
      const double val[], int dir, double eps);
/* perform dual ratio test */

void glp_analyze_bound(glp_prob *P, int k, double *value1, int *var1,
      double *value2, int *var2);
/* analyze active bound of non-basic variable */

void glp_analyze_coef(glp_prob *P, int k, double *coef1, int *var1,
      double *value1, double *coef2, int *var2, double *value2);
/* analyze objective coefficient at basic variable */

int glp_ios_reason(glp_tree *T);
/* determine reason for calling the callback routine */

glp_prob *glp_ios_get_prob(glp_tree *T);
/* access the problem object */

void glp_ios_tree_size(glp_tree *T, int *a_cnt, int *n_cnt,
      int *t_cnt);
/* determine size of the branch-and-bound tree */

int glp_ios_curr_node(glp_tree *T);
/* determine current active subproblem */

int glp_ios_next_node(glp_tree *T, int p);
/* determine next active subproblem */

int glp_ios_prev_node(glp_tree *T, int p);
/* determine previous active subproblem */

int glp_ios_up_node(glp_tree *T, int p);
/* determine parent subproblem */

int glp_ios_node_level(glp_tree *T, int p);
/* determine subproblem level */

double glp_ios_node_bound(glp_tree *T, int p);
/* determine subproblem local bound */

int glp_ios_best_node(glp_tree *T);
/* find active subproblem with best local bound */

double glp_ios_mip_gap(glp_tree *T);
/* compute relative MIP gap */

void *glp_ios_node_data(glp_tree *T, int p);
/* access subproblem application-specific data */

void glp_ios_row_attr(glp_tree *T, int i, glp_attr *attr);
/* retrieve additional row attributes */

int glp_ios_pool_size(glp_tree *T);
/* determine current size of the cut pool */

int glp_ios_add_row(glp_tree *T,
      const char *name, int klass, int flags, int len, const int ind[],
      const double val[], int type, double rhs);
/* add row (constraint) to the cut pool */

void glp_ios_del_row(glp_tree *T, int i);
/* remove row (constraint) from the cut pool */

void glp_ios_clear_pool(glp_tree *T);
/* remove all rows (constraints) from the cut pool */

int glp_ios_can_branch(glp_tree *T, int j);
/* check if can branch upon specified variable */

void glp_ios_branch_upon(glp_tree *T, int j, int sel);
/* choose variable to branch upon */

void glp_ios_select_node(glp_tree *T, int p);
/* select subproblem to continue the search */

int glp_ios_heur_sol(glp_tree *T, const double x[]);
/* provide solution found by heuristic */

void glp_ios_terminate(glp_tree *T);
/* terminate the solution process */

void glp_init_mpscp(glp_mpscp *parm);
/* initialize MPS format control parameters */

int glp_read_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
      const char *fname);
/* read problem data in MPS format */

int glp_write_mps(glp_prob *P, int fmt, const glp_mpscp *parm,
      const char *fname);
/* write problem data in MPS format */

void glp_init_cpxcp(glp_cpxcp *parm);
/* initialize CPLEX LP format control parameters */

int glp_read_lp(glp_prob *P, const glp_cpxcp *parm, const char *fname);
/* read problem data in CPLEX LP format */

int glp_write_lp(glp_prob *P, const glp_cpxcp *parm, const char *fname);
/* write problem data in CPLEX LP format */

int glp_read_prob(glp_prob *P, int flags, const char *fname);
/* read problem data in GLPK format */

int glp_write_prob(glp_prob *P, int flags, const char *fname);
/* write problem data in GLPK format */

glp_tran *glp_mpl_alloc_wksp(void);
/* allocate the MathProg translator workspace */

int glp_mpl_read_model(glp_tran *tran, const char *fname, int skip);
/* read and translate model section */

int glp_mpl_read_data(glp_tran *tran, const char *fname);
/* read and translate data section */

int glp_mpl_generate(glp_tran *tran, const char *fname);
/* generate the model */

void glp_mpl_build_prob(glp_tran *tran, glp_prob *prob);
/* build LP/MIP problem instance from the model */

int glp_mpl_postsolve(glp_tran *tran, glp_prob *prob, int sol);
/* postsolve the model */

void glp_mpl_free_wksp(glp_tran *tran);
/* free the MathProg translator workspace */

int glp_main(int argc, const char *argv[]);
/* stand-alone LP/MIP solver */

/**********************************************************************/

#ifndef GLP_LONG_DEFINED
#define GLP_LONG_DEFINED
typedef struct { int lo, hi; } glp_long;
/* long integer data type */
#endif

/* igraph-specific hack;
 * glp_at_error() is not present in this old GLPK version
 * but we need it for proper error reporting */
int glp_at_error(void);

int glp_init_env(void);
/* initialize GLPK environment */

const char *glp_version(void);
/* determine library version */

int glp_free_env(void);
/* free GLPK environment */

void glp_printf(const char *fmt, ...);
/* write formatted output to terminal */

void glp_vprintf(const char *fmt, va_list arg);
/* write formatted output to terminal */

int glp_term_out(int flag);
/* enable/disable terminal output */

void glp_term_hook(int (*func)(void *info, const char *s), void *info);
/* install hook to intercept terminal output */

int glp_open_tee(const char *fname);
/* start copying terminal output to text file */

int glp_close_tee(void);
/* stop copying terminal output to text file */

#ifndef GLP_ERROR_DEFINED
#define GLP_ERROR_DEFINED
typedef void (*_glp_error)(const char *fmt, ...);
#endif

#define glp_error glp_error_(__FILE__, __LINE__)
_glp_error glp_error_(const char *file, int line);
/* display error message and terminate execution */

#define glp_assert(expr) \
      ((void)((expr) || (glp_assert_(#expr, __FILE__, __LINE__), 1)))
void glp_assert_(const char *expr, const char *file, int line);
/* check for logical condition */

void glp_error_hook(void (*func)(void *info), void *info);
/* install hook to intercept abnormal termination */

void *glp_malloc(int size);
/* allocate memory block */

void *glp_calloc(int n, int size);
/* allocate memory block */

void glp_free(void *ptr);
/* free memory block */

void glp_mem_limit(int limit);
/* set memory usage limit */

void glp_mem_usage(int *count, int *cpeak, glp_long *total,
      glp_long *tpeak);
/* get memory usage information */

glp_long glp_time(void);
/* determine current universal time */

double glp_difftime(glp_long t1, glp_long t0);
/* compute difference between two time values */

/**********************************************************************/

#ifndef GLP_DATA_DEFINED
#define GLP_DATA_DEFINED
typedef struct glp_data glp_data;
/* plain data file */
#endif

glp_data *glp_sdf_open_file(const char *fname);
/* open plain data file */

void glp_sdf_set_jump(glp_data *data, void *jump);
/* set up error handling */

void glp_sdf_error(glp_data *data, const char *fmt, ...);
/* print error message */

void glp_sdf_warning(glp_data *data, const char *fmt, ...);
/* print warning message */

int glp_sdf_read_int(glp_data *data);
/* read integer number */

double glp_sdf_read_num(glp_data *data);
/* read floating-point number */

const char *glp_sdf_read_item(glp_data *data);
/* read data item */

const char *glp_sdf_read_text(glp_data *data);
/* read text until end of line */

int glp_sdf_line(glp_data *data);
/* determine current line number */

void glp_sdf_close_file(glp_data *data);
/* close plain data file */

/**********************************************************************/

typedef struct _glp_graph glp_graph;
typedef struct _glp_vertex glp_vertex;
typedef struct _glp_arc glp_arc;

struct _glp_graph
{     /* graph descriptor */
      void *pool; /* DMP *pool; */
      /* memory pool to store graph components */
      char *name;
      /* graph name (1 to 255 chars); NULL means no name is assigned
         to the graph */
      int nv_max;
      /* length of the vertex list (enlarged automatically) */
      int nv;
      /* number of vertices in the graph, 0 <= nv <= nv_max */
      int na;
      /* number of arcs in the graph, na >= 0 */
      glp_vertex **v; /* glp_vertex *v[1+nv_max]; */
      /* v[i], 1 <= i <= nv, is a pointer to i-th vertex */
      void *index; /* AVL *index; */
      /* vertex index to find vertices by their names; NULL means the
         index does not exist */
      int v_size;
      /* size of data associated with each vertex (0 to 256 bytes) */
      int a_size;
      /* size of data associated with each arc (0 to 256 bytes) */
};

struct _glp_vertex
{     /* vertex descriptor */
      int i;
      /* vertex ordinal number, 1 <= i <= nv */
      char *name;
      /* vertex name (1 to 255 chars); NULL means no name is assigned
         to the vertex */
      void *entry; /* AVLNODE *entry; */
      /* pointer to corresponding entry in the vertex index; NULL means
         that either the index does not exist or the vertex has no name
         assigned */
      void *data;
      /* pointer to data associated with the vertex */
      void *temp;
      /* working pointer */
      glp_arc *in;
      /* pointer to the (unordered) list of incoming arcs */
      glp_arc *out;
      /* pointer to the (unordered) list of outgoing arcs */
};

struct _glp_arc
{     /* arc descriptor */
      glp_vertex *tail;
      /* pointer to the tail endpoint */
      glp_vertex *head;
      /* pointer to the head endpoint */
      void *data;
      /* pointer to data associated with the arc */
      void *temp;
      /* working pointer */
      glp_arc *t_prev;
      /* pointer to previous arc having the same tail endpoint */
      glp_arc *t_next;
      /* pointer to next arc having the same tail endpoint */
      glp_arc *h_prev;
      /* pointer to previous arc having the same head endpoint */
      glp_arc *h_next;
      /* pointer to next arc having the same head endpoint */
};

glp_graph *glp_create_graph(int v_size, int a_size);
/* create graph */

void glp_set_graph_name(glp_graph *G, const char *name);
/* assign (change) graph name */

int glp_add_vertices(glp_graph *G, int nadd);
/* add new vertices to graph */

void glp_set_vertex_name(glp_graph *G, int i, const char *name);
/* assign (change) vertex name */

glp_arc *glp_add_arc(glp_graph *G, int i, int j);
/* add new arc to graph */

void glp_del_vertices(glp_graph *G, int ndel, const int num[]);
/* delete vertices from graph */

void glp_del_arc(glp_graph *G, glp_arc *a);
/* delete arc from graph */

void glp_erase_graph(glp_graph *G, int v_size, int a_size);
/* erase graph content */

void glp_delete_graph(glp_graph *G);
/* delete graph */

void glp_create_v_index(glp_graph *G);
/* create vertex name index */

int glp_find_vertex(glp_graph *G, const char *name);
/* find vertex by its name */

void glp_delete_v_index(glp_graph *G);
/* delete vertex name index */

int glp_read_graph(glp_graph *G, const char *fname);
/* read graph from plain text file */

int glp_write_graph(glp_graph *G, const char *fname);
/* write graph to plain text file */

void glp_mincost_lp(glp_prob *P, glp_graph *G, int names, int v_rhs,
      int a_low, int a_cap, int a_cost);
/* convert minimum cost flow problem to LP */

int glp_mincost_okalg(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, double *sol, int a_x, int v_pi);
/* find minimum-cost flow with out-of-kilter algorithm */

void glp_maxflow_lp(glp_prob *P, glp_graph *G, int names, int s,
      int t, int a_cap);
/* convert maximum flow problem to LP */

int glp_maxflow_ffalg(glp_graph *G, int s, int t, int a_cap,
      double *sol, int a_x, int v_cut);
/* find maximal flow with Ford-Fulkerson algorithm */

int glp_check_asnprob(glp_graph *G, int v_set);
/* check correctness of assignment problem data */

/* assignment problem formulation: */
#define GLP_ASN_MIN        1  /* perfect matching (minimization) */
#define GLP_ASN_MAX        2  /* perfect matching (maximization) */
#define GLP_ASN_MMP        3  /* maximum matching */

int glp_asnprob_lp(glp_prob *P, int form, glp_graph *G, int names,
      int v_set, int a_cost);
/* convert assignment problem to LP */

int glp_asnprob_okalg(int form, glp_graph *G, int v_set, int a_cost,
      double *sol, int a_x);
/* solve assignment problem with out-of-kilter algorithm */

int glp_asnprob_hall(glp_graph *G, int v_set, int a_x);
/* find bipartite matching of maximum cardinality */

double glp_cpp(glp_graph *G, int v_t, int v_es, int v_ls);
/* solve critical path problem */

int glp_read_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const char *fname);
/* read min-cost flow problem data in DIMACS format */

int glp_write_mincost(glp_graph *G, int v_rhs, int a_low, int a_cap,
      int a_cost, const char *fname);
/* write min-cost flow problem data in DIMACS format */

int glp_read_maxflow(glp_graph *G, int *s, int *t, int a_cap,
      const char *fname);
/* read maximum flow problem data in DIMACS format */

int glp_write_maxflow(glp_graph *G, int s, int t, int a_cap,
      const char *fname);
/* write maximum flow problem data in DIMACS format */

int glp_read_asnprob(glp_graph *G, int v_set, int a_cost, const char
      *fname);
/* read assignment problem data in DIMACS format */

int glp_write_asnprob(glp_graph *G, int v_set, int a_cost, const char
      *fname);
/* write assignment problem data in DIMACS format */

int glp_read_ccdata(glp_graph *G, int v_wgt, const char *fname);
/* read graph in DIMACS clique/coloring format */

int glp_write_ccdata(glp_graph *G, int v_wgt, const char *fname);
/* write graph in DIMACS clique/coloring format */

int glp_netgen(glp_graph *G, int v_rhs, int a_cap, int a_cost,
      const int parm[1+15]);
/* Klingman's network problem generator */

int glp_gridgen(glp_graph *G, int v_rhs, int a_cap, int a_cost,
      const int parm[1+14]);
/* grid-like network problem generator */

int glp_rmfgen(glp_graph *G, int *s, int *t, int a_cap,
      const int parm[1+5]);
/* Goldfarb's maximum flow problem generator */

int glp_weak_comp(glp_graph *G, int v_num);
/* find all weakly connected components of graph */

int glp_strong_comp(glp_graph *G, int v_num);
/* find all strongly connected components of graph */

int glp_top_sort(glp_graph *G, int v_num);
/* topological sorting of acyclic digraph */

int glp_wclique_exact(glp_graph *G, int v_wgt, double *sol, int v_set);
/* find maximum weight clique with exact algorithm */

/***********************************************************************
*  NOTE: All symbols defined below are obsolete and kept here only for
*        backward compatibility.
***********************************************************************/

#define LPX glp_prob

/* problem class: */
#define LPX_LP          100   /* linear programming (LP) */
#define LPX_MIP         101   /* mixed integer programming (MIP) */

/* type of auxiliary/structural variable: */
#define LPX_FR          110   /* free variable */
#define LPX_LO          111   /* variable with lower bound */
#define LPX_UP          112   /* variable with upper bound */
#define LPX_DB          113   /* double-bounded variable */
#define LPX_FX          114   /* fixed variable */

/* optimization direction flag: */
#define LPX_MIN         120   /* minimization */
#define LPX_MAX         121   /* maximization */

/* status of primal basic solution: */
#define LPX_P_UNDEF     132   /* primal solution is undefined */
#define LPX_P_FEAS      133   /* solution is primal feasible */
#define LPX_P_INFEAS    134   /* solution is primal infeasible */
#define LPX_P_NOFEAS    135   /* no primal feasible solution exists */

/* status of dual basic solution: */
#define LPX_D_UNDEF     136   /* dual solution is undefined */
#define LPX_D_FEAS      137   /* solution is dual feasible */
#define LPX_D_INFEAS    138   /* solution is dual infeasible */
#define LPX_D_NOFEAS    139   /* no dual feasible solution exists */

/* status of auxiliary/structural variable: */
#define LPX_BS          140   /* basic variable */
#define LPX_NL          141   /* non-basic variable on lower bound */
#define LPX_NU          142   /* non-basic variable on upper bound */
#define LPX_NF          143   /* non-basic free variable */
#define LPX_NS          144   /* non-basic fixed variable */

/* status of interior-point solution: */
#define LPX_T_UNDEF     150   /* interior solution is undefined */
#define LPX_T_OPT       151   /* interior solution is optimal */

/* kind of structural variable: */
#define LPX_CV          160   /* continuous variable */
#define LPX_IV          161   /* integer variable */

/* status of integer solution: */
#define LPX_I_UNDEF     170   /* integer solution is undefined */
#define LPX_I_OPT       171   /* integer solution is optimal */
#define LPX_I_FEAS      172   /* integer solution is feasible */
#define LPX_I_NOFEAS    173   /* no integer solution exists */

/* status codes reported by the routine lpx_get_status: */
#define LPX_OPT         180   /* optimal */
#define LPX_FEAS        181   /* feasible */
#define LPX_INFEAS      182   /* infeasible */
#define LPX_NOFEAS      183   /* no feasible */
#define LPX_UNBND       184   /* unbounded */
#define LPX_UNDEF       185   /* undefined */

/* exit codes returned by solver routines: */
#define LPX_E_OK        200   /* success */
#define LPX_E_EMPTY     201   /* empty problem */
#define LPX_E_BADB      202   /* invalid initial basis */
#define LPX_E_INFEAS    203   /* infeasible initial solution */
#define LPX_E_FAULT     204   /* unable to start the search */
#define LPX_E_OBJLL     205   /* objective lower limit reached */
#define LPX_E_OBJUL     206   /* objective upper limit reached */
#define LPX_E_ITLIM     207   /* iterations limit exhausted */
#define LPX_E_TMLIM     208   /* time limit exhausted */
#define LPX_E_NOFEAS    209   /* no feasible solution */
#define LPX_E_INSTAB    210   /* numerical instability */
#define LPX_E_SING      211   /* problems with basis matrix */
#define LPX_E_NOCONV    212   /* no convergence (interior) */
#define LPX_E_NOPFS     213   /* no primal feas. sol. (LP presolver) */
#define LPX_E_NODFS     214   /* no dual feas. sol. (LP presolver) */
#define LPX_E_MIPGAP    215   /* relative mip gap tolerance reached */

/* control parameter identifiers: */
#define LPX_K_MSGLEV    300   /* lp->msg_lev */
#define LPX_K_SCALE     301   /* lp->scale */
#define LPX_K_DUAL      302   /* lp->dual */
#define LPX_K_PRICE     303   /* lp->price */
#define LPX_K_RELAX     304   /* lp->relax */
#define LPX_K_TOLBND    305   /* lp->tol_bnd */
#define LPX_K_TOLDJ     306   /* lp->tol_dj */
#define LPX_K_TOLPIV    307   /* lp->tol_piv */
#define LPX_K_ROUND     308   /* lp->round */
#define LPX_K_OBJLL     309   /* lp->obj_ll */
#define LPX_K_OBJUL     310   /* lp->obj_ul */
#define LPX_K_ITLIM     311   /* lp->it_lim */
#define LPX_K_ITCNT     312   /* lp->it_cnt */
#define LPX_K_TMLIM     313   /* lp->tm_lim */
#define LPX_K_OUTFRQ    314   /* lp->out_frq */
#define LPX_K_OUTDLY    315   /* lp->out_dly */
#define LPX_K_BRANCH    316   /* lp->branch */
#define LPX_K_BTRACK    317   /* lp->btrack */
#define LPX_K_TOLINT    318   /* lp->tol_int */
#define LPX_K_TOLOBJ    319   /* lp->tol_obj */
#define LPX_K_MPSINFO   320   /* lp->mps_info */
#define LPX_K_MPSOBJ    321   /* lp->mps_obj */
#define LPX_K_MPSORIG   322   /* lp->mps_orig */
#define LPX_K_MPSWIDE   323   /* lp->mps_wide */
#define LPX_K_MPSFREE   324   /* lp->mps_free */
#define LPX_K_MPSSKIP   325   /* lp->mps_skip */
#define LPX_K_LPTORIG   326   /* lp->lpt_orig */
#define LPX_K_PRESOL    327   /* lp->presol */
#define LPX_K_BINARIZE  328   /* lp->binarize */
#define LPX_K_USECUTS   329   /* lp->use_cuts */
#define LPX_K_BFTYPE    330   /* lp->bfcp->type */
#define LPX_K_MIPGAP    331   /* lp->mip_gap */

#define LPX_C_COVER     0x01  /* mixed cover cuts */
#define LPX_C_CLIQUE    0x02  /* clique cuts */
#define LPX_C_GOMORY    0x04  /* Gomory's mixed integer cuts */
#define LPX_C_MIR       0x08  /* mixed integer rounding cuts */
#define LPX_C_ALL       0xFF  /* all cuts */

typedef struct
{     /* this structure contains results reported by the routines which
         checks Karush-Kuhn-Tucker conditions (for details see comments
         to those routines) */
      /*--------------------------------------------------------------*/
      /* xR - A * xS = 0 (KKT.PE) */
      double pe_ae_max;
      /* largest absolute error */
      int    pe_ae_row;
      /* number of row with largest absolute error */
      double pe_re_max;
      /* largest relative error */
      int    pe_re_row;
      /* number of row with largest relative error */
      int    pe_quality;
      /* quality of primal solution:
         'H' - high
         'M' - medium
         'L' - low
         '?' - primal solution is wrong */
      /*--------------------------------------------------------------*/
      /* l[k] <= x[k] <= u[k] (KKT.PB) */
      double pb_ae_max;
      /* largest absolute error */
      int    pb_ae_ind;
      /* number of variable with largest absolute error */
      double pb_re_max;
      /* largest relative error */
      int    pb_re_ind;
      /* number of variable with largest relative error */
      int    pb_quality;
      /* quality of primal feasibility:
         'H' - high
         'M' - medium
         'L' - low
         '?' - primal solution is infeasible */
      /*--------------------------------------------------------------*/
      /* A' * (dR - cR) + (dS - cS) = 0 (KKT.DE) */
      double de_ae_max;
      /* largest absolute error */
      int    de_ae_col;
      /* number of column with largest absolute error */
      double de_re_max;
      /* largest relative error */
      int    de_re_col;
      /* number of column with largest relative error */
      int    de_quality;
      /* quality of dual solution:
         'H' - high
         'M' - medium
         'L' - low
         '?' - dual solution is wrong */
      /*--------------------------------------------------------------*/
      /* d[k] >= 0 or d[k] <= 0 (KKT.DB) */
      double db_ae_max;
      /* largest absolute error */
      int    db_ae_ind;
      /* number of variable with largest absolute error */
      double db_re_max;
      /* largest relative error */
      int    db_re_ind;
      /* number of variable with largest relative error */
      int    db_quality;
      /* quality of dual feasibility:
         'H' - high
         'M' - medium
         'L' - low
         '?' - dual solution is infeasible */
      /*--------------------------------------------------------------*/
      /* (x[k] - bound of x[k]) * d[k] = 0 (KKT.CS) */
      double cs_ae_max;
      /* largest absolute error */
      int    cs_ae_ind;
      /* number of variable with largest absolute error */
      double cs_re_max;
      /* largest relative error */
      int    cs_re_ind;
      /* number of variable with largest relative error */
      int    cs_quality;
      /* quality of complementary slackness:
         'H' - high
         'M' - medium
         'L' - low
         '?' - primal and dual solutions are not complementary */
} LPXKKT;

#define lpx_create_prob _glp_lpx_create_prob
LPX *lpx_create_prob(void);
/* create problem object */

#define lpx_set_prob_name _glp_lpx_set_prob_name
void lpx_set_prob_name(LPX *lp, const char *name);
/* assign (change) problem name */

#define lpx_set_obj_name _glp_lpx_set_obj_name
void lpx_set_obj_name(LPX *lp, const char *name);
/* assign (change) objective function name */

#define lpx_set_obj_dir _glp_lpx_set_obj_dir
void lpx_set_obj_dir(LPX *lp, int dir);
/* set (change) optimization direction flag */

#define lpx_add_rows _glp_lpx_add_rows
int lpx_add_rows(LPX *lp, int nrs);
/* add new rows to problem object */

#define lpx_add_cols _glp_lpx_add_cols
int lpx_add_cols(LPX *lp, int ncs);
/* add new columns to problem object */

#define lpx_set_row_name _glp_lpx_set_row_name
void lpx_set_row_name(LPX *lp, int i, const char *name);
/* assign (change) row name */

#define lpx_set_col_name _glp_lpx_set_col_name
void lpx_set_col_name(LPX *lp, int j, const char *name);
/* assign (change) column name */

#define lpx_set_row_bnds _glp_lpx_set_row_bnds
void lpx_set_row_bnds(LPX *lp, int i, int type, double lb, double ub);
/* set (change) row bounds */

#define lpx_set_col_bnds _glp_lpx_set_col_bnds
void lpx_set_col_bnds(LPX *lp, int j, int type, double lb, double ub);
/* set (change) column bounds */

#define lpx_set_obj_coef _glp_lpx_set_obj_coef
void lpx_set_obj_coef(glp_prob *lp, int j, double coef);
/* set (change) obj. coefficient or constant term */

#define lpx_set_mat_row _glp_lpx_set_mat_row
void lpx_set_mat_row(LPX *lp, int i, int len, const int ind[],
      const double val[]);
/* set (replace) row of the constraint matrix */

#define lpx_set_mat_col _glp_lpx_set_mat_col
void lpx_set_mat_col(LPX *lp, int j, int len, const int ind[],
      const double val[]);
/* set (replace) column of the constraint matrix */

#define lpx_load_matrix _glp_lpx_load_matrix
void lpx_load_matrix(LPX *lp, int ne, const int ia[], const int ja[],
      const double ar[]);
/* load (replace) the whole constraint matrix */

#define lpx_del_rows _glp_lpx_del_rows
void lpx_del_rows(LPX *lp, int nrs, const int num[]);
/* delete specified rows from problem object */

#define lpx_del_cols _glp_lpx_del_cols
void lpx_del_cols(LPX *lp, int ncs, const int num[]);
/* delete specified columns from problem object */

#define lpx_delete_prob _glp_lpx_delete_prob
void lpx_delete_prob(LPX *lp);
/* delete problem object */

#define lpx_get_prob_name _glp_lpx_get_prob_name
const char *lpx_get_prob_name(LPX *lp);
/* retrieve problem name */

#define lpx_get_obj_name _glp_lpx_get_obj_name
const char *lpx_get_obj_name(LPX *lp);
/* retrieve objective function name */

#define lpx_get_obj_dir _glp_lpx_get_obj_dir
int lpx_get_obj_dir(LPX *lp);
/* retrieve optimization direction flag */

#define lpx_get_num_rows _glp_lpx_get_num_rows
int lpx_get_num_rows(LPX *lp);
/* retrieve number of rows */

#define lpx_get_num_cols _glp_lpx_get_num_cols
int lpx_get_num_cols(LPX *lp);
/* retrieve number of columns */

#define lpx_get_row_name _glp_lpx_get_row_name
const char *lpx_get_row_name(LPX *lp, int i);
/* retrieve row name */

#define lpx_get_col_name _glp_lpx_get_col_name
const char *lpx_get_col_name(LPX *lp, int j);
/* retrieve column name */

#define lpx_get_row_type _glp_lpx_get_row_type
int lpx_get_row_type(LPX *lp, int i);
/* retrieve row type */

#define lpx_get_row_lb _glp_lpx_get_row_lb
double lpx_get_row_lb(LPX *lp, int i);
/* retrieve row lower bound */

#define lpx_get_row_ub _glp_lpx_get_row_ub
double lpx_get_row_ub(LPX *lp, int i);
/* retrieve row upper bound */

#define lpx_get_row_bnds _glp_lpx_get_row_bnds
void lpx_get_row_bnds(LPX *lp, int i, int *typx, double *lb,
      double *ub);
/* retrieve row bounds */

#define lpx_get_col_type _glp_lpx_get_col_type
int lpx_get_col_type(LPX *lp, int j);
/* retrieve column type */

#define lpx_get_col_lb _glp_lpx_get_col_lb
double lpx_get_col_lb(LPX *lp, int j);
/* retrieve column lower bound */

#define lpx_get_col_ub _glp_lpx_get_col_ub
double lpx_get_col_ub(LPX *lp, int j);
/* retrieve column upper bound */

#define lpx_get_col_bnds _glp_lpx_get_col_bnds
void lpx_get_col_bnds(LPX *lp, int j, int *typx, double *lb,
      double *ub);
/* retrieve column bounds */

#define lpx_get_obj_coef _glp_lpx_get_obj_coef
double lpx_get_obj_coef(LPX *lp, int j);
/* retrieve obj. coefficient or constant term */

#define lpx_get_num_nz _glp_lpx_get_num_nz
int lpx_get_num_nz(LPX *lp);
/* retrieve number of constraint coefficients */

#define lpx_get_mat_row _glp_lpx_get_mat_row
int lpx_get_mat_row(LPX *lp, int i, int ind[], double val[]);
/* retrieve row of the constraint matrix */

#define lpx_get_mat_col _glp_lpx_get_mat_col
int lpx_get_mat_col(LPX *lp, int j, int ind[], double val[]);
/* retrieve column of the constraint matrix */

#define lpx_create_index _glp_lpx_create_index
void lpx_create_index(LPX *lp);
/* create the name index */

#define lpx_find_row _glp_lpx_find_row
int lpx_find_row(LPX *lp, const char *name);
/* find row by its name */

#define lpx_find_col _glp_lpx_find_col
int lpx_find_col(LPX *lp, const char *name);
/* find column by its name */

#define lpx_delete_index _glp_lpx_delete_index
void lpx_delete_index(LPX *lp);
/* delete the name index */

#define lpx_scale_prob _glp_lpx_scale_prob
void lpx_scale_prob(LPX *lp);
/* scale problem data */

#define lpx_unscale_prob _glp_lpx_unscale_prob
void lpx_unscale_prob(LPX *lp);
/* unscale problem data */

#define lpx_set_row_stat _glp_lpx_set_row_stat
void lpx_set_row_stat(LPX *lp, int i, int stat);
/* set (change) row status */

#define lpx_set_col_stat _glp_lpx_set_col_stat
void lpx_set_col_stat(LPX *lp, int j, int stat);
/* set (change) column status */

#define lpx_std_basis _glp_lpx_std_basis
void lpx_std_basis(LPX *lp);
/* construct standard initial LP basis */

#define lpx_adv_basis _glp_lpx_adv_basis
void lpx_adv_basis(LPX *lp);
/* construct advanced initial LP basis */

#define lpx_cpx_basis _glp_lpx_cpx_basis
void lpx_cpx_basis(LPX *lp);
/* construct Bixby's initial LP basis */

#define lpx_simplex _glp_lpx_simplex
int lpx_simplex(LPX *lp);
/* easy-to-use driver to the simplex method */

#define lpx_exact _glp_lpx_exact
int lpx_exact(LPX *lp);
/* easy-to-use driver to the exact simplex method */

#define lpx_get_status _glp_lpx_get_status
int lpx_get_status(LPX *lp);
/* retrieve generic status of basic solution */

#define lpx_get_prim_stat _glp_lpx_get_prim_stat
int lpx_get_prim_stat(LPX *lp);
/* retrieve primal status of basic solution */

#define lpx_get_dual_stat _glp_lpx_get_dual_stat
int lpx_get_dual_stat(LPX *lp);
/* retrieve dual status of basic solution */

#define lpx_get_obj_val _glp_lpx_get_obj_val
double lpx_get_obj_val(LPX *lp);
/* retrieve objective value (basic solution) */

#define lpx_get_row_stat _glp_lpx_get_row_stat
int lpx_get_row_stat(LPX *lp, int i);
/* retrieve row status (basic solution) */

#define lpx_get_row_prim _glp_lpx_get_row_prim
double lpx_get_row_prim(LPX *lp, int i);
/* retrieve row primal value (basic solution) */

#define lpx_get_row_dual _glp_lpx_get_row_dual
double lpx_get_row_dual(LPX *lp, int i);
/* retrieve row dual value (basic solution) */

#define lpx_get_row_info _glp_lpx_get_row_info
void lpx_get_row_info(LPX *lp, int i, int *tagx, double *vx,
      double *dx);
/* obtain row solution information */

#define lpx_get_col_stat _glp_lpx_get_col_stat
int lpx_get_col_stat(LPX *lp, int j);
/* retrieve column status (basic solution) */

#define lpx_get_col_prim _glp_lpx_get_col_prim
double lpx_get_col_prim(LPX *lp, int j);
/* retrieve column primal value (basic solution) */

#define lpx_get_col_dual _glp_lpx_get_col_dual
double lpx_get_col_dual(glp_prob *lp, int j);
/* retrieve column dual value (basic solution) */

#define lpx_get_col_info _glp_lpx_get_col_info
void lpx_get_col_info(LPX *lp, int j, int *tagx, double *vx,
      double *dx);
/* obtain column solution information (obsolete) */

#define lpx_get_ray_info _glp_lpx_get_ray_info
int lpx_get_ray_info(LPX *lp);
/* determine what causes primal unboundness */

#define lpx_check_kkt _glp_lpx_check_kkt
void lpx_check_kkt(LPX *lp, int scaled, LPXKKT *kkt);
/* check Karush-Kuhn-Tucker conditions */

#define lpx_warm_up _glp_lpx_warm_up
int lpx_warm_up(LPX *lp);
/* "warm up" LP basis */

#define lpx_eval_tab_row _glp_lpx_eval_tab_row
int lpx_eval_tab_row(LPX *lp, int k, int ind[], double val[]);
/* compute row of the simplex table */

#define lpx_eval_tab_col _glp_lpx_eval_tab_col
int lpx_eval_tab_col(LPX *lp, int k, int ind[], double val[]);
/* compute column of the simplex table */

#define lpx_transform_row _glp_lpx_transform_row
int lpx_transform_row(LPX *lp, int len, int ind[], double val[]);
/* transform explicitly specified row */

#define lpx_transform_col _glp_lpx_transform_col
int lpx_transform_col(LPX *lp, int len, int ind[], double val[]);
/* transform explicitly specified column */

#define lpx_prim_ratio_test _glp_lpx_prim_ratio_test
int lpx_prim_ratio_test(LPX *lp, int len, const int ind[],
      const double val[], int how, double tol);
/* perform primal ratio test */

#define lpx_dual_ratio_test _glp_lpx_dual_ratio_test
int lpx_dual_ratio_test(LPX *lp, int len, const int ind[],
      const double val[], int how, double tol);
/* perform dual ratio test */

#define lpx_interior _glp_lpx_interior
int lpx_interior(LPX *lp);
/* easy-to-use driver to the interior point method */

#define lpx_ipt_status _glp_lpx_ipt_status
int lpx_ipt_status(LPX *lp);
/* retrieve status of interior-point solution */

#define lpx_ipt_obj_val _glp_lpx_ipt_obj_val
double lpx_ipt_obj_val(LPX *lp);
/* retrieve objective value (interior point) */

#define lpx_ipt_row_prim _glp_lpx_ipt_row_prim
double lpx_ipt_row_prim(LPX *lp, int i);
/* retrieve row primal value (interior point) */

#define lpx_ipt_row_dual _glp_lpx_ipt_row_dual
double lpx_ipt_row_dual(LPX *lp, int i);
/* retrieve row dual value (interior point) */

#define lpx_ipt_col_prim _glp_lpx_ipt_col_prim
double lpx_ipt_col_prim(LPX *lp, int j);
/* retrieve column primal value (interior point) */

#define lpx_ipt_col_dual _glp_lpx_ipt_col_dual
double lpx_ipt_col_dual(LPX *lp, int j);
/* retrieve column dual value (interior point) */

#define lpx_set_class _glp_lpx_set_class
void lpx_set_class(LPX *lp, int klass);
/* set problem class */

#define lpx_get_class _glp_lpx_get_class
int lpx_get_class(LPX *lp);
/* determine problem klass */

#define lpx_set_col_kind _glp_lpx_set_col_kind
void lpx_set_col_kind(LPX *lp, int j, int kind);
/* set (change) column kind */

#define lpx_get_col_kind _glp_lpx_get_col_kind
int lpx_get_col_kind(LPX *lp, int j);
/* retrieve column kind */

#define lpx_get_num_int _glp_lpx_get_num_int
int lpx_get_num_int(LPX *lp);
/* retrieve number of integer columns */

#define lpx_get_num_bin _glp_lpx_get_num_bin
int lpx_get_num_bin(LPX *lp);
/* retrieve number of binary columns */

#define lpx_integer _glp_lpx_integer
int lpx_integer(LPX *lp);
/* easy-to-use driver to the branch-and-bound method */

#define lpx_intopt _glp_lpx_intopt
int lpx_intopt(LPX *lp);
/* easy-to-use driver to the branch-and-bound method */

#define lpx_mip_status _glp_lpx_mip_status
int lpx_mip_status(LPX *lp);
/* retrieve status of MIP solution */

#define lpx_mip_obj_val _glp_lpx_mip_obj_val
double lpx_mip_obj_val(LPX *lp);
/* retrieve objective value (MIP solution) */

#define lpx_mip_row_val _glp_lpx_mip_row_val
double lpx_mip_row_val(LPX *lp, int i);
/* retrieve row value (MIP solution) */

#define lpx_mip_col_val _glp_lpx_mip_col_val
double lpx_mip_col_val(LPX *lp, int j);
/* retrieve column value (MIP solution) */

#define lpx_check_int _glp_lpx_check_int
void lpx_check_int(LPX *lp, LPXKKT *kkt);
/* check integer feasibility conditions */

#define lpx_reset_parms _glp_lpx_reset_parms
void lpx_reset_parms(LPX *lp);
/* reset control parameters to default values */

#define lpx_set_int_parm _glp_lpx_set_int_parm
void lpx_set_int_parm(LPX *lp, int parm, int val);
/* set (change) integer control parameter */

#define lpx_get_int_parm _glp_lpx_get_int_parm
int lpx_get_int_parm(LPX *lp, int parm);
/* query integer control parameter */

#define lpx_set_real_parm _glp_lpx_set_real_parm
void lpx_set_real_parm(LPX *lp, int parm, double val);
/* set (change) real control parameter */

#define lpx_get_real_parm _glp_lpx_get_real_parm
double lpx_get_real_parm(LPX *lp, int parm);
/* query real control parameter */

#define lpx_read_mps _glp_lpx_read_mps
LPX *lpx_read_mps(const char *fname);
/* read problem data in fixed MPS format */

#define lpx_write_mps _glp_lpx_write_mps
int lpx_write_mps(LPX *lp, const char *fname);
/* write problem data in fixed MPS format */

#define lpx_read_bas _glp_lpx_read_bas
int lpx_read_bas(LPX *lp, const char *fname);
/* read LP basis in fixed MPS format */

#define lpx_write_bas _glp_lpx_write_bas
int lpx_write_bas(LPX *lp, const char *fname);
/* write LP basis in fixed MPS format */

#define lpx_read_freemps _glp_lpx_read_freemps
LPX *lpx_read_freemps(const char *fname);
/* read problem data in free MPS format */

#define lpx_write_freemps _glp_lpx_write_freemps
int lpx_write_freemps(LPX *lp, const char *fname);
/* write problem data in free MPS format */

#define lpx_read_cpxlp _glp_lpx_read_cpxlp
LPX *lpx_read_cpxlp(const char *fname);
/* read problem data in CPLEX LP format */

#define lpx_write_cpxlp _glp_lpx_write_cpxlp
int lpx_write_cpxlp(LPX *lp, const char *fname);
/* write problem data in CPLEX LP format */

#define lpx_read_model _glp_lpx_read_model
LPX *lpx_read_model(const char *model, const char *data,
      const char *output);
/* read LP/MIP model written in GNU MathProg language */

#define lpx_print_prob _glp_lpx_print_prob
int lpx_print_prob(LPX *lp, const char *fname);
/* write problem data in plain text format */

#define lpx_print_sol _glp_lpx_print_sol
int lpx_print_sol(LPX *lp, const char *fname);
/* write LP problem solution in printable format */

#define lpx_print_sens_bnds _glp_lpx_print_sens_bnds
int lpx_print_sens_bnds(LPX *lp, const char *fname);
/* write bounds sensitivity information */

#define lpx_print_ips _glp_lpx_print_ips
int lpx_print_ips(LPX *lp, const char *fname);
/* write interior point solution in printable format */

#define lpx_print_mip _glp_lpx_print_mip
int lpx_print_mip(LPX *lp, const char *fname);
/* write MIP problem solution in printable format */

#define lpx_is_b_avail _glp_lpx_is_b_avail
int lpx_is_b_avail(LPX *lp);
/* check if LP basis is available */

#define lpx_write_pb _glp_lpx_write_pb
int lpx_write_pb(LPX *lp, const char *fname, int normalized,
      int binarize);
/* write problem data in (normalized) OPB format */

#define lpx_main _glp_lpx_main
int lpx_main(int argc, const char *argv[]);
/* stand-alone LP/MIP solver */

#ifdef __cplusplus
}
#endif

#endif

/* eof */
