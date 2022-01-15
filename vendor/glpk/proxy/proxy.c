/* proxy.c (proximity search heuristic algorithm) */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2013, 2016 Free Software Foundation, Inc.
*  Written by Giorgio Sartor <0gioker0@gmail.com>.
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
*
************************************************************************
*
* THIS CODE IS AN IMPLEMENTATION OF THE ALGORITHM PROPOSED IN
*
* M. Fischetti, M. Monaci,
* "Proximity Search for 0-1 Mixed-Integer Convex Programming"
* Technical Report DEI, University of Padua, March 2013.
*
* AVAILABLE AT
*       http://www.dei.unipd.it/~fisch/papers/proximity_search.pdf
*
* THE CODE HAS BEEN WRITTEN BY GIORGIO SARTOR, " 0gioker0@gmail.com "
*
* BASIC IDEA:
*
* The initial feasible solution x_tilde is defined. This initial
* solution can be found by an ad-hoc heuristic and proxy can be used to
* refine it by exploiting an underlying MIP model whose solution from
* scratch turned out to be problematic. Otherwise, x_tilde can be found
* by running the GLPK mip solver until a first feasible solution is
* found, setting a conservative time limit of 10 minutes (by default).
* Time limit can be modified passing variable tlim [ms].
*
* Then the cutoff tolerance "delta" is defined. The default tolerance
* is 1% of the last feasible solution obj value--rounded to integer if
* all the variables and obj coefficients are integer.
*
* Next, the objective function c' x is replaced by the Hamming distance
* between x (the actual obj coefficients) and x_tilde (the given
* solution). Distance is only computed wrt the binary variables.
*
* The GLPK solver is then invoked to hopefully find a new incumbent
* x_star with cost c' x_star <= c' x_tilde - delta. A crucial property
* here is that the root-node solution of the LP relaxation is expected
* to be not too different from x_tilde, as this latter solution would
* be optimal without the cutoff constraint, that for a small delta can
* typically be fulfilled with just local adjustments.
*
* If no new solution x_star is found within the time limit the
* algorithm stops. Of course, if the MIP solver proved infeasibility
* for the given delta, we have that c' x_tilde - delta is a valid lower
* bound (in case of minimazation) on the optimal value of the original
* MIP.
*
* The new solution x_star, if any, is possibly improved by solving a
* simplified problem (refinement) where all binary variables have been
* fixed to their value in x_star so as to find the best solution within
* the neighborhood.
*
* Finally, the approach is reapplied on x_star (that replaces x_tilde)
* so as to recenter the distance Hamming function and by modifying the
* cutoff tolerance delta.
*
* In this way, there will be a series of hopefully not-too-difficult
* sub-MIPs to solve, each leading to an improvement of the incumbent.
* More aggressive policies on the definition of tolerance delta can
* lead to a better performance, but would require an ad-hoc tuning.
*
************************************************************************
*
* int proxy(glp_prob *lp, double *zstar, double *xstar,
*           const double[] initsol, double rel_impr, int tlim,
*           int verbose)
*
* lp       : GLPK problem pointer to a MIP with binary variables
*
* zstar    : the value of objective function of the best solution found
*
* xstar    : best solution with components xstar[1],...,xstar[ncols]
*
* initsol  : pointer to a initial feasible solution, see
*            glp_ios_heur_sol
*            If initsol = NULL, the procedure finds the first solution
*            by itself.
*
* rel_impr : minimum relative obj improvement to be achieved at each
*            internal step; if <= 0.0 a default value of 0.01 (1%) is
*            used; for some problems (e.g., set covering with small
*            integer costs) a more-conservative choice of 0.001 (0.1%)
*            can lead to a better final solution; values larger than
*            0.05 (5%) are typically too aggressive and do not work
*            well.
*
* tlim     : time limit to find a new solution, in ms.
*            If tlim = 0, it is set to its default value, 600000 ms
*
* verbose  : if 1 the output is activated. If 0 only errors are
*            displayed
*
* The procedure returns -1 if an error occurred, 0 otherwise (possibly,
* time limit)
*
***********************************************************************/

/**********************************************************************/
/* 1. INCLUDE                                                         */
/**********************************************************************/

#include "glpk.h"
#include "env.h"
#include "proxy.h"

/**********************************************************************/
/* 2. PARAMETERS AND CONSTANTS                                        */
/**********************************************************************/

#define TDAY            86400.0
#define TRUE                1
#define FALSE               0
#define EPS              1e-6
#define RINF             1e38
#define MAXVAL           1e20
#define MINVAL          -1e20
#if 0 /* by gioker */
    #define PROXY_DEBUG
#endif

/**********************************************************************/
/* 3. GLOBAL VARIABLES                                                */
/**********************************************************************/

struct csa {

int integer_obj;        /* TRUE if each feasible solution has an
                           integral cost */
int b_vars_exist;       /* TRUE if there is at least one binary
                           variable in the problem */
int i_vars_exist;       /* TRUE if there is at least one general
                           integer variable in the problem */
const double *startsol; /* Pointer to the initial solution */

int *ckind;             /* Store the kind of the structural variables
                           of the problem */
double *clb;            /* Store the lower bound on the structural
                           variables of the problem */
double *cub;            /* Store the upper bound on the structural
                           variables of the problem */
double *true_obj;       /* Store the obj coefficients of the problem */

int dir;                /* Minimization or maximization problem */
int ncols;              /* Number of structural variables of the
                           problem */

time_t GLOtstart;       /* starting time of the algorithm */

glp_prob *lp_ref;       /* glp problem for refining only*/

};

/**********************************************************************/
/* 4. FUNCTIONS PROTOTYPES                                            */
/**********************************************************************/

static void callback(glp_tree *tree, void *info);
static void get_info(struct csa *csa, glp_prob *lp);
static int is_integer(struct csa *csa);
static void check_integrality(struct csa *csa);
static int check_ref(struct csa *csa, glp_prob *lp, double *xref);
static double second(void);
static int add_cutoff(struct csa *csa, glp_prob *lp);
static void get_sol(struct csa *csa, glp_prob *lp, double *xstar);
static double elapsed_time(struct csa *csa);
static void redefine_obj(glp_prob *lp, double *xtilde, int ncols,
                         int *ckind, double *clb, double *cub);
static double update_cutoff(struct csa *csa, glp_prob *lp,
                            double zstar, int index, double rel_impr);
static double compute_delta(struct csa *csa, double z,
                            double rel_impr);
static double objval(int ncols, double *x, double *true_obj);
static void array_copy(int begin, int end, double *source,
                       double *destination);
static int do_refine(struct csa *csa, glp_prob *lp_ref, int ncols,
                     int *ckind, double *xref, int *tlim, int tref_lim,
                     int verbose);
static void deallocate(struct csa *csa, int refine);

/**********************************************************************/
/* 5. FUNCTIONS                                                       */
/**********************************************************************/

int proxy(glp_prob *lp, double *zfinal, double *xfinal,
          const double initsol[], double rel_impr, int tlim,
          int verbose)

{   struct csa csa_, *csa = &csa_;
    glp_iocp parm;
    glp_smcp parm_lp;
    size_t tpeak;
    int refine, tref_lim, err, cutoff_row, niter, status, i, tout;
    double *xref, *xstar, zstar, tela, cutoff, zz;

    memset(csa, 0, sizeof(struct csa));


    /**********                         **********/
    /********** RETRIEVING PROBLEM INFO **********/
    /**********                         **********/

    /* getting problem direction (min or max) */
    csa->dir = glp_get_obj_dir(lp);

    /* getting number of variables */
    csa->ncols = glp_get_num_cols(lp);

    /* getting kind, bounds and obj coefficient of each variable
     information is stored in ckind, cub, clb, true_obj */
    get_info(csa, lp);

    /* checking if the objective function is always integral */
    check_integrality(csa);

    /* Proximity search cannot be used if there are no binary
       variables */
    if (csa->b_vars_exist == FALSE) {
        if (verbose) {
            xprintf("The problem has not binary variables. Proximity se"
                    "arch cannot be used.\n");
        }
        tfree(csa->ckind);
        tfree(csa->clb);
        tfree(csa->cub);
        tfree(csa->true_obj);
        return -1;
    }

    /* checking if the problem needs refinement, i.e., not all
       variables are binary. If so, the routine creates a copy of the
       lp problem named lp_ref and initializes the solution xref to
       zero. */
    xref = talloc(csa->ncols+1, double);
#if 0 /* by mao */
    memset(xref, 0, sizeof(double)*(csa->ncols+1));
#endif
    refine = check_ref(csa, lp, xref);
#ifdef PROXY_DEBUG
    xprintf("REFINE = %d\n",refine);
#endif

    /* Initializing the solution */
    xstar = talloc(csa->ncols+1, double);
#if 0 /* by mao */
    memset(xstar, 0, sizeof(double)*(csa->ncols+1));
#endif

    /**********                         **********/
    /********** FINDING FIRST SOLUTION  **********/
    /**********                         **********/

    if (verbose) {
        xprintf("Applying PROXY heuristic...\n");
    }

    /* get the initial time */
    csa->GLOtstart = second();

    /* setting the optimization parameters */
    glp_init_iocp(&parm);
    glp_init_smcp(&parm_lp);
#if 0 /* by gioker */
    /* Preprocessing should be disabled because the mip passed
     to proxy is already preprocessed */
    parm.presolve = GLP_ON;
#endif
#if 1 /* by mao */
    /* best projection backtracking seems to be more efficient to find
       any integer feasible solution */
    parm.bt_tech = GLP_BT_BPH;
#endif

    /* Setting the default value of the minimum relative improvement
       to 1% */
    if ( rel_impr <= 0.0 ) {
        rel_impr = 0.01;
    }

    /* Setting the default value of time limit to 10 minutes */
    if (tlim <= 0) {
        tlim = INT_MAX;
    }
    if (verbose) {
        xprintf("Proxy's time limit set to %d seconds.\n",tlim/1000);
        xprintf("Proxy's relative improvement "
                "set to %2.2lf %c.\n",rel_impr*100,37);
    }

    parm_lp.tm_lim = tlim;

    parm.mip_gap = 9999999.9; /* to stop the optimization at the first
                                 feasible solution found */

    /* finding the first solution */
    if (verbose) {
        xprintf("Searching for a feasible solution...\n");
    }

    /* verifying the existence of an input starting solution */
    if (initsol != NULL) {
        csa->startsol = initsol;
        parm.cb_func = callback;
        parm.cb_info = csa;
        if (verbose) {
            xprintf("Input solution found.\n");
        }
    }

    tout = glp_term_out(GLP_OFF);
    err = glp_simplex(lp,&parm_lp);
    glp_term_out(tout);

    status = glp_get_status(lp);

    if (status != GLP_OPT) {
        if (verbose) {
            xprintf("Proxy heuristic terminated.\n");
        }
#ifdef  PROXY_DEBUG
        /* For debug only */
        xprintf("GLP_SIMPLEX status = %d\n",status);
        xprintf("GLP_SIMPLEX error code = %d\n",err);
#endif
        tfree(xref);
        tfree(xstar);
        deallocate(csa, refine);
        return -1;
    }

    tela = elapsed_time(csa);
    if (tlim-tela*1000 <= 0) {
        if (verbose) {
            xprintf("Time limit exceeded. Proxy could not "
                    "find optimal solution to LP relaxation.\n");
            xprintf("Proxy heuristic aborted.\n");
        }
        tfree(xref);
        tfree(xstar);
        deallocate(csa, refine);
        return -1;
    }

    parm.tm_lim = tlim - tela*1000;
    tref_lim = (tlim - tela *1000) / 20;

    tout = glp_term_out(GLP_OFF);
    err = glp_intopt(lp, &parm);
    glp_term_out(tout);

    status = glp_mip_status(lp);

    /***** If no solution was found *****/

    if (status == GLP_NOFEAS || status == GLP_UNDEF) {
        if (err == GLP_ETMLIM) {
            if (verbose) {
                xprintf("Time limit exceeded. Proxy could not "
                        "find an initial integer feasible solution.\n");
                xprintf("Proxy heuristic aborted.\n");
            }
        }
        else {
            if (verbose) {
                xprintf("Proxy could not "
                        "find an initial integer feasible solution.\n");
                xprintf("Proxy heuristic aborted.\n");
            }
        }
        tfree(xref);
        tfree(xstar);
        deallocate(csa, refine);
        return -1;
    }

    /* getting the first solution and its value */
    get_sol(csa, lp,xstar);
    zstar = glp_mip_obj_val(lp);

    if (verbose) {
        xprintf(">>>>> first solution = %e;\n", zstar);
    }

    /* If a feasible solution was found but the time limit is
       exceeded */
    if (err == GLP_ETMLIM) {
        if (verbose) {
          xprintf("Time limit exceeded. Proxy heuristic terminated.\n");
        }
        goto done;
    }

    tela = elapsed_time(csa);
    tpeak = 0;
    glp_mem_usage(NULL, NULL, NULL, &tpeak);
    if (verbose) {
        xprintf("Time used: %3.1lf secs.  Memory used: %2.1lf Mb\n",
                tela,(double)tpeak/1048576);
        xprintf("Starting proximity search...\n");
    }

    /**********                                 **********/
    /********** PREPARING THE PROBLEM FOR PROXY **********/
    /**********                                 **********/

    /* adding a dummy cutoff constraint */
    cutoff_row = add_cutoff(csa, lp);

    /* proximity search needs minimization direction
       even if the problem is a maximization one */
    if (csa->dir == GLP_MAX) {
        glp_set_obj_dir(lp, GLP_MIN);
    }

    /**********                           **********/
    /********** STARTING PROXIMITY SEARCH **********/
    /**********                           **********/


    niter = 0;

    while (TRUE) {
        niter++;

        /********** CHANGING THE OBJ FUNCTION **********/

        redefine_obj(lp,xstar, csa->ncols, csa->ckind, csa->clb,
                     csa->cub);

        /********** UPDATING THE CUTOFF CONSTRAINT **********/

        cutoff = update_cutoff(csa, lp,zstar, cutoff_row, rel_impr);

#ifdef PROXY_DEBUG
        xprintf("TRUE_OBJ[0] = %f\n",csa->true_obj[0]);
        xprintf("ZSTAR  = %f\n",zstar);
        xprintf("CUTOFF = %f\n",cutoff);
#endif

        /********** SEARCHING FOR A BETTER SOLUTION **********/

        tela = elapsed_time(csa);
        if (tlim-tela*1000 <= 0) {
            if (verbose) {
                xprintf("Time limit exceeded. Proxy heuristic "
                        "terminated.\n");
            }
            goto done;
        }
#ifdef PROXY_DEBUG
        xprintf("TELA = %3.1lf\n",tela*1000);
        xprintf("TLIM = %3.1lf\n",tlim - tela*1000);
#endif
        parm_lp.tm_lim = tlim -tela*1000;

        tout = glp_term_out(GLP_OFF);
        err = glp_simplex(lp,&parm_lp);
        glp_term_out(tout);

        status = glp_get_status(lp);

        if (status != GLP_OPT) {
            if (status == GLP_NOFEAS) {
                if (verbose) {
                    xprintf("Bound exceeded = %f. ",cutoff);
                }
            }
            if (verbose) {
                xprintf("Proxy heuristic terminated.\n");
            }
#ifdef PROXY_DEBUG
            xprintf("GLP_SIMPLEX status = %d\n",status);
            xprintf("GLP_SIMPLEX error code = %d\n",err);
#endif
            goto done;
        }

        tela = elapsed_time(csa);
        if (tlim-tela*1000 <= 0) {
            if (verbose) {
                xprintf("Time limit exceeded. Proxy heuristic "
                        "terminated.\n");
            }
            goto done;
        }
        parm.tm_lim = tlim - tela*1000;
        parm.cb_func = NULL;
#if 0 /* by gioker */
        /* Preprocessing should be disabled because the mip passed
         to proxy is already preprocessed */
        parm.presolve = GLP_ON;
#endif
        tout = glp_term_out(GLP_OFF);
        err = glp_intopt(lp, &parm);
        glp_term_out(tout);

        /********** MANAGEMENT OF THE SOLUTION **********/

        status = glp_mip_status(lp);

        /***** No feasible solutions *****/

        if (status == GLP_NOFEAS) {
            if (verbose) {
                xprintf("Bound exceeded = %f. Proxy heuristic "
                        "terminated.\n",cutoff);
            }
            goto done;
        }

        /***** Undefined solution *****/

        if (status == GLP_UNDEF) {
            if (err == GLP_ETMLIM) {
                if (verbose) {
                    xprintf("Time limit exceeded. Proxy heuristic "
                            "terminated.\n");
                }
            }
            else {
                if (verbose) {
                    xprintf("Proxy terminated unexpectedly.\n");
#ifdef PROXY_DEBUG
                    xprintf("GLP_INTOPT error code = %d\n",err);
#endif
                }
            }
            goto done;
        }

        /***** Feasible solution *****/

        if ((status == GLP_FEAS) || (status == GLP_OPT)) {

            /* getting the solution and computing its value */
            get_sol(csa, lp,xstar);
            zz = objval(csa->ncols, xstar, csa->true_obj);

            /* Comparing the incumbent solution with the current best
               one */
#ifdef PROXY_DEBUG
            xprintf("ZZ = %f\n",zz);
            xprintf("ZSTAR = %f\n",zstar);
            xprintf("REFINE = %d\n",refine);
#endif
            if (((zz<zstar) && (csa->dir == GLP_MIN)) ||
                ((zz>zstar) && (csa->dir == GLP_MAX))) {

                /* refining (possibly) the solution */
                if (refine) {

                    /* copying the incumbent solution in the refinement
                       one */
                    array_copy(1, csa->ncols +1, xstar, xref);
                    err = do_refine(csa, csa->lp_ref, csa->ncols,
                          csa->ckind, xref, &tlim, tref_lim, verbose);
                    if (!err) {
                        double zref = objval(csa->ncols, xref,
                                             csa->true_obj);
                        if (((zref<zz) && (csa->dir == GLP_MIN)) ||
                            ((zref>zz) && (csa->dir == GLP_MAX))) {
                            zz = zref;
                            /* copying the refinement solution in the
                               incumbent one */
                            array_copy(1, csa->ncols +1, xref, xstar);
                        }
                    }
                }
                zstar = zz;
                tela = elapsed_time(csa);
                if (verbose) {
                    xprintf(">>>>> it: %3d:   mip = %e;   elapsed time "
                            "%3.1lf sec.s\n", niter,zstar,tela);
                }
            }
        }
    }

done:
    tela = elapsed_time(csa);
    glp_mem_usage(NULL, NULL, NULL, &tpeak);
    if (verbose) {
        xprintf("Time used: %3.1lf.  Memory used: %2.1lf Mb\n",
                tela,(double)tpeak/1048576);
    }


    /* Exporting solution and obj val */
    *zfinal = zstar;

    for (i=1; i < (csa->ncols + 1); i++) {
        xfinal[i]=xstar[i];
    }

    /* Freeing allocated memory */
    tfree(xref);
    tfree(xstar);
    deallocate(csa, refine);

    return 0;
}

/**********************************************************************/
static void callback(glp_tree *tree, void *info){
/**********************************************************************/
    struct csa *csa = info;
    switch(glp_ios_reason(tree)) {
        case GLP_IHEUR:
            glp_ios_heur_sol(tree, csa->startsol);
            break;
        default: break;
    }
}

/**********************************************************************/
static void get_info(struct csa *csa, glp_prob *lp)
/**********************************************************************/
{
    int i;

    /*  Storing helpful info of the problem  */

    csa->ckind = talloc(csa->ncols+1, int);
#if 0 /* by mao */
    memset(csa->ckind, 0, sizeof(int)*(csa->ncols+1));
#endif
    csa->clb = talloc(csa->ncols+1, double);
#if 0 /* by mao */
    memset(csa->clb, 0, sizeof(double)*(csa->ncols+1));
#endif
    csa->cub = talloc(csa->ncols+1, double);
#if 0 /* by mao */
    memset(csa->cub, 0, sizeof(double)*(csa->ncols+1));
#endif
    csa->true_obj = talloc(csa->ncols+1, double);
#if 0 /* by mao */
    memset(csa->true_obj, 0, sizeof(double)*(csa->ncols+1));
#endif
        for( i = 1 ; i < (csa->ncols + 1); i++ ) {
            csa->ckind[i] = glp_get_col_kind(lp, i);
            csa->clb[i] = glp_get_col_lb(lp, i);
            csa->cub[i] = glp_get_col_ub(lp, i);
            csa->true_obj[i] = glp_get_obj_coef(lp, i);
        }
    csa->true_obj[0] = glp_get_obj_coef(lp, 0);
}

/**********************************************************************/
static int is_integer(struct csa *csa)
/**********************************************************************/
{
    int i;
    csa->integer_obj = TRUE;
    for ( i = 1; i < (csa->ncols + 1); i++ ) {
        if (fabs(csa->true_obj[i]) > INT_MAX ) {
            csa->integer_obj = FALSE;
        }
        if (fabs(csa->true_obj[i]) <= INT_MAX) {
            double tmp, rem;
            if (fabs(csa->true_obj[i]) - floor(fabs(csa->true_obj[i]))
                < 0.5) {
                tmp = floor(fabs(csa->true_obj[i]));
            }
            else {
                tmp = ceil(fabs(csa->true_obj[i]));
            }
            rem = fabs(csa->true_obj[i]) - tmp;
            rem = fabs(rem);
            if (rem > EPS) {
                csa->integer_obj = FALSE;
            }

        }
    }
    return csa->integer_obj;
}

/**********************************************************************/
static void check_integrality(struct csa *csa)
/**********************************************************************/
{
    /*
     Checking if the problem has binary, integer or continuos variables.
     integer_obj is TRUE if the problem has no continuous variables
     and all the obj coefficients are integer (and < INT_MAX).
     */

    int i;
    csa->integer_obj = is_integer(csa);
    csa->b_vars_exist = FALSE;
    csa->i_vars_exist = FALSE;
    for ( i = 1; i < (csa->ncols + 1); i++ ) {
        if ( csa->ckind[i] == GLP_IV ){
            csa->i_vars_exist = TRUE;
            continue;
        }
        if ( csa->ckind[i] == GLP_BV ){
            csa->b_vars_exist =TRUE;
            continue;
        }
        csa->integer_obj = FALSE;
    }
}

/**********************************************************************/
static int check_ref(struct csa *csa, glp_prob *lp, double *xref)
/**********************************************************************/
{
    /*
     checking if the problem has continuos or integer variables. If so,
     refinement is prepared.
     */
    int refine = FALSE;
    int i;
    for ( i = 1; i < (csa->ncols + 1); i++ ) {
        if ( csa->ckind[i] != GLP_BV) {
            refine = TRUE;
            break;
        }
    }

    /* possibly creating a mip clone for refinement only */
    if ( refine ) {
        csa->lp_ref = glp_create_prob();
        glp_copy_prob(csa->lp_ref, lp, GLP_ON);
    }

    return refine;
}

/**********************************************************************/
static double second(void)
/**********************************************************************/
{
#if 0 /* by mao */
    return ((double)clock()/(double)CLOCKS_PER_SEC);
#else
    return xtime() / 1000.0;
#endif
}

/**********************************************************************/
static int add_cutoff(struct csa *csa, glp_prob *lp)
/**********************************************************************/
{
    /*
     Adding a cutoff constraint to set an upper bound (in case of
     minimaztion) on the obj value of the next solution, i.e., the next
     value of the true obj function that we would like to find
     */

    /* store non-zero coefficients in the objective function */
    int *obj_index = talloc(csa->ncols+1, int);
#if 0 /* by mao */
    memset(obj_index, 0, sizeof(int)*(csa->ncols+1));
#endif
    double *obj_value = talloc(csa->ncols+1, double);
#if 0 /* by mao */
    memset(obj_value, 0, sizeof(double)*(csa->ncols+1));
#endif
    int obj_nzcnt = 0;
    int i, irow;
    const char *rowname;
    for ( i = 1; i < (csa->ncols + 1); i++ ) {
        if ( fabs(csa->true_obj[i]) > EPS ) {
            obj_nzcnt++;
            obj_index[obj_nzcnt] = i;
            obj_value[obj_nzcnt] = csa->true_obj[i];
        }
    }

    irow = glp_add_rows(lp, 1);
    rowname = "Cutoff";
    glp_set_row_name(lp, irow, rowname);
    if (csa->dir == GLP_MIN) {
        /* minimization problem */
        glp_set_row_bnds(lp, irow, GLP_UP, MAXVAL, MAXVAL);
    }
    else {
        /* maximization problem */
        glp_set_row_bnds(lp, irow, GLP_LO, MINVAL, MINVAL);
    }

    glp_set_mat_row(lp, irow, obj_nzcnt, obj_index, obj_value);

    tfree(obj_index);
    tfree(obj_value);

    return irow;
}

/**********************************************************************/
static void get_sol(struct csa *csa, glp_prob *lp, double *xstar)
/**********************************************************************/
{
    /* Retrieving and storing the coefficients of the solution */

    int i;
    for (i = 1; i < (csa->ncols +1); i++) {
        xstar[i] = glp_mip_col_val(lp, i);
    }
}

/**********************************************************************/
static double elapsed_time(struct csa *csa)
/**********************************************************************/
{
    double tela = second() - csa->GLOtstart;
    if ( tela < 0 ) tela += TDAY;
    return(tela);
}

/**********************************************************************/
static void redefine_obj(glp_prob *lp, double *xtilde, int ncols,
                         int *ckind, double *clb, double *cub)
/**********************************************************************/

/*
 Redefine the lp objective function obj as the distance-to-integrality
 (Hamming distance) from xtilde (the incumbent feasible solution), wrt
 to binary vars only
 */

{
    int j;
    double *delta = talloc(ncols+1, double);
#if 0 /* by mao */
    memset(delta, 0, sizeof(double)*(ncols+1));
#endif

    for ( j = 1; j < (ncols +1); j++ ) {
        delta[j] = 0.0;
        /* skip continuous variables */
        if ( ckind[j] == GLP_CV ) continue;

        /* skip integer variables that have been fixed */
        if ( cub[j]-clb[j] < 0.5 ) continue;

        /* binary variable */
        if ( ckind[j] == GLP_BV ) {
            if ( xtilde[j] > 0.5 ) {
                delta[j] = -1.0;
            }
            else {
                delta[j] = 1.0;
            }
        }
    }

    /* changing the obj coeff. for all variables, including continuous
       ones */
    for ( j = 1; j < (ncols +1); j++ ) {
        glp_set_obj_coef(lp, j, delta[j]);
    }
    glp_set_obj_coef(lp, 0, 0.0);

    tfree(delta);
}

/**********************************************************************/
static double update_cutoff(struct csa *csa, glp_prob *lp,
                            double zstar, int cutoff_row,
                            double rel_impr)
/**********************************************************************/
{
    /*
     Updating the cutoff constraint with the value we would like to
     find during the next optimization
     */
    double cutoff;
    zstar -= csa->true_obj[0];
    if (csa->dir == GLP_MIN) {
        cutoff = zstar - compute_delta(csa, zstar, rel_impr);
        glp_set_row_bnds(lp, cutoff_row, GLP_UP, cutoff, cutoff);
    }
    else {
        cutoff = zstar + compute_delta(csa, zstar, rel_impr);
        glp_set_row_bnds(lp, cutoff_row, GLP_LO, cutoff, cutoff);
    }

    return cutoff;
}

/**********************************************************************/
static double compute_delta(struct csa *csa, double z, double rel_impr)
/**********************************************************************/
{
    /* Computing the offset for the next best solution */

    double delta = rel_impr * fabs(z);
    if ( csa->integer_obj ) delta = ceil(delta);

    return(delta);
}

/**********************************************************************/
static double objval(int ncols, double *x, double *true_obj)
/**********************************************************************/
{
    /* Computing the true cost of x (using the original obj coeff.s) */

    int j;
    double z = 0.0;
    for ( j = 1; j < (ncols +1); j++ ) {
        z += x[j] * true_obj[j];
    }
    return z + true_obj[0];
}

/**********************************************************************/
static void array_copy(int begin, int end, double *source,
                       double *destination)
/**********************************************************************/
{
    int i;
    for (i = begin; i < end; i++) {
        destination[i] = source[i];
    }
}
/**********************************************************************/
static int do_refine(struct csa *csa, glp_prob *lp_ref, int ncols,
                     int *ckind, double *xref, int *tlim, int tref_lim,
                     int verbose)
/**********************************************************************/
{
    /*
     Refinement is applied when the variables of the problem are not
     all binary. Binary variables are fixed to their value and
     remaining ones are optimized. If there are only continuos
     variables (in addition to those binary) the problem becomes just
     an LP. Otherwise, it remains a MIP but of smaller size.
     */

    int j, tout;
    double refineStart = second();
    double val, tela, tlimit;

    if ( glp_get_num_cols(lp_ref) != ncols ) {
        if (verbose) {
            xprintf("Error in Proxy refinement: ");
            xprintf("wrong number of columns (%d vs %d).\n",
                    ncols, glp_get_num_cols(lp_ref));
        }
        return 1;
    }

    val = -1.0;

    /* fixing all binary variables to their current value in xref */
    for ( j = 1; j < (ncols + 1); j++ ) {
        if ( ckind[j] == GLP_BV ) {
            val = 0.0;
            if ( xref[j] > 0.5 ) val = 1.0;
            glp_set_col_bnds(lp_ref, j, GLP_FX, val, val);
        }
    }

    /* re-optimizing (refining) if some bound has been changed */
    if ( val > -1.0 ) {
        glp_iocp parm_ref;
        glp_smcp parm_ref_lp;
        int err, status;

        glp_init_iocp(&parm_ref);
        parm_ref.presolve = GLP_ON;
        glp_init_smcp(&parm_ref_lp);
        /*
         If there are no general integer variable the problem becomes
         an LP (after fixing the binary variables) and can be solved
         quickly. Otherwise the problem is still a MIP problem and a
         timelimit has to be set.
         */
        parm_ref.tm_lim = tref_lim;
        if (parm_ref.tm_lim > *tlim) {
            parm_ref.tm_lim = *tlim;
        }
        parm_ref_lp.tm_lim = parm_ref.tm_lim;
#ifdef PROXY_DEBUG
        xprintf("***** REFINING *****\n");
#endif
        tout = glp_term_out(GLP_OFF);
        if (csa->i_vars_exist == TRUE) {
            err = glp_intopt(lp_ref, &parm_ref);
        }
        else {
            err = glp_simplex(lp_ref, &parm_ref_lp);
        }
        glp_term_out(tout);

        if (csa->i_vars_exist == TRUE) {
            status = glp_mip_status(lp_ref);
        }
        else {
            status = glp_get_status(lp_ref);
        }

#if 1 /* 29/II-2016 by mao as reported by Chris */
      switch (status)
      {  case GLP_OPT:
         case GLP_FEAS:
            break;
         default:
            status = GLP_UNDEF;
            break;
      }
#endif

#ifdef PROXY_DEBUG
        xprintf("STATUS REFINING = %d\n",status);
#endif
        if (status == GLP_UNDEF) {
            if (err == GLP_ETMLIM) {
#ifdef PROXY_DEBUG
                    xprintf("Time limit exceeded on Proxy refining.\n");
#endif
                return 1;
            }
        }
        for( j = 1 ; j < (ncols + 1); j++ ){
            if (ckind[j] != GLP_BV) {
                if (csa->i_vars_exist == TRUE) {
                    xref[j] = glp_mip_col_val(lp_ref, j);
                }
                else{
                    xref[j] = glp_get_col_prim(lp_ref, j);
                }
            }
        }
    }
    tela = second() - refineStart;
#ifdef PROXY_DEBUG
    xprintf("REFINE TELA = %3.1lf\n",tela*1000);
#endif
    return 0;
}
/**********************************************************************/
static void deallocate(struct csa *csa, int refine)
/**********************************************************************/
{
    /* Deallocating routine */

    if (refine) {
        glp_delete_prob(csa->lp_ref);
    }

    tfree(csa->ckind);
    tfree(csa->clb);
    tfree(csa->cub);
    tfree(csa->true_obj);

}

/* eof */
