
#include "igraph_lsap.h"

#include "igraph_error.h"
#include "igraph_memory.h"

/* #include <stdio.h> */
#include <stdlib.h>
#include <float.h>      /* DBL_MAX */

/* constants used for improving readability of code */

typedef enum covered_t {
    COVERED = 1,
    UNCOVERED = 0
} covered_t;

typedef enum assigned_t {
    ASSIGNED = 1,
    UNASSIGNED = 0
} assigned_t;

typedef enum marked_t {
    MARKED = 1,
    UNMARKED = 0
} marked_t;

typedef enum reduce_t {
    REDUCE = 1,
    NOREDUCE = 0
} reduce_t;

typedef struct {
    igraph_integer_t    n;            /* order of problem             */
    double            **C;            /* cost matrix          */
    double            **c;            /* reduced cost matrix      */
    igraph_integer_t   *s;            /* assignment                   */
    igraph_integer_t   *f;            /* column i is assigned to f[i] */
    igraph_integer_t   na;            /* number of assigned items;    */
    igraph_integer_t runs;            /* number of iterations     */
    double           cost;            /* minimum cost         */
} AP;

/* public interface */

/* constructors and destructor */
static igraph_error_t ap_create_problem(AP **problem, const double *t, const igraph_integer_t n);
/* static AP     *ap_create_problem_from_matrix(double **t, int n); */
/* static AP     *ap_read_problem(char *file); */
static void    ap_free(AP *p);

static igraph_integer_t ap_get_result(AP *p, igraph_integer_t *res);
/* static int     ap_costmatrix(AP *p, double **m); */
/* static int     ap_datamatrix(AP *p, double **m); */
/* static int     ap_iterations(AP *p); */
static igraph_error_t ap_hungarian(AP *p);
/* static double  ap_mincost(AP *p); */
/* static void    ap_print_solution(AP *p); */
/* static void    ap_show_data(AP *p); */
/* static int     ap_size(AP *p); */
/* static int     ap_time(AP *p); */

/* error reporting */
/* static void ap_error(char *message); */

/* private functions */
static void    preprocess(AP *p);
static igraph_error_t preassign(AP *p);
static igraph_error_t cover(AP *p, covered_t *ri, covered_t *ci, reduce_t *res);
static void    reduce(AP *p, const covered_t *ri, const covered_t *ci);

/* Error message used on memory allocation failure. */
static const char *memerr = "Insufficient memory for LSAP.";

igraph_error_t ap_hungarian(AP *p) {
    covered_t *ri;  /* covered rows    */
    covered_t *ci;  /* covered columns */

    const igraph_integer_t n = p->n; /* size of problem */
    p->runs = 0;

    /* allocate memory */
    /* Note: p is already on the finally stack here. */
    p->s = IGRAPH_CALLOC(1 + n, igraph_integer_t);
    IGRAPH_CHECK_OOM(p->s, memerr);

    p->f = IGRAPH_CALLOC(1 + n, igraph_integer_t);
    IGRAPH_CHECK_OOM(p->f, memerr);

    ri = IGRAPH_CALLOC(1 + n, covered_t);
    IGRAPH_CHECK_OOM(ri, memerr);
    IGRAPH_FINALLY(igraph_free, ri);

    ci = IGRAPH_CALLOC(1 + n, covered_t);
    IGRAPH_CHECK_OOM(ci, memerr);
    IGRAPH_FINALLY(igraph_free, ci);

    preprocess(p);
    IGRAPH_CHECK(preassign(p));

    while (p->na < n) {
        reduce_t res;
        IGRAPH_CHECK(cover(p, ri, ci, &res));
        if (REDUCE == res) {
            reduce(p, ri, ci);
        }
        ++p->runs;
    }

    /* check if assignment is a permutation of (1..n) */
    for (igraph_integer_t i = 1; i <= n; i++) {
        igraph_integer_t ok = 0;
        for (igraph_integer_t j = 1; j <= n; j++)
            if (p->s[j] == i) {
                ++ok;
            }
        if (ok != 1)
            IGRAPH_ERROR("ap_hungarian: error in assignment, is not a permutation",
                         IGRAPH_EINVAL);
    }

    /* calculate cost of assignment */
    p->cost = 0;
    for (igraph_integer_t i = 1; i <= n; i++) {
        p->cost += p->C[i][p->s[i]];
    }

    /* reset result back to base-0 indexing */
    for (igraph_integer_t i = 1; i <= n; i++) {
        p->s[i - 1] = p->s[i] - 1;
    }

    /* free memory */

    IGRAPH_FREE(ri);
    IGRAPH_FREE(ci);
    IGRAPH_FINALLY_CLEAN(2);

    return IGRAPH_SUCCESS;
}

/* abbreviated interface */
igraph_integer_t ap_get_result(AP *p, igraph_integer_t *res) {
    for (igraph_integer_t i = 0; i < p->n; i++) {
        res[i] = p->s[i];
    }

    return p->n;
}


/*******************************************************************/
/* constructors                                                    */
/* read data from file                                             */
/*******************************************************************/

#if 0
AP *ap_read_problem(char *file) {
    FILE *f;
    int i, j, c;
    int m, n;
    double x;
    double **t;
    int nrow, ncol;
    AP *p;

    f = fopen(file, "r");
    if (f == NULL) {
        return NULL;
    }

    t = (double **)malloc(sizeof(double*));

    m = 0;
    n = 0;

    nrow = 0;
    ncol = 0;

    while (EOF != (i = fscanf(f, "%lf", &x))) {
        if (i == 1) {
            if (n == 0) {
                t = (double **) realloc(t, (m + 1) * sizeof(double *));
                t[m] = (double *) malloc(sizeof(double));
            } else {
                t[m] = (double *) realloc(t[m], (n + 1) * sizeof(double));
            }

            t[m][n++] = x;

            ncol = (ncol < n) ? n : ncol;
            c = fgetc(f);
            if (c == '\n') {
                n = 0;
                ++m;
                nrow = (nrow < m) ? m : nrow;
            }
        }
    }
    fclose(f);

    /* prepare data */

    if (nrow != ncol) {
        /*
          fprintf(stderr,"ap_read_problem: problem not quadratic\nrows =%d, cols = %d\n",nrow,ncol);
        */
        IGRAPH_WARNINGF("ap_read_problem: problem not quadratic; rows = %d, cols = %d.", nrow, ncol);
        return NULL;
    }

    p = (AP*) malloc(sizeof(AP));
    p->n = ncol;

    p->C  = (double **) malloc((1 + nrow) * sizeof(double *));
    p->c  = (double **) malloc((1 + nrow) * sizeof(double *));
    if (p->C == NULL || p->c == NULL) {
        return NULL;
    }

    for (i = 1; i <= nrow; i++) {
        p->C[i] = (double *) calloc(ncol + 1, sizeof(double));
        p->c[i] = (double *) calloc(ncol + 1, sizeof(double));
        if (p->C[i] == NULL || p->c[i] == NULL) {
            return NULL;
        }
    }

    for (i = 1; i <= nrow; i++)
        for ( j = 1; j <= ncol; j++) {
            p->C[i][j] = t[i - 1][j - 1];
            p->c[i][j] = t[i - 1][j - 1];
        }

    for (i = 0; i < nrow; i++) {
        free(t[i]);
    }
    free(t);

    p->cost = 0;
    p->s = NULL;
    p->f = NULL;
    return p;
}
#endif

#if 0
AP     *ap_create_problem_from_matrix(double **t, int n) {
    int i, j;
    AP *p;

    p = (AP*) malloc(sizeof(AP));
    if (p == NULL) {
        return NULL;
    }

    p->n = n;

    p->C  = (double **) malloc((n + 1) * sizeof(double *));
    p->c  = (double **) malloc((n + 1) * sizeof(double *));
    if (p->C == NULL || p->c == NULL) {
        return NULL;
    }

    for (i = 1; i <= n; i++) {
        p->C[i] = (double *) calloc(n + 1, sizeof(double));
        p->c[i] = (double *) calloc(n + 1, sizeof(double));
        if (p->C[i] == NULL || p->c[i] == NULL) {
            return NULL;
        }
    }


    for (i = 1; i <= n; i++)
        for ( j = 1; j <= n; j++) {
            p->C[i][j] = t[i - 1][j - 1];
            p->c[i][j] = t[i - 1][j - 1];
        }
    p->cost = 0;
    p->s = NULL;
    p->f = NULL;
    return p;
}
#endif

/* read data from vector */
igraph_error_t ap_create_problem(AP **problem, const double *t, const igraph_integer_t n) {
    *problem = IGRAPH_CALLOC(1, AP);
    IGRAPH_CHECK_OOM(*problem, memerr);
    IGRAPH_FINALLY(ap_free, *problem);

    AP *p = *problem;

    p->n = n;

    p->C = IGRAPH_CALLOC(n+1, double *);
    IGRAPH_CHECK_OOM(p->C, memerr);

    p->c = IGRAPH_CALLOC(n+1, double *);
    IGRAPH_CHECK_OOM(p->c, memerr);

    for (igraph_integer_t i = 1; i <= n; i++) {
        p->C[i] = IGRAPH_CALLOC(n+1, double);
        IGRAPH_CHECK_OOM(p->C[i], memerr);
        p->c[i] = IGRAPH_CALLOC(n+1, double);
        IGRAPH_CHECK_OOM(p->c[i], memerr);
    }

    for (igraph_integer_t i = 1; i <= n; i++) {
        for (igraph_integer_t j = 1; j <= n; j++) {
            p->C[i][j] = t[n * (j - 1) + i - 1];
            p->c[i][j] = t[n * (j - 1) + i - 1];
        }
    }
    p->cost = 0;
    p->s = NULL;
    p->f = NULL;

    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* destructor */
void ap_free(AP *p) {
    IGRAPH_FREE(p->s);
    IGRAPH_FREE(p->f);

    for (igraph_integer_t i = 1; i <= p->n; i++) {
        IGRAPH_FREE(p->C[i]);
        IGRAPH_FREE(p->c[i]);
    }

    IGRAPH_FREE(p->C);
    IGRAPH_FREE(p->c);
    IGRAPH_FREE(p);
}

/* set + get functions */

/*
void ap_show_data(AP *p)
{
    igraph_integer_t i, j;

    for(i = 1; i <= p->n; i++){
    for(j = 1; j <= p->n; j++)
        printf("%6.2f ", p->c[i][j]);
    printf("\n");
    }
}

double ap_mincost(AP *p) {
    if (p->s == NULL) {
        ap_hungarian(p);
    }

    return p->cost;
}

igraph_integer_t ap_size(AP *p) {
    return p->n;
}

int ap_time(AP *p) {
    return (int) p->rtime;
}

int ap_iterations(AP *p) {
    return p->runs;
}

void ap_print_solution(AP *p)
{
    igraph_integer_t i;

    printf("%d itertations, %d secs.\n",p->runs, (int)p->rtime);
    printf("Min Cost: %10.4f\n",p->cost);

    for(i = 0; i < p->n; i++)
    printf("%4d",p->s[i]);
    printf("\n");
}

int ap_costmatrix(AP *p, double **m) {
    igraph_integer_t i, j;

    for (i = 0; i < p->n; i++)
        for (j = 0; j < p->n; j++) {
            m[i][j] = p->C[i + 1][j + 1];
        }

    return p->n;
}

int ap_datamatrix(AP *p, double **m) {
    igraph_integer_t i, j;

    for (i = 0; i < p->n; i++)
        for (j = 0; j < p->n; j++) {
            m[i][j] = p->c[i + 1][j + 1];
        }

    return p->n;
}
*/

/* error reporting */

/*
void ap_error(char *message)
{
    fprintf(stderr,"%s\n",message);
    exit(1);
}
*/

/*************************************************************/
/* these functions are used internally                       */
/* by ap_hungarian                                           */
/*************************************************************/

igraph_error_t cover(AP *p, covered_t *ri, covered_t *ci, reduce_t *res) {
    marked_t *mr;
    igraph_integer_t r;

    const igraph_integer_t n = p->n;

    mr = IGRAPH_CALLOC(1 + p->n, marked_t);
    IGRAPH_CHECK_OOM(mr, memerr);

    /* reset cover indices */
    for (igraph_integer_t i = 1; i <= n; i++) {
        if (p->s[i] == UNASSIGNED) {
            ri[i] = UNCOVERED;
            mr[i] = MARKED;
        } else {
            ri[i] = COVERED;
        }
        ci[i] = UNCOVERED;
    }

    while (true) {
        /* find marked row */
        r = 0;
        for (igraph_integer_t i = 1; i <= n; i++)
            if (mr[i] == MARKED) {
                r = i;
                break;
            }

        if (r == 0) {
            break;
        }
        for (igraph_integer_t i = 1; i <= n; i++)
            if (p->c[r][i] == 0 && ci[i] == UNCOVERED) {
                if (p->f[i]) {
                    ri[p->f[i]] = UNCOVERED;
                    mr[p->f[i]] = MARKED;
                    ci[i] = COVERED;
                } else {
                    if (p->s[r] == UNASSIGNED) {
                        ++p->na;
                    }

                    p->f[p->s[r]] = 0;
                    p->f[i] = r;
                    p->s[r] = i;

                    IGRAPH_FREE(mr);
                    *res = NOREDUCE;
                    return IGRAPH_SUCCESS;
                }
            }
        mr[r] = UNMARKED;
    }

    IGRAPH_FREE(mr);
    *res = REDUCE;
    return IGRAPH_SUCCESS;
}

void reduce(AP *p, const covered_t *ri, const covered_t *ci) {
    double min;
    const igraph_integer_t n = p->n;

    /* find minimum in uncovered c-matrix */
    min = DBL_MAX;
    for (igraph_integer_t i = 1; i <= n; i++)
        for (igraph_integer_t j = 1; j <= n; j++)
            if (ri[i] == UNCOVERED && ci[j] == UNCOVERED) {
                if (p->c[i][j] < min) {
                    min = p->c[i][j];
                }
            }

    /* subtract min from each uncovered element and add it to each element */
    /* which is covered twice                                              */
    for (igraph_integer_t i = 1; i <= n; i++)
        for (igraph_integer_t j = 1; j <= n; j++) {
            if (ri[i] == UNCOVERED && ci[j] == UNCOVERED) {
                p->c[i][j] -= min;
            }
            if (ri[i] == COVERED && ci[j] == COVERED) {
                p->c[i][j] += min;
            }
        }
}

igraph_error_t preassign(AP *p) {
    igraph_integer_t min, r, c, n, count;
    assigned_t *ri, *ci;
    igraph_integer_t *rz, *cz;

    n = p->n;
    p->na = 0;

    /* row and column markers */
    ri = IGRAPH_CALLOC(1 + n, assigned_t);
    IGRAPH_CHECK_OOM(ri, memerr);
    IGRAPH_FINALLY(igraph_free, ri);

    ci = IGRAPH_CALLOC(1 + n, assigned_t);
    IGRAPH_CHECK_OOM(ci, memerr);
    IGRAPH_FINALLY(igraph_free, ci);

    /* row and column counts of zeroes */
    rz = IGRAPH_CALLOC(1 + n, igraph_integer_t);
    IGRAPH_CHECK_OOM(rz, memerr);
    IGRAPH_FINALLY(igraph_free, rz);

    cz = IGRAPH_CALLOC(1 + n, igraph_integer_t);
    IGRAPH_CHECK_OOM(cz, memerr);
    IGRAPH_FINALLY(igraph_free, cz);

    for (igraph_integer_t i = 1; i <= n; i++) {
        count = 0;
        for (igraph_integer_t j = 1; j <= n; j++)
            if (p->c[i][j] == 0) {
                ++count;
            }
        rz[i] = count;
    }

    for (igraph_integer_t i = 1; i <= n; i++) {
        count = 0;
        for (igraph_integer_t j = 1; j <= n; j++)
            if (p->c[j][i] == 0) {
                ++count;
            }
        cz[i] = count;
    }

    while (true) {
        /* find unassigned row with least number of zeroes > 0 */
        min = IGRAPH_INTEGER_MAX;
        r = 0;
        for (igraph_integer_t i = 1; i <= n; i++)
            if (rz[i] > 0 && rz[i] < min && ri[i] == UNASSIGNED) {
                min = rz[i];
                r = i;
            }
        /* check if we are done */
        if (r == 0) {
            break;
        }

        /* find unassigned column in row r with least number of zeroes */
        c = 0;
        min = IGRAPH_INTEGER_MAX;
        for (igraph_integer_t i = 1; i <= n; i++)
            if (p->c[r][i] == 0 && cz[i] < min && ci[i] == UNASSIGNED) {
                min = cz[i];
                c = i;
            }

        if (c) {
            ++p->na;
            p->s[r] = c;
            p->f[c] = r;

            ri[r] = ASSIGNED;
            ci[c] = ASSIGNED;

            /* adjust zero counts */
            cz[c] = 0;
            for (igraph_integer_t i = 1; i <= n; i++)
                if (p->c[i][c] == 0) {
                    --rz[i];
                }
        }
    }

    /* free memory */
    IGRAPH_FREE(ri);
    IGRAPH_FREE(ci);
    IGRAPH_FREE(rz);
    IGRAPH_FREE(cz);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}

void preprocess(AP *p) {
    double min;
    const igraph_integer_t n = p->n;

    /* subtract column minima in each row */
    for (igraph_integer_t i = 1; i <= n; i++) {
        min = p->c[i][1];
        for (igraph_integer_t j = 2; j <= n; j++)
            if (p->c[i][j] < min) {
                min = p->c[i][j];
            }
        for (igraph_integer_t j = 1; j <= n; j++) {
            p->c[i][j] -= min;
        }
    }

    /* subtract row minima in each column */
    for (igraph_integer_t i = 1; i <= n; i++) {
        min = p->c[1][i];
        for (igraph_integer_t j = 2; j <= n; j++)
            if (p->c[j][i] < min) {
                min = p->c[j][i];
            }
        for (igraph_integer_t j = 1; j <= n; j++) {
            p->c[j][i] -= min;
        }
    }
}

/**
 * \function igraph_solve_lsap
 * \brief Solve a balanced linear assignment problem.
 *
 * This functions solves a linear assignment problem using the Hungarian
 * method. A number of tasks, an equal number of agents, and the cost
 * of each agent to perform the tasks is given. This function then
 * assigns one task to each agent in such a way that the total cost is
 * minimized.
 *
 * </param><param>
 * If the cost should be maximized instead of minimized, the cost matrix
 * should be negated.
 *
 * </param><param>
 * To solve an unbalanced assignment problem, where the number of agents
 * is greater than the number of tasks, extra tasks with zero costs
 * should be added.
 *
 * \param c The assignment problem, where the number of rows is the
 *          number of agents, the number of columns is the number of
 *          tasks, and each element is the cost of an agent to perform
 *          the task.
 * \param n The number of rows and columns of \p c.
 * \param p An initialized vector which will store the result. The nth
 *          entry gives the task the nth agent is assigned to minimize
 *          the total cost.
 * \return Error code.
 *
 * Time complexity: O(n^3), where n is the number of agents.
 */
igraph_error_t igraph_solve_lsap(const igraph_matrix_t *c, igraph_integer_t n,
                      igraph_vector_int_t *p) {
    AP *ap;

    if (n != igraph_matrix_nrow(c)) {
        IGRAPH_ERRORF("n (%" IGRAPH_PRId ") "
                      "not equal to number of agents (%" IGRAPH_PRId ").", IGRAPH_EINVAL,
                      n, igraph_matrix_nrow(c));
    }
    if (n != igraph_matrix_ncol(c)) {
        IGRAPH_ERRORF("n (%" IGRAPH_PRId ") "
                      "not equal to number of tasks (%" IGRAPH_PRId ").", IGRAPH_EINVAL,
                      n, igraph_matrix_ncol(c));
    }
    IGRAPH_CHECK(igraph_vector_int_resize(p, n));
    igraph_vector_int_null(p);

    IGRAPH_CHECK(ap_create_problem(&ap, &MATRIX(*c, 0, 0), n));
    IGRAPH_FINALLY(ap_free, ap);
    IGRAPH_CHECK(ap_hungarian(ap));
    ap_get_result(ap, VECTOR(*p));
    ap_free(ap);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
