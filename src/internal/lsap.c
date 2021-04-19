
#include "igraph_lsap.h"
#include "igraph_error.h"

/* #include <stdio.h> */
#include <stdlib.h>
#include <limits.h>     /* INT_MAX */
#include <float.h>      /* DBL_MAX */
#include <time.h>

/* constants used for improving readability of code */

#define COVERED       1
#define UNCOVERED     0
#define ASSIGNED      1
#define UNASSIGNED    0
#define TRUE          1
#define FALSE         0

#define MARKED        1
#define UNMARKED      0

#define REDUCE        1
#define NOREDUCE      0

typedef struct {
    int        n;            /* order of problem             */
    double   **C;            /* cost matrix          */
    double   **c;            /* reduced cost matrix      */
    int       *s;            /* assignment                   */
    int       *f;            /* column i is assigned to f[i] */
    int       na;            /* number of assigned items;    */
    int     runs;            /* number of iterations     */
    double  cost;            /* minimum cost         */
    time_t rtime;            /* time                         */
} AP;

/* public interface */

/* constructors and destructor */
static AP     *ap_create_problem(double *t, int n);
/* static AP     *ap_create_problem_from_matrix(double **t, int n); */
/* static AP     *ap_read_problem(char *file); */
static void    ap_free(AP *p);

static int     ap_assignment(AP *p, int *res);
/* static int     ap_costmatrix(AP *p, double **m); */
/* static int     ap_datamatrix(AP *p, double **m); */
/* static int     ap_iterations(AP *p); */
static int     ap_hungarian(AP *p);
/* static double  ap_mincost(AP *p); */
/* static void    ap_print_solution(AP *p); */
/* static void    ap_show_data(AP *p); */
/* static int     ap_size(AP *p); */
/* static int     ap_time(AP *p); */

/* error reporting */
/* static void ap_error(char *message); */

/* private functions */
static void    preprocess(AP *p);
static void    preassign(AP *p);
static int     cover(AP *p, int *ri, int *ci);
static void    reduce(AP *p, int *ri, int *ci);

int ap_hungarian(AP *p) {
    int      n;            /* size of problem */
    int    *ri;            /* covered rows    */
    int    *ci;            /* covered columns */
    time_t start, end;     /* timer           */
    int i, j, ok;

    start = time(0);

    n = p->n;
    p->runs = 0;

    /* allocate memory */
    p->s = calloc(1 + n, sizeof(int));
    p->f = calloc(1 + n, sizeof(int));

    ri = calloc(1 + n, sizeof(int));
    ci = calloc(1 + n, sizeof(int));

    if (ri == NULL || ci == NULL || p->s == NULL || p->f == NULL) {
        IGRAPH_ERROR("ap_hungarian: could not allocate memory", IGRAPH_ENOMEM);
    }

    preprocess(p);
    preassign(p);

    while (p->na < n) {
        if (REDUCE == cover(p, ri, ci)) {
            reduce(p, ri, ci);
        }
        ++p->runs;
    }

    end = time(0);

    p->rtime = end - start;

    /* check if assignment is a permutation of (1..n) */
    for (i = 1; i <= n; i++) {
        ok = 0;
        for (j = 1; j <= n; j++)
            if (p->s[j] == i) {
                ++ok;
            }
        if (ok != 1)
            IGRAPH_ERROR("ap_hungarian: error in assignment, is not a permutation",
                         IGRAPH_EINVAL);
    }

    /* calculate cost of assignment */
    p->cost = 0;
    for (i = 1; i <= n; i++) {
        p->cost += p->C[i][p->s[i]];
    }

    /* reset result back to base-0 indexing */
    for (i = 1; i <= n; i++) {
        p->s[i - 1] = p->s[i] - 1;
    }

    /* free memory */

    free(ri);
    free(ci);

    return 0;
}

/* abbreviated interface */
int ap_assignment(AP *p, int *res) {
    int i;

    if (p->s == NULL) {
        ap_hungarian(p);
    }

    for (i = 0; i < p->n; i++) {
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
AP *ap_create_problem(double *t, int n) {
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
            p->C[i][j] = t[n * (j - 1) + i - 1];
            p->c[i][j] = t[n * (j - 1) + i - 1];
        }
    p->cost = 0;
    p->s = NULL;
    p->f = NULL;
    return p;
}

/* destructor */
void ap_free(AP *p) {
    int i;

    free(p->s);
    free(p->f);

    for (i = 1; i <= p->n; i++) {
        free(p->C[i]);
        free(p->c[i]);
    }

    free(p->C);
    free(p->c);
    free(p);
}

/* set + get functions */

/*
void ap_show_data(AP *p)
{
    int i, j;

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

int ap_size(AP *p) {
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
    int i;

    printf("%d itertations, %d secs.\n",p->runs, (int)p->rtime);
    printf("Min Cost: %10.4f\n",p->cost);

    for(i = 0; i < p->n; i++)
    printf("%4d",p->s[i]);
    printf("\n");
}

int ap_costmatrix(AP *p, double **m) {
    int i, j;

    for (i = 0; i < p->n; i++)
        for (j = 0; j < p->n; j++) {
            m[i][j] = p->C[i + 1][j + 1];
        }

    return p->n;
}

int ap_datamatrix(AP *p, double **m) {
    int i, j;

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

int cover(AP *p, int *ri, int *ci) {
    int *mr, i, r;
    int n;

    n = p->n;
    mr = calloc(1 + p->n, sizeof(int));

    /* reset cover indices */
    for (i = 1; i <= n; i++) {
        if (p->s[i] == UNASSIGNED) {
            ri[i] = UNCOVERED;
            mr[i] = MARKED;
        } else {
            ri[i] = COVERED;
        }
        ci[i] = UNCOVERED;
    }

    while (TRUE) {
        /* find marked row */
        r = 0;
        for (i = 1; i <= n; i++)
            if (mr[i] == MARKED) {
                r = i;
                break;
            }

        if (r == 0) {
            break;
        }
        for (i = 1; i <= n; i++)
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

                    free(mr);
                    return NOREDUCE;
                }
            }
        mr[r] = UNMARKED;
    }
    free(mr);
    return REDUCE;
}

void reduce(AP *p, int *ri, int *ci) {
    int i, j, n;
    double min;

    n = p->n;

    /* find minimum in uncovered c-matrix */
    min = DBL_MAX;
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++)
            if (ri[i] == UNCOVERED && ci[j] == UNCOVERED) {
                if (p->c[i][j] < min) {
                    min = p->c[i][j];
                }
            }

    /* subtract min from each uncovered element and add it to each element */
    /* which is covered twice                                              */
    for (i = 1; i <= n; i++)
        for (j = 1; j <= n; j++) {
            if (ri[i] == UNCOVERED && ci[j] == UNCOVERED) {
                p->c[i][j] -= min;
            }
            if (ri[i] == COVERED && ci[j] == COVERED) {
                p->c[i][j] += min;
            }
        }
}

void preassign(AP *p) {
    int i, j, min, r, c, n, count;
    int *ri, *ci, *rz, *cz;

    n = p->n;
    p->na = 0;

    /* row and column markers */
    ri = calloc(1 + n, sizeof(int));
    ci = calloc(1 + n, sizeof(int));

    /* row and column counts of zeroes */
    rz = calloc(1 + n, sizeof(int));
    cz = calloc(1 + n, sizeof(int));

    for (i = 1; i <= n; i++) {
        count = 0;
        for (j = 1; j <= n; j++)
            if (p->c[i][j] == 0) {
                ++count;
            }
        rz[i] = count;
    }

    for (i = 1; i <= n; i++) {
        count = 0;
        for (j = 1; j <= n; j++)
            if (p->c[j][i] == 0) {
                ++count;
            }
        cz[i] = count;
    }

    while (TRUE) {
        /* find unassigned row with least number of zeroes > 0 */
        min = INT_MAX;
        r = 0;
        for (i = 1; i <= n; i++)
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
        min = INT_MAX;
        for (i = 1; i <= n; i++)
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
            for (i = 1; i <= n; i++)
                if (p->c[i][c] == 0) {
                    --rz[i];
                }
        }
    }

    /* free memory */
    free(ri);
    free(ci);
    free(rz);
    free(cz);
}

void preprocess(AP *p) {
    int i, j, n;
    double min;

    n = p->n;

    /* subtract column minima in each row */
    for (i = 1; i <= n; i++) {
        min = p->c[i][1];
        for (j = 2; j <= n; j++)
            if (p->c[i][j] < min) {
                min = p->c[i][j];
            }
        for (j = 1; j <= n; j++) {
            p->c[i][j] -= min;
        }
    }

    /* subtract row minima in each column */
    for (i = 1; i <= n; i++) {
        min = p->c[1][i];
        for (j = 2; j <= n; j++)
            if (p->c[j][i] < min) {
                min = p->c[j][i];
            }
        for (j = 1; j <= n; j++) {
            p->c[j][i] -= min;
        }
    }
}

/**
 * \function igraph_solve_lsap
 * \brief Solve a balanced linear assignment problem.
 *
 * This functions solves a linear assinment problem using the Hungarian
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
 * is greater than the number of tasks, an extra task with zero cost
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
int igraph_solve_lsap(igraph_matrix_t *c, igraph_integer_t n,
                      igraph_vector_int_t *p) {
    AP *ap;

    if(n != igraph_matrix_nrow(c)) {
        IGRAPH_ERRORF("n (%" IGRAPH_PRId ") "
                      "not equal to number of agents (%ld).", IGRAPH_EINVAL,
                      n, igraph_matrix_nrow(c));
    }
    if(n != igraph_matrix_ncol(c)) {
        IGRAPH_ERRORF("n (%" IGRAPH_PRId ") "
                      "not equal to number of tasks (%ld).", IGRAPH_EINVAL,
                      n, igraph_matrix_ncol(c));
    }
    IGRAPH_CHECK(igraph_vector_int_resize(p, n));
    igraph_vector_int_null(p);

    ap = ap_create_problem(&MATRIX(*c, 0, 0), n);
    ap_hungarian(ap);
    ap_assignment(ap, VECTOR(*p));
    ap_free(ap);

    return IGRAPH_SUCCESS;
}
