
#include "igraph_lsap.h"
#include "igraph_error.h"

/* #include <stdio.h> */
#include <stdlib.h>
#include <math.h>
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
    igraph_long_t        n;            /* order of problem             */
    double   **C;            /* cost matrix          */
    double   **c;            /* reduced cost matrix      */
    igraph_long_t       *s;            /* assignment                   */
    igraph_long_t       *f;            /* column i is assigned to f[i] */
    igraph_long_t       na;            /* number of assigned items;    */
    igraph_long_t     runs;            /* number of iterations     */
    double  cost;            /* minimum cost         */
    time_t rtime;            /* time                         */
} AP;

/* public interface */

/* constructors and destructor */
static AP     *ap_create_problem(double *t, igraph_long_t n);
/* static AP     *ap_create_problem_from_matrix(double **t, igraph_long_t n); */
/* static AP     *ap_read_problem(char *file); */
static void    ap_free(AP *p);

static igraph_long_t     ap_assignment(AP *p, igraph_long_t *res);
/* static igraph_long_t     ap_costmatrix(AP *p, double **m); */
/* static igraph_long_t     ap_datamatrix(AP *p, double **m); */
/* static igraph_long_t     ap_iterations(AP *p); */
static igraph_long_t     ap_hungarian(AP *p);
/* static double  ap_mincost(AP *p); */
/* static void    ap_print_solution(AP *p); */
/* static void    ap_show_data(AP *p); */
/* static igraph_long_t     ap_size(AP *p); */
/* static igraph_long_t     ap_time(AP *p); */

/* error reporting */
/* static void ap_error(char *message); */

/* private functions */
static void    preprocess(AP *p);
static void    preassign(AP *p);
static igraph_long_t     cover(AP *p, igraph_long_t *ri, igraph_long_t *ci);
static void    reduce(AP *p, igraph_long_t *ri, igraph_long_t *ci);

igraph_long_t ap_hungarian(AP *p) {
    igraph_long_t      n;            /* size of problem */
    igraph_long_t    *ri;            /* covered rows    */
    igraph_long_t    *ci;            /* covered columns */
    time_t start, end;     /* timer           */
    igraph_long_t i, j, ok;

    start = time(0);

    n = p->n;
    p->runs = 0;

    /* allocate memory */
    p->s = calloc(1 + n, sizeof(igraph_long_t));
    p->f = calloc(1 + n, sizeof(igraph_long_t));

    ri = calloc(1 + n, sizeof(igraph_long_t));
    ci = calloc(1 + n, sizeof(igraph_long_t));

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
igraph_long_t ap_assignment(AP *p, igraph_long_t *res) {
    igraph_long_t i;

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

AP *ap_read_problem(char *file) {
    FILE *f;
    igraph_long_t i, j, c;
    igraph_long_t m, n;
    double x;
    double **t;
    igraph_long_t nrow, ncol;
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
        igraph_warningf("ap_read_problem: problem not quadratic\nrows = %d, cols = %d\n",
                        IGRAPH_FILE_BASENAME, __LINE__, -1, nrow, ncol);
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

AP     *ap_create_problem_from_matrix(double **t, igraph_long_t n) {
    igraph_long_t i, j;
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

/* read data from vector */
AP *ap_create_problem(double *t, igraph_long_t n) {
    igraph_long_t i, j;
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
    igraph_long_t i;

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
    igraph_long_t i, j;

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

igraph_long_t ap_size(AP *p) {
    return p->n;
}

igraph_long_t ap_time(AP *p) {
    return (igraph_long_t) p->rtime;
}

igraph_long_t ap_iterations(AP *p) {
    return p->runs;
}

void ap_print_solution(AP *p)
{
    igraph_long_t i;

    printf("%d itertations, %d secs.\n",p->runs, (igraph_long_t)p->rtime);
    printf("Min Cost: %10.4f\n",p->cost);

    for(i = 0; i < p->n; i++)
    printf("%4d",p->s[i]);
    printf("\n");
}

igraph_long_t ap_costmatrix(AP *p, double **m) {
    igraph_long_t i, j;

    for (i = 0; i < p->n; i++)
        for (j = 0; j < p->n; j++) {
            m[i][j] = p->C[i + 1][j + 1];
        }

    return p->n;
}

igraph_long_t ap_datamatrix(AP *p, double **m) {
    igraph_long_t i, j;

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

igraph_long_t cover(AP *p, igraph_long_t *ri, igraph_long_t *ci) {
    igraph_long_t *mr, i, r;
    igraph_long_t n;

    n = p->n;
    mr = calloc(1 + p->n, sizeof(igraph_long_t));

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

void reduce(AP *p, igraph_long_t *ri, igraph_long_t *ci) {
    igraph_long_t i, j, n;
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
    igraph_long_t i, j, min, r, c, n, count;
    igraph_long_t *ri, *ci, *rz, *cz;

    n = p->n;
    p->na = 0;

    /* row and column markers */
    ri = calloc(1 + n, sizeof(igraph_long_t));
    ci = calloc(1 + n, sizeof(igraph_long_t));

    /* row and column counts of zeroes */
    rz = calloc(1 + n, sizeof(igraph_long_t));
    cz = calloc(1 + n, sizeof(igraph_long_t));

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
    igraph_long_t i, j, n;
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

igraph_long_t igraph_solve_lsap(igraph_matrix_t *c, igraph_long_t n,
                      igraph_vector_long_t *p) {
    AP *ap;

    IGRAPH_CHECK(igraph_vector_long_resize(p, n));
    igraph_vector_long_null(p);

    ap = ap_create_problem(&MATRIX(*c, 0, 0), n);
    ap_hungarian(ap);
    ap_assignment(ap, VECTOR(*p));
    ap_free(ap);

    return 0;
}
