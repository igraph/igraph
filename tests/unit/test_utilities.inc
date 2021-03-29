#ifndef TEST_UTILITIES_INC
#define TEST_UTILITIES_INC

/*
 * This file contains functions that are useful when writing tests.
 * Include it in the test program using #include "test_utilities.inc"
 */

#include <igraph.h>
#include <stdio.h>
#include <string.h>

/* Print an igraph_real_t value. Be consistent in printing NaN/Inf across platforms. */
void print_real(FILE *f, igraph_real_t x, const char *format) {
    igraph_bool_t g8 = !strcmp(format, "%8g");
    if (igraph_finite(x)) {
        fprintf(f, format, x);
    } else if (igraph_is_nan(x)) {
        fprintf(f, g8 ? "     NaN" : "NaN");
    } else if (igraph_is_posinf(x)) {
        fprintf(f, g8 ? "     Inf" : "Inf");
    } else if (igraph_is_neginf(x)) {
        fprintf(f, g8 ? "    -Inf" : "-Inf");
    }
}

void print_vector_format(const igraph_vector_t *v, FILE *f, const char *format) {
    long i, n = igraph_vector_size(v);
    fprintf(f, "(");
    for (i=0; i < n; i++) {
        fprintf(f, " ");
        print_real(f, VECTOR(*v)[i], format);
    }
    fprintf(f, " )\n");
}

/* Print elements of a vector. Use parentheses to make it clear when a vector has size zero. */
void print_vector(const igraph_vector_t *v) {
    print_vector_format(v, stdout, "%g");
}

/* Round elements of a vector to integers and print them. */
/* This is meant to be used when the elements of a vector are integer values. */
void print_vector_round(const igraph_vector_t *v) {
    print_vector_format(v, stdout, "%.f");
}


/* Print elements of an integer vector */
void print_vector_int(const igraph_vector_int_t *v) {
    long i, n = igraph_vector_int_size(v);
    printf("(");
    for (i=0; i < n; i++) {
        printf(" %d", VECTOR(*v)[i]);
    }
    printf(" )\n");
}


/* Print elements of a matrix. Use brackets to make it clear when a vector has size zero. */
void print_matrix_format(const igraph_matrix_t *m, FILE *f, const char *format) {
    long i, j, nrow = igraph_matrix_nrow(m), ncol = igraph_matrix_ncol(m);
    for (i = 0; i < nrow; i++) {
        fprintf(f, i == 0 ? "[" : " ");
        for (j = 0; j < ncol; j++) {
            fprintf(f, " ");
            print_real(f, MATRIX(*m, i, j), format);
        }
        fprintf(f, i == nrow-1 ? " ]\n" : "\n");
    }
}

void print_matrix(const igraph_matrix_t *m) {
    print_matrix_format(m, stdout, "%8g");
}

/* Round elements of a matrix to integers and print them. */
/* This is meant to be used when the elements of a matrix are integer values. */
void print_matrix_round(const igraph_matrix_t *m) {
    print_matrix_format(m, stdout, "%4.f");
}


/* Print an adjacency list. Use brackets around each vector and also use
 * brackets around the entire adjacency list to make it clear when the list
 * is empty.
 */
void print_adjlist(const igraph_adjlist_t *adjlist) {
    long vcount = igraph_adjlist_size(adjlist);
    long i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %ld: ", i);
        print_vector_int(igraph_adjlist_get(adjlist, i));
    }
    printf("}\n");
}

/* Print a graph. Use brackets to make it obvious when the edge list is empty. */
void print_graph(const igraph_t *graph) {
    long ecount = igraph_ecount(graph);
    long vcount = igraph_vcount(graph);
    long i;

    printf("directed: %s\n", igraph_is_directed(graph) ? "true" : "false");
    printf("vcount: %ld\n", vcount);
    printf("edges: {\n");
    for (i=0; i < ecount; ++i)
        printf("%" IGRAPH_PRId " %" IGRAPH_PRId "\n", IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
    printf("}\n");
}

/* Print an incidence list. Use brackets around each vector and also use
 * brackets around the entire incidence list to make it clear when the list
 * is empty.
 */
void print_inclist(const igraph_inclist_t *inclist) {
    long vcount = igraph_inclist_size(inclist);
    long i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %ld: ", i);
        print_vector_int(igraph_inclist_get(inclist, i));
    }
    printf("}\n");
}

/* Print a lazy adjacency list. Use brackets around each vector and also use
 * brackets around the entire lazy adjacency list to make it clear when the list
 * is empty.
 */
void print_lazy_adjlist(igraph_lazy_adjlist_t *adjlist) {
    long vcount = igraph_lazy_adjlist_size(adjlist);
    long i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %ld: ", i);
        print_vector_int(igraph_lazy_adjlist_get(adjlist, i));
    }
    printf("}\n");
}

/* Print a lazy incidence list. Use brackets around each vector and also use
 * brackets around the entire incidence list to make it clear when the list
 * is empty.
 */
void print_lazy_inclist(igraph_lazy_inclist_t *inclist) {
    long vcount = igraph_lazy_inclist_size(inclist);
    long i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %ld: ", i);
        print_vector_int(igraph_lazy_inclist_get(inclist, i));
    }
    printf("}\n");
}

/* Edge comparisong function used for sorting in print_graph_canon(). */
int edge_compare(const void *e1, const void *e2) {
    const igraph_real_t *edge1 = (igraph_real_t *) e1, *edge2 = (igraph_real_t *) e2;
    if (edge1[0] < edge2[0]) {
        return -1;
    } else if (edge1[0] > edge2[0]) {
        return 1;
    } else if (edge1[1] < edge2[1]) {
        return -1;
    } else if (edge1[1] > edge2[1]) {
        return 1;
    } else {
        return 0;
    }
}

/* Print a graph using a sorted edge list. Other than sorting (i.e. canonicalizing) the edge list,
 * this function is identical to print_graph(). */
void print_graph_canon(const igraph_t *graph) {
    long ecount = igraph_ecount(graph);
    long vcount = igraph_vcount(graph);
    long i;
    igraph_vector_t edges;

    printf("directed: %s\n", igraph_is_directed(graph) ? "true" : "false");
    printf("vcount: %ld\n", vcount);
    printf("edges: {\n");

    igraph_vector_init(&edges, 0);
    igraph_get_edgelist(graph, &edges, 0);

    /* If the graph is undirected, we make sure that the first vertex of undirected edges
     * is always the one with the lower ID. */
    if (! igraph_is_directed(graph)) {
        for (i=0; i < ecount; ++i) {
            if (VECTOR(edges)[2*i] > VECTOR(edges)[2*i+1]) {
                igraph_real_t tmp = VECTOR(edges)[2*i];
                VECTOR(edges)[2*i] = VECTOR(edges)[2*i+1];
                VECTOR(edges)[2*i+1] = tmp;
            }
        }
    }

    /* Sort the edge list */
    igraph_qsort(&VECTOR(edges)[0], ecount, 2*sizeof(igraph_real_t), &edge_compare);

    for (i=0; i < ecount; ++i) {
        printf("%ld %ld\n", (long) VECTOR(edges)[2*i], (long) VECTOR(edges)[2*i+1]);
    }

    igraph_vector_destroy(&edges);

    printf("}\n");
}

/* Print a vector, ensuring that the first nonzero element is positive. */
void print_vector_first_nonzero_element_positive(const igraph_vector_t *vector, const char* format) {
    igraph_vector_t copy;
    long i, n;

    igraph_vector_copy(&copy, vector);

    n = igraph_vector_size(&copy);

    for (i = 0; i < n; i++) {
        if (VECTOR(copy)[i] < 0) {
            for (; i < n; i++) {
                if (VECTOR(copy)[i] != 0) {
                    VECTOR(copy)[i] *= -1;
                }
            }
            break;
        } else if (VECTOR(copy)[i] > 0) {
            break;
        }
    }

    igraph_vector_printf(&copy, format);
    igraph_vector_destroy(&copy);
}

/* Print a complex vector, ensuring that the first element with nonzero real
 * part has a positive real part. */
void print_vector_complex_first_nonzero_real_part_positive(const igraph_vector_complex_t *vector) {
    igraph_vector_complex_t copy;
    long i, n;

    igraph_vector_complex_copy(&copy, vector);

    n = igraph_vector_complex_size(&copy);

    for (i = 0; i < n; i++) {
        if (IGRAPH_REAL(VECTOR(copy)[i]) < 0) {
            for (; i < n; i++) {
                if (IGRAPH_REAL(VECTOR(copy)[i]) != 0) {
                    IGRAPH_REAL(VECTOR(copy)[i]) *= -1;
                }
                if (IGRAPH_IMAG(VECTOR(copy)[i]) != 0) {
                    IGRAPH_IMAG(VECTOR(copy)[i]) *= -1;
                }
            }
            break;
        } else if (IGRAPH_REAL(VECTOR(copy)[i]) > 0) {
            break;
        }
    }

    igraph_vector_complex_print(&copy);
    igraph_vector_complex_destroy(&copy);
}

/* Print a matrix, ensuring that the first nonzero element in each column is
 * positive. */
void print_matrix_first_row_positive(const igraph_matrix_t *matrix, const char* format) {
    igraph_matrix_t copy;
    long i, j, nrow, ncol;

    igraph_matrix_copy(&copy, matrix);

    nrow = igraph_matrix_nrow(&copy);
    ncol = igraph_matrix_ncol(&copy);

    for (i = 0; i < ncol; i++) {
        for (j = 0; j < nrow; j++) {
            if (MATRIX(copy, j, i) < 0) {
                for (; j < nrow; j++) {
                    if (MATRIX(copy, j, i) != 0) {
                        MATRIX(copy, j, i) *= -1;
                    }
                }
                break;
            } else if (MATRIX(copy, j, i) > 0) {
                break;
            }
        }
    }

    igraph_matrix_printf(&copy, format);
    igraph_matrix_destroy(&copy);
}

/* Print a complex matrix, ensuring that the first element with nonzero real
 * part in each column has a positive real part. */
void print_matrix_complex_first_row_positive(const igraph_matrix_complex_t *matrix) {
    igraph_matrix_complex_t copy;
    long i, j, nrow, ncol;
    igraph_complex_t z;
    char buf[256];
    size_t len;

    igraph_matrix_complex_copy(&copy, matrix);

    nrow = igraph_matrix_complex_nrow(&copy);
    ncol = igraph_matrix_complex_ncol(&copy);

    for (i = 0; i < ncol; i++) {
        for (j = 0; j < nrow; j++) {
            if (IGRAPH_REAL(MATRIX(copy, j, i)) < 0) {
                for (; j < nrow; j++) {
                    if (IGRAPH_REAL(MATRIX(copy, j, i)) != 0) {
                        IGRAPH_REAL(MATRIX(copy, j, i)) *= -1;
                    }
                    if (IGRAPH_IMAG(MATRIX(copy, j, i)) != 0) {
                        IGRAPH_IMAG(MATRIX(copy, j, i)) *= -1;
                    }
                }
                break;
            } else if (IGRAPH_REAL(MATRIX(copy, j, i)) > 0) {
                break;
            }
        }
    }

    for (i = 0; i < nrow; i++) {
        for (j = 0; j < ncol; j++) {
            z = MATRIX(copy, i, j);
            if (j != 0) {
                putchar(' ');
            }

            snprintf(buf, sizeof(buf), "%g%+gi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
            len = strlen(buf);

            /* ensure that we don't print -0 in the imaginary part */
            if (len > 3 && buf[len-3] == '-' && buf[len-2] == '0' && buf[len-1] == 'i') {
              buf[len-3] = '+';
            }

            /* ensure that we don't print -0 in the real part either */
            if (buf[0] == '-' && buf[1] == '0' && (buf[2] == '+' || buf[2] == '-')) {
                printf("%s", buf + 1);
            } else {
                printf("%s", buf);
            }
        }
        printf("\n");
    }

    igraph_matrix_complex_destroy(&copy);
}

void matrix_init_int_row_major(igraph_matrix_t *mat, int nrow, int ncol, int* elem) {
    int c, r;
    int i_elem = 0;
    igraph_matrix_init(mat, nrow, ncol);
    for (r = 0; r < nrow; r++) {
        for (c = 0; c < ncol; c++) {
            MATRIX(*mat, r, c) = elem[i_elem];
            i_elem++;
        }
    }
}

void matrix_init_real_row_major(igraph_matrix_t *mat, int nrow, int ncol, igraph_real_t* elem) {
    int c, r;
    int i_elem = 0;
    igraph_matrix_init(mat, nrow, ncol);
    for (r = 0; r < nrow; r++) {
        for (c = 0; c < ncol; c++) {
            MATRIX(*mat, r, c) = elem[i_elem];
            i_elem++;
        }
    }
}

void matrix_chop(igraph_matrix_t *mat, igraph_real_t cutoff) {
    int i;
    for (i = 0; i < igraph_matrix_nrow(mat) * igraph_matrix_ncol(mat); i++) {
        if (fabs(VECTOR(mat->data)[i]) < cutoff) {
            VECTOR(mat->data)[i] = 0;
        }
    }
}

void print_spmatrix(igraph_spmatrix_t *m) {
    long int i, j;
    for (i = 0; i < igraph_spmatrix_nrow(m); i++) {
        for (j = 0; j < igraph_spmatrix_ncol(m); j++) {
            printf(" %8g", igraph_spmatrix_e(m, i, j));
        }
        printf("\n");
    }
}

#define VERIFY_FINALLY_STACK() \
    if (!IGRAPH_FINALLY_STACK_EMPTY) { \
        printf( \
          "Finally stack is not empty. " \
          "Check that the number in IGRAPH_FINALLY_CLEAN matches the IGRAPH_FINALLY count.\n"); \
        return 100; \
    }

/* Run a test in a separate function; return the return value of the function
 * if it is nonzero. Also verify the FINALLY stack and bail out if it is not
 * empty. Needs an integer variable named 'retval' in the local context. */
#define RUN_TEST(func) \
    { \
        retval = func(); \
        if (retval) { \
            return retval; \
        } \
        VERIFY_FINALLY_STACK(); \
    }

#endif /* TEST_UTILITIES_INC */
