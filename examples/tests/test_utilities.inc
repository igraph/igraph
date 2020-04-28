#ifndef TEST_UTILITIES_INC
#define TEST_UTILITIES_INC

/*
 * This file contains functions that are useful when writing tests.
 * Include it in the test program using #include "test_utilities.inc"
 */

#include <igraph.h>
#include <stdio.h>
#include <string.h>

/* Print elements of a vector. Use parentheses to make it clear when a vector has size zero. */
void print_vector(const igraph_vector_t *v, FILE *f) {
    long i;
    fprintf(f, "(");
    for (i=0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %f", VECTOR(*v)[i]);
    }
    fprintf(f, " )\n");
}


/* Round elements of a vector to integers and print them. */
/* This is meant to be used when the elements of a vector are integer values. */
void print_vector_round(const igraph_vector_t *v, FILE *f) {
    long i;
    fprintf(f, "(");
    for (i=0; i < igraph_vector_size(v); i++) {
        fprintf(f, " %li", (long int) VECTOR(*v)[i]);
    }
    fprintf(f, " )\n");
}


/* Print elements of an integer vector */
void print_vector_int(const igraph_vector_int_t *v, FILE *f) {
    long i;
    fprintf(f, "(");
    for (i=0; i < igraph_vector_int_size(v); i++) {
        fprintf(f, " %d", VECTOR(*v)[i]);
    }
    fprintf(f, " )\n");
}


/* Print a graph. Use brackets to make it obvious when the edge list is empty. */
void print_graph(const igraph_t *graph, FILE *f) {
    long ecount = igraph_ecount(graph);
    long vcount = igraph_vcount(graph);
    long i;

    fprintf(f, "directed: %s\n", igraph_is_directed(graph) ? "true" : "false");
    fprintf(f, "vcount: %ld\n", vcount);
    fprintf(f, "edges: {\n");
    for (i=0; i < ecount; ++i)
        fprintf(f, "%d %d\n", IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i));
    fprintf(f, "}\n");
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

#endif /* TEST_UTILITIES_INC */
