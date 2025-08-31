
/*
 * This file contains functions that are useful when writing tests.
 * Add #include "test_utilities.h" to the test program to use them.
 */

#include "test_utilities.h"
#include <stdio.h>
#include <string.h>

/* Print an igraph_real_t value. Be consistent in printing NaN/Inf across platforms. */
void print_real(FILE *f, igraph_real_t x, const char *format) {
    igraph_bool_t g8 = !strcmp(format, "%8g");
    if (isfinite(x)) {
        if (x == 0 && signbit(x)) {
            /* print negative zeros as positive zeros for sake of consistency */
            x = +0.0;
        }
        fprintf(f, format, x);
    } else if (isnan(x)) {
        fprintf(f, g8 ? "     NaN" : "NaN");
    } else if (isinf(x) && x > 0) {
        fprintf(f, g8 ? "     Inf" : "Inf");
    } else if (isinf(x) && x < 0) {
        fprintf(f, g8 ? "    -Inf" : "-Inf");
    }
}

void print_vector_format(const igraph_vector_t *v, FILE *f, const char *format) {
    igraph_int_t i, n = igraph_vector_size(v);
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
    igraph_int_t i, n = igraph_vector_int_size(v);
    printf("(");
    for (i=0; i < n; i++) {
        printf(" %" IGRAPH_PRId, VECTOR(*v)[i]);
    }
    printf(" )\n");
}

/* Print all vectors in an integer vector list. Use brackets around each vector
 * and also use brackets around the entire list to make it clear when the list
 * is empty */
void print_vector_int_list(const igraph_vector_int_list_t *v) {
    igraph_int_t i, n = igraph_vector_int_list_size(v);
    printf("{\n");
    for (i = 0; i < n; ++i) {
        printf("  %" IGRAPH_PRId ": ", i);
        print_vector_int(igraph_vector_int_list_get_ptr(v, i));
    }
    printf("}\n");
}

void print_matrix_format(const igraph_matrix_t *m, FILE *f, const char *format) {
    print_matrix_format_indent(m, f, format, "");
}

/* Print elements of a matrix. Use brackets to make it clear when a vector has size zero. */
void print_matrix_format_indent(const igraph_matrix_t *m, FILE *f, const char *format, const char *indent) {
    igraph_int_t i, j, nrow = igraph_matrix_nrow(m), ncol = igraph_matrix_ncol(m);
    if (nrow == 0 || ncol == 0) {
        /* When the matrix is empty, output the dimensions */
        fprintf(f, "%s[ %" IGRAPH_PRId "-by-%" IGRAPH_PRId" ]\n", indent, nrow, ncol);
        return;
    }
    for (i = 0; i < nrow; i++) {
        fprintf(f, i == 0 ? "%s[" : "%s ", indent);
        for (j = 0; j < ncol; j++) {
            fprintf(f, " ");
            print_real(f, MATRIX(*m, i, j), format);
        }
        fprintf(f, i == nrow-1 ? " ]\n" : "\n");
    }
}

void print_matrix(const igraph_matrix_t *m) {
    print_matrix_indent(m, "");
}

void print_matrix_indent(const igraph_matrix_t *m, const char *indent) {
    print_matrix_format_indent(m, stdout, "%8g", indent);
}

void print_matrix_int(const igraph_matrix_int_t *m) {
    igraph_int_t i, j, nrow = igraph_matrix_int_nrow(m), ncol = igraph_matrix_int_ncol(m);
    if (nrow == 0 || ncol == 0) {
        /* When the matrix is empty, output the dimensions */
        printf("[ %" IGRAPH_PRId "-by-%" IGRAPH_PRId" ]\n", nrow, ncol);
        return;
    }
    for (i = 0; i < nrow; i++) {
        printf(i == 0 ? "[" : " ");
        for (j = 0; j < ncol; j++) {
            printf(" ");
            printf("%8" IGRAPH_PRId, MATRIX(*m, i, j));
        }
        printf(i == nrow-1 ? " ]\n" : "\n");
    }
}

/* Round elements of a matrix to integers and print them. */
/* This is meant to be used when the elements of a matrix are integer values. */
void print_matrix_round(const igraph_matrix_t *m) {
    print_matrix_format(m, stdout, "%4.f");
}

void print_matrix_complex_round(const igraph_matrix_complex_t *m) {

    igraph_int_t nr = igraph_matrix_complex_nrow(m);
    igraph_int_t nc = igraph_matrix_complex_ncol(m);
    igraph_int_t i, j;
    for (i = 0; i < nr; i++) {
        for (j = 0; j < nc; j++) {
            igraph_complex_t z = MATRIX(*m, i, j);
            if (j != 0) {
                putchar(' ');
            }
            printf("%.f%+.fi", IGRAPH_REAL(z), IGRAPH_IMAG(z));
        }
        printf("\n");
    }
}

void print_matrix_list(const igraph_matrix_list_t *m) {
    igraph_int_t i, n = igraph_matrix_list_size(m);
    printf("{\n");
    for (i = 0; i < n; ++i) {
        igraph_matrix_t *mat = igraph_matrix_list_get_ptr(m, i);
        if (igraph_matrix_nrow(mat) < 2) {
            printf("  %2" IGRAPH_PRId ": ", i);
            print_matrix(mat);
        } else {
            printf("  %2" IGRAPH_PRId ":\n", i);
            print_matrix_indent(mat, "      ");
        }
    }
    printf("}\n");
}

/* Print an adjacency list. Use brackets around each vector and also use
 * brackets around the entire adjacency list to make it clear when the list
 * is empty.
 */
void print_adjlist(const igraph_adjlist_t *adjlist) {
    igraph_int_t vcount = igraph_adjlist_size(adjlist);
    igraph_int_t i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %" IGRAPH_PRId ": ", i);
        print_vector_int(igraph_adjlist_get(adjlist, i));
    }
    printf("}\n");
}

/* Print a graph. Use brackets to make it obvious when the edge list is empty. */
void print_graph(const igraph_t *graph) {
    print_weighted_graph(graph, NULL);
}

/* Print a graph with edge weights from a vector. Use brackets to make it
 * obvious when the edge list is empty. */
void print_weighted_graph(const igraph_t *graph, const igraph_vector_t* weights) {
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t i;

    printf("directed: %s\n", igraph_is_directed(graph) ? "true" : "false");
    printf("vcount: %" IGRAPH_PRId "\n", vcount);
    printf("edges: {\n");
    for (i=0; i < ecount; ++i) {
        printf(
            "%" IGRAPH_PRId " %" IGRAPH_PRId,
            IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i)
        );
        if (weights) {
            printf(": ");
            print_real(stdout, VECTOR(*weights)[i], "%g");
        }
        printf("\n");
    }
    printf("}\n");
}

/* Print a graph with edge weights from an edge attribute. Use brackets to make
 * it obvious when the edge list is empty. */
void print_weighted_graph_attr(const igraph_t *graph, const char* attr) {
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_int_t i;

    printf("directed: %s\n", igraph_is_directed(graph) ? "true" : "false");
    printf("vcount: %" IGRAPH_PRId "\n", vcount);
    printf("edges: {\n");
    for (i=0; i < ecount; ++i)
        printf
            ("%" IGRAPH_PRId " %" IGRAPH_PRId ": %g\n",
            IGRAPH_FROM(graph, i), IGRAPH_TO(graph, i),
            EAN(graph, attr, i)
        );
    printf("}\n");
}

/* Print an incidence list. Use brackets around each vector and also use
 * brackets around the entire incidence list to make it clear when the list
 * is empty.
 */
void print_inclist(const igraph_inclist_t *inclist) {
    igraph_int_t vcount = igraph_inclist_size(inclist);
    igraph_int_t i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %" IGRAPH_PRId ": ", i);
        print_vector_int(igraph_inclist_get(inclist, i));
    }
    printf("}\n");
}

/* Print a lazy adjacency list. Use brackets around each vector and also use
 * brackets around the entire lazy adjacency list to make it clear when the list
 * is empty.
 */
void print_lazy_adjlist(igraph_lazy_adjlist_t *adjlist) {
    igraph_int_t vcount = igraph_lazy_adjlist_size(adjlist);
    igraph_int_t i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %" IGRAPH_PRId ": ", i);
        print_vector_int(igraph_lazy_adjlist_get(adjlist, i));
    }
    printf("}\n");
}

/* Print a lazy incidence list. Use brackets around each vector and also use
 * brackets around the entire incidence list to make it clear when the list
 * is empty.
 */
void print_lazy_inclist(igraph_lazy_inclist_t *inclist) {
    igraph_int_t vcount = igraph_lazy_inclist_size(inclist);
    igraph_int_t i;

    printf("{\n");
    for (i = 0; i < vcount; ++i) {
        printf("  %" IGRAPH_PRId ": ", i);
        print_vector_int(igraph_lazy_inclist_get(inclist, i));
    }
    printf("}\n");
}

/* Edge comparison function used for sorting in print_graph_canon(). */
int edge_compare(void *pedges, const void *pi1, const void *pi2) {
    const igraph_int_t i1 = * (const igraph_int_t *) pi1;
    const igraph_int_t i2 = * (const igraph_int_t *) pi2;
    const igraph_vector_int_t *edges = (const igraph_vector_int_t *) pedges;
    if (VECTOR(*edges)[2*i1] < VECTOR(*edges)[2*i2]) {
        return -1;
    } else if (VECTOR(*edges)[2*i1] > VECTOR(*edges)[2*i2]) {
        return 1;
    } else if (VECTOR(*edges)[2*i1+1] < VECTOR(*edges)[2*i2+1]) {
        return -1;
    } else if (VECTOR(*edges)[2*i1+1] > VECTOR(*edges)[2*i2+1]) {
        return 1;
    } else {
        return 0;
    }
}

/* Print a weighted graph using a sorted edge list. Other than sorting (i.e. canonicalizing)
 * the edge list, this function is identical to print_graph(). */
void print_weighted_graph_canon(const igraph_t *graph, const igraph_vector_t *weights) {
    igraph_int_t ecount = igraph_ecount(graph);
    igraph_int_t vcount = igraph_vcount(graph);
    igraph_vector_int_t edges, idx;

    printf("directed: %s\n", igraph_is_directed(graph) ? "true" : "false");
    printf("vcount: %" IGRAPH_PRId "\n", vcount);
    printf("edges: {\n");

    igraph_vector_int_init(&edges, 0);
    igraph_get_edgelist(graph, &edges, false);

    /* If the graph is undirected, we make sure that the first vertex of undirected edges
     * is always the one with the lower ID. */
    if (! igraph_is_directed(graph)) {
        for (igraph_int_t i=0; i < ecount; i++) {
            if (VECTOR(edges)[2*i] > VECTOR(edges)[2*i+1]) {
                igraph_int_t tmp = VECTOR(edges)[2*i];
                VECTOR(edges)[2*i] = VECTOR(edges)[2*i+1];
                VECTOR(edges)[2*i+1] = tmp;
            }
        }
    }

    igraph_vector_int_init_range(&idx, 0, igraph_ecount(graph));

    /* Sort the edge list */
    igraph_qsort_r(&VECTOR(idx)[0], ecount, sizeof(igraph_int_t), &edges, &edge_compare);

    for (igraph_int_t i=0; i < ecount; i++) {
        const igraph_int_t k = VECTOR(idx)[i];
        printf("%" IGRAPH_PRId " %" IGRAPH_PRId, VECTOR(edges)[2*k], VECTOR(edges)[2*k+1]);
        if (weights) {
            printf(": ");
            print_real(stdout, VECTOR(*weights)[k], "%g");
        }
        printf("\n");
    }

    printf("}\n");

    igraph_vector_int_destroy(&idx);
    igraph_vector_int_destroy(&edges);
}

/* Print a graph using a sorted edge list. Other than sorting (i.e. canonicalizing)
 * the edge list, this function is identical to print_graph(). */
void print_graph_canon(const igraph_t *graph) {
    print_weighted_graph_canon(graph, NULL);
}

/* Print a vector, ensuring that the first nonzero element is positive. */
void print_vector_first_nonzero_element_positive(const igraph_vector_t *vector, const char* format) {
    igraph_vector_t copy;
    igraph_int_t i, n;

    igraph_vector_init_copy(&copy, vector);

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
    igraph_int_t i, n;

    igraph_vector_complex_init_copy(&copy, vector);

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
    igraph_int_t i, j, nrow, ncol;

    igraph_matrix_init_copy(&copy, matrix);

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
    igraph_int_t i, j, nrow, ncol;
    igraph_complex_t z;
    char buf[256];
    size_t len;

    igraph_matrix_complex_init_copy(&copy, matrix);

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

void matrix_init_int_row_major(igraph_matrix_t *mat, igraph_int_t nrow, igraph_int_t ncol, const int *elem) {
    igraph_int_t c, r;
    size_t i_elem = 0;
    igraph_matrix_init(mat, nrow, ncol);
    for (r = 0; r < nrow; r++) {
        for (c = 0; c < ncol; c++) {
            MATRIX(*mat, r, c) = elem[i_elem];
            i_elem++;
        }
    }
}

void matrix_int_init_int_row_major(igraph_matrix_int_t *mat, igraph_int_t nrow, igraph_int_t ncol, const int *elem) {
    igraph_int_t c, r;
    size_t i_elem = 0;
    igraph_matrix_int_init(mat, nrow, ncol);
    for (r = 0; r < nrow; r++) {
        for (c = 0; c < ncol; c++) {
            MATRIX(*mat, r, c) = elem[i_elem];
            i_elem++;
        }
    }
}

void matrix_init_real_row_major(igraph_matrix_t *mat, igraph_int_t nrow, igraph_int_t ncol, const igraph_real_t *elem) {
    igraph_int_t c, r;
    size_t i_elem = 0;
    igraph_matrix_init(mat, nrow, ncol);
    for (r = 0; r < nrow; r++) {
        for (c = 0; c < ncol; c++) {
            MATRIX(*mat, r, c) = elem[i_elem];
            i_elem++;
        }
    }
}

void matrix_chop(igraph_matrix_t *mat, igraph_real_t cutoff) {
    igraph_int_t nelems = igraph_matrix_nrow(mat) * igraph_matrix_ncol(mat);
    for (igraph_int_t i = 0; i < nelems; i++) {
        if (fabs(VECTOR(mat->data)[i]) < cutoff) {
            VECTOR(mat->data)[i] = 0;
        }
    }
}

void vector_chop(igraph_vector_t *vec, igraph_real_t cutoff) {
    igraph_int_t nelems = igraph_vector_size(vec);
    for (igraph_int_t i = 0; i < nelems; i++) {
        if (fabs(VECTOR(*vec)[i]) < cutoff) {
            VECTOR(*vec)[i] = 0;
        }
    }
}

/* print all graph, edge and vertex attributes of a graph */
void print_attributes(const igraph_t *g) {
    igraph_vector_int_t gtypes, vtypes, etypes;
    igraph_strvector_t gnames, vnames, enames;
    igraph_int_t i;

    igraph_int_t j;

    igraph_vector_int_init(&gtypes, 0);
    igraph_vector_int_init(&vtypes, 0);
    igraph_vector_int_init(&etypes, 0);
    igraph_strvector_init(&gnames, 0);
    igraph_strvector_init(&vnames, 0);
    igraph_strvector_init(&enames, 0);

    igraph_cattribute_list(g, &gnames, &gtypes, &vnames, &vtypes,
                           &enames, &etypes);

    /* Graph attributes */
    for (i = 0; i < igraph_strvector_size(&gnames); i++) {
        if (i != 0)
            putchar(' ');
        printf("%s=", igraph_strvector_get(&gnames, i));
        if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
            igraph_real_printf(GAN(g, igraph_strvector_get(&gnames, i)));
        } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_BOOLEAN) {
            printf("%d", GAB(g, igraph_strvector_get(&gnames, i)));
        } else {
            printf("\"%s\"", GAS(g, igraph_strvector_get(&gnames, i)));
        }
    }
    if (igraph_strvector_size(&gnames))
        printf("\n");

    for (i = 0; i < igraph_vcount(g); i++) {
        printf("Vertex %" IGRAPH_PRId ":", i);
        for (j = 0; j < igraph_strvector_size(&vnames); j++) {
            putchar(' ');
            printf("%s=", igraph_strvector_get(&vnames, j));
            if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(VAN(g, igraph_strvector_get(&vnames, j), i));
            } else if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                printf("%d", VAB(g, igraph_strvector_get(&vnames, j), i));
            } else {
                printf("\"%s\"", VAS(g, igraph_strvector_get(&vnames, j), i));
            }
        }
        printf("\n");
    }

    for (i = 0; i < igraph_ecount(g); i++) {
        printf("Edge %" IGRAPH_PRId " (%" IGRAPH_PRId "-%" IGRAPH_PRId "):", i, IGRAPH_FROM(g, i), IGRAPH_TO(g, i));
        for (j = 0; j < igraph_strvector_size(&enames); j++) {
            putchar(' ');
            printf("%s=", igraph_strvector_get(&enames, j));
            if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
                igraph_real_printf(EAN(g, igraph_strvector_get(&enames, j), i));
            } else if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_BOOLEAN) {
                printf("%d", EAB(g, igraph_strvector_get(&enames, j), i));
            } else {
                printf("\"%s\"", EAS(g, igraph_strvector_get(&enames, j), i));
            }
        }
        printf("\n");
    }
    printf("\n");

    igraph_strvector_destroy(&enames);
    igraph_strvector_destroy(&vnames);
    igraph_strvector_destroy(&gnames);
    igraph_vector_int_destroy(&etypes);
    igraph_vector_int_destroy(&vtypes);
    igraph_vector_int_destroy(&gtypes);
}

expect_warning_context_t expect_warning_ctx;

void record_last_warning(const char *reason, const char *file, int line) {
    IGRAPH_UNUSED(file); IGRAPH_UNUSED(line);

    if (expect_warning_ctx.observed) {
        igraph_free(expect_warning_ctx.observed);
    }

    expect_warning_ctx.observed = strdup(reason);
}

void print_bitset(const igraph_bitset_t* bitset) {
    printf("(");
    for (igraph_int_t i = bitset->size - 1; i >= 0; --i) {
        printf(" %d", !!IGRAPH_BIT_TEST(*bitset, i));
    }
    printf(" )\n");
}

void print_bitset_list(const igraph_bitset_list_t *v) {
    igraph_int_t i, n = igraph_bitset_list_size(v);
    printf("{\n");
    for (i = 0; i < n; ++i) {
        printf("  %" IGRAPH_PRId ": ", i);
        print_bitset(igraph_bitset_list_get_ptr(v, i));
    }
    printf("}\n");
}
