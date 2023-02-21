#ifndef TEST_UTILITIES_H
#define TEST_UTILITIES_H

/*
 * This file declares functions that are useful when writing tests.
 */

#include <igraph.h>
#include <stdio.h>
#include <string.h>

/* Structure used by the EXPECT_WARNING macro to test whether a warning was
 * triggered in a section of the code in unit tests */
typedef struct {
    const char* expected;
    char* observed;
} expect_warning_context_t;

extern expect_warning_context_t expect_warning_ctx;

/* Print an igraph_real_t value. Be consistent in printing NaN/Inf across platforms. */
void print_real(FILE *f, igraph_real_t x, const char *format);

void print_vector_format(const igraph_vector_t *v, FILE *f, const char *format);

/* Print elements of a vector. Use parentheses to make it clear when a vector has size zero. */
void print_vector(const igraph_vector_t *v);

/* Round elements of a vector to integers and print them. */
/* This is meant to be used when the elements of a vector are integer values. */
void print_vector_round(const igraph_vector_t *v);

/* Print elements of an integer vector */
void print_vector_int(const igraph_vector_int_t *v);

/* Print all vectors in an integer vector list. Use brackets around each vector
 * and also use brackets around the entire list to make it clear when the list
 * is empty */
void print_vector_int_list(const igraph_vector_int_list_t *v);

/* Print elements of a matrix. Use brackets to make it clear when a vector has size zero. */
void print_matrix_format(const igraph_matrix_t *m, FILE *f, const char *format);

void print_matrix(const igraph_matrix_t *m);

void print_matrix_int(const igraph_matrix_int_t *m);

/* Round elements of a matrix to integers and print them. */
/* This is meant to be used when the elements of a matrix are integer values. */
void print_matrix_round(const igraph_matrix_t *m);

/* Round elements of a complex matrix to integers and print them. */
/* This is meant to be used when the elements of a matrix are integer values. */
void print_matrix_complex_round(const igraph_matrix_complex_t *m);

/* Print an adjacency list. Use brackets around each vector and also use
 * brackets around the entire adjacency list to make it clear when the list
 * is empty.
 */
void print_adjlist(const igraph_adjlist_t *adjlist);

/* Print a graph. Use brackets to make it obvious when the edge list is empty. */
void print_graph(const igraph_t *graph);

/* Print a graph with edge weights from a vector. Use brackets to make it
 * obvious when the edge list is empty. */
void print_weighted_graph(const igraph_t *graph, const igraph_vector_t* weights);

/* Print a graph with edge weights from an edge attribute. Use brackets to make
 * it obvious when the edge list is empty. */
void print_weighted_graph_attr(const igraph_t *graph, const char* attr);

/* Print an incidence list. Use brackets around each vector and also use
 * brackets around the entire incidence list to make it clear when the list
 * is empty.
 */
void print_inclist(const igraph_inclist_t *inclist);

/* Print a lazy adjacency list. Use brackets around each vector and also use
 * brackets around the entire lazy adjacency list to make it clear when the list
 * is empty.
 */
void print_lazy_adjlist(igraph_lazy_adjlist_t *adjlist);

/* Print a lazy incidence list. Use brackets around each vector and also use
 * brackets around the entire incidence list to make it clear when the list
 * is empty.
 */
void print_lazy_inclist(igraph_lazy_inclist_t *inclist);

/* Edge comparison function used for sorting in print_graph_canon(). */
int edge_compare(const void *e1, const void *e2);

/* Print a graph using a sorted edge list. Other than sorting (i.e. canonicalizing) the edge list,
 * this function is identical to print_graph(). */
void print_graph_canon(const igraph_t *graph);

/* Print a vector, ensuring that the first nonzero element is positive. */
void print_vector_first_nonzero_element_positive(const igraph_vector_t *vector, const char* format);

/* Print a complex vector, ensuring that the first element with nonzero real
 * part has a positive real part. */
void print_vector_complex_first_nonzero_real_part_positive(const igraph_vector_complex_t *vector);

/* Print a matrix, ensuring that the first nonzero element in each column is
 * positive. */
void print_matrix_first_row_positive(const igraph_matrix_t *matrix, const char* format);

/* Print a complex matrix, ensuring that the first element with nonzero real
 * part in each column has a positive real part. */
void print_matrix_complex_first_row_positive(const igraph_matrix_complex_t *matrix);

/* print all graph, edge and vertex attributes of a graph */
void print_attributes(const igraph_t *g);

void matrix_init_int_row_major(igraph_matrix_t *mat, igraph_integer_t nrow, igraph_integer_t ncol, const int *elem);
void matrix_int_init_int_row_major(igraph_matrix_int_t *mat, igraph_integer_t nrow, igraph_integer_t ncol, const int *elem);
void matrix_init_real_row_major(igraph_matrix_t *mat, igraph_integer_t nrow, igraph_integer_t ncol, const igraph_real_t *elem);

void matrix_chop(igraph_matrix_t *mat, igraph_real_t cutoff);
void vector_chop(igraph_vector_t *vec, igraph_real_t cutoff);

#define VERIFY_FINALLY_STACK() \
    if (!IGRAPH_FINALLY_STACK_EMPTY) { \
        printf( \
          "%s:%d : " \
          "Finally stack is not empty (stack size is %d). " \
          "Check that the number in IGRAPH_FINALLY_CLEAN matches the IGRAPH_FINALLY count.\n", \
          IGRAPH_FILE_BASENAME, __LINE__, IGRAPH_FINALLY_STACK_SIZE()); \
        abort(); \
    }

/* Run a test in a separate function; return the return value of the function
 * if it is nonzero. Also verify the FINALLY stack and bail out if it is not
 * empty. Needs an integer variable named 'retval' in the local context. */
#define RUN_TEST(func) \
    { \
        int retval = func(); \
        if (retval) { \
            return retval; \
        } \
        VERIFY_FINALLY_STACK(); \
    }

#define CHECK_ERROR(funcall, expected_err) \
    do { \
        igraph_error_handler_t *handler; \
        handler = igraph_set_error_handler(igraph_error_handler_ignore); \
        IGRAPH_ASSERT(funcall == expected_err); \
        igraph_set_error_handler(handler); \
    } while (0)

void record_last_warning(const char *reason, const char *file, int line);

#define EXPECT_WARNING(funcall, expected_warning) \
    do { \
        igraph_warning_handler_t *handler; \
        expect_warning_ctx.expected = expected_warning; \
        expect_warning_ctx.observed = NULL; \
        handler = igraph_set_warning_handler(record_last_warning); \
        IGRAPH_ASSERT(IGRAPH_SUCCESS == funcall); \
        igraph_set_warning_handler(handler); \
        if (expect_warning_ctx.observed == NULL) { \
            printf("Expected this warning but none was raised:\n  %s\n", expected_warning); \
            abort(); \
        } else if (strcmp(expect_warning_ctx.observed, expect_warning_ctx.expected)) { \
            printf("Expected warning:\n  %s\ngot:\n  %s\n", expected_warning, expect_warning_ctx.observed); \
            free(expect_warning_ctx.observed); \
            abort(); \
        } else { \
            free(expect_warning_ctx.observed); \
        } \
    } while (0)
#endif /* TEST_UTILITIES_H */
