
#include <igraph.h>

#include "test_utilities.inc"

/* Matrix-vector multiplication: y = A.x */
void matmul(const igraph_matrix_t *A, const igraph_vector_t *x, igraph_vector_t *y, igraph_real_t beta) {
    long int i, j, nr = igraph_matrix_nrow(A), nc = igraph_matrix_ncol(A);

    IGRAPH_ASSERT(nc == igraph_vector_size(x));
    IGRAPH_ASSERT(nr == igraph_vector_size(y));

    igraph_vector_scale(y, beta);

    for (i=0; i < nr; ++i) {
        for (j=0; j < nc; ++j) {
            VECTOR(*y)[i] += MATRIX(*A, i, j) * VECTOR(*x)[j];
        }
    }
}

int main() {
    igraph_matrix_t A;
    igraph_vector_t x, y1, y2;
    long int i, j;
    const long int nr = 5, nc = 8;

    igraph_rng_seed(igraph_rng_default(), 54632);

    igraph_matrix_init(&A, nr, nc);
    igraph_vector_init(&x, nc);

    /* Fill with arbitrary values. Should be zeroes by beta. */
    igraph_vector_init_seq(&y1, 1, nr);
    igraph_vector_copy(&y2, &y1);

    for (i=0; i < nr; ++i) {
        for (j=0; j < nc; ++j) {
            MATRIX(A, i, j) = (igraph_real_t) RNG_INTEGER(-10, 10);
        }
    }

    for (j=0; j < nc; ++j) {
        VECTOR(x)[j] = (igraph_real_t) RNG_INTEGER(-10, 10);
    }

    printf("Input matrix A:\n");
    print_matrix(&A);

    printf("\nInput vector x:\n");
    print_vector(&x);

    igraph_blas_dgemv(0, 1, &A, &x, 0, &y1);
    matmul(&A, &x, &y2, 0);

    printf("\nResult vector DGEMV:\n");
    print_vector(&y1);

    printf("\nResult vector naive:\n");
    print_vector(&y2);

    /* Results should be exact since all values are integers */
    IGRAPH_ASSERT(igraph_vector_all_e(&y1, &y2));

    printf("\nAdding to previous result with beta=2:\n");

    igraph_blas_dgemv(0, 1, &A, &x, 2, &y1);
    matmul(&A, &x, &y2, 2);

    printf("\nResult vector DGEMV:\n");
    print_vector(&y1);

    printf("\nResult vector naive:\n");
    print_vector(&y2);

    IGRAPH_ASSERT(igraph_vector_all_e(&y1, &y2));

    igraph_vector_destroy(&y2);
    igraph_vector_destroy(&y1);
    igraph_vector_destroy(&x);
    igraph_matrix_destroy(&A);

    VERIFY_FINALLY_STACK();

    return 0;
}
