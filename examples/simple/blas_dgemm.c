#include <igraph.h>

int main(void) {
    igraph_matrix_t a, b, c;

    igraph_matrix_init(&a, 2, 2);
    MATRIX(a, 0, 0) = 1;
    MATRIX(a, 0, 1) = 2;
    MATRIX(a, 1, 0) = 3;
    MATRIX(a, 1, 1) = 4;

    igraph_matrix_init(&b, 2, 2);
    MATRIX(b, 0, 0) = 5;
    MATRIX(b, 0, 1) = 6;
    MATRIX(b, 1, 0) = 7;
    MATRIX(b, 1, 1) = 8;

    igraph_matrix_init(&c, 2, 2);

    igraph_blas_dgemm(1, 1, 0.5, &a, &b, 0, &c);
    igraph_matrix_printf(&c, "%g");

    igraph_matrix_destroy(&a);
    igraph_matrix_destroy(&b);
    igraph_matrix_destroy(&c);
    return 0;
}
