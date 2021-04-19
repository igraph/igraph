
#include <igraph.h>

int main() {
    igraph_matrix_t m;
    igraph_vector_t x, y, z;
    igraph_real_t xz, xx;

    igraph_vector_init_real(&x, 3, 1.0, 2.0, 3.0);
    igraph_vector_init_real(&y, 4, 4.0, 5.0, 6.0, 7.0);
    igraph_vector_init_real(&z, 3, -1.0, 0.0, 0.5);

    igraph_matrix_init(&m, 4, 3);
    MATRIX(m, 0, 0) = 1;
    MATRIX(m, 0, 1) = 2;
    MATRIX(m, 0, 2) = 3;
    MATRIX(m, 1, 0) = 2;
    MATRIX(m, 1, 1) = 3;
    MATRIX(m, 1, 2) = 4;
    MATRIX(m, 2, 0) = 3;
    MATRIX(m, 2, 1) = 4;
    MATRIX(m, 2, 2) = 5;
    MATRIX(m, 3, 0) = 4;
    MATRIX(m, 3, 1) = 5;
    MATRIX(m, 3, 2) = 6;

    /* Compute 2 m.x + 3 y and store it in y. */
    igraph_blas_dgemv(/* transpose= */ 0, /* alpha= */ 2, &m, &x, /* beta= */ 3, &y);
    igraph_vector_print(&y);

    /* Compute the squared norm of x, as well as the dor product of x and z. */
    igraph_blas_ddot(&x, &x, &xx);
    igraph_blas_ddot(&x, &z, &xz);
    printf("x.x = %g, x.z = %g\n", xx, xz);

    igraph_matrix_destroy(&m);
    igraph_vector_destroy(&z);
    igraph_vector_destroy(&y);
    igraph_vector_destroy(&x);

    return 0;
}

