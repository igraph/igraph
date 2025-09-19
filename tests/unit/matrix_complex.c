/* igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include "core/math.h" /* M_PI */
#include "test_utilities.h"

int main(void) {
    igraph_matrix_complex_t c;
    igraph_matrix_t real;
    igraph_matrix_t imag;
    int e[] = {1, 2, 3, 4};
    int e2[] = {5, 6, 7, 8};

    matrix_init_int_row_major(&real, 2, 2, e);
    matrix_init_int_row_major(&imag, 2, 2, e2);

    printf("Complex matrix:\n");
    igraph_matrix_complex_create(&c, &real, &imag);
    igraph_matrix_complex_fprint(&c, stdout);

    printf("Real part:\n");
    igraph_matrix_resize(&real, 0, 0);
    igraph_matrix_complex_real(&c, &real);
    igraph_matrix_print(&real);

    printf("Imaginary part:\n");
    igraph_matrix_resize(&imag, 0, 0);
    igraph_matrix_complex_imag(&c, &imag);
    igraph_matrix_print(&imag);

    printf("Real and imaginary part:\n");
    igraph_matrix_resize(&real, 0, 0);
    igraph_matrix_resize(&imag, 0, 0);
    igraph_matrix_complex_realimag(&c, &real, &imag);
    igraph_matrix_print(&real);
    igraph_matrix_print(&imag);
    igraph_matrix_complex_destroy(&c);

    igraph_matrix_destroy(&real);
    igraph_matrix_destroy(&imag);

    {
        printf("Complex matrix from polar:\n");
        igraph_matrix_complex_t p;
        igraph_real_t r_e[] = {1, 2, 3, 4};
        igraph_real_t theta_e[] = {0, .5 * M_PI, M_PI, 1.5 * M_PI};
        const igraph_matrix_t r = igraph_matrix_view(r_e, 2, 2);
        const igraph_matrix_t theta = igraph_matrix_view(theta_e, 2, 2);

        igraph_matrix_complex_create_polar(&p, &r, &theta);
        print_matrix_complex_round(&p);
        igraph_matrix_complex_destroy(&p);
    }

    VERIFY_FINALLY_STACK();

    {
        igraph_real_t e3[] = {1, 2, 3};
        igraph_real_t e4[] = {5, 6, 7, 8};
        printf("Check if unequal number of imaginary and real rows is handled correctly.\n");
        const igraph_matrix_t real = igraph_matrix_view(e3, 1, 2);
        const igraph_matrix_t imag = igraph_matrix_view(e4, 2, 2);
        CHECK_ERROR(igraph_matrix_complex_create(&c, &real, &imag), IGRAPH_EINVAL);
        CHECK_ERROR(igraph_matrix_complex_create_polar(&c, &real, &imag), IGRAPH_EINVAL);
    }
    {
        igraph_real_t e3[] = {1, 2, 3};
        igraph_real_t e4[] = {5, 6, 7, 8};
        printf("Check if unequal number of imaginary and real columns is handled correctly.\n");
        const igraph_matrix_t real = igraph_matrix_view(e3, 2, 1);
        const igraph_matrix_t imag = igraph_matrix_view(e4, 2, 2);
        CHECK_ERROR(igraph_matrix_complex_create(&c, &real, &imag), IGRAPH_EINVAL);
        CHECK_ERROR(igraph_matrix_complex_create_polar(&c, &real, &imag), IGRAPH_EINVAL);
    }

    VERIFY_FINALLY_STACK();
    return 0;
}
