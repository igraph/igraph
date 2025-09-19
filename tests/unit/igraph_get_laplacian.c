/*
   igraph library.
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
#include "test_utilities.h"

void test_laplacian(igraph_t *g, const igraph_vector_t *w, igraph_bool_t dir, igraph_laplacian_normalization_t normalization) {
    igraph_matrix_t m, m_converted_sparse;
    igraph_sparsemat_t m_sparse;

    igraph_matrix_init(&m, 0, 0);
    igraph_matrix_init(&m_converted_sparse, 0, 0);
    igraph_sparsemat_init(&m_sparse, 0, 0, 0);

    igraph_get_laplacian(g, &m, IGRAPH_OUT, normalization, w);
    igraph_get_laplacian_sparse(g, &m_sparse, IGRAPH_OUT, normalization, w);
    igraph_matrix_print(&m);
    igraph_sparsemat_as_matrix(&m_converted_sparse, &m_sparse);
    IGRAPH_ASSERT(igraph_matrix_all_almost_e(&m, &m_converted_sparse, 1e-15));

    igraph_matrix_destroy(&m);
    igraph_matrix_destroy(&m_converted_sparse);
    igraph_sparsemat_destroy(&m_sparse);
}

int main(void) {
    igraph_t g_e_un, g_e_dir, g_un, g_dir;
    igraph_vector_t weights, weights_e;
    char *n[]  = {"none", "symmetric", "left", "right"};

    igraph_vector_init_int(&weights, 9, 1, 2, 3, 4, 5, 2, 2, 3, 3);
    igraph_vector_init_int(&weights_e, 0);
    igraph_small(&g_e_un, 0, IGRAPH_UNDIRECTED, -1);
    igraph_small(&g_e_dir, 0, IGRAPH_DIRECTED, -1);
    igraph_small(&g_un, 6, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,3, 3,4, 4,0, 1,1, 2,2, 1,2, 3,4, -1);
    igraph_small(&g_dir, 6, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 3,4, 4,0, 1,1, 2,2, 1,2, 3,4, -1);

    for (int normalization = 0; normalization < 4; normalization++) {
        for (int weighted = 0; weighted < 2; weighted++) {
            for (int directed = 0; directed < 2; directed++) {
                printf("=== normalization: %s, %sweighted, %sdirected\n",
                        n[normalization],
                        (weighted ? "" : "un"),
                        (directed ? "" : "un")
                      );

                test_laplacian(directed ? &g_e_dir : &g_e_un, weighted ? &weights_e : NULL, directed, (igraph_laplacian_normalization_t) normalization);
                test_laplacian(directed ? &g_dir : &g_un, weighted ? &weights : NULL, directed, (igraph_laplacian_normalization_t) normalization);
            }
        }
    }

    igraph_vector_destroy(&weights_e);
    igraph_destroy(&g_e_un);
    igraph_destroy(&g_e_dir);
    igraph_vector_destroy(&weights);
    igraph_destroy(&g_un);
    igraph_destroy(&g_dir);

    VERIFY_FINALLY_STACK();

    {
        igraph_t g;
        igraph_matrix_t m;
        igraph_sparsemat_t m_sparse;

        igraph_matrix_init(&m, 0, 0);
        igraph_sparsemat_init(&m_sparse, 0, 0, 0);

        printf("Check errors for zero out/in and non-zero in/out degree.\n");
        /* For IGRAPH_OUT, SYMMETRIC and RIGHT should fail,
           for IGRAPH_IN, SYMMETRIC and LEFT should fail. */

        igraph_small(&g, 2, IGRAPH_DIRECTED, 0,1, -1);
        for (int mode_n = 0; mode_n < 2; mode_n++) {
            for (int lsr = mode_n; lsr < 2 + mode_n; lsr++) {
                igraph_neimode_t mode = mode_n ? IGRAPH_OUT : IGRAPH_IN;
                igraph_laplacian_normalization_t normalization;
                if (lsr == 0) {
                    normalization = IGRAPH_LAPLACIAN_LEFT;
                } else if (lsr == 1) {
                    normalization = IGRAPH_LAPLACIAN_SYMMETRIC;
                } else {
                    normalization = IGRAPH_LAPLACIAN_RIGHT;
                }
                CHECK_ERROR(igraph_get_laplacian(&g, &m, mode,
                            normalization, NULL), IGRAPH_EINVAL);
                CHECK_ERROR(igraph_get_laplacian_sparse(&g, &m_sparse, mode,
                            normalization, NULL), IGRAPH_EINVAL);;
            }
        }
        igraph_destroy(&g);
        igraph_matrix_destroy(&m);
        igraph_sparsemat_destroy(&m_sparse);
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
