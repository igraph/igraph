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
#include "test_utilities.h"

void all_degs(const igraph_t *g, igraph_vector_int_t *res, igraph_neimode_t mode, igraph_loops_t loops) {
    igraph_int_t n = igraph_vcount(g);
    igraph_vector_int_resize(res, n);
    for (igraph_int_t i = 0; i < n; i++) {
        igraph_degree_1(g, &VECTOR(*res)[i], i, mode, loops);
    }
}

int main(void) {
    igraph_vector_int_t v, v2, seq;
    igraph_t g;

    igraph_vector_int_init(&seq, 3);
    igraph_vector_int_init(&v, 0);
    igraph_vector_int_init(&v2, 0);
    VECTOR(seq)[0] = 2;
    VECTOR(seq)[1] = 0;
    VECTOR(seq)[2] = 2;

    igraph_small(&g, 4, IGRAPH_DIRECTED, 0,1, 1,2, 2,3, 2,2, 3,2, 2,3, -1);

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf("out, loops twice: ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS_TWICE);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_OUT, IGRAPH_LOOPS_TWICE);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf(" in, loops twice: ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS_TWICE);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_IN, IGRAPH_LOOPS_TWICE);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf("all, loops twice: ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf("out, loops once:  ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS_ONCE);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_OUT, IGRAPH_LOOPS_ONCE);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf(" in, loops once:  ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_LOOPS_ONCE);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_IN, IGRAPH_LOOPS_ONCE);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf("all, loops once:  ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_ALL, IGRAPH_LOOPS_ONCE);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf("out, no loops:    ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_OUT, IGRAPH_NO_LOOPS);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_OUT, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf(" in, no loops:    ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_IN, IGRAPH_NO_LOOPS);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_IN, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    igraph_vector_int_clear(&v);
    igraph_vector_int_clear(&v2);
    printf("all, no loops:    ");
    igraph_degree(&g, &v, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    print_vector_int(&v);
    all_degs(&g, &v2, IGRAPH_ALL, IGRAPH_NO_LOOPS);
    IGRAPH_ASSERT(igraph_vector_int_all_e(&v, &v2));

    /* Invalid mode */
    CHECK_ERROR(igraph_degree(&g, &v, igraph_vss_vector(&seq), (igraph_neimode_t)0,
                        IGRAPH_LOOPS), IGRAPH_EINVMODE);

    /* Vertex does not exist */
    VECTOR(seq)[0] = 4;
    CHECK_ERROR(igraph_degree(&g, &v, igraph_vss_vector(&seq), IGRAPH_ALL, IGRAPH_LOOPS), IGRAPH_EINVVID);

    igraph_destroy(&g);
    igraph_vector_int_destroy(&v2);
    igraph_vector_int_destroy(&v);
    igraph_vector_int_destroy(&seq);

    VERIFY_FINALLY_STACK();

    return 0;
}
