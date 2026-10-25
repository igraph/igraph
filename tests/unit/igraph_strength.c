/* igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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

/* Computes the strength of 'vids' both without weights and with unit weights,
 * prints the former, and checks that the two agree. */
void check_strength(const igraph_t *g, igraph_vs_t vids,
                    igraph_neimode_t mode, igraph_loops_t loops) {
    igraph_vector_t s, sw, weights;

    igraph_vector_init(&s, 0);
    igraph_vector_init(&sw, 0);
    igraph_vector_init(&weights, igraph_ecount(g));
    igraph_vector_fill(&weights, 1);

    igraph_strength(g, &s, vids, mode, loops, NULL);
    igraph_strength(g, &sw, vids, mode, loops, &weights);
    print_vector(&s);
    IGRAPH_ASSERT(igraph_vector_all_e(&s, &sw));

    igraph_vector_destroy(&weights);
    igraph_vector_destroy(&sw);
    igraph_vector_destroy(&s);
}

int main(void) {
    igraph_t g;
    igraph_vector_t s, weights;
    igraph_vs_t vids;

    igraph_vector_init(&s, 0);

    printf("Null graph:\n");
    igraph_empty(&g, 0, IGRAPH_UNDIRECTED);
    check_strength(&g, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS);
    igraph_destroy(&g);

    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0,1, 1,2, 2,3, 2,2, 3,2, 2,3,
                 -1);

    printf("\nAll vertices, out, loops twice:\n");
    check_strength(&g, igraph_vss_all(), IGRAPH_OUT, IGRAPH_LOOPS_TWICE);

    printf("\nAll vertices, all, no loops:\n");
    check_strength(&g, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

    /* The unweighted case used to return a result of length vcount
     * regardless of 'vids'. */
    printf("\nSingle vertex 2, in, loops once:\n");
    check_strength(&g, igraph_vss_1(2), IGRAPH_IN, IGRAPH_LOOPS_ONCE);

    printf("\nVertices 3, 2, 3, all, loops twice:\n");
    igraph_vs_vector_small(&vids, 3, 2, 3, -1);
    check_strength(&g, vids, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
    igraph_vs_destroy(&vids);

    printf("\nMore vertices than the vertex count:\n");
    igraph_vs_vector_small(&vids, 0, 1, 2, 3, 2, 1, 0, -1);
    check_strength(&g, vids, IGRAPH_ALL, IGRAPH_LOOPS_TWICE);
    igraph_vs_destroy(&vids);

    printf("\nNo vertices:\n");
    check_strength(&g, igraph_vss_none(), IGRAPH_ALL, IGRAPH_LOOPS_TWICE);

    printf("\nWeighted, vertices 2, 0:\n");
    igraph_vector_init_int(&weights, 6, 1, 2, 3, 4, 5, 6);
    igraph_vs_vector_small(&vids, 2, 0, -1);
    igraph_strength(&g, &s, vids, IGRAPH_ALL, IGRAPH_LOOPS_TWICE, &weights);
    print_vector(&s);
    igraph_vs_destroy(&vids);

    printf("\nInvalid vertex:\n");
    CHECK_ERROR(igraph_strength(&g, &s, igraph_vss_1(4), IGRAPH_ALL, IGRAPH_LOOPS, NULL), IGRAPH_EINVVID);
    CHECK_ERROR(igraph_strength(&g, &s, igraph_vss_1(4), IGRAPH_ALL, IGRAPH_LOOPS, &weights), IGRAPH_EINVVID);
    printf("OK\n");

    printf("\nInvalid weight vector length:\n");
    igraph_vector_resize(&weights, 5);
    CHECK_ERROR(igraph_strength(&g, &s, igraph_vss_all(), IGRAPH_ALL, IGRAPH_LOOPS, &weights), IGRAPH_EINVAL);
    printf("OK\n");

    igraph_vector_destroy(&weights);
    igraph_destroy(&g);
    igraph_vector_destroy(&s);

    VERIFY_FINALLY_STACK();

    return 0;
}
