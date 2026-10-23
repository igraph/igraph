/*
   igraph library.
   Copyright (C) 2008-2026  The igraph development team <igraph@igraph.org>

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
#include <math.h>

#include "test_utilities.h"

void check_projection(const igraph_t *graph,
                      const igraph_vector_bool_t *types,
                      const igraph_t *proj1,
                      const igraph_t *proj2) {
    igraph_int_t vcount1, ecount1, vcount2, ecount2;
    igraph_bipartite_projection_size(graph, types, &vcount1, &ecount1,
                                     &vcount2, &ecount2);
    IGRAPH_ASSERT(!proj1 || igraph_vcount(proj1) == vcount1);
    IGRAPH_ASSERT(!proj1 || igraph_ecount(proj1) == ecount1);
    IGRAPH_ASSERT(!proj2 || igraph_vcount(proj2) == vcount2);
    IGRAPH_ASSERT(!proj2 || igraph_ecount(proj2) == ecount2);
}

int main(void) {
    igraph_t g, p1, p2, full, ring;
    igraph_vector_bool_t types;
    igraph_bool_t iso;
    igraph_int_t i, m2 = 0, w, f, t;
    igraph_vector_int_t mult1, mult2;

    /*******************************************************/
    /* Full bipartite graph -> full graphs                 */
    /*******************************************************/

    igraph_vector_bool_init(&types, 0);
    igraph_full_bipartite(&g, &types, 5, 3, IGRAPH_UNDIRECTED, IGRAPH_ALL);

    /* Get both projections */
    igraph_bipartite_projection(&g, &types, &p1, &p2, NULL, NULL, /*probe1=*/ -1);
    check_projection(&g, &types, &p1, &p2);

    /* Check first projection */
    igraph_full(&full, igraph_vcount(&p1), IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_isomorphic_bliss(&p1, &full, NULL, NULL, &iso, NULL, NULL,
                            IGRAPH_BLISS_FM, NULL, NULL);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&full);

    /* Check second projection */
    igraph_full(&full, igraph_vcount(&p2), IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    igraph_isomorphic_bliss(&p2, &full, NULL, NULL, &iso, NULL, NULL,
                            IGRAPH_BLISS_FM, NULL, NULL);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&full);

    igraph_destroy(&p1);
    igraph_destroy(&p2);
    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);

    /*******************************************************/
    /* More sophisticated test                             */
    /*******************************************************/

    igraph_ring(&g, 100, IGRAPH_DIRECTED, /*mutual=*/ true, /*circular=*/ true);
    igraph_vector_bool_init(&types, igraph_vcount(&g));
    for (i = 0; i < igraph_vector_bool_size(&types); i++) {
        VECTOR(types)[i] = (i % 2 == 0);
    }

    /* Get both projections */
    igraph_bipartite_projection(&g, &types, &p1, &p2, NULL, NULL, /*probe1=*/ -1);
    check_projection(&g, &types, &p1, &p2);

    /* Check first projection */
    igraph_ring(&ring, igraph_vcount(&g) / 2, IGRAPH_UNDIRECTED,
                /*mutual=*/ false, /*circular=*/ true);
    igraph_isomorphic_bliss(&p1, &ring, NULL, NULL, &iso, NULL, NULL,
                            IGRAPH_BLISS_FM, NULL, NULL);
    IGRAPH_ASSERT(iso);

    /* Check second projection */
    igraph_isomorphic_bliss(&p2, &ring, NULL, NULL, &iso, NULL, NULL,
                            IGRAPH_BLISS_FM, NULL, NULL);
    IGRAPH_ASSERT(iso);
    igraph_destroy(&ring);

    igraph_destroy(&p1);
    igraph_destroy(&p2);
    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);

    /*******************************************************/
    /* Multiplicity test                                   */
    /*******************************************************/

    igraph_small(&g, 10, IGRAPH_UNDIRECTED,
                 0, 8, 1, 8, 2, 8, 3, 8, 4, 8, 4, 9, 5, 9, 6, 9, 7, 9, 0, 9,
                 -1);
    igraph_vector_bool_init(&types, igraph_vcount(&g));
    igraph_vector_bool_fill(&types, true);
    VECTOR(types)[8] = VECTOR(types)[9] = false;

    igraph_vector_int_init(&mult1, 0);
    igraph_vector_int_init(&mult2, 0);
    igraph_bipartite_projection(&g, &types, &p1, &p2, &mult1, &mult2,
                                /*probe=*/ -1);
    check_projection(&g, &types, &p1, &p2);

    IGRAPH_ASSERT(igraph_vector_int_size(&mult1) == igraph_ecount(&p1));
    IGRAPH_ASSERT(igraph_vector_int_size(&mult2) == igraph_ecount(&p2));
    IGRAPH_ASSERT(VECTOR(mult1)[0] == 2);
    for (i = 0; i < igraph_vector_int_size(&mult2); i++) {
        IGRAPH_ASSERT(VECTOR(mult2)[i] == 1 || VECTOR(mult2)[i] == 2);
        if (VECTOR(mult2)[i] == 2) {
            m2++;
            w = i;
        }
    }
    IGRAPH_ASSERT(m2 == 1);
    f = IGRAPH_FROM(&p2, w);
    t = IGRAPH_TO(&p2, w);
    IGRAPH_ASSERT(fmin(f, t) == 0 && fmax(f, t) == 4);

    igraph_vector_int_destroy(&mult1);
    igraph_vector_int_destroy(&mult2);
    igraph_destroy(&p1);
    igraph_destroy(&p2);
    igraph_destroy(&g);
    igraph_vector_bool_destroy(&types);

    VERIFY_FINALLY_STACK();
    return 0;
}
