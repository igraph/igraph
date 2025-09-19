/*
  igraph library.
  Copyright (C) 2024 The igraph development team <igraph@igraph.org>

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program. If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_bool_t is_clique, is_indep_set;
    igraph_vector_int_t vids;

    igraph_small(&g, 10, IGRAPH_DIRECTED,
                 0, 2, 0, 3, 0, 4, 0, 6, 0, 7, 0, 8, 0, 9, 1, 2, 1, 4, 2, 1, 2, 3, 2, \
                 7, 2, 8, 2, 9, 3, 0, 3, 1, 3, 2, 3, 4, 3, 5, 3, 6, 3, 7, 3, 9, 4, 0, \
                 4, 2, 4, 3, 4, 5, 4, 6, 5, 1, 5, 2, 5, 4, 5, 7, 5, 9, 6, 0, 6, 1, 6, \
                 4, 6, 8, 7, 2, 7, 3, 7, 6, 7, 9, 8, 2, 8, 4, 8, 6, 8, 7, 9, 0, 9, 2, \
                 9, 3, 9, 4, 9, 7, 9, 8,
                 -1);

    /* Empty set */
    igraph_vector_int_init(&vids, 0);
    is_clique = false;
    igraph_is_clique(&g, igraph_vss_vector(&vids), true, &is_clique);
    IGRAPH_ASSERT(is_clique);
    igraph_vector_int_destroy(&vids);

    /* Singleton set */
    is_clique = false;
    igraph_is_clique(&g, igraph_vss_1(0), true, &is_clique);
    IGRAPH_ASSERT(is_clique);

    igraph_vector_int_init_int_end(&vids, -1,
                                   1, 2, -1);
    is_clique = false;
    igraph_is_clique(&g, igraph_vss_vector(&vids), true, &is_clique);
    IGRAPH_ASSERT(is_clique);
    igraph_vector_int_destroy(&vids);

    /* It is a clique: Duplicates handled? */
    igraph_vector_int_init_int_end(&vids, -1,
                                   1, 2, 2, 2, 1, -1);
    is_clique = false;
    igraph_is_clique(&g, igraph_vss_vector(&vids), true, &is_clique);
    IGRAPH_ASSERT(is_clique);
    igraph_vector_int_destroy(&vids);

    igraph_vector_int_init_int_end(&vids, -1,
                                   9, 2, 7, 3, -1);
    is_clique = false;
    igraph_is_clique(&g, igraph_vss_vector(&vids), true, &is_clique);
    IGRAPH_ASSERT(is_clique);
    igraph_vector_int_destroy(&vids);

    /* Not a clique, in either a directed or undirected sense. */
    igraph_vector_int_init_int_end(&vids, -1,
                                   3, 8, 5, -1);
    is_clique = true;
    igraph_is_clique(&g, igraph_vss_vector(&vids), true, &is_clique);
    IGRAPH_ASSERT(! is_clique);
    is_clique = true;
    igraph_is_clique(&g, igraph_vss_vector(&vids), false, &is_clique);
    IGRAPH_ASSERT(! is_clique);
    igraph_vector_int_destroy(&vids);

    /* Not a clique: Duplicates handled? */
    igraph_vector_int_init_int_end(&vids, -1,
                                   5, 3, 4, 3, 4, 5, 5, 5, -1);
    is_clique = true;
    igraph_is_clique(&g, igraph_vss_vector(&vids), true, &is_clique);
    IGRAPH_ASSERT(! is_clique);
    igraph_vector_int_destroy(&vids);

    /* This is only an undirected clique, but not a directed one. */
    igraph_vector_int_init_int_end(&vids, -1,
                                   4, 0, 8, 9, 2, -1);
    is_clique = true;
    igraph_is_clique(&g, igraph_vss_vector(&vids), true, &is_clique);
    IGRAPH_ASSERT(! is_clique);
    igraph_is_clique(&g, igraph_vss_vector(&vids), false, &is_clique);
    IGRAPH_ASSERT(is_clique);
    igraph_vector_int_destroy(&vids);

    /* Complete vertex set */

    is_clique = true;
    igraph_is_clique(&g, igraph_vss_all(), true, &is_clique);
    IGRAPH_ASSERT(! is_clique);

    is_clique = true;
    igraph_is_clique(&g, igraph_vss_all(), false, &is_clique);
    IGRAPH_ASSERT(! is_clique);

    /* Independent sets */

    /* Empty set */
    igraph_vector_int_init(&vids, 0);
    is_indep_set = false;
    igraph_is_independent_vertex_set(&g, igraph_vss_vector(&vids), &is_indep_set);
    IGRAPH_ASSERT(is_indep_set);
    igraph_vector_int_destroy(&vids);

    /* Singleton set */
    is_indep_set = false;
    igraph_is_independent_vertex_set(&g, igraph_vss_1(5), &is_indep_set);
    IGRAPH_ASSERT(is_indep_set);

    igraph_vector_int_init_int_end(&vids, -1,
                                   6, 9, -1);
    is_indep_set = false;
    igraph_is_independent_vertex_set(&g, igraph_vss_vector(&vids), &is_indep_set);
    IGRAPH_ASSERT(is_indep_set);
    igraph_vector_int_destroy(&vids);

    igraph_vector_int_init_int_end(&vids, -1,
                                   5, 6, 9, -1);
    is_indep_set = true;
    igraph_is_independent_vertex_set(&g, igraph_vss_vector(&vids), &is_indep_set);
    IGRAPH_ASSERT(! is_indep_set);
    igraph_vector_int_destroy(&vids);

    is_indep_set = true;
    igraph_is_independent_vertex_set(&g, igraph_vss_all(), &is_indep_set);
    IGRAPH_ASSERT(! is_indep_set);

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
