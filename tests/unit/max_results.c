/* igraph library.
   Copyright (C) 2025  The igraph development team <igraph@igraph.org>

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

igraph_bool_t veclist_is_equal(const igraph_vector_int_list_t *a, const igraph_vector_int_list_t *b) {
    igraph_int_t n = igraph_vector_int_list_size(a);
    if (igraph_vector_int_list_size(b) != n) {
        return false;
    }
    for (igraph_int_t i=0; i < n; i++) {
        if (!igraph_vector_int_is_equal(igraph_vector_int_list_get_ptr(a, i),
                                        igraph_vector_int_list_get_ptr(b, i))) {
            return false;
        }
    }
    return true;
}

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t results_full, results_limited;
    igraph_vector_t vertex_weights;
    igraph_int_t max_results;

    igraph_vector_int_list_init(&results_full, 0);
    igraph_vector_int_list_init(&results_limited, 0);
    igraph_famous(&graph, "Zachary");
    igraph_vector_init_range(&vertex_weights, 1, igraph_vcount(&graph) + 1);

    /* Cliques */

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 0;
    igraph_cliques(&graph, &results_full, 2, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_cliques(&graph, &results_limited, 2, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 3;
    igraph_cliques(&graph, &results_full, 2, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_cliques(&graph, &results_limited, 2, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    /* Weighted cliques */

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 0;
    igraph_weighted_cliques(&graph, &vertex_weights, &results_full, false, 5, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_weighted_cliques(&graph, &vertex_weights, &results_limited, false, 5, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 5;
    igraph_weighted_cliques(&graph, &vertex_weights, &results_full, false, 4, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_weighted_cliques(&graph, &vertex_weights, &results_limited, false, 4, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    /* Maximal cliques */

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 0;
    igraph_maximal_cliques(&graph, &results_full, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_maximal_cliques(&graph, &results_limited, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 1;
    igraph_maximal_cliques(&graph, &results_full, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_maximal_cliques(&graph, &results_limited, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    /* Independent vertex sets */

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 0;
    igraph_independent_vertex_sets(&graph, &results_full, IGRAPH_UNLIMITED, 3, IGRAPH_UNLIMITED);
    igraph_independent_vertex_sets(&graph, &results_limited, IGRAPH_UNLIMITED, 3, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 10;
    igraph_independent_vertex_sets(&graph, &results_full, IGRAPH_UNLIMITED, 4, IGRAPH_UNLIMITED);
    igraph_independent_vertex_sets(&graph, &results_limited, IGRAPH_UNLIMITED, 4, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 7;
    igraph_independent_vertex_sets(&graph, &results_full, 2, 4, IGRAPH_UNLIMITED);
    igraph_independent_vertex_sets(&graph, &results_limited, 2, 4, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    /* Maximal independent vertex sets */

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 0;
    igraph_maximal_independent_vertex_sets(&graph, &results_full, IGRAPH_UNLIMITED, 8, IGRAPH_UNLIMITED);
    igraph_maximal_independent_vertex_sets(&graph, &results_limited, IGRAPH_UNLIMITED, 8, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 5;
    igraph_maximal_independent_vertex_sets(&graph, &results_full, IGRAPH_UNLIMITED, 9, IGRAPH_UNLIMITED);
    igraph_maximal_independent_vertex_sets(&graph, &results_limited, IGRAPH_UNLIMITED, 9, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 4;
    igraph_maximal_independent_vertex_sets(&graph, &results_full, 7, 10, IGRAPH_UNLIMITED);
    igraph_maximal_independent_vertex_sets(&graph, &results_limited, 7, 10, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    /* Simple cycles */

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 0;
    igraph_simple_cycles(&graph, NULL, &results_full, IGRAPH_ALL, IGRAPH_UNLIMITED, 3, IGRAPH_UNLIMITED);
    igraph_simple_cycles(&graph, NULL, &results_limited, IGRAPH_ALL, IGRAPH_UNLIMITED, 3, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 3;
    igraph_simple_cycles(&graph, NULL, &results_full, IGRAPH_ALL, IGRAPH_UNLIMITED, 4, IGRAPH_UNLIMITED);
    igraph_simple_cycles(&graph, NULL, &results_limited, IGRAPH_ALL, IGRAPH_UNLIMITED, 4, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 8;
    igraph_simple_cycles(&graph, NULL, &results_full, IGRAPH_ALL, 5, 5, IGRAPH_UNLIMITED);
    igraph_simple_cycles(&graph, NULL, &results_limited, IGRAPH_ALL, 5, 5, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    /* Simple paths */

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 0;
    igraph_get_all_simple_paths(&graph, &results_full, 0, igraph_vss_1(4), IGRAPH_ALL, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_get_all_simple_paths(&graph, &results_limited, 0, igraph_vss_1(4), IGRAPH_ALL, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 5;
    igraph_get_all_simple_paths(&graph, &results_full, 0, igraph_vss_1(4), IGRAPH_ALL, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_get_all_simple_paths(&graph, &results_limited, 0, igraph_vss_1(4), IGRAPH_ALL, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 4;
    igraph_get_all_simple_paths(&graph, &results_full, 0, igraph_vss_1(4), IGRAPH_ALL, 3, IGRAPH_UNLIMITED, IGRAPH_UNLIMITED);
    igraph_get_all_simple_paths(&graph, &results_limited, 0, igraph_vss_1(4), IGRAPH_ALL, 3, IGRAPH_UNLIMITED, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_int_list_clear(&results_full);
    igraph_vector_int_list_clear(&results_limited);
    max_results = 3;
    igraph_get_all_simple_paths(&graph, &results_full, 0, igraph_vss_1(4), IGRAPH_ALL, 3, 4, IGRAPH_UNLIMITED);
    igraph_get_all_simple_paths(&graph, &results_limited, 0, igraph_vss_1(4), IGRAPH_ALL, 3, 4, max_results);
    IGRAPH_ASSERT(igraph_vector_int_list_size(&results_limited) == max_results);
    igraph_vector_int_list_resize(&results_full, max_results);
    IGRAPH_ASSERT(veclist_is_equal(&results_limited, &results_full));

    igraph_vector_destroy(&vertex_weights);
    igraph_destroy(&graph);
    igraph_vector_int_list_destroy(&results_limited);
    igraph_vector_int_list_destroy(&results_full);

    VERIFY_FINALLY_STACK();

    return 0;
}
