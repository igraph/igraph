/* IGraph library.
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

int main(void) {
    igraph_set_attribute_table(&igraph_cattribute_table);
    igraph_t g;
    igraph_es_t es;

    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, -1);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 3, 2, -1);

    /* error test, no such edge to delete */
    CHECK_ERROR(igraph_delete_edges(&g, es), IGRAPH_EINVAL);
    if (igraph_ecount(&g) != 3) {
        return 3;
    }

    /* error test, invalid vertex ID */
    igraph_es_destroy(&es);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 10, 2, -1);
    CHECK_ERROR(igraph_delete_edges(&g, es), IGRAPH_EINVVID);
    if (igraph_ecount(&g) != 3) {
        return 5;
    }

    /* error test, invalid (odd) length */
    igraph_es_destroy(&es);
    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 0, 1, 2, -1);
    CHECK_ERROR(igraph_delete_edges(&g, es), IGRAPH_EINVAL);
    if (igraph_ecount(&g) != 3) {
        return 7;
    }

    /* no error test, deleting edges with attributes */
    // igraph_vector_int_t values;
    // igraph_vector_int_init(&values, 0);
    // for (size_t i = 0; i< igraph_ecount(&g); i++) {
    //     igraph_vector_int_push_back(&values, i);
    // }
    // igraph_cattribute_EAB_setv(&g, "test", &values);
    // igraph_vector_destroy(&values);
    igraph_cattribute_EAN_set(&g, "test", 2, 9);
    igraph_delete_edges(&g, igraph_ess_1(0)); // CHECK_NO_ERROR
    igraph_attribute_combination_t comb;
    igraph_attribute_combination_init(&comb);
    igraph_simplify(&g, true, true, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_delete_edges(&g, igraph_ess_1(1)); // CHECK_NO_ERROR
    if (igraph_ecount(&g) != 1) {
        return 9;
    }

    igraph_es_destroy(&es);
    igraph_destroy(&g);

    // no error test, deleting edges with attributes
    igraph_t g2;
    igraph_empty(&g2, 0, IGRAPH_UNDIRECTED);
    igraph_add_vertices(&g2, 8, 0);
    igraph_vector_t values;
    igraph_vector_init(&values, 0);
    for (size_t i = 0; i< igraph_vcount(&g2); i++) {
        igraph_vector_push_back(&values, i+1);
    }
    igraph_cattribute_VAN_setv(&g2, "id", &values);
    igraph_vector_destroy(&values);
    igraph_vector_int_t new_edges;
    igraph_vector_int_init(&new_edges, 14);
    igraph_vector_int_set(&new_edges, 0, 0);
    igraph_vector_int_set(&new_edges, 1, 1);
    igraph_vector_int_set(&new_edges, 2, 1);
    igraph_vector_int_set(&new_edges, 3, 2);
    igraph_vector_int_set(&new_edges, 4, 2);
    igraph_vector_int_set(&new_edges, 5, 5);
    igraph_vector_int_set(&new_edges, 6, 5);
    igraph_vector_int_set(&new_edges, 7, 4);
    igraph_vector_int_set(&new_edges, 8, 4);
    igraph_vector_int_set(&new_edges, 9, 6);
    igraph_vector_int_set(&new_edges, 10, 6);
    igraph_vector_int_set(&new_edges, 11, 0);
    igraph_vector_int_set(&new_edges, 12, 6);
    igraph_vector_int_set(&new_edges, 13, 7);
    igraph_add_edges(&g2, &new_edges, 0);

    igraph_vector_init(&values, 0);
    for (size_t i = 0; i< igraph_ecount(&g2); i++) {
        igraph_vector_push_back(&values, i == 2 ? 11 : 1);
    }
    igraph_cattribute_EAN_setv(&g2, "type", &values);
    igraph_vector_destroy(&values);
    igraph_delete_edges(&g2, igraph_ess_1(0)); // CHECK_NO_ERROR
    igraph_attribute_combination_t comb2;
    igraph_attribute_combination_init(&comb2);
    igraph_simplify(&g2, true, true, &comb2);
    igraph_attribute_combination_destroy(&comb2);
    igraph_delete_edges(&g2, igraph_ess_1(1)); // CHECK_NO_ERROR
    if (igraph_ecount(&g2) != 5) {
        return 11;
    }
    igraph_destroy(&g2);


    VERIFY_FINALLY_STACK();

    return 0;
}
