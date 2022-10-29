/*
   IGraph library.
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

/* This test attempts to verify that when adding vertices or edges to a graph fails,
 * the graph, as well as its attribute table, are left in a consistent state.
 * This test operates with the C attribute handler. Since this attribute table is
 * not directly accessible, only limited checks are done to verify its consistency.
 *
 * We use the printignore() error handler insted of ignore() as this eases debugging.
 * Thus, this test has not corresponding expected output file. All checks are done
 * directly in the code.
 */

/* Check that the graph stays in a consistent state upon failure.
 * This macro does not destroy the graph. */
#define CHECK1(funcall) \
    do { \
        igraph_integer_t vcount = igraph_vcount(&graph); \
        igraph_integer_t ecount = igraph_ecount(&graph); \
        igraph_error_handler_t *handler; \
        handler = igraph_set_error_handler(igraph_error_handler_printignore); \
        IGRAPH_ASSERT(funcall != IGRAPH_SUCCESS); \
        igraph_set_error_handler(handler); \
        verify_graph(&graph, vcount, ecount); \
    } while (0)

/* Check that there is no invalid memory access if the graph is on the finally
 * stack when failure occurs. This macro unrolls the finally stack and destroys
 * the graph. */
#define CHECK2(funcall) \
    do { \
        igraph_error_handler_t *handler; \
        IGRAPH_FINALLY(igraph_destroy, &graph); \
        handler = igraph_set_error_handler(igraph_error_handler_printignore); \
        IGRAPH_ASSERT(funcall != IGRAPH_SUCCESS); \
        igraph_set_error_handler(handler); \
    } while (0)

/* Do both checks, destroying the graph and resetting the finally stack in the process. */
#define CHECK_DESTROY(funcall) \
    do { CHECK1(funcall); CHECK2(funcall); } while (0)

void verify_graph(const igraph_t *graph, igraph_integer_t vcount, igraph_integer_t ecount) {
    IGRAPH_ASSERT(igraph_vector_int_size(&graph->from) == ecount);
    IGRAPH_ASSERT(igraph_vector_int_size(&graph->to) == ecount);
    IGRAPH_ASSERT(igraph_vector_int_size(&graph->oi) == ecount);
    IGRAPH_ASSERT(igraph_vector_int_size(&graph->ii) == ecount);
    IGRAPH_ASSERT(igraph_vector_int_size(&graph->os) == vcount + 1);
    IGRAPH_ASSERT(igraph_vector_int_size(&graph->is) == vcount + 1);
}

int main(void) {
    igraph_t graph;

    /* Enable attribute handler */
    igraph_set_attribute_table(&igraph_cattribute_table);

    /* Helper vectors */

    igraph_vector_t nvalues;
    igraph_vector_init(&nvalues, 0);

    igraph_strvector_t svalues;
    igraph_strvector_init(&svalues, 0);

    igraph_vector_bool_t bvalues;
    igraph_vector_bool_init(&bvalues, 0);

    const igraph_vector_int_t edges; /* 2 edges */
    igraph_vector_int_init_int((igraph_vector_int_t *) &edges, 4, 0,1, 1,2);

    /* attr1, numeric, 3 values */

    igraph_vector_ptr_t attr1;
    igraph_vector_ptr_init(&attr1, 1);

    igraph_attribute_record_t a1;
    a1.name = "foo";
    a1.type = IGRAPH_ATTRIBUTE_NUMERIC;

    igraph_vector_t vec_a1;
    igraph_vector_init_int(&vec_a1, 3, 1, 2, 3);
    a1.value = &vec_a1;

    igraph_vs_t vs1;
    igraph_vs_vector_small(&vs1, 0, 1, 2, -1);

    VECTOR(attr1)[0] = &a1;

    /* attr2, string, 2 values */

    igraph_vector_ptr_t attr2;
    igraph_vector_ptr_init(&attr2, 1);

    igraph_attribute_record_t a2;
    a2.name = "bar";
    a2.type = IGRAPH_ATTRIBUTE_STRING;

    igraph_strvector_t vec_a2;
    igraph_strvector_init(&vec_a2, 2);
    igraph_strvector_set(&vec_a2, 0, "one");
    igraph_strvector_set(&vec_a2, 1, "two");
    a2.value = &vec_a2;

    VECTOR(attr2)[0] = &a2;

    /* attr3, boolean, 2 values */

    igraph_vector_ptr_t attr3;
    igraph_vector_ptr_init(&attr3, 1);

    igraph_attribute_record_t a3;
    a3.name = "baz";
    a3.type = IGRAPH_ATTRIBUTE_BOOLEAN;

    igraph_vector_bool_t vec_a3;
    igraph_vector_bool_init_int(&vec_a3, 2, 1, 0);
    a3.value = &vec_a3;

    igraph_vs_t vs3;
    igraph_vs_vector_small(&vs3, 3, 4, -1);

    VECTOR(attr3)[0] = &a3;

    /* Perform tests */

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    CHECK_DESTROY(igraph_add_vertices(&graph, 2, &attr1)); /* wrong number of vertices */

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_add_vertices(&graph, 3, &attr1);
    CHECK_DESTROY(igraph_add_edges(&graph, &edges, &attr1)); /* wrong number of edges */

    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    igraph_add_vertices(&graph, 3, &attr1);
    igraph_add_edges(&graph, &edges, &attr2);
    igraph_add_vertices(&graph, 2, &attr3);

    CHECK1(igraph_add_vertices(&graph, 2, &attr1)); /* wrong number of vertices */
    /* Attribute handling still works and graph has correct attributes after failure to add vertices? */
    IGRAPH_ASSERT(igraph_cattribute_VANV(&graph, a1.name, vs1, &nvalues) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vector_all_e(&vec_a1, &nvalues));
    IGRAPH_ASSERT(igraph_cattribute_EASV(&graph, a2.name, igraph_ess_all(IGRAPH_EDGEORDER_ID), &svalues) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_strvector_size(&svalues) == igraph_strvector_size(&vec_a2));
    for (igraph_integer_t i=0; i < igraph_strvector_size(&svalues); ++i) {
        IGRAPH_ASSERT(strcmp(STR(vec_a2, i), STR(svalues, i)) == 0);
    }
    IGRAPH_ASSERT(igraph_cattribute_VABV(&graph, a3.name, vs3, &bvalues) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vector_bool_all_e(&vec_a3, &bvalues));

    CHECK1(igraph_add_edges(&graph, &edges, &attr1)); /* wrong number of edges */
    /* Attribute handling still works and graph has correct attributes after failure to add vertices? */
    IGRAPH_ASSERT(igraph_cattribute_VANV(&graph, a1.name, vs1, &nvalues) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vector_all_e(&vec_a1, &nvalues));
    IGRAPH_ASSERT(igraph_cattribute_EASV(&graph, a2.name, igraph_ess_all(IGRAPH_EDGEORDER_ID), &svalues) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_strvector_size(&svalues) == igraph_strvector_size(&vec_a2));
    for (igraph_integer_t i=0; i < igraph_strvector_size(&svalues); ++i) {
        IGRAPH_ASSERT(strcmp(STR(vec_a2, i), STR(svalues, i)) == 0);
    }
    IGRAPH_ASSERT(igraph_cattribute_VABV(&graph, a3.name, vs3, &bvalues) == IGRAPH_SUCCESS);
    IGRAPH_ASSERT(igraph_vector_bool_all_e(&vec_a3, &bvalues));

    igraph_destroy(&graph);

    /* Release attr1 */

    igraph_vector_destroy(&vec_a1);
    igraph_vector_ptr_destroy(&attr1);
    igraph_vs_destroy(&vs1);

    /* Release attr2 */

    igraph_strvector_destroy(&vec_a2);
    igraph_vector_ptr_destroy(&attr2);

    /* Release attr3 */

    igraph_vector_bool_destroy(&vec_a3);
    igraph_vector_ptr_destroy(&attr3);
    igraph_vs_destroy(&vs3);

    /* Release helpers */

    igraph_vector_int_destroy((igraph_vector_int_t *) &edges);
    igraph_strvector_destroy(&svalues);
    igraph_vector_destroy(&nvalues);
    igraph_vector_bool_destroy(&bvalues);

    VERIFY_FINALLY_STACK();

    return 0;
}
