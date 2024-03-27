/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include "igraph_adjlist.h"
#include "igraph_community.h"
#include "igraph_components.h"
#include "igraph_interface.h"
#include "igraph_reachability.h"

void _print_vector_int(const igraph_vector_int_t *v) {
    igraph_integer_t i, n = igraph_vector_int_size(v);
    printf("(");
    for (i=0; i < n; i++) {
        printf(" %" IGRAPH_PRId, VECTOR(*v)[i]);
    }
    printf(" )\n");
}

igraph_integer_t igraph_bitset_popcount(igraph_vector_int_t *bitset)
{
    igraph_integer_t i, n, count = 0;
    n = igraph_vector_int_size(bitset);
    for (i = 0; i < n; ++i)
    {
        // TODO: Check if __cpuid claims the operation is supported by CPU
        // https://learn.microsoft.com/en-us/cpp/intrinsics/popcnt16-popcnt-popcnt64?view=msvc-170
        #ifdef _MSC_VER
        count += __popcnt64(VECTOR(*bitset)[i]);
        #else
        count += __builtin_popcountll(VECTOR(*bitset)[i]);
        //printf("Current count at %ld: %ld\n", i, (igraph_integer_t)__builtin_popcountll(VECTOR(*bitset)[i]));
        #endif
    }
    return count;
}

igraph_error_t igraph_reachability_directed(
    const igraph_t *graph,
    igraph_vector_int_t *membership,
    igraph_vector_int_t *csize,
    igraph_integer_t *no_of_components,
    igraph_vector_int_list_t *reach)
{
    igraph_integer_t no_of_nodes;
    igraph_integer_t i, j, k, n;
    igraph_adjlist_t adjlist, dag;
    igraph_vector_int_t *dag_neighbours, *neighbours, *from_bitset, *to_bitset;

    no_of_nodes = igraph_vcount(graph);

    igraph_connected_components(graph, membership, csize, no_of_components, IGRAPH_STRONG);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_OUT, IGRAPH_LOOPS_ONCE, IGRAPH_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_vector_int_list_init(reach, *no_of_components));
    IGRAPH_FINALLY(igraph_vector_int_list_destroy, reach);

    for (i = 0; i < *no_of_components; ++i)
    {
        IGRAPH_CHECK(igraph_vector_int_resize(igraph_vector_int_list_get_ptr(reach, i), BITNSLOTS(no_of_nodes)));
    }
    for (i = 0; i < no_of_nodes; ++i)
    {
        BITSET(VECTOR(*igraph_vector_int_list_get_ptr(reach, VECTOR(*membership)[i])), i);
    }

    IGRAPH_CHECK(igraph_adjlist_init_empty(&dag, *no_of_components));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &dag);

    for (i = 0; i < no_of_nodes; ++i)
    {
        neighbours = igraph_adjlist_get(&adjlist, i);
        dag_neighbours = igraph_adjlist_get(&dag, VECTOR(*membership)[i]);
        n = igraph_vector_int_size(neighbours);
        for (j = 0; j < n; j++)
        {
            if (VECTOR(*membership)[i] != VECTOR(*membership)[VECTOR(*neighbours)[j]])
            {
                igraph_vector_int_push_back(dag_neighbours, VECTOR(*membership)[VECTOR(*neighbours)[j]]);
            }
        }
    }

    for (i = *no_of_components - 1; i >= 0; --i)
    {
        dag_neighbours = igraph_adjlist_get(&dag, i);
        from_bitset = igraph_vector_int_list_get_ptr(reach, i);
        n = igraph_vector_int_size(dag_neighbours);
        for (j = 0; j < n; j++)
        {
            to_bitset = igraph_vector_int_list_get_ptr(reach, VECTOR(*dag_neighbours)[j]);
            for (k = 0; k < BITNSLOTS(no_of_nodes); ++k)
            {
                VECTOR(*from_bitset)[k] |= VECTOR(*to_bitset)[k];
            }
        }
    }

    igraph_adjlist_destroy(&adjlist);
    igraph_adjlist_destroy(&dag);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}
