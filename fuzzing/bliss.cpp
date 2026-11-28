/*
   igraph library.
   Copyright (C) 2021-2022  The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA
*/

#include <igraph.h>
#include <cstdlib>

// Checks that each generator of the automorphism group found by Bliss
// is indeed an automorphism. 'graph' must not have multi-edges.
void check_automorphisms(const igraph_t *graph, igraph_bliss_sh_t heur) {
    igraph_vector_int_list_t generators;
    igraph_bliss_info_t info;

    igraph_vector_int_list_init(&generators, 0);
    igraph_automorphism_group_bliss(graph, nullptr, &generators, heur, &info);
    igraph_free(info.group_size);
    for (igraph_int_t i=0; i < igraph_vector_int_list_size(&generators); i++) {
        igraph_t mapped;
        igraph_bool_t same;
        igraph_permute_vertices(graph, &mapped, igraph_vector_int_list_get_ptr(&generators, i));
        igraph_is_same_graph(graph, &mapped, &same);
        IGRAPH_ASSERT(same);
        igraph_destroy(&mapped);
    }
    igraph_vector_int_list_destroy(&generators);
}

// Checks that Bliss finds 'graph' isomorphic to a randomly permuted copy
// of itself, and that the returned mapping is an isomorphism. This also checks
// the canonical labelings, as igraph_isomorphic_bliss() compares the canonical
// forms of the two graphs. 'graph' must not have multi-edges.
void check_isomorphism(const igraph_t *graph, igraph_bliss_sh_t heur) {
    igraph_t permuted, mapped;
    igraph_vector_int_t perm, map21;
    igraph_bliss_info_t info1, info2;
    igraph_bool_t iso, same;

    igraph_vector_int_init_range(&perm, 0, igraph_vcount(graph));
    igraph_vector_int_shuffle(&perm);
    igraph_permute_vertices(graph, &permuted, &perm);
    igraph_vector_int_destroy(&perm);

    // map21[i] is the vertex of 'graph' that corresponds to vertex i of 'permuted',
    // which is the convention of igraph_permute_vertices().
    igraph_vector_int_init(&map21, 0);
    igraph_isomorphic_bliss(graph, &permuted, nullptr, nullptr, &iso, nullptr, &map21, heur, &info1, &info2);
    igraph_free(info1.group_size);
    igraph_free(info2.group_size);
    IGRAPH_ASSERT(iso);
    igraph_permute_vertices(graph, &mapped, &map21);
    igraph_is_same_graph(&mapped, &permuted, &same);
    IGRAPH_ASSERT(same);
    igraph_destroy(&mapped);
    igraph_vector_int_destroy(&map21);

    igraph_destroy(&permuted);
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    // Data[0]  - splitting heuristic (Data[0] % 6) and the check to perform (Data[0] / 6)
    // Data[1:] - edges

    igraph_t graph;
    igraph_vector_int_t edges;

    igraph_bliss_sh_t heurs[] = {
        IGRAPH_BLISS_F, IGRAPH_BLISS_FL,
        IGRAPH_BLISS_FS, IGRAPH_BLISS_FM,
        IGRAPH_BLISS_FLM, IGRAPH_BLISS_FSM
    };
    igraph_bliss_sh_t heur;

    const igraph_int_t max_vcount  = 64;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 2 == 0 || Size > 512+1 || Size < 1) {
        return 0;
    }

    // Bliss runs can be slow. To limit the running time per input, only one of
    // the two checks, each needing one or two Bliss runs, is done per input.
    const size_t heur_count = sizeof(heurs) / sizeof(heurs[0]);
    if (Data[0] >= 2 * heur_count) {
        return 0;
    }

    heur = heurs[Data[0] % heur_count];
    const bool isomorphism_mode = Data[0] >= heur_count;

    igraph_vector_int_init(&edges, Size-1);
    for (size_t i=0; i < Size-1; ++i) {
        VECTOR(edges)[i] = Data[i+1];
    }

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Undirected */
    if (igraph_create(&graph, &edges, 0, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS) {
        if (igraph_vcount(&graph) <= max_vcount) {
            igraph_bool_t multi;

            igraph_has_multiple(&graph, &multi);

            /* Bliss does not support multigraphs and the input is currently not checked */
            if (! multi) {
                if (isomorphism_mode) {
                    check_isomorphism(&graph, heur);
                } else {
                    check_automorphisms(&graph, heur);
                }
            }
        }

        igraph_destroy(&graph);
    }

    /* Directed */
    if (igraph_create(&graph, &edges, 0, IGRAPH_DIRECTED) == IGRAPH_SUCCESS) {
        if (igraph_vcount(&graph) <= max_vcount) {
            igraph_bool_t multi;

            igraph_has_multiple(&graph, &multi);

            /* Bliss does not support multigraphs and the input is currently not checked */
            if (! multi) {
                if (isomorphism_mode) {
                    check_isomorphism(&graph, heur);
                } else {
                    check_automorphisms(&graph, heur);
                }
            }
        }

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
