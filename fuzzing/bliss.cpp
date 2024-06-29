/*
   IGraph library.
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

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    // Data[0]  - splitting heuristic
    // Data[1:] - edges

    igraph_t graph;
    igraph_vector_int_t edges;

    igraph_bliss_sh_t heurs[] = {
        IGRAPH_BLISS_F, IGRAPH_BLISS_FL,
        IGRAPH_BLISS_FS, IGRAPH_BLISS_FM,
        IGRAPH_BLISS_FLM, IGRAPH_BLISS_FSM
    };
    igraph_bliss_sh_t heur;

    const igraph_integer_t max_vcount  = 64;

    igraph_set_warning_handler(igraph_warning_handler_ignore);

    if (Size % 2 == 0 || Size > 512+1 || Size < 1) {
        return 0;
    }

    if (Data[0] >= sizeof(heurs) / sizeof(heurs[0])) {
        return 0;
    }

    heur = (igraph_bliss_sh_t) Data[0];

    igraph_vector_int_init(&edges, Size-1);
    for (size_t i=0; i < Size-1; ++i) {
        VECTOR(edges)[i] = Data[i+1];
    }

    /* Undirected */
    if (igraph_create(&graph, &edges, 0, IGRAPH_UNDIRECTED) == IGRAPH_SUCCESS) {
        if (igraph_vcount(&graph) <= max_vcount) {
            igraph_bool_t multi;

            igraph_has_multiple(&graph, &multi);

            /* Bliss does not support multigraphs and the input is currently not checked */
            if (! multi) {
                igraph_bliss_info_t info;
                igraph_vector_int_list_t generators;
                igraph_vector_int_list_init(&generators, 0);
                igraph_automorphism_group(&graph, nullptr, &generators, heur, &info);
                igraph_free(info.group_size);
                igraph_vector_int_list_destroy(&generators);
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
                igraph_bliss_info_t info;
                igraph_vector_int_list_t generators;
                igraph_vector_int_list_init(&generators, 0);
                igraph_automorphism_group(&graph, nullptr, &generators, heur, &info);
                igraph_free(info.group_size);
                igraph_vector_int_list_destroy(&generators);
            }

        }

        igraph_destroy(&graph);
    }

    igraph_vector_int_destroy(&edges);

    IGRAPH_ASSERT(IGRAPH_FINALLY_STACK_EMPTY);

    return 0;  // Non-zero return values are reserved for future use.
}
