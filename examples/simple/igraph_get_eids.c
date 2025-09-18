/*
   igraph library.
   Copyright (C) 2008-2025  The igraph development team <igraph@igraph.org>

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
#include <stdlib.h>

int main(void) {
    igraph_t g;
    igraph_int_t nodes = 100;
    igraph_int_t edges = 1000;
    igraph_real_t p = 3.0 / nodes;
    igraph_int_t runs = 10;
    igraph_int_t r, e, ecount;
    igraph_vector_int_t eids, pairs, path;

    /* Initialize the library. */
    igraph_setup();

    igraph_rng_seed(igraph_rng_default(), 42); /* make tests deterministic */

    igraph_vector_int_init(&pairs, edges * 2);
    igraph_vector_int_init(&path, 0);
    igraph_vector_int_init(&eids, 0);

    for (r = 0; r < runs; r++) {
        igraph_vector_int_resize(&pairs, edges * 2);
        igraph_vector_int_clear(&path);
        igraph_vector_int_clear(&eids);

        igraph_erdos_renyi_game_gnp(&g, nodes, p, IGRAPH_UNDIRECTED, IGRAPH_SIMPLE_SW, IGRAPH_EDGE_UNLABELED);
        ecount = igraph_ecount(&g);
        for (e = 0; e < edges; e++) {
            igraph_int_t edge = RNG_INTEGER(0, ecount - 1);
            VECTOR(pairs)[2 * e] = IGRAPH_FROM(&g, edge);
            VECTOR(pairs)[2 * e + 1] = IGRAPH_TO(&g, edge);
        }
        igraph_get_eids(&g, &eids, &pairs, IGRAPH_UNDIRECTED, /*error=*/ true);
        for (e = 0; e < edges; e++) {
            igraph_int_t edge = VECTOR(eids)[e];
            igraph_int_t from1 = VECTOR(pairs)[2 * e];
            igraph_int_t to1 = VECTOR(pairs)[2 * e + 1];
            igraph_int_t from2 = IGRAPH_FROM(&g, edge);
            igraph_int_t to2 = IGRAPH_TO(&g, edge);
            igraph_int_t min1 = from1 < to1 ? from1 : to1;
            igraph_int_t max1 = from1 < to1 ? to1 : from1;
            igraph_int_t min2 = from2 < to2 ? from2 : to2;
            igraph_int_t max2 = from2 < to2 ? to2 : from2;
            if (min1 != min2 || max1 != max2) {
                return 11;
            }
        }

        igraph_diameter(&g, /*weights=*/ NULL, /*res=*/ NULL, /*from=*/ NULL, /*to=*/ NULL, &path, NULL,
                        IGRAPH_UNDIRECTED, /*unconn=*/ true);
        igraph_vector_int_update(&pairs, &path);
        igraph_expand_path_to_pairs(&pairs);
        igraph_get_eids(&g, &eids, &pairs, /* directed= */ false, /*error=*/ true);
        for (e = 0; e < igraph_vector_int_size(&path) - 1; e++) {
            igraph_int_t edge = VECTOR(eids)[e];
            igraph_int_t from1 = VECTOR(path)[e];
            igraph_int_t to1 = VECTOR(path)[e + 1];
            igraph_int_t from2 = IGRAPH_FROM(&g, edge);
            igraph_int_t to2 = IGRAPH_TO(&g, edge);
            igraph_int_t min1 = from1 < to1 ? from1 : to1;
            igraph_int_t max1 = from1 < to1 ? to1 : from1;
            igraph_int_t min2 = from2 < to2 ? from2 : to2;
            igraph_int_t max2 = from2 < to2 ? to2 : from2;
            if (min1 != min2 || max1 != max2) {
                return 12;
            }
        }

        igraph_destroy(&g);
    }

    igraph_vector_int_destroy(&path);
    igraph_vector_int_destroy(&pairs);
    igraph_vector_int_destroy(&eids);

    return 0;
}
