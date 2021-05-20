/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <stdio.h>
#include <stdlib.h>

#include "test_utilities.inc"

int main() {

    igraph_t g1, g2;
    igraph_t ring1, ring2;
    igraph_vector_int_t color1, color2;
    igraph_vector_t perm;
    igraph_bool_t iso;

    igraph_bliss_sh_t sh_values[] = {
        IGRAPH_BLISS_F, IGRAPH_BLISS_FL,
        IGRAPH_BLISS_FS, IGRAPH_BLISS_FM,
        IGRAPH_BLISS_FLM, IGRAPH_BLISS_FSM
    };

    const char *sh_names[] = {
        "F", "FL", "FS", "FM", "FLM", "FSM"
    };

    size_t i;

    /* necessary because of igraph_vector_shuffle() below */
    igraph_rng_seed(igraph_rng_default(), 137);

    for (i=0; i < sizeof(sh_values) / sizeof(igraph_bliss_sh_t); ++i) {

        igraph_bliss_sh_t sh = sh_values[i];
        printf("Splitting heuristic: %s\n", sh_names[i]);

        igraph_ring(&ring1, 100, /*directed=*/ 0, /*mutual=*/ 0, /*circular=*/1);
        igraph_vector_init_seq(&perm, 0, igraph_vcount(&ring1) - 1);
        igraph_vector_shuffle(&perm);
        igraph_permute_vertices(&ring1, &ring2, &perm);

        /* Without colors */
        printf("Without vertex colors: ");
        iso = 0;
        igraph_isomorphic_bliss(&ring1, &ring2, 0, 0, &iso, 0, 0, sh, 0, 0);
        printf("%s\n", iso ? "YES" : "NO");

        /* Everything has the same colors */
        printf("All vertices having the same color: ");
        igraph_vector_int_init(&color1, igraph_vcount(&ring1));
        igraph_vector_int_init(&color2, igraph_vcount(&ring2));
        igraph_vector_int_fill(&color1, 1);
        igraph_vector_int_fill(&color2, 1);

        iso = 0;
        igraph_isomorphic_bliss(&ring1, &ring2, &color1, &color2, &iso, 0, 0, sh, 0, 0);
        printf("%s\n", iso ? "YES" : "NO");

        /* Try a negative result */
        printf("Non-matching colors 1: ");
        igraph_vector_int_fill(&color1, 0);
        igraph_vector_int_fill(&color2, 0);
        VECTOR(color1)[0] = 1;

        iso = 1;
        igraph_isomorphic_bliss(&ring1, &ring2, &color1, &color2, &iso, 0, 0, sh, 0, 0);
        printf("%s\n", iso ? "YES" : "NO");

        /* Another negative, same color distribution, different topology */
        printf("Non-matching colors 2: ");
        igraph_vector_int_fill(&color1, 0);
        igraph_vector_int_fill(&color2, 0);
        VECTOR(color1)[0] = 1;
        VECTOR(color1)[1] = 1;
        VECTOR(color2)[0] = 1;
        VECTOR(color2)[2] = 1;

        iso = 1;
        igraph_isomorphic_bliss(&ring1, &ring2, &color1, &color2, &iso, 0, 0, sh, 0, 0);
        printf("%s\n", iso ? "YES" : "NO");

        igraph_vector_int_destroy(&color1);
        igraph_vector_int_destroy(&color2);

        igraph_vector_destroy(&perm);
        igraph_destroy(&ring2);
        igraph_destroy(&ring1);

        /* More complicated test with colors */
        printf("Isomorphic colored graphs: ");

        igraph_small(&g1, 8, IGRAPH_DIRECTED,
                     0, 4, 0, 5, 0, 6, 1, 4, 1, 5, 1, 7, 2, 4, 2, 6, 2, 7, 3, 5, 3, 6, 3, 7, -1
                    );
        igraph_small(&g2, 8, IGRAPH_DIRECTED,
                     0, 1, 0, 3, 0, 4, 2, 3, 2, 1, 2, 6, 5, 1, 5, 4, 5, 6, 7, 3, 7, 6, 7, 4, -1
                    );

        igraph_vector_int_init(&color1, 8);
        igraph_vector_int_init(&color2, 8);

        VECTOR(color1)[1] = 1;
        VECTOR(color1)[3] = 1;
        VECTOR(color1)[5] = 1;
        VECTOR(color1)[7] = 1;

        VECTOR(color2)[2] = 1;
        VECTOR(color2)[3] = 1;
        VECTOR(color2)[6] = 1;
        VECTOR(color2)[7] = 1;

        iso = 0;
        igraph_isomorphic_bliss(&g1, &g2, &color1, &color2, &iso, 0, 0, sh, 0, 0);
        printf("%s\n", iso ? "YES" : "NO");

        igraph_vector_int_destroy(&color1);
        igraph_vector_int_destroy(&color2);

        igraph_destroy(&g2);
        igraph_destroy(&g1);

        printf("\n");

        VERIFY_FINALLY_STACK();

    }

    return 0;
}
