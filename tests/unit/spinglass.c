/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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

#include "test_utilities.inc"

int main() {
    igraph_t g;
    igraph_real_t  modularity, temperature;
    igraph_vector_t membership, csize;
    /* long int i; */
    igraph_real_t cohesion, adhesion;
    igraph_integer_t inner_links;
    igraph_integer_t outer_links;

    igraph_rng_seed(igraph_rng_default(), 137);

    /* Two 5-cliques connected by a single edge */
    igraph_small(&g, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9, 0, 5, -1);
    igraph_vector_init(&membership, 0);
    igraph_vector_init(&csize, 0);

    printf("\nOriginal implementation.\n");
    igraph_community_spinglass(&g,
                               NULL, /* no weights */
                               &modularity,
                               &temperature,
                               &membership,
                               &csize,
                               10,   /* no of spins */
                               0,    /* parallel update */
                               1.0,  /* start temperature */
                               0.01, /* stop temperature */
                               0.99, /* cooling factor */
                               IGRAPH_SPINCOMM_UPDATE_CONFIG,
                               1.0, /* gamma */
                               IGRAPH_SPINCOMM_IMP_ORIG,
                               /*gamma_minus =*/ 0);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_size(&csize) == igraph_vector_max(&membership) + 1);

    /* The following depend on the random seed, however, for this graph,
       the result is almost always the same (i.e. two clusters). */
    printf("Modularity: %g\n", modularity);
    print_vector_round(&membership);

    printf("\nOriginal implementation, parallel updating.\n");
    igraph_community_spinglass(&g,
                               NULL, /* no weights */
                               &modularity,
                               &temperature,
                               &membership,
                               &csize,
                               10,   /* no of spins */
                               1,    /* parallel update */
                               1.0,  /* start temperature */
                               0.01, /* stop temperature */
                               0.99, /* cooling factor */
                               IGRAPH_SPINCOMM_UPDATE_CONFIG,
                               1.0, /* gamma */
                               IGRAPH_SPINCOMM_IMP_ORIG,
                               /*gamma_minus =*/ 0);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_size(&csize) == igraph_vector_max(&membership) + 1);

    /* The following depend on the random seed, however, for this graph,
       the result is almost always the same (i.e. two clusters). */
    printf("Modularity: %g\n", modularity);
    print_vector_round(&membership);

    printf("\nNegative implementation.\n");
    igraph_community_spinglass(&g,
                               NULL, /* no weights */
                               &modularity,
                               &temperature,
                               &membership,
                               &csize,
                               10,   /* no of spins */
                               0,    /* parallel update */
                               1.0,  /* start temperature */
                               0.01, /* stop temperature */
                               0.99, /* cooling factor */
                               IGRAPH_SPINCOMM_UPDATE_CONFIG,
                               1.0, /* gamma */
                               IGRAPH_SPINCOMM_IMP_NEG,
                               /*gamma_minus =*/ 0);

    IGRAPH_ASSERT(igraph_vector_size(&membership) == igraph_vcount(&g));
    IGRAPH_ASSERT(igraph_vector_size(&csize) == igraph_vector_max(&membership) + 1);

    /* The following depend on the random seed, however, for this graph,
       the result is almost always the same (i.e. two clusters). */
    printf("Modularity: %g\n", modularity);
    print_vector_round(&membership);

    /* Try to call this as well, we don't check the results currently.... */

    igraph_community_spinglass_single(&g,
                                      /*weights=  */ 0,
                                      /*vertex=   */ 0,
                                      /*community=*/ &membership,
                                      /*cohesion= */ &cohesion,
                                      /*adhesion= */ &adhesion,
                                      /*inner_links= */ &inner_links,
                                      /*outer_links= */ &outer_links,
                                      /*spins=       */ 2,
                                      /*update_rule= */ IGRAPH_SPINCOMM_UPDATE_CONFIG,
                                      /*gamma=       */ 1.0);

    igraph_destroy(&g);
    igraph_vector_destroy(&membership);
    igraph_vector_destroy(&csize);

    VERIFY_FINALLY_STACK();
    return 0;
}
