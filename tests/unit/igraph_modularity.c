/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_vector_t weights;
    igraph_vector_int_t membership;
    igraph_real_t modularity, resolution;
    igraph_attribute_combination_t comb;

    /* Turn on attribute handling */
    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_attribute_combination(&comb,
                                 "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 IGRAPH_NO_MORE_ATTRIBUTES);

    /* Set default seed to get reproducible results */
    igraph_rng_seed(igraph_rng_default(), 0);

    /* Null graph */
    igraph_vector_int_init(&membership, 0);
    igraph_small(&graph, 0, IGRAPH_UNDIRECTED, -1);
    igraph_modularity(&graph, &membership, 0, /* resolution */ 1, IGRAPH_UNDIRECTED, &modularity);
    IGRAPH_ASSERT(isnan(modularity));

    igraph_destroy(&graph);
    igraph_small(&graph, 0, IGRAPH_DIRECTED, -1);
    igraph_modularity(&graph, &membership, 0, /* resolution */ 1, IGRAPH_UNDIRECTED, &modularity);
    IGRAPH_ASSERT(isnan(modularity));

    /* Should not crash if we omit 'modularity' */
    igraph_modularity(&graph, &membership, 0, /* resolution */ 1, IGRAPH_UNDIRECTED, /* modularity = */ NULL);
    igraph_destroy(&graph);
    igraph_vector_int_destroy(&membership);

    /* Simple unweighted graph */
    igraph_small(&graph, 10, IGRAPH_UNDIRECTED,
                 0, 1, 0, 2, 0, 3, 0, 4, 1, 2, 1, 3, 1, 4, 2, 3, 2, 4, 3, 4,
                 5, 6, 5, 7, 5, 8, 5, 9, 6, 7, 6, 8, 6, 9, 7, 8, 7, 9, 8, 9,
                 0, 5, -1);

    /* Set weights */
    igraph_vector_init(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1.0);
    SETEANV(&graph, "weight", &weights);

    /* Set membership */
    igraph_vector_int_init_int(&membership, igraph_vcount(&graph),
                               0, 0, 0, 0, 0, 1, 1, 1, 1, 1);

    /* Calculate modularity */
    for (resolution = 0.5; resolution <= 1.5; resolution += 0.5) {
        igraph_modularity(&graph, &membership, &weights,
                        /* resolution */ resolution,
                        IGRAPH_DIRECTED, &modularity);
        printf("Modularity (resolution %.2f) is %g.\n", resolution, modularity);
    }

    igraph_to_directed(&graph, IGRAPH_TO_DIRECTED_MUTUAL);
    igraph_vector_resize(&weights, igraph_ecount(&graph));
    igraph_vector_fill(&weights, 1.0);
    for (resolution = 0.5; resolution <= 1.5; resolution += 0.5) {
        igraph_modularity(&graph, &membership, &weights,
                        /* resolution */ resolution,
                        IGRAPH_DIRECTED, &modularity);
        printf("Modularity (resolution %.2f) is %g on directed graph.\n", resolution, modularity);
    }

    /* Recalculate modularity on contracted graph */
    igraph_contract_vertices(&graph, &membership, NULL);
    igraph_vector_int_destroy(&membership);

    igraph_simplify(&graph, /* multiple */ true, /* loops */ false, &comb);

    igraph_vector_int_init_range(&membership, 0, igraph_vcount(&graph));
    EANV(&graph, "weight", &weights);
    for (resolution = 0.5; resolution <= 1.5; resolution += 0.5) {
        igraph_modularity(&graph, &membership, &weights,
                        /* resolution */ resolution,
                        IGRAPH_DIRECTED, &modularity);
        printf("Modularity (resolution %.2f) is %g after aggregation.\n", resolution, modularity);
    }

    igraph_vector_int_destroy(&membership);
    igraph_vector_destroy(&weights);
    igraph_destroy(&graph);
    igraph_attribute_combination_destroy(&comb);

    VERIFY_FINALLY_STACK();

    return 0;
}
