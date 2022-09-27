/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

igraph_error_t mf(const igraph_vector_t *input, igraph_real_t *output) {
    *output = 0.0;
    return IGRAPH_SUCCESS;
}

static void simplify_write_destroy(igraph_t *g, igraph_attribute_combination_t *comb) {
    igraph_simplify(g, /*multiple=*/ true, /*loops=*/ true, comb);
    igraph_write_graph_graphml(g, stdout, /*prefixattr=*/ true);
    igraph_attribute_combination_destroy(comb);
    igraph_destroy(g);
}

static void weight_test(igraph_t *g, igraph_attribute_combination_type_t weight_attr) {
    igraph_t g2;
    igraph_attribute_combination_t comb;

    igraph_copy(&g2, g);
    igraph_attribute_combination(&comb,
                                 "weight", weight_attr,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    simplify_write_destroy(&g2, &comb);
}

int main(void) {

    igraph_t g, g2;
    igraph_vector_t weight;
    igraph_attribute_combination_t comb;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0, 1, 0, 1, 0, 1,
                 1, 2, 2, 3,
                 -1);

    igraph_vector_init_range(&weight, 1, igraph_ecount(&g) + 1);
    SETEANV(&g, "weight", &weight);
    igraph_vector_destroy(&weight);

    weight_test(&g, IGRAPH_ATTRIBUTE_COMBINE_SUM);
    weight_test(&g, IGRAPH_ATTRIBUTE_COMBINE_PROD);
    weight_test(&g, IGRAPH_ATTRIBUTE_COMBINE_MIN);
    weight_test(&g, IGRAPH_ATTRIBUTE_COMBINE_MAX);
    weight_test(&g, IGRAPH_ATTRIBUTE_COMBINE_FIRST);
    weight_test(&g, IGRAPH_ATTRIBUTE_COMBINE_LAST);
    weight_test(&g, IGRAPH_ATTRIBUTE_COMBINE_MEAN);

    /* ****************************************************** */

    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "weight", IGRAPH_ATTRIBUTE_COMBINE_FUNCTION, mf,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    simplify_write_destroy(&g2, &comb);

    /* ****************************************************** */

    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_MEAN,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    simplify_write_destroy(&g2, &comb);

    /* ****************************************************** */

    igraph_destroy(&g);

    return 0;
}
