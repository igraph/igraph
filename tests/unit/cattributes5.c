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

#include "test_utilities.inc"

int main() {

    igraph_t g, g2;
    igraph_attribute_combination_t comb;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0, 1, 0, 1, 0, 1,
                 1, 2, 2, 3,
                 -1);

    SETEAB(&g, "type", 0, 1);
    SETEAB(&g, "type", 1, 1);
    SETEAB(&g, "type", 2, 0);
    SETEAB(&g, "type", 3, 0);
    SETEAB(&g, "type", 4, 1);

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_FIRST,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_LAST,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_LAST,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_PROD,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_MIN,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_MAX,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_MEAN,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_MEDIAN,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    igraph_destroy(&g);

    VERIFY_FINALLY_STACK();

    return 0;
}
