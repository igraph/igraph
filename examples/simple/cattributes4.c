/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2021  The igraph development team

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

int main() {

    igraph_t g, g2;
    igraph_attribute_combination_t comb;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_small(&g, 4, IGRAPH_DIRECTED,
                 0, 1, 0, 1, 0, 1,
                 1, 2, 2, 3,
                 -1);

    SETEAS(&g, "color", 0, "green");
    SETEAS(&g, "color", 1, "red");
    SETEAS(&g, "color", 2, "blue");
    SETEAS(&g, "color", 3, "white");
    SETEAS(&g, "color", 4, "black");

    /* ****************************************************** */
    igraph_copy(&g2, &g);
    igraph_attribute_combination(&comb,
                                 "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 "color",  IGRAPH_ATTRIBUTE_COMBINE_FIRST,
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
                                 "color",  IGRAPH_ATTRIBUTE_COMBINE_CONCAT,
                                 IGRAPH_NO_MORE_ATTRIBUTES);
    igraph_simplify(&g2, /*multiple=*/ 1, /*loops=*/ 1, &comb);
    igraph_attribute_combination_destroy(&comb);
    igraph_write_graph_graphml(&g2, stdout, /*prefixattr=*/ 1);
    igraph_destroy(&g2);
    /* ****************************************************** */

    igraph_destroy(&g);

    return 0;
}
