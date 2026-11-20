/*
   igraph library.
   Copyright (C) 2026  The igraph development team <igraph@igraph.org>

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
#include <string.h>

#include "test_utilities.h"

/* Combination function for string attributes: returns the smallest string. */
igraph_error_t smallest_string(const igraph_strvector_t *input, char **output) {
    igraph_int_t n = igraph_strvector_size(input);
    const char *smallest = "";
    for (igraph_int_t i = 0; i < n; i++) {
        const char *str = igraph_strvector_get(input, i);
        if (i == 0 || strcmp(str, smallest) < 0) {
            smallest = str;
        }
    }
    *output = IGRAPH_CALLOC(strlen(smallest) + 1, char);
    IGRAPH_ASSERT(*output != NULL);
    strcpy(*output, smallest);
    return IGRAPH_SUCCESS;
}

/* Contracts vertices 0-3 into vertex 1, 4-5 into vertex 0, and 6 into vertex 2,
 * then checks the combined values of the "name" attribute. Vertex 1 is formed
 * from more vertices than the contracted graph has. */
void contract_and_check(igraph_attribute_combination_type_t type, igraph_function_pointer_t func,
                        const char *expected0, const char *expected1, const char *expected2) {
    igraph_t g;
    igraph_attribute_combination_t comb;
    igraph_vector_int_t mapping;
    igraph_int_t mappings[7] = { 1, 1, 1, 1, 0, 0, 2 };

    igraph_small(&g, 7, IGRAPH_UNDIRECTED, -1);
    SETVAS(&g, "name", 0, "x");
    SETVAS(&g, "name", 1, "x");
    SETVAS(&g, "name", 2, "x");
    SETVAS(&g, "name", 3, "x");
    SETVAS(&g, "name", 4, "y");
    SETVAS(&g, "name", 5, "y");
    SETVAS(&g, "name", 6, "z");

    igraph_attribute_combination_init(&comb);
    igraph_attribute_combination_add(&comb, "name", type, func);
    mapping = igraph_vector_int_view(mappings, 7);
    igraph_contract_vertices(&g, &mapping, &comb);

    IGRAPH_ASSERT(igraph_vcount(&g) == 3);
    IGRAPH_ASSERT(strcmp(VAS(&g, "name", 0), expected0) == 0);
    IGRAPH_ASSERT(strcmp(VAS(&g, "name", 1), expected1) == 0);
    IGRAPH_ASSERT(strcmp(VAS(&g, "name", 2), expected2) == 0);

    igraph_attribute_combination_destroy(&comb);
    igraph_destroy(&g);
}

int main(void) {
    igraph_set_attribute_table(&igraph_cattribute_table);

    contract_and_check(IGRAPH_ATTRIBUTE_COMBINE_CONCAT, NULL, "yy", "xxxx", "z");
    contract_and_check(IGRAPH_ATTRIBUTE_COMBINE_RANDOM, NULL, "y", "x", "z");
    contract_and_check(IGRAPH_ATTRIBUTE_COMBINE_FUNCTION,
                       (igraph_function_pointer_t) smallest_string, "y", "x", "z");

    VERIFY_FINALLY_STACK();

    return 0;
}
