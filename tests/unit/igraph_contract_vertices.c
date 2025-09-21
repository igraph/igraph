/* igraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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
#include "test_utilities.h"

int main(void) {
    igraph_t g;
    igraph_vector_int_t mapping;
    igraph_attribute_combination_t comb;

    igraph_set_attribute_table(&igraph_cattribute_table);

    igraph_attribute_combination(&comb,
                                 "an", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 "ab", IGRAPH_ATTRIBUTE_COMBINE_FIRST,
                                 "as", IGRAPH_ATTRIBUTE_COMBINE_LAST,
                                 "ean", IGRAPH_ATTRIBUTE_COMBINE_PROD,
                                 IGRAPH_NO_MORE_ATTRIBUTES);

    printf("Graph with no vertices:\n");
    {
        igraph_small(&g, 0, IGRAPH_DIRECTED, -1);
        mapping = igraph_vector_int_view(NULL, 0);
        igraph_contract_vertices(&g, &mapping, &comb);
        igraph_destroy(&g);
    }

    printf("Graph with loops and multiple edges:\n");
    {
        igraph_small(&g, 6, IGRAPH_DIRECTED, 0,1, 0,2, 1,1, 1,3, 2,0, 2,3, 3,4, 3,4, -1);
        igraph_bool_t booleans[6] = {0, 1, 0, 1, 0, 1};
        igraph_real_t real_numbers[6] = {0, 1, 2, 3, 4, 5};
        igraph_real_t real_numbers_edges[8] = {0, 1, 2, 3, 4, 5, 6, 7};
        igraph_vector_t reals;
        igraph_vector_t reals_edges;
        igraph_vector_bool_t bools;

        reals = igraph_vector_view(real_numbers, 6);
        reals_edges = igraph_vector_view(real_numbers_edges, 8);
        bools = igraph_vector_bool_view(booleans, 6);

        SETVABV(&g, "ab", &bools);
        SETVANV(&g, "an", &reals);
        SETVAS(&g, "as", 0, "zero");
        SETVAS(&g, "as", 1, "one");
        SETVAS(&g, "as", 2, "two");
        SETVAS(&g, "as", 3, "three");
        SETVAS(&g, "as", 4, "four");
        SETVAS(&g, "as", 5, "five");
        SETEANV(&g, "ean", &reals_edges);

        {
            printf("No contraction:\n");
            igraph_int_t mappings[6] = {0, 1, 2, 3, 4, 5};
            mapping = igraph_vector_int_view(mappings, 6);
            igraph_contract_vertices(&g, &mapping, &comb);
            print_attributes(&g);
        }
        {
            printf("Contract 5 into 4:\n");
            igraph_int_t mappings[6] = {0, 1, 2, 3, 4, 4};
            mapping = igraph_vector_int_view(mappings, 6);
            igraph_contract_vertices(&g, &mapping, &comb);
            print_attributes(&g);
        }
        {
            printf("Contract 4 into 3:\n");
            igraph_int_t mappings[5] = {0, 1, 2, 3, 3};
            mapping = igraph_vector_int_view(mappings, 5);
            igraph_contract_vertices(&g, &mapping, &comb);
            print_attributes(&g);
        }
        {
            printf("Contract 3 into 2:\n");
            igraph_int_t mappings[4] = {0, 1, 2, 2};
            mapping = igraph_vector_int_view(mappings, 4);
            igraph_contract_vertices(&g, &mapping, &comb);
            print_attributes(&g);
        }
        {
            printf("Contract 2 into 0:\n");
            igraph_int_t mappings[3] = {0, 1, 0};
            mapping = igraph_vector_int_view(mappings, 3);
            igraph_contract_vertices(&g, &mapping, &comb);
            print_attributes(&g);
        }
    }

    VERIFY_FINALLY_STACK();

    {
        printf("Check incorrect mapping size error.\n");
        igraph_int_t mappings[3] = {0, 1, 0};
        mapping = igraph_vector_int_view(mappings, 3);
        CHECK_ERROR(igraph_contract_vertices(&g, &mapping, &comb), IGRAPH_EINVAL);
    }

    igraph_attribute_combination_destroy(&comb);
    igraph_destroy(&g);
    VERIFY_FINALLY_STACK();
    return 0;
}
