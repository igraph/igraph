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

void print_comb(igraph_attribute_combination_t *comb) {
    int i;
    igraph_attribute_combination_record_t *r;
    for (i = 0; i < igraph_vector_ptr_size(&comb->list); i++) {
        r = VECTOR(comb->list)[i];
        if (r->name) {
            printf("name: %s", r->name);
        } else {
            printf("name: NULL");
        }
        printf(", type: %d\n", r->type);
    }
    printf("\n");
}

int main() {
    igraph_attribute_combination_t comb;

    igraph_attribute_combination(&comb,
                                 "weight", IGRAPH_ATTRIBUTE_COMBINE_SUM,
                                 "type",   IGRAPH_ATTRIBUTE_COMBINE_FIRST,
                                 "",       IGRAPH_ATTRIBUTE_COMBINE_IGNORE,
                                 IGRAPH_NO_MORE_ATTRIBUTES);

    printf("starting combinations:\n");
    print_comb(&comb);

    printf("Removing nonexistent combination:\n");
    igraph_attribute_combination_remove(&comb, "nonexistent_name");
    print_comb(&comb);

    printf("Removing weight:\n");
    igraph_attribute_combination_remove(&comb, "weight");
    print_comb(&comb);

    printf("Removing type and NULL:\n");
    igraph_attribute_combination_remove(&comb, "type");
    igraph_attribute_combination_remove(&comb, NULL);
    print_comb(&comb);

    printf("Removing nonexistent combination again:\n");
    igraph_attribute_combination_remove(&comb, "nonexistent_name");
    igraph_attribute_combination_remove(&comb, NULL);
    print_comb(&comb);

    igraph_attribute_combination_destroy(&comb);

    VERIFY_FINALLY_STACK();
    return 0;
}
