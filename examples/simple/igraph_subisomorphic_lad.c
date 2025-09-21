/*
   igraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

void print_maps(igraph_vector_int_t *map, igraph_vector_int_list_t *maps) {
    igraph_int_t n, i;
    igraph_vector_int_print(map);
    n = igraph_vector_int_list_size(maps);
    for (i = 0; i < n; i++) {
        igraph_vector_int_t *v = igraph_vector_int_list_get_ptr(maps, i);
        igraph_vector_int_print(v);
    }
    igraph_vector_int_list_clear(maps);
}

int main(void) {
    igraph_t target, pattern;
    igraph_bool_t iso;
    igraph_vector_int_t map;
    igraph_vector_int_list_t maps;
    igraph_int_t i;
    int domainsvec[] = { 0, 2, 8, -1,
                         4, 5, 6, 7, -1,
                         1, 3, 5, 6, 7, 8, -1,
                         0, 2, 8, -1,
                         1, 3, 7, 8, -1, -2
                       };
    igraph_vector_int_list_t domains;
    igraph_vector_int_t v;

    /* Initialize the library. */
    igraph_setup();

    igraph_small(&target, 9, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4, 0, 6,
                 1, 4, 1, 2,
                 2, 3,
                 3, 4, 3, 5, 3, 7, 3, 8,
                 4, 5, 4, 6,
                 5, 6, 5, 8,
                 7, 8,
                 -1);

    igraph_small(&pattern, 5, IGRAPH_UNDIRECTED,
                 0, 1, 0, 4,
                 1, 4, 1, 2,
                 2, 3,
                 3, 4,
                 -1);

    igraph_vector_int_init(&map, 0);
    igraph_vector_int_list_init(&maps, 0);

    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ NULL, &iso, &map,
                             &maps, /*induced=*/ false);

    if (!iso) {
        return 1;
    }
    print_maps(&map, &maps);

    printf("---------\n");

    igraph_subisomorphic_lad(&pattern, &target, /*domains=*/ NULL, &iso, &map,
                             &maps, /*induced=*/ true);

    if (!iso) {
        return 2;
    }
    print_maps(&map, &maps);

    printf("---------\n");

    igraph_vector_int_list_init(&domains, 0);
    i = 0;
    igraph_vector_int_init(&v, 0);
    while (1) {
        if (domainsvec[i] == -2) {
            break;
        } else if (domainsvec[i] == -1) {
            igraph_vector_int_list_push_back_copy(&domains, &v);
            igraph_vector_int_clear(&v);
        } else {
            igraph_vector_int_push_back(&v, domainsvec[i]);
        }
        i++;
    }
    igraph_vector_int_destroy(&v);

    igraph_subisomorphic_lad(&pattern, &target, &domains, &iso, &map, &maps,
                             /*induced=*/ false);

    if (!iso) {
        return 3;
    }
    print_maps(&map, &maps);

    igraph_vector_int_list_destroy(&domains);
    igraph_vector_int_destroy(&map);
    igraph_vector_int_list_destroy(&maps);

    igraph_destroy(&pattern);
    igraph_destroy(&target);


    return 0;
}
