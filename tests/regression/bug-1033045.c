/*
   igraph library.
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA
*/

#include <igraph.h>

#include "../unit/test_utilities.h"

int main(void) {
    igraph_t graph;
    igraph_vector_int_list_t separators;

    igraph_small(&graph, 0, /*directed=*/ 0,
                 0, 1, 0, 2, 1, 3, 1, 4, 2, 3, 2, 5, 3, 4, 3, 5, 4, 6, 5, 6, -1);
    igraph_vector_int_list_init(&separators, 0);

    igraph_all_minimal_st_separators(&graph, &separators);

    print_vector_int_list(&separators);
    igraph_vector_int_list_destroy(&separators);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
