/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "../unit/test_utilities.inc"

int main() {

    igraph_t graph;
    igraph_vector_t mod;

    igraph_empty(&graph, 25, IGRAPH_UNDIRECTED);
    igraph_vector_init(&mod, 0);
    igraph_community_multilevel(&graph, /*weights=*/ 0, /*resolution=*/ 1,
                                /*membership=*/ 0, /*memberships=*/ 0, &mod);

    if (igraph_vector_size(&mod) != 1 ||
        !igraph_is_nan(VECTOR(mod)[0])) {
        return 1;
    }

    igraph_vector_destroy(&mod);
    igraph_destroy(&graph);

    VERIFY_FINALLY_STACK();

    return 0;
}
