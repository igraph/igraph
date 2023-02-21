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

int main(void) {

    igraph_t g;
    igraph_es_t es;

    igraph_small(&g, 4, IGRAPH_UNDIRECTED, 0,1, 1,2, 2,2, 2,3, -1);

    igraph_es_pairs_small(&es, IGRAPH_DIRECTED, 3, 2, -1);
    igraph_delete_edges(&g, es);
    if (igraph_ecount(&g) != 3) {
        return 1;
    }

    igraph_es_destroy(&es);
    igraph_destroy(&g);

    return 0;
}
