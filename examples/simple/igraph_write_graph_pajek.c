/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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

    igraph_t g;
    igraph_strvector_t names;

    igraph_set_attribute_table(&igraph_cattribute_table);

    /* save a simple ring graph */
    igraph_ring(&g, 10, IGRAPH_DIRECTED, 0 /* mutual */, 1 /* circular */);
    igraph_write_graph_pajek(&g, stdout);

    /* add some vertex attributes */
    igraph_strvector_init(&names, 0);
    igraph_strvector_add(&names, "A");
    igraph_strvector_add(&names, "B");
    igraph_strvector_add(&names, "C");
    igraph_strvector_add(&names, "D");
    igraph_strvector_add(&names, "E");
    igraph_strvector_add(&names, "F");
    igraph_strvector_add(&names, "G");
    igraph_strvector_add(&names, "H");
    igraph_strvector_add(&names, "I");
    igraph_strvector_add(&names, "J");
    SETVASV(&g, "id", &names);
    igraph_strvector_destroy(&names);

    /* save the graph with vertex names */
    igraph_write_graph_pajek(&g, stdout);

    igraph_strvector_init(&names, 0);
    igraph_strvector_add(&names, "square");
    igraph_strvector_add(&names, "square");
    igraph_strvector_add(&names, "square");
    igraph_strvector_add(&names, "square");
    igraph_strvector_add(&names, "escaping spaces");
    igraph_strvector_add(&names, "square");
    igraph_strvector_add(&names, "square");
    igraph_strvector_add(&names, "escaping \\backslashes\\");
    igraph_strvector_add(&names, "square");
    igraph_strvector_add(&names, "escaping \"quotes\"");
    SETVASV(&g, "shape", &names);
    igraph_strvector_destroy(&names);

    /* save the graph with escaped shapes */
    igraph_write_graph_pajek(&g, stdout);

    /* destroy the graph */
    igraph_destroy(&g);
    return 0;
}
