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

int igraph_rewire_core(igraph_t *graph, igraph_integer_t n, igraph_rewiring_t mode, igraph_bool_t use_adjlist);

static void check_rewiring(igraph_tree_mode_t tree_mode, igraph_bool_t use_adjlist, const char* description) {

    igraph_t g;
    igraph_vector_t indegree_before, outdegree_before, indegree_after, outdegree_after;

    igraph_tree(&g, 10, 3, tree_mode);

    igraph_vector_init(&indegree_before, 0);
    igraph_vector_init(&outdegree_before, 0);
    igraph_degree(&g, &indegree_before, igraph_vss_all(), IGRAPH_IN, 0);
    igraph_degree(&g, &outdegree_before, igraph_vss_all(), IGRAPH_OUT, 0);

    igraph_rewire_core(&g, 1000, IGRAPH_REWIRING_SIMPLE, use_adjlist);

    igraph_vector_init(&indegree_after, 0);
    igraph_vector_init(&outdegree_after, 0);
    igraph_degree(&g, &indegree_after, igraph_vss_all(), IGRAPH_IN, 0);
    igraph_degree(&g, &outdegree_after, igraph_vss_all(), IGRAPH_OUT, 0);

    if ((!igraph_vector_all_e(&indegree_before, &indegree_after)) ||
        (!igraph_vector_all_e(&outdegree_before, &outdegree_after))) {

        fprintf(stderr, "%s graph degrees changed\n", description);
        exit(1);

    }

    igraph_destroy(&g);
    igraph_vector_destroy(&indegree_before);
    igraph_vector_destroy(&outdegree_before);
    igraph_vector_destroy(&indegree_after);
    igraph_vector_destroy(&outdegree_after);

}

int main() {

    check_rewiring(IGRAPH_TREE_OUT, 0, "Directed, standard-method");
    check_rewiring(IGRAPH_TREE_OUT, 1, "Directed, adjlist-method");
    check_rewiring(IGRAPH_TREE_UNDIRECTED, 0, "Undirected, standard-method");
    check_rewiring(IGRAPH_TREE_UNDIRECTED, 1, "Undirected, adjlist-method");

    return 0;

}
