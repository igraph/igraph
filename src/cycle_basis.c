/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   Copyright (C) 2019 the igraph team.
   
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

#include "igraph_structural.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
#include "igraph_interrupt_internal.h"
#include "igraph_memory.h"
#include "igraph_types_internal.h"
#include "igraph_visitor.h"

/**
 * \function igraph_fundamental_cycle_basis
 * \brief Finds a cycle basis for an unweighted, undirected graph
 * 
 * A cycle basis is a set of cycles such that the set difference
 * between them (i.e. the set of edges that are counted an odd
 * number of times by the cycles) is a Eulerian path across the
 * graph, which means it touches all vertices once and never
 * walks the same edge more than once.
 *
 * A fundamental cycle basis is found by calculating a spanning tree
 * and closing each pair of nonadjacent vertices with an external edge
 * if that edge is present in the graph.
 *
 * \param graph The input graph, it might be directed, but edge
 *    direction is ignored.
 * \param basis A pointer to a vector of vectors. Each fundamental cycle
 *    will be a member as a list of edge ids.
 * \return Error code.
 * 
 * Time complexity: O(|V| + |E|) (spanning tree) + ??.
 * 
 */
int igraph_fundamental_cycle_basis(const igraph_t *graph,
    igraph_vector_ptr_t *basis) {

    long int no_of_nodes=igraph_vcount(graph);
    long int no_of_edges=igraph_ecount(graph);
    long int no_of_cycles=no_of_edges - (no_of_nodes - 1);
    igraph_vector_char_t added_nodes, added_edges;
    igraph_vector_t parent_nodes, parent_edges, dist, inc_edges;
    igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;

    long int i, j, k;

    // Given a spanning tree (forest), fundamental cycles are given by
    // unions of:
    // - a path within the tree
    // - a single edge connecting end to start nodes outside the tree
    // The number of fundamental cycles is #edges - (#vertices -1),
    // because the quantity in parenthesis is the number of edges in a
    // tree (no cycles).
    
    // The key is that every non-tree edge generates a fundamental cycle
    // so we run a BFS and whenever we encounter a vertex which is not
    // new we are:
    // - either trying to walk back the last edge we just walked
    // - or using a non-tree edge (since both vertices are known already)
    // Note: every edge must be walked once, so even after the tree
    // is complete we still have to finish walking the remaining edges.
    // Note: we cannot use the igraph_bfs function because we have to
    // inspect outgoing edges that point to vertices that have already
    // been added

    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&added_nodes, no_of_nodes);
    IGRAPH_VECTOR_CHAR_INIT_FINALLY(&added_edges, no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&parent_nodes, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&parent_edges, no_of_edges);
    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    IGRAPH_VECTOR_INIT_FINALLY(&inc_edges, 0);

    for (i = 0; i < no_of_nodes; i++) {
        /* Start of a new connected component including the first one */
        if (added_nodes[i] > 0) { continue; }

        IGRAPH_ALLOW_INTERRUPTION();

	/* This is the root, mark it visited */
        added_nodes[i] = 1;

	/* The queue stores the next vertex, the edge to it, and its distance to root */
        IGRAPH_CHECK(igraph_dqueue_push(&q, i));
        IGRAPH_CHECK(igraph_dqueue_push(&q, -1));
        IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
        while (! igraph_dqueue_empty(&q)) {
          /* current node, incoming edge, and distance from tree root */
          long int actnode=(long int) igraph_dqueue_pop(&q);
	  long int actedge=(long int) igraph_dqueue_pop(&q);
          long int actdist=(long int) igraph_dqueue_pop(&q);

	  /* check all edges connected to this vertex */
          IGRAPH_CHECK(igraph_incident(graph, &inc_edges, (igraph_integer_t) actnode,
				       IGRAPH_ALL));

          for (j=0; j<igraph_vector_size(&inc_edges); j++) {
              long int edge=(long int) VECTOR(inc_edges)[j];

	      /* visited edges have been given their cycles already */
              if (added_edges[edge] == 0) {

		/* find the neighbour identity */
                igraph_integer_t from, to;
                igraph_edge(graph, (igraph_integer_t) edge, &from, &to);
		
		/* we only know that either one is actnode, swap if needed */
                if (to==actnode) {
		  to = from;
		  from = to;
		}

		/* if the edge goes to a new vertex, expand the spanning tree */
                if (added_nodes[(long int) to] == 0) {
                    added_nodes[(long int) to] = 1;
                    added_edges[edge] = 1;
		    parent_nodes[(long int) to] = actnode;
		    dist[(long int) to] = actdist + 1;
		    parent_edges[edge] = actedge;

		    /* and continue the search from there  */
                    IGRAPH_CHECK(igraph_dqueue_push(&q, to));
                    IGRAPH_CHECK(igraph_dqueue_push(&q, edge));
                    IGRAPH_CHECK(igraph_dqueue_push(&q, actdist + 1));
                }

               /* if it's a new edge to a known vertex, it's outside the tree
                * so it generates a cycle */
                else {
                    igraph_i_fundamental_cycle_basis_add(graph,
				    basis,
				    &parent_nodes, &parent_edges,
				    &added_edges, dist,
				    actedge, edge,
				    (long int) from, (long int) to);
                }
              }
          }
        }
    }
    
    igraph_vector_destroy(&inc_edges);
    igraph_dqueue_destroy(&q);
    igraph_vector_destroy(&dist);
    igraph_vector_destroy(&parent_edges);
    igraph_vector_destroy(&parent_nodes);
    igraph_vector_char_destroy(&added_edges);
    igraph_vector_char_destroy(&added_nodes);
    IGRAPH_FINALLY_CLEAN(7);

    return IGRAPH_SUCCESS;
}

/* Note: this function also works for self-edges */
int igraph_i_fundamental_cycle_basis_add(const igraph_t *graph,
	igraph_vector_ptr_t *basis,
	const igraph_vector_t *parent_nodes, const igraph_vector_t *parent_edges,
	const igraph_vector_chat_t *added_edges, const igraph_vector_t *dist,
	long int edge_from, long int edge_to,
	long int from, long int to) {

    /* 'from' is the last node visited. 'to' is the neighbor of 'from' that
     * was already in the tree. They are connected by 'edge_to'. We will need
     * to backtrack both 'from' and 'to' to their common parent in the tree
     * so we use two edge variables for that:
     * - 'edge_from' backtracks 'from' (since it starts with the correct parent)
     * - 'edge_to' backtracks 'to' (after it's added to the cycle)
     * */

    igraph_vector_t *cycle=igraph_Calloc(1, igraph_vector_t);
    igraph_vector_t to_backtrack, to_edges;
    long int n, i;
    
    IGRAPH_VECTOR_INIT_FINALLY(&to_backtrack, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&to_edges, 0);

    /* start by adding the non-tree edge */
    igraph_vector_init(cycle, 1);
    cycle[0] = edge_to;

    /* find the edge connecting 'to' to its parent in the tree */
    IGRAPH_CHECK(igraph_incident(graph, &to_edges, (igraph_integer_t) to,
	IGRAPH_ALL));
    for (i = 0; i < igraph_vector_size(&to_edges); i++) {
        igraph_integer_t to_from, to_to;

        edge_to = (long int) VECTOR(to_edges)[i];

        /* only a visited edge can be in the tree */
        if (added_edges[edge]==0) {
	    continue
        }

        /* find the vertices of this edge */
        IGRAPH_CHECK(igraph_edge(graph,
				 (igraph_integer_t) edge_to,
				 &to_from, &to_to);)

	/* we only know that either one is to, swap if needed */
        if (to_to == to) {
	  to_to = to_from;
	  to_from = to_to;
	}

	if (to_to == parent_nodes[to]) {
	    break;
	}
    }
    igraph_vector_destroy(&to_edges);
    IGRAPH_FINALLY_CLEAN(1);

    /* construct the cycle by backtracing the farthest one
     * first, and then both together */
    while (from != to) {
        if (dist[from] > dist[to]) {
	    /* backtrack 'from' */
            IGRAPH_CHECK(igraph_vector_push_back(cycle, edge_from));
	    from = parent_nodes[from];
	    edge_from = parent_edges[edge_from];
	}
	else if (dist[to] > dist[from]) {
	    /* backtrack 'to' */
            IGRAPH_CHECK(igraph_vector_push_back(to_backtrack, edge_to));
	    to = parent_nodes[to];
	    edge_to = parent_edges[edge_to];
	}
	else {
	    /* same distance but different, backtrack both until they merge */
            IGRAPH_CHECK(igraph_vector_push_back(cycle, edge_from));
            IGRAPH_CHECK(igraph_vector_push_back(to_backtrack, edge_to));
	    from = parent_nodes[from];
	    to = parent_nodes[to];
	    edge_from = parent_edges[edge_from];
	    edge_to = parent_edges[edge_to];
	}
    }

    /* add the backtracked edges from 'to' in reverse order*/
    n = igraph_vector_size(to_backtrack);
    for(i=0; i<n; i++) {
        IGRAPH_CHECK(igraph_vector_push_back(cycle, to_backtrack[n - 1 - i]));
    }
    igraph_vector_destroy(&to_backtrack);
    IGRAPH_FINALLY_CLEAN(1);

    /* add the cycle to the basis */
    IGRAPH_CHECK(igraph_vector_ptr_push_back(basis, cycle));	

    return IGRAPH_SUCCESS;
}


int igraph_minimum_cycle_basis(const igraph_t *graph,
                igraph_vector_ptr_t *basis,
                const igraph_vector_t *weights) {




    /* TODO: implement */
    return IGRAPH_SUCCESS;
}


