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
 * \brief Find a fundamental cycle basis for an unweighted, undirected graph
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
 * Time complexity: O(|V| + |E|) (spanning tree) + ?? (backtracking).
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


/**
 * \function igraph_minimum_weight_cycle_basis
 * \brief Find a minimum weight cycle basis for a weighted, undirected graph
 * 
 * A cycle basis is a set of cycles such that the set difference
 * between them (i.e. the set of edges that are counted an odd
 * number of times by the cycles) is a Eulerian path across the
 * graph, which means it touches all vertices once and never
 * walks the same edge more than once.
 *
 * A minimum weight cycle basis is the equivalent for weighted graphs
 * of a fundamental cycle basis for unweighted graphs.
 *
 * \param graph The input graph, it might be directed, but edge
 *    direction is ignored.
 * \param basis A pointer to a vector of vectors. Each fundamental cycle
 *    will be a member as a list of edge ids.
 * \param weights A const pointer to a vector of weights.
 * \return Error code.
 * 
 * Time complexity: O(|V| |E|) (Horton cycles) ??.
 *
 * A review:
 * https://www.sciencedirect.com/science/article/pii/S1574013709000483?casa_token=B7FmnsCBElwAAAAA:bfOYZDhzHpP24ym_2S2n5VuEU-dYsdAQJdmH3nvk8efkFrVmJkEY-KhHXawLZfoMFAKPIdHE
 *
 * The algorithm implemented here:
 * https://dl.acm.org/doi/abs/10.1145/1644015.1644023
 * 
 */


int igraph_minimum_cycle_basis(const igraph_t *graph,
                const igraph_vector_t *weights,
                igraph_vector_ptr_t *basis,
		) {

    /* General plan:
     *
     * 0. Find a feedback vertex set Z: all cycles pass through a vertex
     *    in Z. Worst case Z := {all vertices}. Basically, 2-degree nodes
     *    are excluded because there's no decision: the path goes in one
     *    way and out another. Of course, 0-degree are also not needed.
     * 1. Construct all shortest path trees 
     * 2. Construct candidate set of cycle bases, sorted in nondecreasing
     *    order of weight
     * 3. for i = 1 to N do
     *     3.1 compute vector Si s.t. <Cj, Si> = 0 for each 1 <= j < i
     *         (orthogonal to the current cycles)
     *     3.2 for all trees, update vertex labels based on Si
     *     3.3 for each candidate cycle, compute <C, Si> using the labels
     *     3.4 append the least-weight cycle non-orthogonal to Si
     *         (so we can normalize it into <C, Si> = 1
     *
     * */
    /* TODO: implement the actual functions */

    long int i;
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t si;
    igraph_vector_int_t nonorth_cycles;

    igraph_vector_t fvs, basis_weight;
    /* list of vertices = tree -> list of trees */
    igraph_vector_ptr_t trees_z;
    /* list of vertices = vector -> list of vectors = basis -> list of bases */
    igraph_vector_ptr_t candidate_cycles;
    /* total weight of each basis (?) */
    igraph_vector_t candidate_weights;

    /* 0. Feedback vertex set */
    IGRAPH_CHECK(igraph_i_feedback_vertex_set_approx(
			    graph, weights, &fvs));

    /* 1. Construct all shortest path trees */
    IGRAPH_CHECK(igraph_i_shortest_path_trees(
			    graph, weights, &fvs, &trees_z));

    /* 2. Construct candidate list */
    IGRAPH_CHECK(igraph_i_candidate_bases(
			    graph, weights,
			    &candidate_cycles, &candidate_weights));

    /* 3. Construct basis one by one */
    for (i = 0; i < no_of_nodes; i++) {
        /* 3.1 compute Si */
	igraph_i_compute_Si(graph, weights, basis, i, &si);

	/* 3.2 update tree vertex labels from Si */
	igraph_i_update_tree_z_vertex_labels(
			&trees_z, &si);

	/* 3.3 find non-orthogonal candidates */
	igraph_i_nonorthogonal_candidates(
			candidate_cycles, &si, &nonorth_cycles);

	/* 3.4 get the shortest nonorthogonal cycle */
	igraph_i_shortest_nonorthogonal_cycle(
			candidate_cycles, &candidate_weights, &nonorth_cycles, basis);

    }

    return IGRAPH_SUCCESS;
}


int igraph_i_feedback_vertex_set_approx(const igraph_t *graph,
		const igraph_vector_t *weights,
		igraph_vector_t *vertices,
		) {

    long int i;
    long int no_of_nodes = igraph_vcount(graph);

    /* TODO: implement */
    /* For now, return all vertices */
    igraph_vector_clear(vertices);
    IGRAPH_CHECK(igraph_vector_reserve(vertices, no_of_nodes));

    igraph_vs_t vertices_all = igraph_vss_all();
    for (i = 0; i < no_of_nodes; i++) {
	    vertices[i] = vertices_all[i];
    }

    return IGRAPH_SUCCESS;
}


int igraph_i_shortest_path_trees(const igraph_t *graph,
		const igraph_vector_t *weights,
		const igraph_vector_t *vertices,
		igraph_vector_ptr_t *trees_z,
		) {

    /* TODO: implement */
    return IGRAPH_SUCCESS;
}

int igraph_i_candidate_bases(const igraph_t *graph,
		const igraph_vector_t *weights,
		igraph_vector_ptr_t *candidate_cycles,
		igraph_vector_t *cadidate_weights,
		){

    /* TODO: implement */
    return IGRAPH_SUCCESS;
}

int igraph_i_compte_Si(const igraph_t *graph,
		const igraph_vector_t *weights,
		const igraph_vector_ptr_t *basis,
		long int i,
		igraph_vector_t *si,
		) {

    /* TODO: implement */
    return IGRAPH_SUCCESS;

}


int igraph_i_update_tree_z_vertex_labels(igraph_vector_ptr_t *trees_z,
		const igraph_vector_t *si,
		) {

    /* TODO: implement */
    return IGRAPH_SUCCESS;
}


int igraph_i_nonorthogonal_candidates(const igraph_vector_ptr_t *candidate_cycles,
		const igraph_vector_t *si,
		igraph_vector_int_t *nonorth_cycles,
		) {

    /* TODO: implement */
    return IGRAPH_SUCCESS;
}


int igraph_i_shortest_nonorthogonal_cycle(const igraph_vector_ptr_t *candidate_cycles,
		const igraph_vector_t *cadidate_weights,
		const igraph_vector_int_t *nonorth_cycles,
		igraph_vector_ptr_t *basis,
		) {

    /* TODO: implement */
    return IGRAPH_SUCCESS;
}

