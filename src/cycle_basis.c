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
     * 2. Construct candidate set of cycles, sorted in nondecreasing
     *    order of weight. Each cycle is attached to a tree (but more than
     *    one cycle can be attached to the same tree).
     * 3. for i = 1 to N do
     *     3.1 compute vector Si s.t. <Cj, Si> = 0 for each 1 <= j < i
     *         (orthogonal to the current cycles)
     *     3.2 for all trees, update vertex labels based on Si
     *     3.3 for each candidate cycle, compute <C, Si> using the labels
     *     3.4 append the least-weight cycle non-orthogonal to Si
     *         (so we can normalize it into <C, Si> = 1
     *
     * */
    /* TODO: finish implementation */

    long int i, n;
    long int no_of_nodes = igraph_vcount(graph);
    igraph_vector_t si;
    igraph_vector_int_t nonorth_cycles;

    igraph_vector_t fvs, basis_weight;
    /* list of vertices = tree -> list of trees */
    igraph_vector_ptr_t trees;
    /* NOTE: [tree_id, edge_id, tree_id2, edge_id2, ...] */
    igraph_vector_int_t candidate_cycles;
    /* total weight of each cycle */
    igraph_vector_t candidate_weights;

    /* 0. Feedback vertex set */
    IGRAPH_CHECK(igraph_i_feedback_vertex_set(
			    graph, &fvs, weights));
    n = igraph_vector_size(&fvs);

    /* 1. Construct all shortest path trees */
    IGRAPH_CHECK(igraph_i_shortest_path_trees(
			    graph, weights, &fvs, &trees));

    /* 2. Construct candidate list */
    IGRAPH_CHECK(igraph_i_candidate_cycles(
			    graph, weights,
			    &candidate_cycles, &candidate_weights));

    /* 3. Construct basis one cycle at a time
     *    similar to a standard Gram Schmidt orthogonalization */
    for (i = 0; i < no_of_nodes; i++) {
        /* 3.1 compute Si, the nonnormalized, minimal vector
	 *     sticking out from the current space */
	IGRAPH_CHECK(igraph_i_compute_Si(graph, weights, basis, i, &si));

	/* 3.2 update tree vertex labels from Si */
	IGRAPH_CHECK(igraph_i_update_trees_vertex_labels(
			&trees, &si));

	/* 3.3 find non-orthogonal candidates */
	IGRAPH_CHECK(igraph_i_nonorthogonal_candidates(
			candidate_cycles, &si, &nonorth_cycles));

	/* 3.4 get the shortest nonorthogonal cycle */
	IGRAPH_CHECK(igraph_i_shortest_nonorthogonal_cycle(
			candidate_cycles, &candidate_weights, &nonorth_cycles, basis));

    }


    /* Clean */
    igraph_vector_destroy(&fvs);
    igraph_vector_ptr_destroy_all(trees);
    IGRAPH_FINALLY_CLEAN(1 + n + 1); // FIXME: Check this number

    return IGRAPH_SUCCESS;
}


/* NOTE: Adaptated from igraph_i_minimum_spanning_tree_prim */
int igraph_i_shortest_path_tree_rooted(const igraph_t *graph,
		const igraph_vector_t *weights,
		const igraph_int_t root,
		igraph_vector_t *res,
		) {

    long int no_of_nodes = igraph_vcount(graph);
    long int no_of_edges = igraph_ecount(graph);
    char *already_added;
    char *added_edges;

    igraph_d_indheap_t heap = IGRAPH_D_INDHEAP_NULL;
    igraph_integer_t mode = IGRAPH_ALL;

    igraph_vector_t adj;
    igraph_vector_int_t node_order;

    long int i, j, ir;

    igraph_vector_clear(res);

    IGRAPH_FINALLY(igraph_free, already_added);
    IGRAPH_CHECK(igraph_d_indheap_init(&heap, 0));
    IGRAPH_FINALLY(igraph_d_indheap_destroy, &heap);
    IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);

    /* use two pointers: ir goes from 0 to n-1, i starts from the root and follows
     * e.g. if root = 2
     *
     * ir   i
     * 0 -> 2
     * 1 -> 0
     * 2 -> 1
     * 3 -> 3
     * 4 -> 4
     *
     * et cetera
     * */
    for (ir = 0; ir < no_of_nodes; ir++) {
        if (ir == 0) {
	    i = root;
	} else if (ir <= root) {
	    i = ir - 1;
	} else {
	    i = ir;
	}

        if (already_added[i] > 0) {
            continue;
        }
        IGRAPH_ALLOW_INTERRUPTION();

        already_added[i] = 1;
        /* add all edges of the first vertex to the heap */
        igraph_incident(graph, &adj, (igraph_integer_t) i, (igraph_neimode_t) mode);
        for (j = 0; j < igraph_vector_size(&adj); j++) {
            long int edgeno = (long int) VECTOR(adj)[j];
            igraph_integer_t edgefrom, edgeto;
            long int neighbor;
            igraph_edge(graph, (igraph_integer_t) edgeno, &edgefrom, &edgeto);
            neighbor = edgefrom != i ? edgefrom : edgeto;
            if (already_added[neighbor] == 0) {
                IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], i,
                                                   edgeno));
            }
        }

        while (! igraph_d_indheap_empty(&heap)) {
            /* Get minimal edge among the heaped ones */
            long int from, edge;
            igraph_integer_t tmp, to;
            igraph_d_indheap_max_index(&heap, &from, &edge);
            igraph_edge(graph, (igraph_integer_t) edge, &tmp, &to);

            /* Erase it */
            igraph_d_indheap_delete_max(&heap);

            /* Is this edge already included? */
            if (added_edges[edge] == 0) {
                if (from == to) {
                    to = tmp;
                }
                /* Does it point to a visited node? */
                if (already_added[(long int)to] == 0) {
                    already_added[(long int)to] = 1;
                    added_edges[edge] = 1;
                    IGRAPH_CHECK(igraph_vector_push_back(res, edge));
                    /* add all outgoing edges */
                    igraph_incident(graph, &adj, to, (igraph_neimode_t) mode);
                    for (j = 0; j < igraph_vector_size(&adj); j++) {
                        long int edgeno = (long int) VECTOR(adj)[j];
                        igraph_integer_t edgefrom, edgeto;
                        long int neighbor;
                        igraph_edge(graph, (igraph_integer_t) edgeno, &edgefrom, &edgeto);
                        neighbor = edgefrom != to ? edgefrom : edgeto;
                        if (already_added[neighbor] == 0) {
                            IGRAPH_CHECK(igraph_d_indheap_push(&heap, -VECTOR(*weights)[edgeno], to,
                                                               edgeno));
                        }
                    }
                } /* for */
            } /* if !already_added */
        } /* while in the same component */
    } /* for all nodes */

    igraph_d_indheap_destroy(&heap);
    igraph_Free(already_added);
    igraph_vector_destroy(&adj);
    igraph_Free(added_edges);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;

}


int igraph_i_shortest_path_trees(const igraph_t *graph,
		const igraph_vector_t *weights,
		const igraph_vector_t *fvs,
		igraph_vector_ptr_t *trees,
		) {

    long int i;
    long int n = igraph_vector_size(fvs);

    igraph_vector_ptr_init(trees, n);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, trees);

    for(i = 0; i < n; i++) {
	/* each tree is a list of edges */
        IGRAPH_CHECK(igraph_i_shortest_path_tree_rooted(graph,
				weights,
				VECTOR(fvs)[i],
				VECTOR(trees)[i]));
    }

    return IGRAPH_SUCCESS;
}

int igraph_i_candidate_cycles(const igraph_t *graph,
		const igraph_vector_t *weights,
		const igraph_vector_ptr_t *trees,
		igraph_vector_int_t *candidate_cycles,
		igraph_vector_t *cadidate_weights,
		){

    /* These can be Horton cycles or a subset, which is more efficient
     * Let's see how far we get with the implementation */

    long int i, ie;
    igraph_vector_int_t subtree_label, cycle_order;
    long int no_of_nodes=igraph_vcount(graph);
    long int no_of_edges=igraph_ecount(graph);
    igraph_vector_es_t edges;
    igraph_es_all(&es, IGRAPH_EDGEORDER_ID);
    
    IGRAPH_CHECK(igraph_vector_int_init(subtree_label, no_of_nodes));

    /* For each tree... */
    for(i = 0; i < igraph_vector_size(trees); i++) {
        /* find edges for which the two vertices only meet in the tree at the root */
	/* how do we find them? well we
	 * 1. get the children of root and create subtrees
	 * 2. we test if any edge has vertices in different subtrees
	 * whether we approach by vertices or edges depends on how sparse the
	 * graph is. For now let's assign every vertex to a subtree, then go over
	 * the edges and check */
         igraph_inclist_t inclist;
         igraph_vector_int_t *root_edges;
         long int j, nchild;
         igraph_vector_t children;

	 igraph_vector_t tree = VECTOR(trees)[i];
	 igraph_integer_t root = VECTOR(tree)[0];

	 /* FIXME: each tree is not a graph but a list of edges, fix this */
         IGRAPH_CHECK(igraph_inclist_init(tree, &inclist, IGRAPH_ALL));
         root_edges = igraph_inclist_get(il, root);
	 nchild = igraph_vector_int_size(edges);
         IGRAPH_CHECK(igraph_vector_init(children, nchild));

	 /* root is not in any subtree, so put it as -1 */
	 VECTOR(subtree_label)[0] = -1;
	 for(j = 0; j < nchild; j++) {
	     VECTOR(children)[j] = IGRAPH_OTHER(tree, VECTOR(edges)[j], root);
	     /* label each child into its own subtree */
	     VECTOR(subtree_label)[VECTOR(children)[j]] = j
	 }

	 /* TODO: label the rest of the tree recursively based on the children */

	 /* Look through the edges for the ones spanning different subtrees */
	 for(j = 0; j < no_of_edges; j++) {
             double weight = 0;
             igraph_integer_t from, igraph_integer_t to;
	     igraph_integer_t edge = VECTOR(edges)[j];
	     igraph_edge(graph, edge, &from, &to);

	     from = VECTOR(subtree_label)[from];
	     to = VECTOR(subtree_label)[to];

	     if ((from != -1) && (to != -1) && (from != to)) {
	         /* since we store the trees anyway, we could make it such as the vector
	          * of candidates contains
	          *
	          * [tree_id, key_edge_id, tree_id2, key_edge_id2]
	          *
	          * and separately the cycle total weight (since it's a float). The only
	          * problem is that we need to sort by increasing weight, so we'll do
	          * argsort and carefully resort this vector too.
	          * */
                 IGRAPH_CHECK(igraph_vector_push_back(candidate_cycles, i));
                 IGRAPH_CHECK(igraph_vector_push_back(candidate_cycles, edge));

		 /* TODO: compute the total weight of this cycle */
                 IGRAPH_CHECK(igraph_vector_push_back(candidate_weights, weight));
	     }
	 }
    
    }

    /* sort candidate list by weight */
    /* TODO implement these functions? or better use structs? */
    igraph_i_argsort(candidate_weights, &cycle_order);
    igraph_i_sort_from_argsort(candidate_cycles, candidate_weights, &cycle_order);

    return IGRAPH_SUCCESS;
}

int igraph_i_compute_Si(const igraph_t *graph,
		const igraph_vector_t *weights,
		const igraph_vector_ptr_t *basis,
		long int i,
		igraph_vector_t *si,
		) {

    /* TODO: implement */
    return IGRAPH_SUCCESS;

}


int igraph_i_update_trees_vertex_labels(igraph_vector_ptr_t *trees,
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

