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
 * \function igraph_cycle_basis_unweighted_undirected
 * Finds a cycle basis for an unweighted, undirected graph
 * 
 * A cycle basis is a set of cycles such that the set difference
 * between them (i.e. the set of edges that are counted an odd
 * number of times by the cycles) is a Eulerian path across the
 * graph, which means it touches all vertices once and never
 * walks the same edge more than once.
 *
 * For unweighted, undirected graphs a cycle basis is always found
 * by calculating a spanning tree and closing it with edges (??)
 *
 * \param graph The input graph, it might be directed, but edge
 *    direction is ignored.
 * \param cycles An initialized vector of graphs. Each cycle will be
 *    a member of this vector (as a graph).
 * \return Error code.
 * 
 * Time complexity: ???.
 * 
 * \sa \ref igraph_maximum_cardinality_search().
 */
int igraph_cycle_basis_unweighted_undirected(const igraph_t* graph,
    igraph_vector_t* cycles) {

    // FIXME: deal with unconnected graphs
    igraph_t mst;
    igraph_vector_t edges_tree=IGRAPH_VECTOR_NULL;
    igraph_vector_t edges_outside=IGRAPH_VECTOR_NULL;

    long int no_of_nodes=igraph_vcount(graph);
    long int no_of_edges=igraph_ecount(graph);
    long int no_of_cycles=no_of_edges - (no_of_nodes - 1);
    char *added_nodes;
    char *added_edges;
    
    igraph_dqueue_t q=IGRAPH_DQUEUE_NULL;
    igraph_vector_t tmp=IGRAPH_VECTOR_NULL;
    long int i, j, k;
    igraph_vector_t parent;

    igraph_vector_clear(cycles);

    // If the graph is a tree, there are no cycles
    if (no_of_edges == no_of_nodes - 1)
	return IGRAPH_SUCCESS;

    // Given a spanning tree, fundamental cycles are given by unions of:
    // - a path within the tree
    // - a single edge connecting end to start nodes outside the tree
    // The number of fundamental cycles is #edges - (#vertices -1),
    // because the quantity in parenthesis is the number of edges in a
    // tree (no cycles). So the way to go is to:
    // 1. find a tree and mark those edges
    // 2. mark edges outside the tree. Each will form a fundamental cycle
    // 3. for each edge, find the unique tree path that runs parallel to it
    // 4. form that cycle by union of the edge and the parallel path

    // 1. find a tree and mark those edges
    // this is a simple BFS
    IGRAPH_VECTOR_INIT_FINALLY(&edges_tree, no_of_nodes-1);
    added_edges=igraph_Calloc(no_of_edges, char);
    if (added_edges==0) {
      IGRAPH_ERROR("unweighted spanning tree failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added_edges);
    added_nodes=igraph_Calloc(no_of_nodes, char);
    if (added_nodes==0) {
      IGRAPH_ERROR("unweighted spanning tree failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_DQUEUE_INIT_FINALLY(&q, 100);
    
    for (i=0; i<no_of_nodes; i++) {
      if (added_nodes[i]>0) { continue; }

      // We are here once for every connected component.
      // Within that component, the following deque performs the BFS
      // that traverses it fully

      IGRAPH_ALLOW_INTERRUPTION();

      added_nodes[i]=1;
      IGRAPH_CHECK(igraph_dqueue_push(&q, i));
      while (! igraph_dqueue_empty(&q)) {
        long int act_node=(long int) igraph_dqueue_pop(&q);
        IGRAPH_CHECK(igraph_incident(graph, &tmp, (igraph_integer_t) act_node,
          			   IGRAPH_ALL));
        for (j=0; j<igraph_vector_size(&tmp); j++) {
          long int edge=(long int) VECTOR(tmp)[j];
          if (added_edges[edge]==0) {
            igraph_integer_t from, to;
            igraph_edge(graph, (igraph_integer_t) edge, &from, &to);
            if (act_node==to) { to=from; }
            if (added_nodes[(long int) to]==0) {
              added_nodes[(long int) to]=1;
              added_edges[edge]=1;
              IGRAPH_CHECK(igraph_vector_push_back(&edges_tree, edge));
              IGRAPH_CHECK(igraph_dqueue_push(&q, to));
            }
          }
        }
      }
    }
    
    igraph_Free(added_nodes);
    igraph_Free(added_edges);
    igraph_dqueue_destroy(&q);
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(4);

    // actually make the mst
    // FIXME: do the vids in the mst differ from the same in the original graph??
    IGRAPH_CHECK(igraph_subgraph_edges(graph,
	&mst, igraph_ess_vector(&edges_tree),
	/* delete_vertices = */ 0));

    // 2. mark edges outside the tree. Each will form a fundamental cycle
    for(j=0; j<no_of_edges; j++) {
        IGRAPH_ALLOW_INTERRUPTION();
	if (!added_edges[j]) {
              IGRAPH_CHECK(igraph_vector_push_back(&edges_outside, j));
	}
    }

    // 3. for each outside edge, find the tree path that runs parallel to it
    IGRAPH_CHECK(igraph_vector_init(&parent, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_destroy, &parent);
    for(k=0; k<igraph_vector_size(&edges_outside); k++) {
    
        igraph_vector_t edges_cycle=IGRAPH_VECTOR_NULL;
	igraph_integer_t root, tgt, actnode;

        IGRAPH_ALLOW_INTERRUPTION();

	// find the vertices of this outside edge
        igraph_edge(graph, (igraph_integer_t) k, &root, &tgt);
	
	// DFS that sets the parents from the tgt to the root within the tree
	// TODO: check that the vid in the mst and original graph match
	// NOTE: they probably don't
	IGRAPH_CHECK(igraph_dfs(mst, root, IGRAPH_ALL,
				/*unreachable=*/ 0,
				/*order=*/ 0,
				/*order_out=*/ 0,
				/*parent */ &parent,
				/*dist=*/ 0,
				/*in_callback=*/
				igraph_i_cycle_basis_unweighted_undirected_dfs_incb,
				/*out_callback=*/ 0,
				/*extra=*/ &tgt));

	// add the edges inside the tree
	actnode = tgt;
        while (actnode != -1) {
	    IGRAPH_CHECK(igraph_vector_push_back(&edges_cycle, actnode));
	    actnode = VECTOR(parent)[(long int) actnode];
	}
	// no need to add the outside edge to the cycle, it goes without saying
	// that there is an edge the cycle was identified

	IGRAPH_CHECK(igraph_vector_push_back(cycles, &edges_cycle));
    }
    
    // clean up
    // TODO: clean up the mst
    igraph_vector_destroy(&edges_tree);
    igraph_vector_destroy(&edges_outside);
    igraph_vector_destroy(&parent);
    IGRAPH_FINALLY_CLEAN(3);

    return IGRAPH_SUCCESS;
}

igraph_bool_t igraph_i_cycle_basis_unweighted_undirected_dfs_incb(const igraph_t *graph,
						    igraph_integer_t vid,
						    igraph_integer_t dist,
						    void *extra) {

    igraph_integer_t *tgt=extra;
    IGRAPH_UNUSED(graph); IGRAPH_UNUSED(dist);
    return *tgt == vid;
}
