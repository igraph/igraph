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

#include "igraph_graphlets.h"
#include "igraph_memory.h"
#include "igraph_constructors.h"

typedef struct {
  igraph_vector_int_t *newidvectors;
  igraph_t *newgraphs;
  igraph_vector_t *newweights;
  int nc;
} igraph_i_subclique_next_free_t;

void igraph_i_subclique_next_free(void *ptr) {
  igraph_i_subclique_next_free_t *data=ptr;
  if (data->newidvectors) {
    /* TODO */
  }
  if (data->newgraphs) {
    /* TODO */
  }
  if (data->newweights) {
    /* TODO */
  }
}

/**
 * \function igraph_subclique_next
 *
 * \param graph Input graph.
 * \param weight Edge weights.
 * \param ids The ids of the vertices in the input graph.
 * \param cliques A list of vectors, vertex ids for cliques.
 * \param result The result is stored here, a list of graphs is stored
 *        here.
 * \param resultids The ids of the vertices in the result graphs is
 *        stored here.
 * \param clique_thr The thresholds for the cliques are stored here,
 *        if not a null pointer.
 * \param next_thr The next thresholds for the cliques are stored
 *        here, if not a null pointer.
 */

int igraph_subclique_next(const igraph_t *graph,
                          const igraph_vector_t *weights,
                          const igraph_vector_int_t *ids,
                          const igraph_vector_ptr_t *cliques,
                          igraph_vector_ptr_t *result,
                          igraph_vector_ptr_t *resultweights,
                          igraph_vector_ptr_t *resultids,
                          igraph_vector_t *clique_thr,
                          igraph_vector_t *next_thr) {

  igraph_vector_int_t mark, map;
  igraph_vector_int_t edges;
  igraph_vector_t neis, newedges;
  igraph_integer_t c, nc=igraph_vector_ptr_size(cliques);
  igraph_integer_t no_of_nodes=igraph_vcount(graph);
  igraph_integer_t no_of_edges=igraph_ecount(graph);
  igraph_vector_int_t *newidvectors=0;
  igraph_t *newgraphs=0;
  igraph_vector_t *newweights;
  igraph_i_subclique_next_free_t freedata={ newidvectors, newgraphs,
                                            newweights, nc };

  if (igraph_vector_size(weights) != no_of_edges) {
    IGRAPH_ERROR("Invalid length of weight vector", IGRAPH_EINVAL);
  }

  if (igraph_vector_int_size(ids) != no_of_nodes) {
    IGRAPH_ERROR("Invalid length of ID vector", IGRAPH_EINVAL);
  }

  if (igraph_vector_ptr_size(result) != nc) {
    IGRAPH_ERROR("Invalid graph list size", IGRAPH_EINVAL);
  }

  if (igraph_vector_ptr_size(resultweights) != nc) {
    IGRAPH_ERROR("Invalid weight list size", IGRAPH_EINVAL);
  }

  if (igraph_vector_ptr_size(resultids) != nc) {
    IGRAPH_ERROR("Invalid id vector size", IGRAPH_EINVAL);
  }

  IGRAPH_FINALLY(igraph_i_subclique_next_free, &freedata);
  newidvectors=igraph_Calloc(nc, igraph_vector_int_t);
  if (!newidvectors) {
    IGRAPH_ERROR("Cannot calculate next cliques", IGRAPH_ENOMEM);
  }
  freedata.newidvectors = newidvectors;
  newweights=igraph_Calloc(nc, igraph_vector_t);
  if (!newweights) {
    IGRAPH_ERROR("Cannot calculate next cliques", IGRAPH_ENOMEM);
  }
  freedata.newweights = newweights;
  newgraphs=igraph_Calloc(nc, igraph_t);
  if (!newgraphs) {
    IGRAPH_ERROR("Cannot calculate next cliques", IGRAPH_ENOMEM);
  }
  freedata.newgraphs = newgraphs;

  igraph_vector_init(&newedges, 100);
  IGRAPH_FINALLY(igraph_vector_destroy, &newedges);
  igraph_vector_int_init(&mark, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_destroy, &mark);
  igraph_vector_int_init(&map, no_of_nodes);
  IGRAPH_FINALLY(igraph_vector_destroy, &map);
  igraph_vector_int_init(&edges, 100);
  IGRAPH_FINALLY(igraph_vector_int_destroy, &edges);
  igraph_vector_init(&neis, 10);
  IGRAPH_FINALLY(igraph_vector_destroy, &neis);

  if (clique_thr) { igraph_vector_resize(clique_thr, nc); }
  if (next_thr)   { igraph_vector_resize(next_thr,   nc); }

  /* Iterate over all cliques. We will create graphs for all
     subgraphs defined by the cliques. */

  for (c=0; c<nc; c++) {
    igraph_vector_t *clique=VECTOR(*cliques)[c];
    igraph_real_t minweight=IGRAPH_INFINITY, nextweight=IGRAPH_INFINITY;
    igraph_integer_t e, v, clsize=igraph_vector_size(clique);
    igraph_integer_t noe, nov=0;
    igraph_vector_int_t *newids=newidvectors+c;
    igraph_vector_t *neww=newweights+c;
    igraph_t *newgraph=newgraphs+c;
    igraph_vector_int_clear(&edges);
    igraph_vector_clear(&newedges);

    /* --------------------------------------------------- */

    /* Iterate over the vertices of a clique and find the
       edges within the clique, put them in a list.
       At the same time, search for the minimum edge weight within
       the clique and the next edge weight if any. */

    for (v=0; v<clsize; v++) {
      igraph_integer_t i, neilen, node=VECTOR(*clique)[v];
      igraph_incident(graph, &neis, node, IGRAPH_ALL);
      neilen=igraph_vector_size(&neis);
      VECTOR(mark)[node] = c+1;
      for (i=0; i<neilen; i++) {
        igraph_integer_t edge=VECTOR(neis)[i];
        igraph_integer_t nei=IGRAPH_OTHER(graph, edge, node);
        if (VECTOR(mark)[nei] == c+1) {
          igraph_real_t w=VECTOR(*weights)[edge];
          igraph_vector_int_push_back(&edges, edge);
          if (w < minweight) {
            nextweight=minweight;
            minweight=w;
          } else if (w > minweight && w < nextweight) {
            nextweight=w;
          }
        }
      }
    } /* v < clsize */

    /* --------------------------------------------------- */

    /* OK, we have stored the edges and found the weight of
       the clique and the next weight to consider */

    if (clique_thr) { VECTOR(*clique_thr)[c] = minweight;  }
    if (next_thr)   { VECTOR(*next_thr  )[c] = nextweight; }

    /* --------------------------------------------------- */

    /* Now we create the subgraph from the edges above the next
       threshold, and their incident vertices. */

    igraph_vector_int_init(newids, 0);
    VECTOR(*resultids)[c] = newids;
    igraph_vector_init(neww, 0);
    VECTOR(*resultweights)[c] = neww;

    /* We use mark[] to denote the vertices already mapped to
       the new graph. If this is -(c+1), then the vertex was
       mapped, otherwise it was not. The mapping itself is in
       map[]. */

    noe=igraph_vector_int_size(&edges);
    for (e=0; e<noe; e++) {
      igraph_integer_t edge=VECTOR(edges)[e];
      igraph_integer_t from, to;
      igraph_real_t w=VECTOR(*weights)[edge];
      igraph_edge(graph, edge, &from, &to);
      if (w >= nextweight) {
        if (VECTOR(mark)[from] == c+1) {
          VECTOR(map)[from] = nov++;
          VECTOR(mark)[from] = -(c+1);
          igraph_vector_int_push_back(newids, VECTOR(*ids)[from]);
        }
        if (VECTOR(mark)[to] == c+1) {
          VECTOR(map)[to] = nov++;
          VECTOR(mark)[to] = -(c+1);
          igraph_vector_int_push_back(newids, VECTOR(*ids)[to]);
        }
        igraph_vector_push_back(neww, w);
        igraph_vector_push_back(&newedges, VECTOR(map)[from]);
        igraph_vector_push_back(&newedges, VECTOR(map)[to]);
      }
    }

    igraph_create(newgraph, &newedges, nov, IGRAPH_UNDIRECTED);
    VECTOR(*result)[c] = newgraph;

    /* --------------------------------------------------- */

  } /* c < nc */

  igraph_vector_destroy(&neis);
  igraph_vector_int_destroy(&edges);
  igraph_vector_int_destroy(&mark);
  igraph_vector_int_destroy(&map);
  igraph_vector_destroy(&newedges);
  IGRAPH_FINALLY_CLEAN(6);      /* +1 for the result */

  return 0;
}
