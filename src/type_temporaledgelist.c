/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_datatype.h"

#ifdef IGRAPH_DATA_TYPE_TEMPORAL_EDGE_LIST

#include "igraph_interface.h"
#include "igraph_attributes.h"
#include "igraph_memory.h"
#include <string.h>		/* memset & co. */
#include "config.h"

/* Internal functions */

int igraph_i_create_start(igraph_vector_t *res, igraph_vector_t *el,
                          igraph_vector_t *index, 
			  igraph_integer_t nodes);

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

int igraph_empty(igraph_t *graph, igraph_integer_t n,
                 igraph_bool_t directed) {
  return igraph_empty_attrs(graph, n, directed, 0);
}

int igraph_empty_attrs(igraph_t *graph, igraph_integer_t n,
                       igraph_bool_t directed, void *attr) {
  if (n<0) {
    IGRAPH_ERROR("cannot create empty graph with negative "
                 "number of vertices", IGRAPH_EINVAL);
  }
  
  if (!IGRAPH_FINITE(n)) {
    IGRAPH_ERROR("number of vertices is not finite (NA, NaN or Inf)",
                 IGRAPH_EINVAL);
  }

  graph->n=0;
  graph->directed=directed;
  IGRAPH_VECTOR_INIT_FINALLY(&graph->from, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->to, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->oi, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->ii, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->os, 1);
  IGRAPH_VECTOR_INIT_FINALLY(&graph->is, 1);

  VECTOR(graph->os)[0]=0;
  VECTOR(graph->is)[0]=0;

  /* time labels */
  igraph_vector_int_set_null(&graph->vb);
  igraph_vector_int_set_null(&graph->eb);
  igraph_vector_int_set_null(&graph->vd);
  igraph_vector_int_set_null(&graph->ed);

  /* init attributes */
  graph->attr=0;
  IGRAPH_CHECK(igraph_i_attribute_init(graph, attr));

  /* add the vertices */
  IGRAPH_CHECK(igraph_add_vertices(graph, n, 0));
  
  IGRAPH_FINALLY_CLEAN(6);
  return 0;
}

int igraph_destroy(igraph_t *graph) {

  IGRAPH_I_ATTRIBUTE_DESTROY(graph);

  igraph_vector_destroy(&graph->from);
  igraph_vector_destroy(&graph->to);
  igraph_vector_destroy(&graph->oi);
  igraph_vector_destroy(&graph->ii);
  igraph_vector_destroy(&graph->os);
  igraph_vector_destroy(&graph->is);

  /* time labels */
  igraph_vector_int_destroy(&graph->vb);
  igraph_vector_int_destroy(&graph->eb);
  igraph_vector_int_destroy(&graph->vd);
  igraph_vector_int_destroy(&graph->ed);
  

  return 0;
}

int igraph_copy(igraph_t *to, const igraph_t *from) {
  
  to->n=from->n;
  to->directed=from->directed;
  IGRAPH_CHECK(igraph_vector_copy(&to->from, &from->from));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->from);
  IGRAPH_CHECK(igraph_vector_copy(&to->to, &from->to));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->to);
  IGRAPH_CHECK(igraph_vector_copy(&to->oi, &from->oi));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->oi);
  IGRAPH_CHECK(igraph_vector_copy(&to->ii, &from->ii));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->ii);
  IGRAPH_CHECK(igraph_vector_copy(&to->os, &from->os));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->os);
  IGRAPH_CHECK(igraph_vector_copy(&to->is, &from->is));
  IGRAPH_FINALLY(igraph_vector_destroy, &to->is);

  IGRAPH_I_ATTRIBUTE_COPY(to, from, 1,1,1); /* does IGRAPH_CHECK */

  /* time labels */

  IGRAPH_CHECK(igraph_vector_int_copy(&to->vb, &from->vb));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &to->vb);
  IGRAPH_CHECK(igraph_vector_int_copy(&to->eb, &from->eb));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &to->eb);
  IGRAPH_CHECK(igraph_vector_int_copy(&to->vd, &from->vd));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &to->vd);
  IGRAPH_CHECK(igraph_vector_int_copy(&to->ed, &from->ed));
  IGRAPH_FINALLY(igraph_vector_int_destroy, &to->ed);

  IGRAPH_FINALLY_CLEAN(10);
  return 0;
}

int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges,
		     void *attr) {
    long int no_of_edges=igraph_vector_size(&graph->from);
  long int edges_to_add=igraph_vector_size(edges)/2;
  long int i=0;
  igraph_error_handler_t *oldhandler;
  int ret1, ret2;
  igraph_vector_t newoi, newii;
  igraph_bool_t directed=igraph_is_directed(graph);

  if (igraph_vector_size(edges) % 2 != 0) {
    IGRAPH_ERROR("invalid (odd) length of edges vector",
                 IGRAPH_EINVEVECTOR);
  }
  if (!igraph_vector_isininterval(edges, 0, igraph_vcount(graph)-1)) {
    IGRAPH_ERROR("cannot add edges", IGRAPH_EINVVID);
  }

  /* from & to */
  IGRAPH_CHECK(igraph_vector_reserve(&graph->from,
                                     no_of_edges+edges_to_add));
  IGRAPH_CHECK(igraph_vector_reserve(&graph->to,
                                     no_of_edges+edges_to_add));

  while (i<edges_to_add*2) {
    if (directed || VECTOR(*edges)[i] > VECTOR(*edges)[i+1]) {
      igraph_vector_push_back(&graph->from,
                              VECTOR(*edges)[i++]); /* reserved */
      igraph_vector_push_back(&graph->to,
                              VECTOR(*edges)[i++]); /* reserved */
    } else {
      igraph_vector_push_back(&graph->to,
                              VECTOR(*edges)[i++]); /* reserved */
      igraph_vector_push_back(&graph->from,
                              VECTOR(*edges)[i++]); /* reserved */
    }      
  }

  /* disable the error handler temporarily */
  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);
    
  /* oi & ii */
  ret1=igraph_vector_init(&newoi, no_of_edges);
  ret2=igraph_vector_init(&newii, no_of_edges);
  if (ret1 != 0 || ret2 != 0) {
    igraph_vector_resize(&graph->from, no_of_edges); /* gets smaller */
    igraph_vector_resize(&graph->to, no_of_edges);   /* gets smaller */
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
  }  
  ret1=igraph_vector_order(&graph->from, &graph->to, &newoi, graph->n);
  ret2=igraph_vector_order(&graph->to  , &graph->from, &newii,
                           graph->n);
  if (ret1 != 0 || ret2 != 0) {
    igraph_vector_resize(&graph->from, no_of_edges);
    igraph_vector_resize(&graph->to, no_of_edges);
    igraph_vector_destroy(&newoi);
    igraph_vector_destroy(&newii);
    igraph_set_error_handler(oldhandler);
    IGRAPH_ERROR("cannot add edges", IGRAPH_ERROR_SELECT_2(ret1, ret2));
  }  

  /* Attributes */
  if (graph->attr) { 
    igraph_set_error_handler(oldhandler);
    ret1=igraph_i_attribute_add_edges(graph, edges, attr);
    igraph_set_error_handler(igraph_error_handler_ignore);
    if (ret1 != 0) {
      igraph_vector_resize(&graph->from, no_of_edges);
      igraph_vector_resize(&graph->to, no_of_edges);
      igraph_vector_destroy(&newoi);
      igraph_vector_destroy(&newii);
      igraph_set_error_handler(oldhandler);
      IGRAPH_ERROR("cannot add edges", ret1);
    }  
  }
  
  /* os & is, its length does not change, error safe */
  igraph_i_create_start(&graph->os, &graph->from, &newoi, graph->n);
  igraph_i_create_start(&graph->is, &graph->to  , &newii, graph->n);

  /* TODO: time labels */

  /* everything went fine  */
  igraph_vector_destroy(&graph->oi);
  igraph_vector_destroy(&graph->ii);
  graph->oi=newoi;
  graph->ii=newii;
  igraph_set_error_handler(oldhandler);
  
  return 0;
}

int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, 
			void *attr) {
  long int ec=igraph_ecount(graph);
  long int i;

  if (nv < 0) {
    IGRAPH_ERROR("cannot add negative number of vertices",
                 IGRAPH_EINVAL);
  }

  IGRAPH_CHECK(igraph_vector_reserve(&graph->os, graph->n+nv+1));
  IGRAPH_CHECK(igraph_vector_reserve(&graph->is, graph->n+nv+1));
  
  igraph_vector_resize(&graph->os, graph->n+nv+1); /* reserved */
  igraph_vector_resize(&graph->is, graph->n+nv+1); /* reserved */
  for (i=graph->n+1; i<graph->n+nv+1; i++) {
    VECTOR(graph->os)[i]=ec;
    VECTOR(graph->is)[i]=ec;
  }
  
  graph->n += nv;   
  
  if (graph->attr) {
    IGRAPH_CHECK(igraph_i_attribute_add_vertices(graph, nv, attr));
  }

  /* TODO: time labels */

  return 0;
}

int igraph_delete_edges(igraph_t *graph, igraph_es_t edges) {
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_nodes=igraph_vcount(graph);
  long int edges_to_remove=0;
  long int remaining_edges;
  igraph_eit_t eit;
  
  igraph_vector_t newfrom, newto, newoi;

  int *mark;
  long int i, j;
  
  mark=igraph_Calloc(no_of_edges, int);
  if (mark==0) {
    IGRAPH_ERROR("Cannot delete edges", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, mark);

  IGRAPH_CHECK(igraph_eit_create(graph, edges, &eit));
  IGRAPH_FINALLY(igraph_eit_destroy, &eit);

  for (IGRAPH_EIT_RESET(eit); !IGRAPH_EIT_END(eit);
       IGRAPH_EIT_NEXT(eit)) {
    long int e=IGRAPH_EIT_GET(eit);
    if (mark[e]==0) {
      edges_to_remove++;
      mark[e]++;
    }
  }
  remaining_edges=no_of_edges-edges_to_remove;

  /* We don't need the iterator any more */
  igraph_eit_destroy(&eit);
  IGRAPH_FINALLY_CLEAN(1);

  IGRAPH_VECTOR_INIT_FINALLY(&newfrom, remaining_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&newto, remaining_edges);
  
  /* Actually remove the edges, move from pos i to pos j in 
     newfrom/newto */
  for (i=0,j=0; j<remaining_edges; i++) {
    if (mark[i]==0) {
      VECTOR(newfrom)[j] = VECTOR(graph->from)[i];
      VECTOR(newto)[j] = VECTOR(graph->to)[i];
      j++;
    }
  }

  /* Create index, this might require additional memory */
  IGRAPH_VECTOR_INIT_FINALLY(&newoi, remaining_edges);
  IGRAPH_CHECK(igraph_vector_order(&newfrom, &newto, &newoi,
                                   no_of_nodes));
  IGRAPH_CHECK(igraph_vector_order(&newto, &newfrom, &graph->ii,
                                   no_of_nodes));

  /* Edge attributes, we need an index that gives the ids of the 
     original edges for every new edge. 
  */
  if (graph->attr) {    
    igraph_vector_t idx;
    IGRAPH_VECTOR_INIT_FINALLY(&idx, remaining_edges);
    for (i=0, j=0; i<no_of_edges; i++) {
      if (mark[i] == 0) {
	VECTOR(idx)[j++] = i;
      }
    }
    IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, graph, &idx));
    igraph_vector_destroy(&idx);
    IGRAPH_FINALLY_CLEAN(1);
  }

  /* Ok, we've all memory needed, free the old structure  */
  igraph_vector_destroy(&graph->from);
  igraph_vector_destroy(&graph->to);
  igraph_vector_destroy(&graph->oi);
  graph->from=newfrom;
  graph->to=newto;
  graph->oi=newoi;
  IGRAPH_FINALLY_CLEAN(3);

  igraph_Free(mark);
  IGRAPH_FINALLY_CLEAN(1);
  
  /* Create start vectors, no memory is needed for this */
  igraph_i_create_start(&graph->os, &graph->from, &graph->oi, 
			(igraph_integer_t) no_of_nodes);
  igraph_i_create_start(&graph->is, &graph->to,   &graph->ii, 
			(igraph_integer_t) no_of_nodes);
  
  /* TODO: time labels */

  /* Nothing to deallocate... */
  return 0;
}

int igraph_delete_vertices(igraph_t *graph,
                           const igraph_vs_t vertices) {
  return igraph_delete_vertices_idx(graph, vertices, /* idx= */ 0, 
				    /* invidx= */ 0);
}

int igraph_delete_vertices_idx(igraph_t *graph,
                               const igraph_vs_t vertices, 
			       igraph_vector_t *idx, 
			       igraph_vector_t *invidx) {
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t edge_recoding, vertex_recoding;
  igraph_vector_t *my_vertex_recoding=&vertex_recoding;
  igraph_vit_t vit;
  igraph_t newgraph = *graph;
  long int i, j;
  long int remaining_vertices, remaining_edges;

  if (idx) {
    my_vertex_recoding=idx;
    IGRAPH_CHECK(igraph_vector_resize(idx, no_of_nodes));
    igraph_vector_null(idx);
  } else { 
    IGRAPH_VECTOR_INIT_FINALLY(&vertex_recoding, no_of_nodes);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edge_recoding, no_of_edges);
 
  IGRAPH_CHECK(igraph_vit_create(graph, vertices, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);
  
  /* mark the vertices to delete */
  for (; !IGRAPH_VIT_END(vit); IGRAPH_VIT_NEXT(vit) ) {
    long int vertex=IGRAPH_VIT_GET(vit);
    if (vertex < 0 || vertex >= no_of_nodes) {
      IGRAPH_ERROR("Cannot delete vertices", IGRAPH_EINVVID);
    }
    VECTOR(*my_vertex_recoding)[vertex]=1;
  }
  /* create vertex recoding vector */
  for (remaining_vertices=0, i=0; i<no_of_nodes; i++) {
    if (VECTOR(*my_vertex_recoding)[i]==0) {
      VECTOR(*my_vertex_recoding)[i]=remaining_vertices+1;
      remaining_vertices++;
    } else { 
      VECTOR(*my_vertex_recoding)[i]=0;
    }
  }
  /* create edge recoding vector */
  for (remaining_edges=0, i=0; i<no_of_edges; i++) {
    long int from=(long int) VECTOR(graph->from)[i];
    long int to=(long int) VECTOR(graph->to)[i];
    if (VECTOR(*my_vertex_recoding)[from] != 0 &&
	VECTOR(*my_vertex_recoding)[to  ] != 0) {
      VECTOR(edge_recoding)[i]=remaining_edges+1;
      remaining_edges++;
    } 
  }

  /* start creating the graph, directed was copied already */
  newgraph.n=(igraph_integer_t) remaining_vertices;

  /* allocate vectors */
  IGRAPH_VECTOR_INIT_FINALLY(&newgraph.from, remaining_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&newgraph.to, remaining_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&newgraph.oi, remaining_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&newgraph.ii, remaining_edges);
  IGRAPH_VECTOR_INIT_FINALLY(&newgraph.os, remaining_vertices+1);
  IGRAPH_VECTOR_INIT_FINALLY(&newgraph.is, remaining_vertices+1);
  
  /* Add the edges */
  for (i=0, j=0; j<remaining_edges; i++) {
    if (VECTOR(edge_recoding)[i]>0) {
      long int from=(long int) VECTOR(graph->from)[i];
      long int to=(long int) VECTOR(graph->to  )[i];
      VECTOR(newgraph.from)[j]=VECTOR(*my_vertex_recoding)[from]-1;
      VECTOR(newgraph.to  )[j]=VECTOR(*my_vertex_recoding)[to]-1;
      j++;
    }
  }
  /* update oi & ii */
  IGRAPH_CHECK(igraph_vector_order(&newgraph.from, &newgraph.to,
                                   &newgraph.oi, 
				   remaining_vertices));
  IGRAPH_CHECK(igraph_vector_order(&newgraph.to, &newgraph.from,
                                   &newgraph.ii, 
				   remaining_vertices));  

  IGRAPH_CHECK(igraph_i_create_start(&newgraph.os, &newgraph.from, 
				     &newgraph.oi, (igraph_integer_t) 
				     remaining_vertices));
  IGRAPH_CHECK(igraph_i_create_start(&newgraph.is, &newgraph.to,
				     &newgraph.ii, (igraph_integer_t) 
				     remaining_vertices));
  
  /* attributes */
  IGRAPH_I_ATTRIBUTE_COPY(&newgraph, graph, 
			  /*graph=*/ 1, /*vertex=*/0, /*edge=*/0);
  IGRAPH_FINALLY_CLEAN(6);
  IGRAPH_FINALLY(igraph_destroy, &newgraph);

  if (newgraph.attr) {
    igraph_vector_t iidx;
    IGRAPH_VECTOR_INIT_FINALLY(&iidx, remaining_vertices);
    for (i=0; i<no_of_nodes; i++) {
      long int jj=(long int) VECTOR(*my_vertex_recoding)[i];
      if (jj != 0) {
	VECTOR(iidx)[ jj-1 ] = i;
      }
    }
    IGRAPH_CHECK(igraph_i_attribute_permute_vertices(graph,
						     &newgraph,
						     &iidx));
    IGRAPH_CHECK(igraph_vector_resize(&iidx, remaining_edges));
    for (i=0; i<no_of_edges; i++) {
      long int jj=(long int) VECTOR(edge_recoding)[i];
      if (jj != 0) {
	VECTOR(iidx)[ jj-1 ] = i;
      }
    }
    IGRAPH_CHECK(igraph_i_attribute_permute_edges(graph, &newgraph,
                                                  &iidx));
    igraph_vector_destroy(&iidx);
    IGRAPH_FINALLY_CLEAN(1);
  }
	       
  igraph_vit_destroy(&vit);
  igraph_vector_destroy(&edge_recoding);
  igraph_destroy(graph);
  *graph=newgraph;

  IGRAPH_FINALLY_CLEAN(3);  

  /* TODO: this is duplicate */
  if (invidx) {
    IGRAPH_CHECK(igraph_vector_resize(invidx, remaining_vertices));
    for (i=0; i<no_of_nodes; i++) {
      long int newid=(long int) VECTOR(*my_vertex_recoding)[i];
      if (newid != 0) {
	VECTOR(*invidx)[newid-1] = i;
      }
    }
  }

  /* TODO: time labels */

  if (!idx) {
    igraph_vector_destroy(my_vertex_recoding);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

igraph_integer_t igraph_vcount(const igraph_t *graph) {
  return graph->n;
}

igraph_integer_t igraph_ecount(const igraph_t *graph) {
  return (igraph_integer_t) igraph_vector_size(&graph->from);
}

int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis,
                     igraph_integer_t pnode, 
		     igraph_neimode_t mode) {
  long int length=0, idx=0;   
  long int i, j;

  long int node=pnode;

  if (node<0 || node>igraph_vcount(graph)-1) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVMODE);
  }

  if (! graph->directed) {
    mode=IGRAPH_ALL;
  }

  /* Calculate needed space first & allocate it*/

  if (mode & IGRAPH_OUT) {
    length += (VECTOR(graph->os)[node+1] - VECTOR(graph->os)[node]);
  }
  if (mode & IGRAPH_IN) {
    length += (VECTOR(graph->is)[node+1] - VECTOR(graph->is)[node]);
  }
  
  IGRAPH_CHECK(igraph_vector_resize(neis, length));
  
  if (!igraph_is_directed(graph) || mode != IGRAPH_ALL) {

    if (mode & IGRAPH_OUT) {
      j=(long int) VECTOR(graph->os)[node+1];
      for (i=(long int) VECTOR(graph->os)[node]; i<j; i++) {
	VECTOR(*neis)[idx++] = 
	  VECTOR(graph->to)[ (long int)VECTOR(graph->oi)[i] ];
      }
    }
    if (mode & IGRAPH_IN) {
      j=(long int) VECTOR(graph->is)[node+1];
      for (i=(long int) VECTOR(graph->is)[node]; i<j; i++) {
	VECTOR(*neis)[idx++] =
	  VECTOR(graph->from)[ (long int)VECTOR(graph->ii)[i] ];
      }
    }
  } else {
    /* both in- and out- neighbors in a directed graph,
       we need to merge the two 'vectors' */
    long int jj1=(long int) VECTOR(graph->os)[node+1];
    long int j2=(long int) VECTOR(graph->is)[node+1];
    long int i1=(long int) VECTOR(graph->os)[node];
    long int i2=(long int) VECTOR(graph->is)[node];
    while (i1 < jj1 && i2 < j2) {
      long int n1=(long int) VECTOR(graph->to)[ 
		       (long int)VECTOR(graph->oi)[i1] ];
      long int n2=(long int) VECTOR(graph->from)[
		       (long int)VECTOR(graph->ii)[i2] ];
      if (n1<n2) {
	VECTOR(*neis)[idx++]=n1;
	i1++;
      } else if (n1>n2) {
	VECTOR(*neis)[idx++]=n2;
	i2++;
      } else {
	VECTOR(*neis)[idx++]=n1;
	VECTOR(*neis)[idx++]=n2;
	i1++;
	i2++;
      }
    }
    while (i1 < jj1) {
      long int n1=(long int) VECTOR(graph->to)[ 
		       (long int)VECTOR(graph->oi)[i1] ];
      VECTOR(*neis)[idx++]=n1;
      i1++;
    }
    while (i2 < j2) {
      long int n2=(long int) VECTOR(graph->from)[ 
		       (long int)VECTOR(graph->ii)[i2] ];
      VECTOR(*neis)[idx++]=n2;
      i2++;
    }
  }

  return 0;
}

igraph_bool_t igraph_is_directed(const igraph_t *graph) {
  return graph->directed;
}

int igraph_degree(const igraph_t *graph, igraph_vector_t *res, 
		  const igraph_vs_t vids, igraph_neimode_t mode, 
		  igraph_bool_t loops)  {
  long int nodes_to_calc;
  long int i, j;
  igraph_vit_t vit;

  IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
  IGRAPH_FINALLY(igraph_vit_destroy, &vit);

  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("degree calculation failed", IGRAPH_EINVMODE);
  }
  
  nodes_to_calc=IGRAPH_VIT_SIZE(vit);
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }

  IGRAPH_CHECK(igraph_vector_resize(res, nodes_to_calc));
  igraph_vector_null(res);

  if (loops) {
    if (mode & IGRAPH_OUT) {
      for (IGRAPH_VIT_RESET(vit), i=0; 
	   !IGRAPH_VIT_END(vit); 
	   IGRAPH_VIT_NEXT(vit), i++) {
	long int vid=IGRAPH_VIT_GET(vit);
	VECTOR(*res)[i] +=
          (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
      }
    }
    if (mode & IGRAPH_IN) {
      for (IGRAPH_VIT_RESET(vit), i=0; 
	   !IGRAPH_VIT_END(vit); 
	   IGRAPH_VIT_NEXT(vit), i++) {
	long int vid=IGRAPH_VIT_GET(vit);
	VECTOR(*res)[i] +=
          (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
      }
    }
  } else { /* no loops */
    if (mode & IGRAPH_OUT) {
      for (IGRAPH_VIT_RESET(vit), i=0; 
	   !IGRAPH_VIT_END(vit); 
	   IGRAPH_VIT_NEXT(vit), i++) {
	long int vid=IGRAPH_VIT_GET(vit);
	VECTOR(*res)[i] +=
          (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
	for (j=(long int) VECTOR(graph->os)[vid]; 
	     j<VECTOR(graph->os)[vid+1]; j++) {
	  if (VECTOR(graph->to)[ (long int)VECTOR(graph->oi)[j] ] == 
              vid) {
	    VECTOR(*res)[i] -= 1;
	  }
	}
      }
    }
    if (mode & IGRAPH_IN) {
      for (IGRAPH_VIT_RESET(vit), i=0; 
	   !IGRAPH_VIT_END(vit);
	   IGRAPH_VIT_NEXT(vit), i++) {
	long int vid=IGRAPH_VIT_GET(vit);
	VECTOR(*res)[i] +=
          (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
	for (j=(long int) VECTOR(graph->is)[vid]; 
	     j<VECTOR(graph->is)[vid+1]; j++) {
	  if (VECTOR(graph->from)[ (long int)VECTOR(graph->ii)[j] ] ==
              vid) {
	    VECTOR(*res)[i] -= 1;
	  }
	}
      }
    }
  }  /* loops */

  igraph_vit_destroy(&vit);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_edge(const igraph_t *graph, igraph_integer_t eid, 
		igraph_integer_t *from, igraph_integer_t *to) {
  *from = (igraph_integer_t) VECTOR(graph->from)[(long int)eid];
  *to   = (igraph_integer_t) VECTOR(graph->to  )[(long int)eid];
  
  if (! igraph_is_directed(graph) && *from > *to) {
    igraph_integer_t tmp=*from;
    *from=*to;
    *to=tmp;
  }

  return 0;
}

int igraph_edges(const igraph_t *graph, igraph_es_t eids,
		 igraph_vector_t *edges) {
    
  igraph_eit_t eit;
  long int n, ptr=0;

  IGRAPH_CHECK(igraph_eit_create(graph, eids, &eit));
  IGRAPH_FINALLY(igraph_eit_destroy, &eit);
  n=IGRAPH_EIT_SIZE(eit);
  IGRAPH_CHECK(igraph_vector_resize(edges, n*2));
  for (; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit)) {
    long int e=IGRAPH_EIT_GET(eit);
    VECTOR(*edges)[ptr++]=IGRAPH_FROM(graph, e);
    VECTOR(*edges)[ptr++]=IGRAPH_TO(graph, e);
  }
  
  igraph_eit_destroy(&eit);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

#define BINSEARCH(start,end,value,iindex,edgelist,N,pos)    \
  do {                                                      \
  while ((start) < (end)) {                                 \
    long int mid=(start)+((end)-(start))/2;                 \
    long int e=(long int) VECTOR((iindex))[mid];	    \
    if (VECTOR((edgelist))[e] < (value)) {                  \
      (start)=mid+1;                                        \
    } else {                                                \
      (end)=mid;                                            \
    }                                                       \
  }                                                         \
  if ((start)<(N)) {                                        \
    long int e=(long int) VECTOR((iindex))[(start)];        \
    if (VECTOR((edgelist))[e] == (value)) {                 \
      *(pos)=(igraph_integer_t) e;			  \
    }                                                       \
  } } while(0)

#define FIND_DIRECTED_EDGE(graph,xfrom,xto,eid)                     \
  do {                                                              \
    long int start=(long int) VECTOR(graph->os)[xfrom];	     \
    long int end=(long int) VECTOR(graph->os)[xfrom+1];	     \
    long int N=end;                                                 \
    long int start2=(long int) VECTOR(graph->is)[xto];	      \
    long int end2=(long int) VECTOR(graph->is)[xto+1];	      \
    long int N2=end2;                                               \
    if (end-start<end2-start2) {                                    \
      BINSEARCH(start,end,xto,graph->oi,graph->to,N,eid);           \
    } else {                                                        \
      BINSEARCH(start2,end2,xfrom,graph->ii,graph->from,N2,eid);    \
    }                                                               \
  } while (0)

#define FIND_UNDIRECTED_EDGE(graph,from,to,eid)                     \
  do {                                                              \
    long int xfrom1= from > to ? from : to;                         \
    long int xto1= from > to ? to : from;                           \
    FIND_DIRECTED_EDGE(graph,xfrom1,xto1,eid);                      \
  } while (0)

int igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid,
		   igraph_integer_t pfrom, igraph_integer_t pto,
		   igraph_bool_t directed, igraph_bool_t error) {
  long int from=pfrom, to=pto;
  long int nov=igraph_vcount(graph);

  if (from < 0 || to < 0 || from > nov-1 || to > nov-1) {
    IGRAPH_ERROR("cannot get edge id", IGRAPH_EINVVID);
  }

  *eid=-1;
  if (igraph_is_directed(graph)) {

    /* Directed graph */
    FIND_DIRECTED_EDGE(graph,from,to,eid);
    if (!directed && *eid < 0) {
      FIND_DIRECTED_EDGE(graph,to,from,eid);
    }
    
  } else {

    /* Undirected graph, they only have one mode */
    FIND_UNDIRECTED_EDGE(graph,from,to,eid);

  }

  if (*eid < 0) {
    if (error) {
      IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
    }
  }
  
  return IGRAPH_SUCCESS;  
}

int igraph_get_eids_pairs(const igraph_t *graph, igraph_vector_t *eids,
			  const igraph_vector_t *pairs, 
			  igraph_bool_t directed, igraph_bool_t error);

int igraph_get_eids_path(const igraph_t *graph, igraph_vector_t *eids,
			 const igraph_vector_t *path, 
			 igraph_bool_t directed, igraph_bool_t error);

int igraph_get_eids_pairs(const igraph_t *graph, igraph_vector_t *eids,
			  const igraph_vector_t *pairs, 
			  igraph_bool_t directed, igraph_bool_t error) {
  long int n=igraph_vector_size(pairs);
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_integer_t eid=-1;
    
  if (n % 2 != 0) {
    IGRAPH_ERROR("Cannot get edge ids, invalid length of edge ids",
		 IGRAPH_EINVAL);
  }
  if (!igraph_vector_isininterval(pairs, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("Cannot get edge ids, invalid vertex id",
                 IGRAPH_EINVVID);
  }

  IGRAPH_CHECK(igraph_vector_resize(eids, n/2));
  
  if (igraph_is_directed(graph)) {
    for (i=0; i<n/2; i++) {
      long int from=(long int) VECTOR(*pairs)[2*i];
      long int to=(long int) VECTOR(*pairs)[2*i+1];

      eid=-1;
      FIND_DIRECTED_EDGE(graph,from,to,&eid);
      if (!directed && eid < 0) {
	FIND_DIRECTED_EDGE(graph,to,from,&eid);
      }
      
      VECTOR(*eids)[i]=eid;
      if (eid < 0 && error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  } else {
    for (i=0; i<n/2; i++) {
      long int from=(long int) VECTOR(*pairs)[2*i];
      long int to=(long int) VECTOR(*pairs)[2*i+1];

      eid=-1;
      FIND_UNDIRECTED_EDGE(graph,from,to,&eid);
      VECTOR(*eids)[i]=eid;
      if (eid < 0 && error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  }
  
  return 0;
}

int igraph_get_eids_path(const igraph_t *graph, igraph_vector_t *eids,
			 const igraph_vector_t *path, 
			 igraph_bool_t directed, igraph_bool_t error) {

  long int n=igraph_vector_size(path);
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  igraph_integer_t eid=-1;

  if (!igraph_vector_isininterval(path, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("Cannot get edge ids, invalid vertex id",
                 IGRAPH_EINVVID);
  }
  
  IGRAPH_CHECK(igraph_vector_resize(eids, n==0 ? 0 : n-1));
  
  if (igraph_is_directed(graph)) {
    for (i=0; i<n-1; i++) {
      long int from=(long int) VECTOR(*path)[i];
      long int to=(long int) VECTOR(*path)[i+1];
      
      eid=-1;
      FIND_DIRECTED_EDGE(graph, from, to, &eid);
      if (!directed && eid < 0) {
	FIND_DIRECTED_EDGE(graph, to, from, &eid);
      }
      
      VECTOR(*eids)[i]=eid;
      if (eid < 0 && error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  } else {
    for (i=0; i<n-1; i++) {
      long int from=(long int) VECTOR(*path)[i];
      long int to=(long int) VECTOR(*path)[i+1];
      
      eid=-1;
      FIND_UNDIRECTED_EDGE(graph, from, to, &eid);
      VECTOR(*eids)[i]=eid;
      if (eid < 0 && error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  }
  
  return 0;
}

int igraph_get_eids(const igraph_t *graph, igraph_vector_t *eids,
                    const igraph_vector_t *pairs,
		    const igraph_vector_t *path,
	    igraph_bool_t directed, igraph_bool_t error) {

  if (!pairs && !path) {
    igraph_vector_clear(eids);
    return 0;
  } else if (pairs && !path) {
    return igraph_get_eids_pairs(graph, eids, pairs, directed, error);
  } else if (!pairs && path) { 
    return igraph_get_eids_path(graph, eids, path, directed, error);
  } else {
    /* both */
    igraph_vector_t tmp;
    IGRAPH_VECTOR_INIT_FINALLY(&tmp, 0);
    IGRAPH_CHECK(igraph_get_eids_pairs(graph, eids, pairs, directed,
                                       error));
    IGRAPH_CHECK(igraph_get_eids_path(graph, &tmp, path, directed,
                                      error));
    IGRAPH_CHECK(igraph_vector_append(eids, &tmp));
    igraph_vector_destroy(&tmp);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
  }
}

#undef BINSEARCH
#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE

#define BINSEARCH(start,end,value,iindex,edgelist,N,pos,seen)    \
  do {                                                           \
  while ((start) < (end)) {                                      \
    long int mid=(start)+((end)-(start))/2;                      \
    long int e=(long int) VECTOR((iindex))[mid];	         \
    if (VECTOR((edgelist))[e] < (value)) {                       \
      (start)=mid+1;                                             \
    } else {                                                     \
      (end)=mid;                                                 \
    }                                                            \
  }                                                              \
  if ((start)<(N)) {                                             \
    long int e=(long int) VECTOR((iindex))[(start)];	     \
    while ((start)<(N) && seen[e] &&                             \
           VECTOR(edgelist)[e] == (value)) {	             \
      (start)++;					         \
      e=(long int) VECTOR(iindex)[(start)];		      \
    }				                            \
    if ((start)<(N) && !(seen[e]) &&                             \
        VECTOR(edgelist)[e] == (value)) {	                \
      *(pos)=(igraph_integer_t) e;			       \
    }                                                            \
  } } while(0)

#define FIND_DIRECTED_EDGE(graph,xfrom,xto,eid,seen)	          \
  do {                                                                \
    long int start=(long int) VECTOR(graph->os)[xfrom];               \
    long int end=(long int) VECTOR(graph->os)[xfrom+1];               \
    long int N=end;                                                   \
    long int start2=(long int) VECTOR(graph->is)[xto];	        \
    long int end2=(long int) VECTOR(graph->is)[xto+1];	        \
    long int N2=end2;                                                 \
    if (end-start<end2-start2) {                                      \
      BINSEARCH(start,end,xto,graph->oi,graph->to,N,eid,seen);	\
    } else {                                                          \
      BINSEARCH(start2,end2,xfrom,graph->ii,graph->from,N2,eid,seen); \
    }                                                                 \
  } while (0)

#define FIND_UNDIRECTED_EDGE(graph,from,to,eid,seen)		  \
  do {                                                                \
    long int xfrom1= from > to ? from : to;                           \
    long int xto1= from > to ? to : from;                             \
    FIND_DIRECTED_EDGE(graph,xfrom1,xto1,eid,seen);		   \
  } while (0)


int igraph_get_eids_multipairs(const igraph_t *graph,
                               igraph_vector_t *eids,
			       const igraph_vector_t *pairs, 
			       igraph_bool_t directed,
                               igraph_bool_t error);

int igraph_get_eids_multipath(const igraph_t *graph,
                              igraph_vector_t *eids,
			      const igraph_vector_t *path, 
			      igraph_bool_t directed,
                              igraph_bool_t error);

int igraph_get_eids_multipairs(const igraph_t *graph,
                               igraph_vector_t *eids,
			       const igraph_vector_t *pairs, 
			       igraph_bool_t directed, 
                               igraph_bool_t error) {

  long int n=igraph_vector_size(pairs);
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_bool_t *seen;
  long int i;
  igraph_integer_t eid=-1;
    
  if (n % 2 != 0) {
    IGRAPH_ERROR("Cannot get edge ids, invalid length of edge ids",
		 IGRAPH_EINVAL);
  }
  if (!igraph_vector_isininterval(pairs, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("Cannot get edge ids, invalid vertex id",
                 IGRAPH_EINVVID);
  }

  seen=igraph_Calloc(no_of_edges, igraph_bool_t);
  if (seen==0) {
    IGRAPH_ERROR("Cannot get edge ids", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, seen);
  IGRAPH_CHECK(igraph_vector_resize(eids, n/2));
  
  if (igraph_is_directed(graph)) {
    for (i=0; i<n/2; i++) {
      long int from=(long int) VECTOR(*pairs)[2*i];
      long int to=(long int) VECTOR(*pairs)[2*i+1];

      eid=-1;
      FIND_DIRECTED_EDGE(graph,from,to,&eid,seen);
      if (!directed && eid < 0) {
	FIND_DIRECTED_EDGE(graph,to,from,&eid,seen);
      }
      
      VECTOR(*eids)[i]=eid;
      if (eid >= 0) {
	seen[(long int)(eid)]=1;
      } else if (error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  } else {
    for (i=0; i<n/2; i++) {
      long int from=(long int) VECTOR(*pairs)[2*i];
      long int to=(long int) VECTOR(*pairs)[2*i+1];

      eid=-1;
      FIND_UNDIRECTED_EDGE(graph,from,to,&eid,seen);
      VECTOR(*eids)[i]=eid;
      if (eid >= 0) {
	seen[(long int)(eid)]=1;
      } else if (error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  }
  
  igraph_Free(seen);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_get_eids_multipath(const igraph_t *graph, 
                              igraph_vector_t *eids,
			      const igraph_vector_t *path, 
			      igraph_bool_t directed,
                              igraph_bool_t error) {
  
  long int n=igraph_vector_size(path);
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_bool_t *seen;
  long int i;
  igraph_integer_t eid=-1;

  if (!igraph_vector_isininterval(path, 0, no_of_nodes-1)) {
    IGRAPH_ERROR("Cannot get edge ids, invalid vertex id",
                 IGRAPH_EINVVID);
  }
  
  seen=igraph_Calloc(no_of_edges, igraph_bool_t);
  if (!seen) {
    IGRAPH_ERROR("Cannot get edge ids", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, seen);
  IGRAPH_CHECK(igraph_vector_resize(eids, n==0 ? 0 : n-1));
  
  if (igraph_is_directed(graph)) {
    for (i=0; i<n-1; i++) {
      long int from=(long int) VECTOR(*path)[i];
      long int to=(long int) VECTOR(*path)[i+1];
      
      eid=-1;
      FIND_DIRECTED_EDGE(graph, from, to, &eid, seen);
      if (!directed && eid < 0) {
	FIND_DIRECTED_EDGE(graph, to, from, &eid, seen);
      }
      
      VECTOR(*eids)[i]=eid;
      if (eid >= 0) {
	seen[(long int)(eid)]=1;
      } else if (error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  } else {
    for (i=0; i<n-1; i++) {
      long int from=(long int) VECTOR(*path)[i];
      long int to=(long int) VECTOR(*path)[i+1];
      
      eid=-1;
      FIND_UNDIRECTED_EDGE(graph, from, to, &eid, seen);
      VECTOR(*eids)[i]=eid;
      if (eid >= 0) {
	seen[(long int)(eid)]=1;
      } else if (error) {
	IGRAPH_ERROR("Cannot get edge id, no such edge", IGRAPH_EINVAL);
      }
    }
  }
  
  igraph_Free(seen);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

#undef BINSEARCH
#undef FIND_DIRECTED_EDGE
#undef FIND_UNDIRECTED_EDGE

int igraph_get_eids_multi(const igraph_t *graph, igraph_vector_t *eids,
			  const igraph_vector_t *pairs, 
			  const igraph_vector_t *path,
			  igraph_bool_t directed, igraph_bool_t error) {
  if (!pairs && !path) {
    igraph_vector_clear(eids);
    return 0;
  } else if (pairs && !path) {
    return igraph_get_eids_multipairs(graph, eids, pairs, directed,
                                      error);
  } else if (!pairs && path) {
    return igraph_get_eids_multipath(graph, eids, path, directed,
                                     error);
  } else { /* both */
    IGRAPH_ERROR("Give `pairs' or `path' but not both", IGRAPH_EINVAL);
  }
}

int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids,
                    igraph_integer_t pnode,
		    igraph_neimode_t mode) {
  IGRAPH_WARNING("igraph_adjacent is deprecated, use igraph_incident");
  return igraph_incident(graph, eids, pnode, mode);  
}

int igraph_incident(const igraph_t *graph, igraph_vector_t *eids,
                    igraph_integer_t pnode,
		    igraph_neimode_t mode) {
  long int length=0, idx=0;   
  long int i, j;

  long int node=pnode;

  if (node<0 || node>igraph_vcount(graph)-1) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVVID);
  }
  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVMODE);
  }

  if (! graph->directed) {
    mode=IGRAPH_ALL;
  }

  /* Calculate needed space first & allocate it*/

  if (mode & IGRAPH_OUT) {
    length += (VECTOR(graph->os)[node+1] - VECTOR(graph->os)[node]);
  }
  if (mode & IGRAPH_IN) {
    length += (VECTOR(graph->is)[node+1] - VECTOR(graph->is)[node]);
  }
  
  IGRAPH_CHECK(igraph_vector_resize(eids, length));
  
  if (mode & IGRAPH_OUT) {
    j=(long int) VECTOR(graph->os)[node+1];
    for (i=(long int) VECTOR(graph->os)[node]; i<j; i++) {
      VECTOR(*eids)[idx++] = VECTOR(graph->oi)[i];
    }
  }
  if (mode & IGRAPH_IN) {
    j=(long int) VECTOR(graph->is)[node+1];
    for (i=(long int) VECTOR(graph->is)[node]; i<j; i++) {
      VECTOR(*eids)[idx++] = VECTOR(graph->ii)[i];
    }
  }

  return 0;
}
                    
/* -------------------------------------------------- */
/* Interface, temporal functions                      */
/* -------------------------------------------------- */

int igraph_empty_at(igraph_t *graph, igraph_integer_t n,
                igraph_bool_t directed,
                const igraph_vector_time_t *v_active,
                const igraph_vector_time_t *v_inactive,
                void *attrs) {
  /* TODO */
}

igraph_integer_t igraph_vcount_at(const igraph_t *graph,
                igraph_time_t at) {
  /* TODO */
}

igraph_integer_t igraph_ecount_at(const igraph_t *graph,
                igraph_time_t at) {
  /* TODO */
}

int igraph_neighbors_at(const igraph_t *graph, igraph_vector_t *neis,
                igraph_integer_t vid, igraph_neimode_t mode,
                igraph_time_t at) {
  /* TODO */
}

int igraph_incident_at(const igraph_t *graph, igraph_vector_t *eids,
                igraph_integer_t vid, igraph_neimode_t mode,
                igraph_time_t at) {
  /* TODO */
}
                      
int igraph_degree_at(const igraph_t *graph, igraph_vector_t *res, 
		const igraph_vs_t vids, igraph_neimode_t mode, 
		igraph_bool_t loops, igraph_time_t at) {
  /* TODO */
}

int igraph_get_eid_at(const igraph_t *graph, igraph_integer_t *eid,
		igraph_integer_t from, igraph_integer_t to,
		igraph_bool_t directed, igraph_bool_t error,
                igraph_time_t at) {
  /* TODO */
}

int igraph_get_eids_at(const igraph_t *graph, igraph_vector_t *eids,
                const igraph_vector_t *pairs,
		const igraph_vector_t *path,
		igraph_bool_t directed, igraph_bool_t error,
                igraph_time_t at) {
  /* TODO */
}

int igraph_get_eids_multi_at(const igraph_t *graph,
                igraph_vector_t *eids, const igraph_vector_t *pairs, 
	        const igraph_vector_t *path, igraph_bool_t directed,
                igraph_bool_t error, igraph_time_t at) {
  /* TODO */
}

int igraph_add_edges_at(igraph_t *graph, const igraph_vector_t *edges,
                const igraph_vector_int_t *e_active,
                const igraph_vector_int_t *e_inactive, void *attr) {
  /* TODO */
}

int igraph_add_vertices_at(igraph_t *graph, igraph_integer_t nv, 
		const igraph_vector_int_t *v_active,
                const igraph_vector_int_t *v_inacive, void *attr) {
  /* TODO */
}

int igraph_i_create_start(igraph_vector_t *res, igraph_vector_t *el,
                          igraph_vector_t *iindex, 
                          igraph_integer_t nodes) {
  
# define EDGE(i) (VECTOR(*el)[ (long int) VECTOR(*iindex)[(i)] ])
  
  long int no_of_nodes;
  long int no_of_edges;
  long int i, j, idx;
  
  no_of_nodes=nodes;
  no_of_edges=igraph_vector_size(el);
  
  /* result */
  
  IGRAPH_CHECK(igraph_vector_resize(res, nodes+1));
  
  /* create the index */

  if (igraph_vector_size(el)==0) {
    /* empty graph */
    igraph_vector_null(res);
  } else {
    idx=-1;
    for (i=0; i<=EDGE(0); i++) {
      idx++; VECTOR(*res)[idx]=0;
    }
    for (i=1; i<no_of_edges; i++) {
      long int n=(long int) (EDGE(i) - 
        EDGE((long int)VECTOR(*res)[idx]));
      for (j=0; j<n; j++) {
	idx++; VECTOR(*res)[idx]=i;
      }
    }
    j=(long int) EDGE((long int)VECTOR(*res)[idx]);
    for (i=0; i<no_of_nodes-j; i++) {
      idx++; VECTOR(*res)[idx]=no_of_edges;
    }
  }

  /* clean */

# undef EDGE
  return 0;
}

#endif /* IGRAPH_TEMPORAL_EDGE_LIST */
