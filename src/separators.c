/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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

#include "igraph_separators.h"
#include "igraph_memory.h"
#include "igraph_adjlist.h"
#include "igraph_dqueue.h"
#include "igraph_vector.h"
#include "igraph_interface.h"

int igraph_i_is_separator(const igraph_t *graph,
			  const igraph_vector_long_t *candidate,
			  igraph_bool_t *res, 
			  igraph_vector_bool_t *removed,
			  igraph_dqueue_t *Q,
			  igraph_vector_t *neis,
			  long int no_of_nodes) {
  
  /* Remove the given vertices from the graph, do a breadth-first
     search and check the number of components */  
  long int i, clen=igraph_vector_long_size(candidate);
  long int start=0;
  
  for (i=0; i<clen; i++) {
    long int v=VECTOR(*candidate)[i];
    VECTOR(*removed)[v] = 1;
  }

  /* Look for the first node that is not removed */
  while (start < no_of_nodes && VECTOR(*removed)[start]) start++;
  
  if (start==no_of_nodes) { 
    IGRAPH_ERROR("All vertices are included in the separator", 
		 IGRAPH_EINVAL);
  }
  
  igraph_dqueue_push(Q, start);
  VECTOR(*removed)[start]=1;
  while (!igraph_dqueue_empty(Q)) {
    long int node=igraph_dqueue_pop(Q);
    long int j, n;
    igraph_neighbors(graph, neis, node, IGRAPH_ALL);
    n=igraph_vector_size(neis);
    for (j=0; j<n; j++) {
      long int nei=VECTOR(*neis)[j];
      if (!VECTOR(*removed)[nei]) {
	IGRAPH_CHECK(igraph_dqueue_push(Q, nei));
	VECTOR(*removed)[nei]=1;
      }
    }
  }
  
  /* Look for the next node that was neighter removed, not visited */
  while (start < no_of_nodes && VECTOR(*removed)[start]) start++;
  
  /* If there is another component, then we have a separator */
  *res = (start < no_of_nodes);
    
  return 0;
}

int igraph_is_separator(const igraph_t *graph, 
			const igraph_vector_long_t *candidate,
			igraph_bool_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_bool_t removed;
  igraph_dqueue_t Q;
  igraph_vector_t neis;

  IGRAPH_CHECK(igraph_vector_bool_init(&removed, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &removed);
  IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  IGRAPH_CHECK(igraph_i_is_separator(graph, candidate, res, &removed, 
				     &Q, &neis, no_of_nodes));

  igraph_vector_destroy(&neis);
  igraph_dqueue_destroy(&Q);
  igraph_vector_bool_destroy(&removed);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

int igraph_is_minimal_separator(const igraph_t *graph,
				const igraph_vector_long_t *candidate, 
				igraph_bool_t *res) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_bool_t removed;
  igraph_dqueue_t Q;
  igraph_vector_t neis;
  igraph_vector_long_t candcopy;
  long int candsize=igraph_vector_long_size(candidate);

  IGRAPH_CHECK(igraph_vector_bool_init(&removed, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &removed);
  IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  /* Is it a separator at all? */
  IGRAPH_CHECK(igraph_i_is_separator(graph, candidate, res, &removed, 
				     &Q, &neis, no_of_nodes));
  if (!(*res)) {
    /* Not a separator at all, nothing to do, *res is already set */
  } else if (candsize == 0) {
    /* Nothing to do, minimal, *res is already set */
  } else if (candsize == 1) {
    /* Check whether the graph was connected at all, by checking
       whether the empty set is a separator. If yes, then 'candidate'
       is not minimal, otherwise it is */
    IGRAPH_CHECK(igraph_vector_long_init(&candcopy, 0));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &candcopy);
    igraph_vector_bool_null(&removed);
    IGRAPH_CHECK(igraph_i_is_separator(graph, &candcopy, res, &removed, 
				       &Q, &neis, no_of_nodes));
    (*res) = (*res) ? 0 : 1;	/* opposite */
    igraph_vector_long_destroy(&candcopy);
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    /* General case, we need to remove each vertex from 'candidate'
     * and check whether the remainder is a separator. If this is
     * alse for all vertices, then 'candidate' is a minimal
     * separator.
     * 
     * Trick: We copy 'candidate', and "mask out" vertices one by one,
     * by putting candcopy[0] over it. The only drawback is that we
     * need to handle the removal of candcopy[0] separately.
     */
    long int i;
    IGRAPH_CHECK(igraph_vector_long_copy(&candcopy, candidate));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &candcopy);

    /* "Remove" candcopy[0] */
    VECTOR(candcopy)[0] = VECTOR(candcopy)[1];
    igraph_vector_bool_null(&removed);
    IGRAPH_CHECK(igraph_i_is_separator(graph, &candcopy, res, &removed, 
				       &Q, &neis, no_of_nodes));
    VECTOR(candcopy)[0] = VECTOR(*candidate)[0];

    /* And the rest of them one by one */
    for (i=1; i<candsize && (!*res); i++) {
      VECTOR(candcopy)[i] = VECTOR(*candidate)[0];
      igraph_vector_bool_null(&removed);
      IGRAPH_CHECK(igraph_i_is_separator(graph, &candcopy, res, &removed, 
					 &Q, &neis, no_of_nodes));    
      VECTOR(candcopy)[i] = VECTOR(*candidate)[i];
    }
    (*res) = (*res) ? 0 : 1;	/* opposite */
    igraph_vector_long_destroy(&candcopy);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  igraph_vector_destroy(&neis);
  igraph_dqueue_destroy(&Q);
  igraph_vector_bool_destroy(&removed);
  IGRAPH_FINALLY_CLEAN(3);

  return 0;
}

/* --------------------------------------------------------------------*/

#define UPDATEMARK() do {                              \
    (*mark)++;					       \
    if (!(*mark)) {				       \
      igraph_vector_long_null(leaveout);	       \
      (*mark)=1;				       \
    }                                                  \
  } while (0)

int igraph_i_clusters_leaveout(const igraph_adjlist_t *adjlist, 
			       igraph_vector_long_t *components, 
			       igraph_vector_long_t *leaveout, 
			       unsigned long int *mark,
			       igraph_dqueue_t *Q) {

  /* Another trick: we use the same 'leaveout' vector to mark the
   * vertices that were already found in the BFS 
   */

  long int i, no_of_nodes=igraph_adjlist_size(adjlist);
  
  igraph_dqueue_clear(Q);
  igraph_vector_long_clear(components);

  for (i=0; i<no_of_nodes; i++) {

    if (VECTOR(*leaveout)[i] == *mark) continue;

    VECTOR(*leaveout)[i]= *mark;
    igraph_dqueue_push(Q, i);
    igraph_vector_long_push_back(components, i);
    
    while (!igraph_dqueue_empty(Q)) {
      long int act_node=igraph_dqueue_pop(Q);
      igraph_vector_t *neis=igraph_adjlist_get(adjlist, act_node);
      long int j, n=igraph_vector_size(neis);
      for (j=0; j<n; j++) {
	long int nei=VECTOR(*neis)[j];
	if (VECTOR(*leaveout)[nei]== *mark) continue;
	IGRAPH_CHECK(igraph_dqueue_push(Q, nei));
	VECTOR(*leaveout)[nei]= *mark;
	igraph_vector_long_push_back(components, nei);
      }
    }
    
    igraph_vector_long_push_back(components, -1);
  }
  
  UPDATEMARK();
  
  return 0;
}

igraph_bool_t igraph_i_separators_newsep(const igraph_vector_ptr_t *comps,
					 const igraph_vector_long_t *newc) {

  long int co, nocomps=igraph_vector_ptr_size(comps);
  
  for (co=0; co<nocomps; co++) {
    igraph_vector_long_t *act=VECTOR(*comps)[co];
    if (igraph_vector_long_is_equal(act, newc)) return 0;
  }

  /* If not found, then it is new */
  return 1;
}

int igraph_i_separators_store(igraph_vector_ptr_t *separators, 
			      const igraph_adjlist_t *adjlist,
			      igraph_vector_long_t *components, 
			      igraph_vector_long_t *leaveout, 
			      unsigned long int *mark, 
			      igraph_vector_long_t *sorter) {
  
  /* We need to stote N(C), the neighborhood of C, but only if it is 
   * not already stored among the separators.
   */
  
  long int cptr=0, next, complen=igraph_vector_long_size(components);

  while (cptr < complen) {
    long int saved=cptr;
    igraph_vector_long_clear(sorter);

    /* Calculate N(C) for the next C */

    while ( (next=VECTOR(*components)[cptr++]) != -1) {
      VECTOR(*leaveout)[next] = *mark;
    }
    cptr=saved;

    while ( (next=VECTOR(*components)[cptr++]) != -1) {
      igraph_vector_t *neis=igraph_adjlist_get(adjlist, next);
      long int j, nn=igraph_vector_size(neis);
      for (j=0; j<nn; j++) {
	long int nei=VECTOR(*neis)[j];
	if (VECTOR(*leaveout)[nei] != *mark) {
	  igraph_vector_long_push_back(sorter, nei);
	  VECTOR(*leaveout)[nei] = *mark;
	}
      }    
    }
    igraph_vector_long_sort(sorter);

    UPDATEMARK();

    /* Add it to the list of separators, if it is new */

    if (igraph_i_separators_newsep(separators, sorter)) {
      igraph_vector_long_t *newc=igraph_Calloc(1, igraph_vector_long_t);
      if (!newc) {
	IGRAPH_ERROR("Cannot calculate minimal separators", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, newc);
      igraph_vector_long_copy(newc, sorter);
      IGRAPH_FINALLY(igraph_vector_long_destroy, newc);
      IGRAPH_CHECK(igraph_vector_ptr_push_back(separators, newc));
      IGRAPH_FINALLY_CLEAN(2);      
    }
  } /* while cptr < complen */

  return 0;
}

void igraph_i_separators_free(igraph_vector_ptr_t *separators) {
  long int i, n=igraph_vector_ptr_size(separators);
  for (i=0; i<n; i++) {
    igraph_vector_t *vec=VECTOR(*separators)[i];
    if (vec) {   
      igraph_vector_destroy(vec);
      igraph_Free(vec);
    }
  }
}

int igraph_minimal_separators_berry(const igraph_t *graph, 
				    igraph_vector_ptr_t *separators) {

  /* 
   * Some notes about the tricks used here. For finding the components
   * of the graph after removing some vertices, we do the
   * following. First we mark the vertices with the actual mark stamp
   * (mark), then run breadth-first search on the graph, but not
   * considering the marked vertices. Then we increase the mark. If
   * there is integer overflow here, then we zero out the mark and set
   * it to one. (We might as well just always zero it out.)
   * 
   * For each separator the vertices are stored in vertex id order. 
   * This facilitates the comparison of the separators when we find a
   * potential new candidate. 
   * 
   * To keep track of which separator we already used as a basis, we
   * keep a boolean vector (already_tried). The try_next pointer show
   * the next separator to try as a basis.
   */

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_long_t leaveout;
  igraph_vector_bool_t already_tried;
  long int try_next=0;
  unsigned long int mark=1;
  long int v;
  
  igraph_adjlist_t adjlist;
  igraph_vector_long_t components;
  igraph_dqueue_t Q;
  igraph_vector_long_t sorter;

  igraph_vector_ptr_clear(separators);
  IGRAPH_FINALLY(igraph_i_separators_free, separators);

  IGRAPH_CHECK(igraph_vector_long_init(&leaveout, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &leaveout);
  IGRAPH_CHECK(igraph_vector_bool_init(&already_tried, 0));
  IGRAPH_FINALLY(igraph_vector_bool_destroy, &already_tried);
  IGRAPH_CHECK(igraph_vector_long_init(&components, 0));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &components);
  IGRAPH_CHECK(igraph_vector_long_reserve(&components, no_of_nodes*2));
  IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL));
  IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
  IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);
  IGRAPH_CHECK(igraph_vector_long_init(&sorter, 0));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &sorter);
  IGRAPH_CHECK(igraph_vector_long_reserve(&sorter, no_of_nodes));
  
  /* --------------------------------------------------------------- 
   * INITIALIZATION, we check whether the neighborhoods of the 
   * vertices separate the graph. The ones that do will form the 
   * initial basis.
   */
  
  for (v=0; v<no_of_nodes; v++) {

    /* Mark v and its neighbors */
    igraph_vector_t *neis=igraph_adjlist_get(&adjlist, v);
    long int i, n=igraph_vector_size(neis);
    VECTOR(leaveout)[v]=mark;
    for (i=0; i<n; i++) {
      long int nei=VECTOR(*neis)[i];
      VECTOR(leaveout)[nei]=mark;
    }

    /* Find the components */
    IGRAPH_CHECK(igraph_i_clusters_leaveout(&adjlist, &components, &leaveout, 
					    &mark, &Q));

    /* Store the corresponding separators, N(C) for each component C */
    IGRAPH_CHECK(igraph_i_separators_store(separators, &adjlist, &components, 
					   &leaveout, &mark, &sorter));
  }

  /* ---------------------------------------------------------------
   * GENERATION, we need to use all already found separators as
   * basis and see if they generate more separators
   */
  
  while (try_next < igraph_vector_ptr_size(separators)) {
    igraph_vector_long_t *basis=VECTOR(*separators)[try_next];
    long int x, basislen=igraph_vector_long_size(basis);
    for (x=0; x<basislen; x++) {

      /* Remove N(x) U basis */
      igraph_vector_t *neis=igraph_adjlist_get(&adjlist, x);
      long int i, n=igraph_vector_size(neis);
      for (i=0; i<basislen; i++) {
	long int sn=VECTOR(*basis)[i];
	VECTOR(leaveout)[sn]=mark;
      }
      for (i=0; i<n; i++) {
	long int nei=VECTOR(*neis)[i];
	VECTOR(leaveout)[nei]=mark;
      }
      
      /* Find the components */
      IGRAPH_CHECK(igraph_i_clusters_leaveout(&adjlist, &components, 
					      &leaveout, &mark, &Q));
      
      /* Store the corresponding separators, N(C) for each component C */
      IGRAPH_CHECK(igraph_i_separators_store(separators, &adjlist, 
					     &components, &leaveout, &mark, 
					     &sorter));
    }
    
    try_next++;
  }
  
  /* --------------------------------------------------------------- */

  igraph_vector_long_destroy(&sorter);
  igraph_dqueue_destroy(&Q);
  igraph_adjlist_destroy(&adjlist);
  igraph_vector_long_destroy(&components);
  igraph_vector_bool_destroy(&already_tried);
  igraph_vector_long_destroy(&leaveout);
  IGRAPH_FINALLY_CLEAN(7);	/* +1 for separators */

  return 0;
}

#undef UPDATEMARK

/**
 * \function igraph_minimal_separators
 * Find all minimal separators in a graph
 * 
 * 
 */

int igraph_minimal_separators(const igraph_t *graph, 
			      igraph_vector_ptr_t *separators, 
			      igraph_separator_algorithm_t algo) {

  if (algo == IGRAPH_SEPARATOR_ALGORITHM_BERRY) {
    return igraph_minimal_separators_berry(graph, separators);
  } else {
    IGRAPH_ERROR("Unknown minimal separator algorithm", IGRAPH_EINVAL);
  }
  
  return IGRAPH_FAILURE;	/* never happens */
}
