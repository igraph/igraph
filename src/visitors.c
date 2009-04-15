/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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

#include "igraph.h"
#include "memory.h"
#include "config.h"

int igraph_bfs(const igraph_t *graph, 
	       igraph_integer_t root, igraph_neimode_t mode,
	       igraph_vector_t *previsit, igraph_vector_t *postvisit,
	       igraph_vector_t *prerank, igraph_vector_t *postrank,
	       igraph_vector_t *pred, igraph_vector_t *succ,
	       igraph_vector_t *dist, igraph_bfshandler_t *callback,
	       void *extra) {
  
  igraph_dqueue_t Q;
  long int no_of_nodes=igraph_vcount(graph);
  long int i, actroot;
  igraph_vector_char_t added;
  igraph_lazy_adjlist_t adjlist;
  
  long int act_prerank=0, act_postrank=0;
  long int pred_vec=-1;
  
  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }
  
  IGRAPH_CHECK(igraph_vector_char_init(&added, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_char_destroy, &added);
  IGRAPH_CHECK(igraph_dqueue_init(&Q, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &Q);

  IGRAPH_CHECK(igraph_lazy_adjlist_init(graph, &adjlist, mode, /*simplify=*/ 0));
  IGRAPH_FINALLY(igraph_lazy_adjlist_destroy, &adjlist);

  /* Resize result vectors */
  if (previsit)  { igraph_vector_resize(previsit,  no_of_nodes); }
  if (postvisit) { igraph_vector_resize(postvisit, no_of_nodes); }
  if (prerank)   { igraph_vector_resize(prerank,   no_of_nodes); }
  if (postrank)  { igraph_vector_resize(postrank,  no_of_nodes); }
  if (pred)      { igraph_vector_resize(pred,      no_of_nodes); }
  if (succ)      { igraph_vector_resize(succ,      no_of_nodes); }
  if (dist)      { igraph_vector_resize(dist,      no_of_nodes); }

  IGRAPH_CHECK(igraph_dqueue_push(&Q, root));
  IGRAPH_CHECK(igraph_dqueue_push(&Q, 0));
  VECTOR(added)[(long int) root] = 1;

  if (prerank) { VECTOR(*prerank) [(long int)root] = act_prerank; }
  if (previsit) { VECTOR(*previsit)[act_prerank++] = root; }

  for (actroot=0; actroot<no_of_nodes; actroot++) {

    /* 'root' first, then all other vertices */
    if (igraph_dqueue_empty(&Q)) {
      if (VECTOR(added)[actroot]) { continue; }
      IGRAPH_CHECK(igraph_dqueue_push(&Q, actroot));
      IGRAPH_CHECK(igraph_dqueue_push(&Q, 0));
      VECTOR(added)[actroot] = 1;

      if (prerank) { VECTOR(*prerank) [(long int)actroot] = act_prerank; }
      if (previsit) { VECTOR(*previsit)[act_prerank++] = actroot; }
    }
      
    while (!igraph_dqueue_empty(&Q)) {
      long int actvect=igraph_dqueue_pop(&Q);
      long int actdist=igraph_dqueue_pop(&Q);
      long int succ_vec;
      igraph_vector_t *neis=igraph_lazy_adjlist_get(&adjlist, actvect);
      long int i, n=igraph_vector_size(neis);    
      
      if (pred) { VECTOR(*pred)[actvect] = pred_vec; }
      if (postrank) { VECTOR(*postrank) [actvect] = act_postrank; }
      if (postvisit) { VECTOR(*postvisit)[act_postrank++] = actvect; }
      if (dist) { VECTOR(*dist)[actvect] = actdist; }
      
      for (i=0; i<n; i++) {
	long int nei=VECTOR(*neis)[i];
	if (! VECTOR(added)[nei]) {
	  VECTOR(added)[nei] = 1;
	  IGRAPH_CHECK(igraph_dqueue_push(&Q, nei));
	  IGRAPH_CHECK(igraph_dqueue_push(&Q, actdist+1));
	}
      }

      succ_vec = igraph_dqueue_empty(&Q) ? -1 : igraph_dqueue_head(&Q);
      if (callback) {
	callback(graph, actvect, pred_vec, succ_vec, act_postrank, actdist, extra);
      }

      if (succ) { VECTOR(*succ)[actvect] = succ_vec; }
      pred_vec=actvect;

    } /* while Q !empty */

  } /* for actroot < no_of_nodes */

  igraph_lazy_adjlist_destroy(&adjlist);
  igraph_dqueue_destroy(&Q);
  igraph_vector_char_destroy(&added);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}

/**
 * \function igraph_i_bfs
 * \ingroup internal
 * 
 * Added in version 0.2.
 * 
 * TODO
 */

int igraph_i_bfs(igraph_t *graph, igraph_integer_t vid, igraph_neimode_t mode,
		 igraph_vector_t *vids, igraph_vector_t *layers,
		 igraph_vector_t *parents) {   

  igraph_dqueue_t q;
  long int vidspos=0;
  igraph_vector_t neis;
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  char *added;
  long int lastlayer=-1;
  
  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  if (mode != IGRAPH_OUT && mode != IGRAPH_IN && 
      mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Invalid mode argument", IGRAPH_EINVMODE);
  }
  
  /* temporary storage */
  added=igraph_Calloc(no_of_nodes, char);
  if (added==0) {
    IGRAPH_ERROR("Cannot calculate BFS", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  IGRAPH_CHECK(igraph_dqueue_init(&q, 100));
  IGRAPH_FINALLY(igraph_dqueue_destroy, &q);

  /* results */
  IGRAPH_CHECK(igraph_vector_resize(vids, no_of_nodes));
  igraph_vector_clear(layers);
  IGRAPH_CHECK(igraph_vector_resize(parents, no_of_nodes));
  
  /* ok start with vid */
  IGRAPH_CHECK(igraph_dqueue_push(&q, vid));
  IGRAPH_CHECK(igraph_dqueue_push(&q, 0));
  IGRAPH_CHECK(igraph_vector_push_back(layers, vidspos)); 
  VECTOR(*vids)[vidspos++]=vid; 
  VECTOR(*parents)[(long int)vid]=vid;
  added[(long int)vid]=1;
  
  while (!igraph_dqueue_empty(&q)) {
    long int actvect=igraph_dqueue_pop(&q);
    long int actdist=igraph_dqueue_pop(&q);
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, actvect, mode));
    for (i=0; i<igraph_vector_size(&neis); i++) {
      long int neighbor=VECTOR(neis)[i];
      if (added[neighbor]==0) {
	added[neighbor]=1;
	VECTOR(*parents)[neighbor]=actvect;
	IGRAPH_CHECK(igraph_dqueue_push(&q, neighbor));
	IGRAPH_CHECK(igraph_dqueue_push(&q, actdist+1));
	if (lastlayer != actdist+1) { 
	  IGRAPH_CHECK(igraph_vector_push_back(layers, vidspos));
	}
	VECTOR(*vids)[vidspos++]=neighbor;
	lastlayer=actdist+1;
      }
    } /* for i in neis */
  } /* while ! dqueue_empty */
  IGRAPH_CHECK(igraph_vector_push_back(layers, vidspos));
  
  igraph_vector_destroy(&neis);
  igraph_dqueue_destroy(&q);
  igraph_Free(added);
  IGRAPH_FINALLY_CLEAN(3);
		 
  return 0;
}
