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

/**
 * \function igraph_bfs
 * \ingroup internal
 * 
 * Added in version 0.2.
 * 
 * TODO
 */

int igraph_bfs(igraph_t *graph, igraph_integer_t vid, igraph_neimode_t mode,
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
