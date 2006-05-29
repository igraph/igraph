/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

int igraph_i_adjlist_init(const igraph_t *graph, igraph_i_adjlist_t *al, 
			  igraph_neimode_t mode) {
  long int i;

  if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_EINVMODE);
  }

  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  al->length=igraph_vcount(graph);
  al->adjs=Calloc(al->length, igraph_vector_t);
  if (al->adjs == 0) {
    IGRAPH_ERROR("Cannot create adjlist view", IGRAPH_ENOMEM);
  }

  IGRAPH_FINALLY(igraph_i_adjlist_destroy, al);
  for (i=0; i<al->length; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_vector_init(&al->adjs[i], 0));
    IGRAPH_CHECK(igraph_neighbors(graph, &al->adjs[i], i, mode));
  }

  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

void igraph_i_adjlist_destroy(igraph_i_adjlist_t *al) {
  long int i;
  for (i=0; i<al->length; i++) {
    igraph_vector_destroy(&al->adjs[i]);
  }
  Free(al->adjs);
}

/* igraph_vector_t *igraph_i_adjlist_get(igraph_i_adjlist_t *al, igraph_integer_t no) { */
/*   return &al->adjs[(long int)no]; */
/* } */

int igraph_i_adjedgelist_init(const igraph_t *graph, 
			      igraph_i_adjedgelist_t *ael, 
			      igraph_neimode_t mode) {
  long int i;

  if (mode != IGRAPH_IN && mode != IGRAPH_OUT && mode != IGRAPH_ALL) {
    IGRAPH_ERROR("Cannot create adjedgelist view", IGRAPH_EINVMODE);
  }

  if (!igraph_is_directed(graph)) { mode=IGRAPH_ALL; }

  ael->length=igraph_vcount(graph);
  ael->adjs=Calloc(ael->length, igraph_vector_t);
  if (ael->adjs == 0) {
    IGRAPH_ERROR("Cannot create adjedgelist view", IGRAPH_ENOMEM);
  }

  IGRAPH_FINALLY(igraph_i_adjlist_destroy, ael);  
  for (i=0; i<ael->length; i++) {
    IGRAPH_ALLOW_INTERRUPTION();
    IGRAPH_CHECK(igraph_vector_init(&ael->adjs[i], 0));
    IGRAPH_CHECK(igraph_adjacent(graph, &ael->adjs[i], i, mode));
  }
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

void igraph_i_adjedgelist_destroy(igraph_i_adjedgelist_t *ael) {
  long int i;
  for (i=0; i<ael->length; i++) {
    /* This works if some igraph_vector_t's are 0, because igraph_vector_destroy can
       handle this. */
    igraph_vector_destroy(&ael->adjs[i]);
  }
  Free(ael->adjs);
}
