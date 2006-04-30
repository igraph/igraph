/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
#include "random.h"

#include <string.h>
#include <stdarg.h>

int igraph_vs_all(igraph_vs_t *vs) {
  vs->type=IGRAPH_VS_ALL;
  return 0;
}

igraph_vs_t igraph_vss_all() {
  igraph_vs_t allvs;
  allvs.type=IGRAPH_VS_ALL;
  return allvs;  
}

int igraph_vs_adj(igraph_vs_t *vs, 
		  igraph_integer_t vid, igraph_neimode_t mode) {
  vs->type=IGRAPH_VS_ADJ;
  vs->data.adj.vid=vid;
  vs->data.adj.mode=mode;
  return 0;
}

int igraph_vs_none(igraph_vs_t *vs) {
  vs->type=IGRAPH_VS_NONE;
  return 0;
}

igraph_vs_t igraph_vss_none() {
  igraph_vs_t nonevs;
  nonevs.type=IGRAPH_VS_NONE;
  return nonevs;
}

int igraph_vs_1(igraph_vs_t *vs, igraph_integer_t vid) {
  vs->type=IGRAPH_VS_1;
  vs->data.vid=vid;
  return 0;
}

igraph_vs_t igraph_vss_1(igraph_integer_t vid) {
  igraph_vs_t onevs;
  onevs.type=IGRAPH_VS_1;
  onevs.data.vid=vid;
  return onevs;
}

int igraph_vs_vector(igraph_vs_t *vs,
		     const igraph_vector_t *v) {
  vs->type=IGRAPH_VS_VECTORPTR;
  vs->data.vecptr=v;
  return 0;
}

igraph_vs_t igraph_vss_vector(const igraph_vector_t *v) {
  igraph_vs_t vecvs;
  vecvs.type=IGRAPH_VS_VECTORPTR;
  vecvs.data.vecptr=v;
  return vecvs;
}

int igraph_vs_vector_small(igraph_vs_t *vs, ...) {
  va_list ap;
  long int i, n=0;
  vs->type=IGRAPH_VS_VECTOR;
  vs->data.vecptr=Calloc(1, igraph_vector_t);
  if (vs->data.vecptr==0) {
    IGRAPH_ERROR("Cannot create vertex selector", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, (igraph_vector_t*)vs->data.vecptr);
  
  va_start(ap, vs);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)vs->data.vecptr, n);
  
  va_start(ap, vs);
  for (i=0; i<n; i++) {
    VECTOR(*vs->data.vecptr)[i]=(igraph_real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;  
}

void igraph_vs_destroy(igraph_vs_t *vs) {
  switch (vs->type) {
  case IGRAPH_VS_ALL:
  case IGRAPH_VS_ADJ:
  case IGRAPH_VS_NONE:
  case IGRAPH_VS_1:
  case IGRAPH_VS_VECTORPTR:
    break;
  case IGRAPH_VS_VECTOR:
    igraph_vector_destroy((igraph_vector_t*)vs->data.vecptr);
    Free(vs->data.vecptr);
    break;
  default:
    break;
  }
}

/***************************************************/

int igraph_vit_create(const igraph_t *graph, 
		      igraph_vs_t vs, igraph_vit_t *vit) {
  switch (vs.type) {
  case IGRAPH_VS_ALL:
    vit->type=IGRAPH_VIT_SEQ;
    vit->pos=0;
    vit->start=0;
    vit->end=igraph_vcount(graph);
    break;
  case IGRAPH_VS_ADJ:
    vit->type=IGRAPH_VIT_VECTOR;
    vit->pos=0;
    vit->start=0;
    vit->vec=Calloc(1, igraph_vector_t);
    if (vit->vec != 0) {
      IGRAPH_ERROR("Cannot create iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) vit->vec);
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)vit->vec, 0);
    IGRAPH_CHECK(igraph_neighbors(graph, (igraph_vector_t*)vit->vec, 
				  vs.data.adj.vid, vs.data.adj.mode));
    vit->end=igraph_vector_size(vit->vec);
    IGRAPH_FINALLY_CLEAN(2);
    break;
  case IGRAPH_VS_NONE:
    vit->type=IGRAPH_VIT_SEQ;
    vit->pos=0;
    vit->start=0;
    vit->end=0;
    break;
  case IGRAPH_VS_1:
    vit->type=IGRAPH_VIT_SEQ;
    vit->pos=vs.data.vid;
    vit->start=vs.data.vid;
    vit->end=vs.data.vid+1;
    if (vit->pos >= igraph_vcount(graph)) {
      IGRAPH_ERROR("Cannot create iterator, invalid vertex id",IGRAPH_EINVVID);
    }
    break;
  case IGRAPH_VS_VECTORPTR:
  case IGRAPH_VS_VECTOR:
    vit->type=IGRAPH_VIT_VECTORPTR;
    vit->pos=0;
    vit->start=0;
    vit->vec=vs.data.vecptr;
    vit->end=igraph_vector_size(vit->vec);
    if (!igraph_vector_isininterval(vit->vec, 0, igraph_vcount(graph)-1)) {
      IGRAPH_ERROR("Cannot create iterator, invalid vertex id",IGRAPH_EINVVID);
    }
    break;
  default:
    IGRAPH_ERROR("Cannot create iterator, invalid selector", IGRAPH_EINVAL);
    break;
  }
  return 0;
}

void igraph_vit_destroy(const igraph_vit_t *vit) {
  switch (vit->type) {
  case IGRAPH_VIT_SEQ:
  case IGRAPH_VIT_VECTORPTR:
    break;
  case IGRAPH_VIT_VECTOR:
    igraph_vector_destroy((igraph_vector_t*)vit->vec);
    igraph_free((igraph_vector_t*)vit->vec);
    break;
  default:
/*     IGRAPH_ERROR("Cannot destroy iterator, unknown type", IGRAPH_EINVAL); */
    break;
  }
}

/*******************************************************/

int igraph_es_all(igraph_es_t *es, 
		  igraph_edgeorder_type_t order) {
  switch (order) {
  case IGRAPH_EDGEORDER_ID:
    es->type=IGRAPH_ES_ALL;
    break;
  case IGRAPH_EDGEORDER_FROM:
    es->type=IGRAPH_ES_ALLFROM;
    break;
  case IGRAPH_EDGEORDER_TO:
    es->type=IGRAPH_ES_ALLTO;
    break;
  default:
    IGRAPH_ERROR("Invalid edge order, cannot create selector", IGRAPH_EINVAL);
    break;
  }
  return 0;
}

igraph_es_t igraph_ess_all(igraph_edgeorder_type_t order) {
  igraph_es_t es;
  igraph_es_all(&es, order); /* cannot fail */
  return es;  
}

int igraph_es_adj(igraph_es_t *es, 
		  igraph_integer_t vid, igraph_neimode_t mode) {
  es->type=IGRAPH_ES_ADJ;
  es->data.adj.vid=vid;
  es->data.adj.mode=mode;
  return 0;
}

int igraph_es_none(igraph_es_t *es) {
  es->type=IGRAPH_ES_NONE;
  return 0;
}

igraph_es_t igraph_ess_none() {
  igraph_es_t es;
  es.type=IGRAPH_ES_NONE;
  return es;
}

int igraph_es_1(igraph_es_t *es, igraph_integer_t eid) {
  es->type=IGRAPH_ES_1;
  es->data.eid=eid;
  return 0;
}

igraph_es_t igraph_ess_1(igraph_integer_t eid) {
  igraph_es_t es;
  es.type=IGRAPH_ES_NONE;
  es.data.eid=eid;
  return es;
}

int igraph_es_vector(igraph_es_t *es,
		     const igraph_vector_t *v) {
  es->type=IGRAPH_ES_VECTORPTR;
  es->data.vecptr=v;
  return 0;
}

igraph_es_t igraph_ess_vector(const igraph_vector_t *v) {
  igraph_es_t es;
  es.type=IGRAPH_ES_VECTORPTR;
  es.data.vecptr=v;
  return es;
}

int igraph_es_fromto(igraph_es_t *es,
		     igraph_vs_t from, igraph_vs_t to) {
  /* TODO */
}

void igraph_es_destroy(igraph_es_t *es) {
  switch (es->type) { 
  case IGRAPH_ES_ALL:
  case IGRAPH_ES_ALLFROM:
  case IGRAPH_ES_ALLTO:
  case IGRAPH_ES_ADJ:
  case IGRAPH_ES_NONE:
  case IGRAPH_ES_1:
  case IGRAPH_ES_VECTORPTR:
    break;
  default:
    break;
  }
}

/**************************************************/

int igraph_i_eit_create_allfromto(const igraph_t *graph,
				  igraph_es_t es, igraph_eit_t *eit, 
				  igraph_neimode_t mode) {
  igraph_vector_t *vec;
  long int no_of_nodes=igraph_vcount(graph);
  long int i;
  
  vec=Calloc(1, igraph_vector_t);
  if (vec==0) {
    IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, vec);
  IGRAPH_VECTOR_INIT_FINALLY(vec, 0);
  IGRAPH_CHECK(igraph_vector_reserve(vec, igraph_ecount(graph)));
  
  if (igraph_is_directed(graph)) {
    igraph_vector_t adj;
    IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
    for (i=0; i<no_of_nodes; i++) {
      igraph_adjacent(graph, &adj, i, mode);
      igraph_vector_append(vec, &adj);
    }
    igraph_vector_destroy(&adj);
    IGRAPH_FINALLY_CLEAN(1);

  } else {

    igraph_vector_t adj;
    igraph_bool_t *added;
    long int j;
    IGRAPH_VECTOR_INIT_FINALLY(&adj, 0);
    added=Calloc(igraph_ecount(graph), igraph_bool_t);
    if (added==0) {
      IGRAPH_ERROR("Cannot create edge iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, added);      
    for (i=0; i<no_of_nodes; i++) {
      igraph_adjacent(graph, &adj, i, IGRAPH_ALL);
      for (j=0; j<igraph_vector_size(&adj); j++) {
	if (!added[ (long int)VECTOR(adj)[j] ]) {
	  igraph_vector_push_back(vec, VECTOR(adj)[j]);
	  added[ (long int)VECTOR(adj)[j] ]+=1;
	}
      }
    }
    igraph_vector_destroy(&adj);
    Free(added);
    IGRAPH_FINALLY_CLEAN(2);
  }

  eit->type=IGRAPH_EIT_VECTORPTR;
  eit->pos=0;
  eit->start=0;
  eit->vec=vec;
  eit->end=igraph_vector_size(eit->vec);

  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_eit_create(const igraph_t *graph, 
		      igraph_es_t es, igraph_eit_t *eit) {
  switch (es.type) {
  case IGRAPH_ES_ALL:
    eit->type=IGRAPH_EIT_SEQ;
    eit->pos=0;
    eit->start=0;
    eit->end=igraph_ecount(graph);
    break;
  case IGRAPH_ES_ALLFROM:
    IGRAPH_CHECK(igraph_i_eit_create_allfromto(graph, es, eit, IGRAPH_OUT));
    break;
  case IGRAPH_ES_ALLTO:
    IGRAPH_CHECK(igraph_i_eit_create_allfromto(graph, es, eit, IGRAPH_IN));
    break;
  case IGRAPH_ES_ADJ:
    eit->type=IGRAPH_EIT_VECTOR;
    eit->pos=0;
    eit->start=0;
    eit->vec=Calloc(1, igraph_vector_t);
    if (eit->vec != 0) {
      IGRAPH_ERROR("Cannot create iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (igraph_vector_t*) eit->vec);
    IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)eit->vec, 0);
    IGRAPH_CHECK(igraph_adjacent(graph, (igraph_vector_t*)eit->vec, 
				 es.data.adj.vid, es.data.adj.mode));
    eit->end=igraph_vector_size(eit->vec);
    IGRAPH_FINALLY_CLEAN(2);
    break;
  case IGRAPH_ES_NONE:
    eit->type=IGRAPH_EIT_SEQ;
    eit->pos=0;
    eit->start=0;
    eit->end=0;
    break;
  case IGRAPH_ES_1:
    eit->type=IGRAPH_EIT_SEQ;
    eit->pos=es.data.eid;
    eit->start=es.data.eid;
    eit->end=es.data.eid+1;
    if (eit->pos >= igraph_ecount(graph)) {
      IGRAPH_ERROR("Cannot create iterator, invalid edge id", IGRAPH_EINVVID);
    }
    break;
  case IGRAPH_ES_VECTORPTR:
    eit->type=IGRAPH_EIT_VECTORPTR;
    eit->pos=0;
    eit->start=0;
    eit->vec=es.data.vecptr;
    eit->end=igraph_vector_size(eit->vec);
    if (!igraph_vector_isininterval(eit->vec, 0, igraph_ecount(graph)-1)) {
      IGRAPH_ERROR("Cannot create iterator, invalid edge id",IGRAPH_EINVVID);
    }
    break;
  default:
    IGRAPH_ERROR("Cannot create iterator, invalid selector", IGRAPH_EINVAL);
    break;
  }
  return 0;
}

void igraph_eit_destroy(const igraph_eit_t *eit) {
  switch (eit->type) {
  case IGRAPH_EIT_SEQ:
  case IGRAPH_EIT_VECTORPTR:
    break;
  case IGRAPH_EIT_VECTOR:
    igraph_vector_destroy((igraph_vector_t*)eit->vec);
    igraph_free((igraph_vector_t*)eit->vec);
    break;
  default:
/*     IGRAPH_ERROR("Cannot destroy iterator, unknown type", IGRAPH_EINVAL); */
    break;
  }
}
