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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "memory.h"

/* Vertex iterators */
#define IGRAPH_ITERATOR_VID      1

/* Edge iterators */
#define IGRAPH_ITERATOR_EID   1001
#define IGRAPH_ITERATOR_ENEIS 1002

/* Constructors */
int igraph_iterator_vid(igraph_t *graph, igraph_iterator_t *it) {  
  it->type=IGRAPH_ITERATOR_VID;
  it->data=Calloc(1, real_t);
  ((real_t*)it->data)[0]=0;
  it->next=igraph_next_vid;
  it->prev=igraph_prev_vid;
  it->end=igraph_end_vid;
  it->getvertex=igraph_get_vertex_vid;
  it->getvertexfrom=it->getvertexto=it->getedge=0;
  it->getvertexnei=0;
  return 0;
}

int igraph_iterator_eid(igraph_t *graph, igraph_iterator_t *it) {  
  it->type=IGRAPH_ITERATOR_EID;
  it->data=Calloc(1, real_t);
  ((real_t*)it->data)[0]=0;
  it->next=igraph_next_eid;
  it->prev=igraph_prev_eid;
  it->end=igraph_end_eid;
  it->getvertex=0;
  it->getvertexfrom=igraph_get_vertex_from_eid;
  it->getvertexto=igraph_get_vertex_to_eid;
  it->getedge=igraph_get_edge_eid;
  it->getvertexnei=0;
  return 0;
}

typedef struct igraph_iterator_eneis_data_t {
  integer_t vid;
  integer_t mode;
  integer_t oidx;
  integer_t iidx;
} igraph_iterator_eneis_data_t;

int igraph_iterator_eneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, integer_t mode) {  

  igraph_iterator_eneis_data_t *data;
  it->type=IGRAPH_ITERATOR_ENEIS;
  it->next=igraph_next_eneis;
  it->prev=0;
  it->end=igraph_end_eneis;
  it->getvertex=0;
  it->getvertexfrom=igraph_get_vertex_from_eneis;
  it->getvertexto=igraph_get_vertex_to_eneis;
  it->getedge=igraph_get_edge_eneis;
  it->getvertexnei=igraph_get_vertex_nei_eneis;

  it->data=Calloc(1, igraph_iterator_eneis_data_t);
  data=it->data;
  data->vid=vid;
  data->mode=mode;
  if ((int) mode & 1) {
    data->oidx=VECTOR(graph->os)[(long int)vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if ((int) mode & 2) {
    data->iidx=VECTOR(graph->is)[(long int)vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }
  
  return 0;
}
  
/* Destructor */
int igraph_iterator_destroy(igraph_t *graph, igraph_iterator_t *it) {
  switch ((long int) it->type) {
  case IGRAPH_ITERATOR_VID:
    Free(it->data);
    break;
  case IGRAPH_ITERATOR_EID:
    free(it->data); 
    break;
  case IGRAPH_ITERATOR_ENEIS:
    Free(it->data);
    break;
  }
  return 0;
}

/* Common */
int igraph_next(igraph_t *graph, igraph_iterator_t *it) {
  it->next(graph, it);
  return 0;
}

int igraph_prev(igraph_t *graph, igraph_iterator_t *it) {
  it->prev(graph, it);
  return 0;
}

bool_t igraph_end(igraph_t *graph, igraph_iterator_t *it) {
  return it->end(graph, it);
}

integer_t igraph_get_vertex_nei(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertexnei(graph, it);
}

/* Vertices */

integer_t igraph_get_vertex(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertex(graph, it);
}

/* Edges */

integer_t igraph_get_vertex_from(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertexfrom(graph, it);
}

integer_t igraph_get_vertex_to(igraph_t *graph, igraph_iterator_t *it) {
  return it->getvertexto(graph, it);
}

integer_t igraph_get_edge(igraph_t *graph, igraph_iterator_t *it) {
  return it->getedge(graph, it);
}

/* Specifics, simple vertex iterator */
int igraph_next_vid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)(it->data))[0] += 1;
  return 0;
}

int igraph_prev_vid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] -= 1;
  return 0;
}

bool_t igraph_end_vid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0] >= igraph_vcount(graph);
}

integer_t igraph_get_vertex_vid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0];
}

/* Specifics, simple edge iterator */
int igraph_next_eid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] += 1;
  return 0;
}

int igraph_prev_eid(igraph_t *graph, igraph_iterator_t *it) {
  ((real_t*)it->data)[0] -= 1;
  return 0;
}

bool_t igraph_end_eid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0] >= igraph_ecount(graph);
}

integer_t igraph_get_vertex_from_eid(igraph_t *graph, igraph_iterator_t *it) {
  return VECTOR(graph->from)[ (long int) ((real_t*)it->data)[0] ];
}

integer_t igraph_get_vertex_to_eid(igraph_t *graph, igraph_iterator_t *it) {
  return VECTOR(graph->to)[ (long int) ((real_t*)it->data)[0] ];
}

integer_t igraph_get_edge_eid(igraph_t *graph, igraph_iterator_t *it) {
  return ((real_t*)it->data)[0];
}

/* Iterates over the edges to and/or from a vertex */
int igraph_iterator_eneis_set(igraph_t *graph, igraph_iterator_t *it, 
			      integer_t vid, integer_t mode) {
  igraph_iterator_eneis_data_t *data=it->data;
  data->vid=vid;
  data->mode=mode;
  if ((int) mode & 1) {
    data->oidx=VECTOR(graph->os)[(long int)vid];
  } else {
    data->oidx=igraph_ecount(graph);
  }
  if ((int) mode & 2) {
    data->iidx=VECTOR(graph->is)[(long int)vid];
  } else {
    data->iidx=igraph_ecount(graph);
  }
  
  return 0;
}

int igraph_next_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  data->oidx ++; 
  if (data->oidx > VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    data->iidx ++;
  }     
}

bool_t igraph_end_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  return (data->oidx >= VECTOR(graph->os)[ (long int) data->vid+1 ] &&
	  data->iidx >= VECTOR(graph->is)[ (long int) data->vid+1 ]);
}

integer_t igraph_get_vertex_from_eneis(igraph_t *graph, 
				       igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data->oidx];
    return VECTOR(graph->from)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data->iidx];
    return VECTOR(graph->from)[idx];
  }  
}

integer_t igraph_get_vertex_to_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data->oidx];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data->iidx];
    return VECTOR(graph->to)[idx];
  }  
}

integer_t igraph_get_vertex_nei_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)data->oidx];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)data->iidx];
    return VECTOR(graph->from)[idx];
  }  
}  

integer_t igraph_get_edge_eneis(igraph_t *graph, igraph_iterator_t *it) {
  igraph_iterator_eneis_data_t* data=it->data;
  if (data->oidx < VECTOR(graph->os)[ (long int)data->vid+1 ]) {
    return VECTOR(graph->oi)[(long int)data->oidx];
  } else {
    return VECTOR(graph->ii)[(long int)data->iidx];
  }
}

/**********************************************************
 * Testing purposes, indexed edgelist type                *
 *********************************************************/

/* int main() { */
  
/*   igraph_t g; */
/*   vector_t edges; */
/*   igraph_iterator_t it; */

/*   vector_init(&edges, 10); */
/*   VECTOR(edges)[0]=0;  VECTOR(edges)[1]=1; */
/*   VECTOR(edges)[2]=0;  VECTOR(edges)[3]=2; */
/*   VECTOR(edges)[4]=3;  VECTOR(edges)[5]=0; */
/*   VECTOR(edges)[6]=0;  VECTOR(edges)[7]=4; */
/*   VECTOR(edges)[8]=3;  VECTOR(edges)[9]=4; */
/*   igraph_create(&g, &edges, 0, 1); */
  
/*   print_igraph(&g); */

/*   igraph_iterator_eneis(&g, &it, 3, 3); */
/*   while (! igraph_end(&g, &it)) { */
/*     printf("%f: %f -> %f\n", igraph_get_edge(&g, &it),  */
/* 	   igraph_get_vertex_from(&g, &it), */
/* 	   igraph_get_vertex_to(&g, &it)); */
/*     igraph_next(&g, &it); */
/*   } */

/*   printf("----------------------\n"); */
/*   igraph_iterator_eneis_set(&g, &it, 0, 3); */
/*   while (! igraph_end(&g, &it)) { */
/*     printf("%f: %f -> %f\n", igraph_get_edge(&g, &it),  */
/* 	   igraph_get_vertex_from(&g, &it), */
/* 	   igraph_get_vertex_to(&g, &it)); */
/*     igraph_next(&g, &it); */
/*   }   */
  
/*   igraph_iterator_destroy(&g, &it); */
/*   igraph_destroy(&g); */
/*   vector_destroy(&edges); */
  
/*   return 0; */
/* } */
	   
