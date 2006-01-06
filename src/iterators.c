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
#include "random.h"

#include <string.h>
#include <stdarg.h>

/* -------------------------------------------------- */
/* Vertex iterator generics                           */
/* -------------------------------------------------- */

void igraph_vs_next(const igraph_t *graph, igraph_vs_t *vs) {
  vs->table->next(graph, (struct igraph_vs_t*)vs);
}

bool_t igraph_vs_end(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->table->end(graph, (struct igraph_vs_t*)vs);
}

void igraph_vs_reset(const igraph_t *graph, igraph_vs_t *vs) {
  vs->table->reset(graph, (struct igraph_vs_t*)vs);
}

integer_t igraph_vs_get(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->table->get(graph, (struct igraph_vs_t*)vs);
}

int igraph_vs_unfold(const igraph_t *graph, const igraph_vs_t *vs, 
		     vector_t *v) {
  return vs->table->unfold(graph, (const struct igraph_vs_t*)vs, v);
}

void igraph_vs_destroy(igraph_vs_t *vs) {
  vs->table->destroy((struct igraph_vs_t*)vs);
}

/* -------------------------------------------------- */
/* Simple vertex iterator                             */
/* -------------------------------------------------- */

void igraph_vs_next_all(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_all(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_all(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_all(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_all(const igraph_t *graph, const igraph_vs_t *vs, 
			 vector_t *v);
void igraph_vs_destroy_all(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_all_table = {
  igraph_vs_next_all, igraph_vs_end_all, igraph_vs_reset_all,
  igraph_vs_get_all, igraph_vs_unfold_all, igraph_vs_destroy_all
};

int igraph_vs_all(const igraph_t *graph, igraph_vs_t *vs) {
  vs->type=IGRAPH_ITERATOR_VS_ALL;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_all_table;
  vs->shorthand=0;
  return 0;
}

const igraph_vs_t *IGRAPH_VS_ALL(const igraph_t *graph) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_vs_all(graph, vs);
  vs->shorthand=1;
  return vs;
}

void igraph_vs_next_all(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_all(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0] >= igraph_vcount(graph);
}

void igraph_vs_reset_all(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=0;
}

integer_t igraph_vs_get_all(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_all(const igraph_t *graph, const igraph_vs_t *vs, 
			 vector_t *v) {
  long int n;
  vector_t v2;
  n=igraph_vcount(graph);
  IGRAPH_CHECK(vector_init_seq(&v2, 0, n-1));
  vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_vs_destroy_all(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Adjacent vertices of a vertex                      */
/* -------------------------------------------------- */

void igraph_vs_next_adj(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_adj(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_adj(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_adj(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_adj(const igraph_t *graph, const igraph_vs_t *vs, 
			 vector_t *v);
void igraph_vs_destroy_adj(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_adj_table = {
  igraph_vs_next_adj, igraph_vs_end_adj, igraph_vs_reset_adj,
  igraph_vs_get_adj, igraph_vs_unfold_adj, igraph_vs_destroy_adj
};

int igraph_vs_adj(const igraph_t *graph, igraph_vs_t *vs,
		   integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->type=IGRAPH_ITERATOR_VS_ADJ;
  vs->table=&igraph_i_vs_adj_table;
  vs->shorthand=0;
  
  vs->stdata[0]=vid;
  vs->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    vs->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    vs->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    vs->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    vs->stdata[3]=igraph_ecount(graph);
  }
  vs->stdata[4]=VECTOR(graph->os)[ (long int) vs->stdata[0]+1 ];
  vs->stdata[5]=VECTOR(graph->is)[ (long int) vs->stdata[0]+1 ];
  
  return 0;
}

void igraph_vs_next_adj(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[2] ++;
  if (vs->stdata[2] > vs->stdata[4]) {
    vs->stdata[3] ++;
  }
}

bool_t igraph_vs_end_adj(const igraph_t *graph, const igraph_vs_t *vs) {
  return (vs->stdata[2] >= vs->stdata[4] &&
	  vs->stdata[3] >= vs->stdata[5]);
}

void igraph_vs_reset_adj(const igraph_t *graph, igraph_vs_t *vs) {
  if ((int)vs->stdata[1] & IGRAPH_OUT) {
    vs->stdata[2]=VECTOR(graph->os)[(long int)vs->stdata[0]];
  } else {
    vs->stdata[2]=igraph_ecount(graph);
  }
  if ((int)vs->stdata[1] & IGRAPH_IN) {
    vs->stdata[3]=VECTOR(graph->is)[(long int)vs->stdata[0]];
  } else {
    vs->stdata[3]=igraph_ecount(graph);
  }
}

integer_t igraph_vs_get_adj(const igraph_t *graph, const igraph_vs_t *vs) {
  if (vs->stdata[2] < vs->stdata[4]) {
    long int idx=VECTOR(graph->oi)[(long int)vs->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)vs->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

int igraph_vs_unfold_adj(const igraph_t *graph, const igraph_vs_t *vs,
			 vector_t *v) {
  IGRAPH_CHECK(igraph_neighbors(graph, v, vs->stdata[0], vs->stdata[1]));
  return 0;
}

void igraph_vs_adj_set(const igraph_t *graph, igraph_vs_t *vs,
			integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->stdata[0]=vid;
  vs->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    vs->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    vs->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    vs->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    vs->stdata[3]=igraph_ecount(graph);
  }
  vs->stdata[4]=VECTOR(graph->os)[ (long int) vs->stdata[0]+1 ];
  vs->stdata[5]=VECTOR(graph->is)[ (long int) vs->stdata[0]+1 ];
}

void igraph_vs_destroy_adj(igraph_vs_t *vs) {
}  

/* -------------------------------------------------- */
/* Random walker                                      */
/* -------------------------------------------------- */

void igraph_vs_next_rw(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_rw(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_rw(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_rw(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_rw(const igraph_t *graph, const igraph_vs_t *vs, 
			vector_t *v);
void igraph_vs_destroy_rw(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_rw_table = {
  igraph_vs_next_rw, igraph_vs_end_rw, igraph_vs_reset_rw,
  igraph_vs_get_rw, igraph_vs_unfold_rw, igraph_vs_destroy_rw
};

int igraph_vs_rw(const igraph_t *graph, igraph_vs_t *vs,
		 integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->type=IGRAPH_ITERATOR_VS_RW;
  vs->shorthand=0;
  
  vs->stdata[0]=vid;		/* actual vertex */
  vs->stdata[1]=vid;		/* start vertex  */
  vs->stdata[2]=mode;		/* mode */
  vs->stdata[3]=0;		/* number of steps so far */

  vs->table=&igraph_i_vs_rw_table;
  return 0;
}

void igraph_vs_next_rw(const igraph_t *graph, igraph_vs_t *vs) {
  long int vid=vs->stdata[0], nvid;
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
  }
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (indegree+outdegree==0) {
    /* TODO: Nowhere to step, isolate vertex. What should we do? */
    return;
  }
  
  RNG_BEGIN();
  nvid=RNG_INTEGER(0, outdegree+indegree-1);
  RNG_END();

  if (nvid < outdegree) {
    long int i=VECTOR(graph->os)[vid]+nvid;
    nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
  } else {
    long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
    nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
  }
  
  vs->stdata[0]=nvid;
  vs->stdata[3] += 1.0;
}

bool_t igraph_vs_end_rw(const igraph_t *graph, const igraph_vs_t *vs) {
  long int vid=vs->stdata[0];
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  return indegree+outdegree == 0;
}

void igraph_vs_reset_rw(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
  vs->stdata[3]=0;
}

integer_t igraph_vs_get_rw(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_rw(const igraph_t *graph, const igraph_vs_t *vs,
			vector_t *v) {
  IGRAPH_ERROR("attempt to unfold random walker", IGRAPH_EUNFOLDINF);
  return 0;			/* return unneccesary */
}

long int igraph_vs_rw_length(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[3];
}

void igraph_vs_destroy_rw(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Random walker with one unit memory                 */
/* -------------------------------------------------- */

void igraph_vs_next_rw1(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_rw1(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_rw1(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_rw1(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_rw1(const igraph_t *graph, const igraph_vs_t *vs, 
			 vector_t *v);
void igraph_vs_destroy_rw1(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_rw1_table = {
  igraph_vs_next_rw1, igraph_vs_end_rw1, igraph_vs_reset_rw1,
  igraph_vs_get_rw1, igraph_vs_unfold_rw1, igraph_vs_destroy_rw1
};

int igraph_vs_rw1(const igraph_t *graph, igraph_vs_t *vs,
		  integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs->type=IGRAPH_ITERATOR_VS_RW1;
  vs->shorthand=0;
  
  vs->stdata[0]=vid;			/* actual vertex */
  vs->stdata[1]=vid;			/* start vertex  */
  vs->stdata[2]=mode;			/* mode */
  vs->stdata[3]=0;			/* number of steps so far */
  vs->stdata[4]=-1;		        /* the previous vertex */

  vs->table=&igraph_i_vs_rw1_table;
  return 0;
}

void igraph_vs_next_rw1(const igraph_t *graph, igraph_vs_t *vs) {

  long int vid=vs->stdata[0], nvid;
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;

  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->os)[vid]);
  }
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  if (indegree+outdegree==0) {
    /* TODO: Nowhere to step, isolate vertex. What should we do? */
    return;
  }
  
  if (indegree+outdegree==1) {
    /* There only one way, we may go back as well */
    nvid=0;
    
    if (nvid < outdegree) {
      long int i=VECTOR(graph->os)[vid]+nvid;
      nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
    } else {
      long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
      nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
    }
    
  } else {
    /* There are other options */
    RNG_BEGIN();
    if (vs->stdata[4] >=0) {
      nvid=RNG_INTEGER(0, outdegree+indegree-2);
    } else {
      nvid=RNG_INTEGER(0, outdegree+indegree-1);
    }
    RNG_END();

    if (nvid < outdegree) {
      long int i=VECTOR(graph->os)[vid]+nvid;
      nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
    } else {
      long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
      nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
    }
    
    /* In case we wanted to step back */
    if (nvid==vs->stdata[4]) {
      nvid=outdegree+indegree-1;

      if (nvid < outdegree) {
	long int i=VECTOR(graph->os)[vid]+nvid;
	nvid=VECTOR(graph->to)[ (long int) VECTOR(graph->oi)[i] ];
      } else {
	long int i=VECTOR(graph->is)[vid]+nvid-outdegree;
	nvid=VECTOR(graph->from)[ (long int) VECTOR(graph->ii)[i] ];
      }
    }
  }

  vs->stdata[4]=vs->stdata[0];
  vs->stdata[0]=nvid;
  vs->stdata[3]+=1.0;
}

bool_t igraph_vs_end_rw1(const igraph_t *graph, const igraph_vs_t *vs) {
  long int vid=vs->stdata[0];
  igraph_neimode_t mode=vs->stdata[2];
  long int indegree=0, outdegree=0;
  
  if (mode & IGRAPH_OUT) {
    outdegree += (VECTOR(graph->os)[vid+1]-VECTOR(graph->is)[vid]);
  }
  if (mode & IGRAPH_IN) {
    indegree += (VECTOR(graph->is)[vid+1]-VECTOR(graph->is)[vid]);
  }

  return indegree+outdegree == 0;
}

integer_t igraph_vs_get_rw1(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

void igraph_vs_reset_rw1(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
  vs->stdata[3]=0;
  vs->stdata[5]=-1;
}

int igraph_vs_unfold_rw1(const igraph_t *graph, const igraph_vs_t *vs,
			 vector_t *v) {
  IGRAPH_ERROR("attempt to unfold random walker", IGRAPH_EUNFOLDINF);
  return 0;
}

long int igraph_vs_rw1_length(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[3];
}

void igraph_vs_destroy_rw1(igraph_vs_t *pvs) {
  if (pvs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Empty vertex set                                   */
/* -------------------------------------------------- */

void igraph_vs_next_none(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_none(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_none(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_none(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_none(const igraph_t *graph, const igraph_vs_t *vs,
			  vector_t *v);
void igraph_vs_destroy_none(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_none_table = {
  igraph_vs_next_none, igraph_vs_end_none, igraph_vs_reset_none,
  igraph_vs_get_none, igraph_vs_unfold_none, igraph_vs_destroy_none
};

int igraph_vs_none(const igraph_t *graph, igraph_vs_t *vs) {
  vs->type=IGRAPH_ITERATOR_VS_NONE;
  vs->table=&igraph_i_vs_none_table;
  vs->shorthand=0;
  return 0;
}

void igraph_vs_next_none(const igraph_t *graph, igraph_vs_t *vs) {
  /* nothing to do */
}

bool_t igraph_vs_end_none(const igraph_t *graph, const igraph_vs_t *vs) {
  return 1;
}

void igraph_vs_reset_none(const igraph_t *graph, igraph_vs_t *vs) {
  /* nothing to do */
}

integer_t igraph_vs_get_none(const igraph_t *graph, const igraph_vs_t *vs) {
  /* ooops this is an error, no way to signal it though... */
  return -1;
}

int igraph_vs_unfold_none(const igraph_t *graph, const igraph_vs_t *vs,
			  vector_t *v) {
  vector_clear(v);
  return 0;
}

void igraph_vs_destroy_none(igraph_vs_t *vs) {
  if (vs->shorthand) {
    Free(vs);
  }
}

/* -------------------------------------------------- */
/* vertex set with single vertex                      */
/* -------------------------------------------------- */

void igraph_vs_next_1(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_1(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_1(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_1(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_1(const igraph_t *graph, const igraph_vs_t *vs, 
		       vector_t *v);
void igraph_vs_destroy_1(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_1_table = {
  igraph_vs_next_1, igraph_vs_end_1, igraph_vs_reset_1,
  igraph_vs_get_1, igraph_vs_unfold_1, igraph_vs_destroy_1
};


int igraph_vs_1(const igraph_t *igraph, igraph_vs_t *vs, integer_t vid) {
  vs->type=IGRAPH_ITERATOR_VS_1;
  vs->stdata[0]=vid;
  vs->stdata[1]=0;		/* write 1 here if end */
  vs->table=&igraph_i_vs_1_table;
  vs->shorthand=0;
  return 0;
}

const igraph_vs_t *IGRAPH_VS_1(const igraph_t *graph, integer_t vid) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_vs_1(graph, vs, vid);
  vs->shorthand=1;
  return vs;
}

void igraph_vs_next_1(const igraph_t *graph, igraph_vs_t *vs) {
  /* signal end */
  vs->stdata[1]=1;
}

bool_t igraph_vs_end_1(const igraph_t *graph, const igraph_vs_t *vs) {
  return (vs->stdata[1]==1);
}

void igraph_vs_reset_1(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[1]=0;
}

integer_t igraph_vs_get_1(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_1(const igraph_t *graph, const igraph_vs_t *vs,
		       vector_t *v) {
  IGRAPH_CHECK(vector_resize(v, 1));
  VECTOR(*v)[0]=vs->stdata[0];
  return 0;
}

void igraph_vs_destroy_1(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Vertex set with sequence of vertices               */
/* -------------------------------------------------- */

void igraph_vs_next_seq(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_seq(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_seq(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_seq(const igraph_t *graph, const igraph_vs_t *vs);
int igraph_vs_unfold_seq(const igraph_t *graph, const igraph_vs_t *vs,
			 vector_t *v);
void igraph_vs_destroy_seq(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_seq_table = {
  igraph_vs_next_seq, igraph_vs_end_seq, igraph_vs_reset_seq,
  igraph_vs_get_seq, igraph_vs_unfold_seq, igraph_vs_destroy_seq
};

int igraph_vs_seq(const igraph_t *igraph, igraph_vs_t *vs, integer_t from,
		  integer_t to) {
  vs->type=IGRAPH_ITERATOR_VS_SEQ;
  vs->stdata[0]=from;
  vs->stdata[1]=from;
  vs->stdata[2]=to;
  vs->table=&igraph_i_vs_seq_table;
  vs->shorthand=0;
  return 0;
}

void igraph_vs_next_seq(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_seq(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0] > vs->stdata[2];
}

void igraph_vs_reset_seq(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
}

integer_t igraph_vs_get_seq(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_seq(const igraph_t *graph, const igraph_vs_t *vs,
			 vector_t *v) {
  vector_t v2;
  IGRAPH_CHECK(vector_init_seq(&v2, vs->stdata[1], vs->stdata[2]));
  vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_vs_destroy_seq(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->shorthand) {
    Free(pvs);
  }
}

/* -------------------------------------------------- */
/* Vertex ids in a vector,                            */
/*   this can be a view or a copy                     */
/* -------------------------------------------------- */

void igraph_vs_next_vectorview(const igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_vectorview(const igraph_t *graph, const igraph_vs_t *vs);
void igraph_vs_reset_vectorview(const igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_vectorview(const igraph_t *graph, 
				   const igraph_vs_t *vs);
int igraph_vs_unfold_vectorview(const igraph_t *graph, const igraph_vs_t *vs,
				vector_t *v);
void igraph_vs_destroy_vectorview(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_vectorview_table = {
  igraph_vs_next_vectorview, igraph_vs_end_vectorview, 
  igraph_vs_reset_vectorview, igraph_vs_get_vectorview, 
  igraph_vs_unfold_vectorview, igraph_vs_destroy_vectorview
};

typedef struct igraph_i_vs_vectorview_pdata_t {
  const vector_t v;
  bool_t destroy;
} igraph_i_vs_vectorview_pdata_t;

int igraph_vs_vectorview(const igraph_t *igraph, igraph_vs_t *vs, 
			 const vector_t *vids) {
  igraph_i_vs_vectorview_pdata_t *data;
  vector_t *fakev;

  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;

  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(vector_t*) &data->v;
  *fakev = *vids;
  data->destroy=0;

  vs->stdata[1]=vector_size(&data->v);

  return 0;
}

int igraph_vs_vector(const igraph_t *igraph, igraph_vs_t *vs,
		     const vector_t *vids) {
  igraph_i_vs_vectorview_pdata_t *data;
  vector_t *fakev;

  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;
  
  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, vs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(vector_t*)&data->v;
  IGRAPH_CHECK(vector_copy(fakev, vids));
  data->destroy=1;

  vs->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_vs_vectorview_it(const igraph_t *graph, const igraph_vs_t *vs,
			    igraph_vs_t *newvs) {
  igraph_i_vs_vectorview_pdata_t *data, *vsdata;
  vector_t *fakev;

  newvs->type=IGRAPH_ITERATOR_VS_VECTOR;
  newvs->stdata[0]=0;
  newvs->table=&igraph_i_vs_vectorview_table;
  newvs->shorthand=0;
  
  newvs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (newvs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, newvs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)newvs->pdata;
  fakev=(vector_t*)&data->v;
  if (vs->type != IGRAPH_ITERATOR_VS_VECTOR) {
    VECTOR_INIT_FINALLY((vector_t*) &data->v, 0);
    IGRAPH_CHECK(igraph_vs_unfold(graph, vs, fakev));
    data->destroy=1;
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    /* If it is a vector we just create a view */
    vsdata=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;    
    *fakev=vsdata->v;
    data->destroy=0;
  }

  newvs->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_vs_vector_small(const igraph_t *igraph, igraph_vs_t *vs, ...) {
  va_list ap;
  igraph_i_vs_vectorview_pdata_t *data;
  vector_t *fakev;
  long int i, n=0;
  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;
  
  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, vs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, vs);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  VECTOR_INIT_FINALLY(fakev, n);
  
  va_start(ap, vs);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  vs->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}  

const igraph_vs_t *IGRAPH_VS_VECTOR(const igraph_t *graph, 
				    const vector_t *vids) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_vs_vectorview(graph, vs, vids);
  vs->shorthand=1;
  return vs;
}

const igraph_vs_t *IGRAPH_VS(const igraph_t *graph, ...) {
  igraph_vs_t *vs=Calloc(1, igraph_vs_t);
  va_list ap;
  igraph_i_vs_vectorview_pdata_t *data;
  vector_t *fakev;
  long int i, n=0;
  int ret;
  if (vs==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }

  vs->type=IGRAPH_ITERATOR_VS_VECTOR;
  vs->stdata[0]=0;
  vs->table=&igraph_i_vs_vectorview_table;
  vs->shorthand=0;
  
  vs->pdata=Calloc(1, igraph_i_vs_vectorview_pdata_t);
  if (vs->pdata==0) {
    igraph_error("Cannot create vector iterator", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;
  }
  IGRAPH_FINALLY(igraph_free, vs->pdata);
  data=(igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  fakev=(vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, graph);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  ret=vector_init(fakev, n);
  if (ret != 0) {
    igraph_error("Cannot create vector for iterator shorthand", __FILE__, 
		 __LINE__, ret);
    return 0;
  }
  IGRAPH_FINALLY(vector_destroy, fakev);
  
  va_start(ap, graph);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  vs->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  vs->shorthand=1;
  return vs;  
}

void igraph_vs_next_vectorview(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_vectorview(const igraph_t *graph, const igraph_vs_t *vs) {
  return vs->stdata[0] >= vs->stdata[1];
}

void igraph_vs_reset_vectorview(const igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=0;
}

integer_t igraph_vs_get_vectorview(const igraph_t *graph, 
				   const igraph_vs_t *vs) {
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  return VECTOR(data->v)[ (long int) (vs->stdata[0]) ]; 
}

int igraph_vs_unfold_vectorview(const igraph_t *graph, const igraph_vs_t *vs,
				vector_t *v) {
  vector_t v2;
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  IGRAPH_CHECK(vector_copy(&v2, &data->v));
  vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_vs_destroy_vectorview(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  if (data->destroy) {
    vector_destroy((vector_t*)&data->v);
  }
  Free(data);
  
  if (vs->shorthand) {
    Free(pvs);
  }
}

const vector_t *igraph_vs_vector_getvector(const igraph_t *graph, 
					   const igraph_vs_t *vs) {
  igraph_i_vs_vectorview_pdata_t *data=
    (igraph_i_vs_vectorview_pdata_t*)vs->pdata;
  return &data->v;
}


/* -------------------------------------------------- */
/* Edge iterator generics                             */
/* -------------------------------------------------- */

void igraph_es_next(const igraph_t *graph, igraph_es_t *es) {
  es->table->next(graph, es);
}

bool_t igraph_es_end(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->end(graph, es);
}

void igraph_es_reset(const igraph_t *graph, igraph_es_t *es) {
  return es->table->reset(graph, es);
}

integer_t igraph_es_from(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->from(graph, es);
}

integer_t igraph_es_to(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->to(graph, es);
}

integer_t igraph_es_get(const igraph_t *graph, const igraph_es_t *es) {
  return es->table->get(graph, es);
}

int igraph_es_unfold(const igraph_t *graph, const igraph_es_t *es,
		     vector_t *v) {
  return es->table->unfold(graph, es, v);
}

void igraph_es_destroy(igraph_es_t *es) {
  es->table->destroy(es);
}

/* -------------------------------------------------- */
/* Simple edge iterator                               */
/* -------------------------------------------------- */

void igraph_es_next_all(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_all(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_all(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_all(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_all(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_all(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_all(const igraph_t *graph, const igraph_es_t *es,
			 vector_t *v);
void igraph_es_destroy_all(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_all_table = {
  igraph_es_next_all, igraph_es_end_all, igraph_es_reset_all,
  igraph_es_get_all, igraph_es_from_all, igraph_es_to_all,
  igraph_es_unfold_all, igraph_es_destroy_all
};

int igraph_es_all(const igraph_t *graph, igraph_es_t *es) {
  es->type=IGRAPH_ITERATOR_ES_ALL;
  es->stdata[0]=0;
  es->table=&igraph_i_es_all_table;
  es->shorthand=0;
  return 0;
}

const igraph_es_t *IGRAPH_ES_ALL(const igraph_t *graph) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_es_all(graph, es);
  es->shorthand=1;
  return es;
}

void igraph_es_next_all(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_all(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] >= igraph_ecount(graph);
}

void igraph_es_reset_all(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_all(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_all(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->from)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_to_all(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->to)[ (long int) es->stdata[0] ];
}

int igraph_es_unfold_all(const igraph_t *graph, const igraph_es_t *es, 
			 vector_t *v) {
  long int n;
  vector_t v2;
  n=igraph_ecount(graph);
  IGRAPH_CHECK(vector_init_seq(&v2, 0, n-1));
  vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_es_destroy_all(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* Ordered edge iterator                              */
/* -------------------------------------------------- */

void igraph_es_next_fromorder(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_fromorder(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_fromorder(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_fromorder(const igraph_t *graph, 
				  const igraph_es_t *es);
integer_t igraph_es_from_fromorder(const igraph_t *graph, 
				   const igraph_es_t *es);
integer_t igraph_es_to_fromorder(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_fromorder(const igraph_t *graph, const igraph_es_t *es,
			       vector_t *v);
void igraph_es_destroy_fromorder(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_fromorder_table = {
  igraph_es_next_fromorder, igraph_es_end_fromorder, igraph_es_reset_fromorder,
  igraph_es_get_fromorder, igraph_es_from_fromorder, igraph_es_to_fromorder,
  igraph_es_unfold_fromorder, igraph_es_destroy_fromorder
};
  
int igraph_es_fromorder(const igraph_t *graph, igraph_es_t *es) {
  es->type=IGRAPH_ITERATOR_ES_FROMORDER;
  es->stdata[0]=0;
  es->table=&igraph_i_es_fromorder_table;
  es->shorthand=0;
  return 0;
}

void igraph_es_next_fromorder(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_fromorder(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] >= igraph_ecount(graph);
}

void igraph_es_reset_fromorder(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_fromorder(const igraph_t *graph, 
				  const igraph_es_t *es) {
  return VECTOR(graph->oi)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_from_fromorder(const igraph_t *graph, 
				   const igraph_es_t *es) {
  long int idx=VECTOR(graph->oi)[ (long int) es->stdata[0] ];
  return VECTOR(graph->from)[ idx ];
}

integer_t igraph_es_to_fromorder(const igraph_t *graph, 
				 const igraph_es_t *es) {
  long int idx=VECTOR(graph->oi)[ (long int) es->stdata[0] ];
  return VECTOR(graph->to)[ idx ];
}

int igraph_es_unfold_fromorder(const igraph_t *graph, const igraph_es_t *es,
			       vector_t *v) {
  long int i=0, n=igraph_ecount(graph);  
  IGRAPH_CHECK(vector_resize(v, n));
  for (i=0; i<n; i++) {
    VECTOR(*v)[i]=VECTOR(graph->oi)[i];
  }
  return 0;
}

void igraph_es_destroy_fromorder(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* Adjacent edges of a vertex                         */
/* -------------------------------------------------- */

void igraph_es_next_adj(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_adj(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_adj(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_adj(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_adj(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_adj(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_adj(const igraph_t *graph, const igraph_es_t *es,
			 vector_t *v);
void igraph_es_destroy_adj(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_adj_table = {
  igraph_es_next_adj, igraph_es_end_adj, igraph_es_reset_adj,
  igraph_es_get_adj, igraph_es_from_adj, igraph_es_to_adj,
  igraph_es_unfold_adj, igraph_es_destroy_adj
};
  
int igraph_es_adj(const igraph_t *graph, igraph_es_t *es,
		  integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  
  es->type=IGRAPH_ITERATOR_ES_ADJ;
  es->shorthand=0;

  es->stdata[0]=vid;
  es->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    es->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    es->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    es->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    es->stdata[3]=igraph_ecount(graph);
  }

  es->table=&igraph_i_es_adj_table;
  return 0;
}

void igraph_es_next_adj(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[2] ++;
  if (es->stdata[2] > VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    es->stdata[3] ++;
  }
}

bool_t igraph_es_end_adj(const igraph_t *graph, const igraph_es_t *es) {
  return (es->stdata[2] >= VECTOR(graph->os)[ (long int) es->stdata[0]+1 ] &&
	  es->stdata[3] >= VECTOR(graph->is)[ (long int) es->stdata[0]+1 ]);
}

void igraph_es_reset_adj(const igraph_t *graph, igraph_es_t *es) {
  if ((int)es->stdata[1] & IGRAPH_OUT) {
    es->stdata[2]=VECTOR(graph->os)[(long int)es->stdata[0]];
  } else {
    es->stdata[2]=igraph_ecount(graph);
  }
  if ((int)es->stdata[1] & IGRAPH_IN) {
    es->stdata[3]=VECTOR(graph->is)[(long int)es->stdata[0]];
  } else {
    es->stdata[3]=igraph_ecount(graph);
  }
}

integer_t igraph_es_get_adj(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    return VECTOR(graph->oi)[(long int)es->stdata[2]];
  } else {
    return VECTOR(graph->ii)[(long int)es->stdata[3]];
  }
}

integer_t igraph_es_from_adj(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->from)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

integer_t igraph_es_to_adj(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->to)[idx];
  }
}

int igraph_es_unfold_adj(const igraph_t *graph, const igraph_es_t *es,
			 vector_t *v) {
  long int length=0, idx=0;
  long int no_of_edges;
  long int i;
  
  long int node=es->stdata[0];
  int mode=es->stdata[1];
  
  no_of_edges=igraph_ecount(graph);

  /* Calculate needed space first & allocate it*/
  
  if (mode & IGRAPH_OUT) {
    length += (VECTOR(graph->os)[node+1] - VECTOR(graph->os)[node]);
  }
  if (mode & IGRAPH_IN) {
    length += (VECTOR(graph->is)[node+1] - VECTOR(graph->is)[node]);
  }
  
  IGRAPH_CHECK(vector_resize(v, length));
  
  if (mode & IGRAPH_OUT) {
    for (i=VECTOR(graph->os)[node]; i<VECTOR(graph->os)[node+1]; i++) {
      VECTOR(*v)[idx++] = VECTOR(graph->oi)[i];
    }
  }
  if (mode & IGRAPH_IN) {
    for (i=VECTOR(graph->is)[node]; i<VECTOR(graph->is)[node+1]; i++) {
      VECTOR(*v)[idx++] = VECTOR(graph->ii)[i];
    }
  }

  return 0;
}

void igraph_es_destroy_adj(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

void igraph_es_adj_set(const igraph_t *graph, igraph_es_t *es,
		       integer_t vid, igraph_neimode_t mode) {
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  es->stdata[0]=vid;
  es->stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    es->stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    es->stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    es->stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    es->stdata[3]=igraph_ecount(graph);
  }
}

integer_t igraph_es_adj_vertex(const igraph_t *graph, const igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

/* -------------------------------------------------- */
/* Empty edge iterator                                */
/* -------------------------------------------------- */

void igraph_es_next_none(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_none(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_none(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_none(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_none(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_none(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_none(const igraph_t *graph, const igraph_es_t *es, 
			  vector_t *v);
void igraph_es_destroy_none(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_none_table = {
  igraph_es_next_none, igraph_es_end_none, igraph_es_reset_none,
  igraph_es_get_none, igraph_es_from_none, igraph_es_to_none,
  igraph_es_unfold_none, igraph_es_destroy_none
};

int igraph_es_none(const igraph_t *graph, igraph_es_t *es) {
  es->type=IGRAPH_ITERATOR_ES_NONE;
  es->shorthand=0;  
  es->table=&igraph_i_es_none_table;
  return 0;
}

void igraph_es_next_none(const igraph_t *graph, igraph_es_t *es) {
  /* nothing to do */
}

bool_t igraph_es_end_none(const igraph_t *graph, const igraph_es_t *es) {
  return 1;
}

void igraph_es_reset_none(const igraph_t *graph, igraph_es_t *es) {
  /* nothing to do */
}

integer_t igraph_es_get_none(const igraph_t *graph, const igraph_es_t *es) {
  /* ooops this is an error, no way to signal it though... */
  return -1;
}

integer_t igraph_es_from_none(const igraph_t *graph, const igraph_es_t *es) {
  /* error */
  return 0;
}

integer_t igraph_es_to_none(const igraph_t *graph, const igraph_es_t *es) {
  /* error */
  return 0;
}

int igraph_es_unfold_none(const igraph_t *graph, const igraph_es_t *es,
			  vector_t *v) {
  vector_clear(v);
  return 0;
}

void igraph_es_destroy_none(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* edge iterator, single edge                         */
/* -------------------------------------------------- */

void igraph_es_next_1(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_1(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_1(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_1(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_1(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_1(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_1(const igraph_t *graph, const igraph_es_t *es, 
		       vector_t *v);
void igraph_es_destroy_1(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_1_table = {
  igraph_es_next_1, igraph_es_end_1, igraph_es_reset_1,
  igraph_es_get_1, igraph_es_from_1, igraph_es_to_1,
  igraph_es_unfold_1, igraph_es_destroy_1
};

int igraph_es_1(const igraph_t *igraph, igraph_es_t *es, integer_t eid) {
  es->type=IGRAPH_ITERATOR_ES_1;
  es->stdata[0]=eid;
  es->stdata[1]=0;		/* write 1 here if end */
  es->table=&igraph_i_es_1_table;
  es->shorthand=0;
  return 0;
}

const igraph_es_t *IGRAPH_ES_1(const igraph_t *graph, integer_t eid) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_es_1(graph, es, eid);
  es->shorthand=1;
  return es;
}

void igraph_es_next_1(const igraph_t *graph, igraph_es_t *es) {
  /* signal end */
  es->stdata[1]=1;
}

bool_t igraph_es_end_1(const igraph_t *graph, const igraph_es_t *es) {
  return (es->stdata[1]==1);
}

void igraph_es_reset_1(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[1]=0;
}

integer_t igraph_es_get_1(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_1(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->from)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_to_1(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->to)[ (long int) es->stdata[0] ];
}

int igraph_es_unfold_1(const igraph_t *graph, const igraph_es_t *es,
		       vector_t *v) {
  IGRAPH_CHECK(vector_resize(v, 1));
  VECTOR(*v)[0]=es->stdata[0];
  return 0;
}

void igraph_es_destroy_1(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->shorthand) {
    Free(es);
  }
}

/* -------------------------------------------------- */
/* Edge set with sequence of vertices               */
/* -------------------------------------------------- */

void igraph_es_next_seq(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_seq(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_seq(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_seq(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_from_seq(const igraph_t *graph, const igraph_es_t *es);
integer_t igraph_es_to_seq(const igraph_t *graph, const igraph_es_t *es);
int igraph_es_unfold_seq(const igraph_t *graph, const igraph_es_t *es,
			 vector_t *v);
void igraph_es_destroy_seq(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_seq_table = {
  igraph_es_next_seq, igraph_es_end_seq, igraph_es_reset_seq,
  igraph_es_get_seq, igraph_es_from_seq, igraph_es_to_seq,
  igraph_es_unfold_seq, igraph_es_destroy_seq
};

int igraph_es_seq(const igraph_t *igraph, igraph_es_t *es, integer_t from,
		  integer_t to) {
  es->type=IGRAPH_ITERATOR_ES_SEQ;
  es->stdata[0]=from;
  es->stdata[1]=from;
  es->stdata[2]=to;
  es->table=&igraph_i_es_seq_table;
  es->shorthand=0;
  return 0;
}

void igraph_es_next_seq(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_seq(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] > es->stdata[2];
}

void igraph_es_reset_seq(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=es->stdata[1];
}

integer_t igraph_es_get_seq(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_seq(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->from)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_to_seq(const igraph_t *graph, const igraph_es_t *es) {
  return VECTOR(graph->to)[ (long int) es->stdata[0] ];
}

int igraph_es_unfold_seq(const igraph_t *graph, const igraph_es_t *es,
			 vector_t *v) {
  vector_t v2;
  IGRAPH_CHECK(vector_init_seq(&v2, es->stdata[1], es->stdata[2]));
  vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_es_destroy_seq(igraph_es_t *pes) {
  igraph_es_t *es=(igraph_es_t*)pes;
  if (es->shorthand) {
    Free(pes);
  }
}

/* -------------------------------------------------- */
/* Edge ids in a vector, this is a view               */
/* -------------------------------------------------- */

void igraph_es_next_vectorview(const igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_vectorview(const igraph_t *graph, const igraph_es_t *es);
void igraph_es_reset_vectorview(const igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_vectorview(const igraph_t *graph, 
				   const igraph_es_t *es);
integer_t igraph_es_from_vectorview(const igraph_t *graph, 
				    const igraph_es_t *es);
integer_t igraph_es_to_vectorview(const igraph_t *graph, 
				  const igraph_es_t *es);
int igraph_es_unfold_vectorview(const igraph_t *graph, const igraph_es_t *es,
				vector_t *v);
void igraph_es_destroy_vectorview(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_vectorview_table = {
  igraph_es_next_vectorview, igraph_es_end_vectorview, 
  igraph_es_reset_vectorview, igraph_es_get_vectorview, 
  igraph_es_from_vectorview, igraph_es_to_vectorview,
  igraph_es_unfold_vectorview, igraph_es_destroy_vectorview
};

typedef igraph_i_vs_vectorview_pdata_t igraph_i_es_vectorview_pdata_t;

int igraph_es_vectorview(const igraph_t *igraph, igraph_es_t *es, 
			 const vector_t *eids) {
  igraph_i_es_vectorview_pdata_t *data;
  vector_t *fakev;

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(vector_t*) &data->v;
  *fakev = *eids;
  data->destroy=0;

  es->stdata[1]=vector_size(&data->v);

  return 0;
}

int igraph_es_vector(const igraph_t *igraph, igraph_es_t *es,
		     const vector_t *eids) {
  igraph_i_es_vectorview_pdata_t *data;
  vector_t *fakev;

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(vector_t*)&data->v;
  IGRAPH_CHECK(vector_copy(fakev, eids));
  data->destroy=1;

  es->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

int igraph_es_vectorview_it(const igraph_t *graph, const igraph_es_t *es,
			    igraph_es_t *newes) {
  igraph_i_es_vectorview_pdata_t *data, *esdata;
  vector_t *fakev;

  newes->type=IGRAPH_ITERATOR_ES_VECTOR;
  newes->stdata[0]=0;
  newes->table=&igraph_i_es_vectorview_table;
  newes->shorthand=0;
  
  newes->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (newes->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, newes->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)newes->pdata;
  VECTOR_INIT_FINALLY((vector_t*) &data->v, 0);
  fakev=(vector_t*)&data->v;
  if (es->type != IGRAPH_ITERATOR_ES_VECTOR) {
    vector_init(fakev, 0);
    IGRAPH_CHECK(igraph_es_unfold(graph, es, fakev));
    data->destroy=1;
  } else {
    /* If it is a vector we just create a view */
    esdata=(igraph_i_es_vectorview_pdata_t*)es->pdata;    
    *fakev=esdata->v;
    data->destroy=0;
  }

  newes->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_es_vector_small(const igraph_t *igraph, igraph_es_t *es, ...) {
  va_list ap;
  igraph_i_es_vectorview_pdata_t *data;
  vector_t *fakev;
  long int i, n=0;
  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    IGRAPH_ERROR("Cannot create vector iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, es);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  VECTOR_INIT_FINALLY(fakev, n);
  
  va_start(ap, es);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  es->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

const igraph_es_t *IGRAPH_ES_VECTOR(const igraph_t *graph, 
				    const vector_t *eids) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }
  igraph_es_vectorview(graph, es, eids);
  es->shorthand=1;
  return es;
}

const igraph_es_t *IGRAPH_ES(const igraph_t *graph, ...) {
  igraph_es_t *es=Calloc(1, igraph_es_t);
  va_list ap;
  igraph_i_es_vectorview_pdata_t *data;
  vector_t *fakev;
  long int i, n=0;
  int ret;
  if (es==0) {
    igraph_error("Cannot create iterator shorthand", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;			/* TODO: how to sign error ??? */
  }

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;
  es->shorthand=0;
  
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata==0) {
    igraph_error("Cannot create vector iterator", __FILE__, __LINE__,
		 IGRAPH_ENOMEM);
    return 0;
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  fakev=(vector_t*)&data->v;
  data->destroy=1;

  va_start(ap, graph);
  while (1) {
    int num = va_arg(ap, int);
    if (num == -1) {
      break;
    }
    n++;
  }
  va_end(ap);

  ret=vector_init(fakev, n);
  if (ret != 0) {
    igraph_error("Cannot create vector for iterator shorthand", __FILE__, 
		 __LINE__, ret);
    return 0;
  }
  IGRAPH_FINALLY(vector_destroy, fakev);
  
  va_start(ap, graph);
  for (i=0; i<n; i++) {
    VECTOR(*fakev)[i]=(real_t) va_arg(ap, int);
  }
  va_end(ap);  
  
  es->stdata[1]=vector_size(&data->v);
  
  IGRAPH_FINALLY_CLEAN(2);
  es->shorthand=1;
  return es;  
}

void igraph_es_next_vectorview(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_vectorview(const igraph_t *graph, const igraph_es_t *es) {
  return es->stdata[0] >= es->stdata[1];
}

void igraph_es_reset_vectorview(const igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_vectorview(const igraph_t *graph, 
				   const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  return VECTOR(data->v)[ (long int) (es->stdata[0]) ]; 
}

integer_t igraph_es_from_vectorview(const igraph_t *graph, 
				    const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  long int id=VECTOR(data->v)[ (long int) (es->stdata[0]) ];
  return VECTOR(graph->from) [id];
}

integer_t igraph_es_to_vectorview(const igraph_t *graph, 
				  const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  long int id=VECTOR(data->v)[ (long int) (es->stdata[0]) ];
  return VECTOR(graph->to) [id];
}

int igraph_es_unfold_vectorview(const igraph_t *graph, const igraph_es_t *es,
				vector_t *v) {
  vector_t v2;
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  IGRAPH_CHECK(vector_copy(&v2, &data->v));
  vector_destroy(v);
  *v=v2;
  return 0;
}

void igraph_es_destroy_vectorview(igraph_es_t *pes) {
  igraph_es_t *es=(igraph_es_t*)pes;
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  if (data->destroy) {
    vector_destroy((vector_t*)&data->v);
  }
  Free(data);
  
  if (es->shorthand) {
    Free(pes);
  }
}

const vector_t *igraph_es_vector_getvector(const igraph_t *graph, 
					   const igraph_es_t *es) {
  igraph_i_es_vectorview_pdata_t *data=
    (igraph_i_es_vectorview_pdata_t*)es->pdata;
  return &data->v;  
}

/* -------------------------------------------------- */
/* edge iterator, all edges between two vertex sets   */
/* -------------------------------------------------- */

int igraph_es_fromto(const igraph_t *graph, igraph_es_t *es, 
		     const igraph_vs_t *from, const igraph_vs_t *to, 
		     bool_t directed) {

  igraph_vs_t myfrom, myto;  
  long int i, j, lfrom;
  igraph_es_t edgeit;
  igraph_neimode_t mode;
  const vector_t *fromvect;
  vector_t tovect;
  
  igraph_i_es_vectorview_pdata_t *data;

  if (directed && igraph_is_directed(graph)) {
    mode=IGRAPH_OUT;
  } else {
    mode=IGRAPH_ALL;
  }

  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, from, &myfrom));
  IGRAPH_FINALLY(igraph_vs_destroy, &myfrom);
  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, to, &myto));
  IGRAPH_FINALLY(igraph_vs_destroy, &myto);  

  fromvect=igraph_vs_vector_getvector(graph, &myfrom);
  IGRAPH_CHECK(vector_copy(&tovect, igraph_vs_vector_getvector(graph, &myto)));
  IGRAPH_FINALLY(vector_destroy, &tovect);

  es->type=IGRAPH_ITERATOR_ES_VECTOR;
  es->shorthand=0;
  es->stdata[0]=0;
  es->table=&igraph_i_es_vectorview_table;    
  es->pdata=Calloc(1, igraph_i_es_vectorview_pdata_t);
  if (es->pdata == 0) {
    IGRAPH_ERROR("Cannot create iterator", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, es->pdata);
  data=(igraph_i_es_vectorview_pdata_t*)es->pdata;
  VECTOR_INIT_FINALLY((vector_t*) &data->v, 0);
  data->destroy=1;

  lfrom=vector_size(fromvect);
  if (lfrom != 0 && vector_size(&tovect) != 0) {
    IGRAPH_CHECK(igraph_es_adj(graph, &edgeit, 0, mode));
    IGRAPH_FINALLY(igraph_es_destroy, &edgeit);
    vector_sort(&tovect);
    for (i=0; i<lfrom; i++) {
      long int vfrom=VECTOR(*fromvect)[i];
      igraph_es_adj_set(graph, &edgeit, vfrom, mode);
      while (!igraph_es_end(graph, &edgeit)) {
	long int vto=igraph_es_adj_vertex(graph, &edgeit);
	if (vector_binsearch(&tovect, vto, 0)) {
	  vector_push_back((vector_t*)&data->v, igraph_es_get(graph, &edgeit));
	}
	igraph_es_next(graph, &edgeit);
      }
    }
    IGRAPH_FINALLY_CLEAN(1);
  }

  es->stdata[1]=vector_size(&data->v);

  /* Clean */
  igraph_vs_destroy(&myfrom);
  igraph_vs_destroy(&myto);
  vector_destroy(&tovect);
  IGRAPH_FINALLY_CLEAN(5);

  if (from->shorthand) { igraph_vs_destroy((igraph_vs_t*) from); }
  if (to->shorthand) { igraph_vs_destroy((igraph_vs_t*) to); }

  return 0;
}

