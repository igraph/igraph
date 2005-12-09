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

/* -------------------------------------------------- */
/* Vertex iterator generics                           */
/* -------------------------------------------------- */

void igraph_vs_next(igraph_t *graph, igraph_vs_t *vs) {
  vs->table->next(graph, (struct igraph_vs_t*)vs);
}

bool_t igraph_vs_end(igraph_t *graph, igraph_vs_t *vs) {
  return vs->table->end(graph, (struct igraph_vs_t*)vs);
}

void igraph_vs_reset(igraph_t *graph, igraph_vs_t *vs) {
  vs->table->reset(graph, (struct igraph_vs_t*)vs);
}

integer_t igraph_vs_get(igraph_t *graph, igraph_vs_t *vs) {
  return vs->table->get(graph, (struct igraph_vs_t*)vs);
}

int igraph_vs_unfold(igraph_t *graph, igraph_vs_t *vs) {
  return vs->table->unfold(graph, (struct igraph_vs_t*)vs);
}

void igraph_vs_destroy(igraph_vs_t *vs) {
  vs->table->destroy((struct igraph_vs_t*)vs);
}

int igraph_vs_create_view_as_vector(igraph_t *graph, igraph_vs_t *vs,
				    igraph_vs_t *newvs) {
  newvs->type = vs->type;
  memcpy(newvs->stdata, vs->stdata, sizeof(real_t)*IGRAPH_I_STDATA_SIZE);
  newvs->pdata = vs->pdata;
  newvs->table = vs->table;
  newvs->v = vs->v;
  newvs->view=1;
  if (vs->v==0) {
    IGRAPH_CHECK(igraph_vs_unfold(graph, newvs));
    newvs->vdestroy=1;
  } else {
    newvs->v=vs->v;
    newvs->vdestroy=0;
  }
  return 0;
}

/* -------------------------------------------------- */
/* Simple vertex iterator                             */
/* -------------------------------------------------- */

void igraph_vs_next_all(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_all(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_all(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_all(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_all(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_all(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_all_table = {
  igraph_vs_next_all, igraph_vs_end_all, igraph_vs_reset_all,
  igraph_vs_get_all, igraph_vs_unfold_all, igraph_vs_destroy_all
};

igraph_vs_t igraph_vs_all(igraph_t *graph) {
  igraph_vs_t vs;
  vs.type=IGRAPH_ITERATOR_VS_ALL;
  vs.stdata[0]=0;
  vs.table=&igraph_i_vs_all_table;
  vs.v=0;
  vs.vdestroy=0;
  vs.view=0;
  return vs;
}

void igraph_vs_next_all(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_all(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0] >= igraph_vcount(graph);
}

void igraph_vs_reset_all(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=0;
}

integer_t igraph_vs_get_all(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_all(igraph_t *graph, igraph_vs_t *vs) {
  long int n;
  if (vs->v==0) {		/* otherwise already unfolded */
    n=igraph_vcount(graph);
    vs->v=Calloc(1, vector_t);
    if (vs->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vs->v); /* free it if cannot allocate vector */
    IGRAPH_CHECK(vector_init_seq(vs->v, 0, n-1));
    IGRAPH_FINALLY_CLEAN(1);
    vs->vdestroy=1;
  }
  return 0;
}

void igraph_vs_destroy_all(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }
}

/* -------------------------------------------------- */
/* Adjacent vertices of a vertex                      */
/* -------------------------------------------------- */

void igraph_vs_next_adj(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_adj(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_adj(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_adj(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_adj(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_adj(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_adj_table = {
  igraph_vs_next_adj, igraph_vs_end_adj, igraph_vs_reset_adj,
  igraph_vs_get_adj, igraph_vs_unfold_adj, igraph_vs_destroy_adj
};

igraph_vs_t igraph_vs_adj(igraph_t *graph,
			  integer_t vid, igraph_neimode_t mode) {
  igraph_vs_t vs;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs.type=IGRAPH_ITERATOR_VS_ADJ;
  vs.table=&igraph_i_vs_adj_table;
  vs.v=0;
  vs.vdestroy=0;
  vs.view=0;
  
  vs.stdata[0]=vid;
  vs.stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    vs.stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    vs.stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    vs.stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    vs.stdata[3]=igraph_ecount(graph);
  }
  
  return vs;
}

void igraph_vs_next_adj(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[2] ++;
  if (vs->stdata[2] > VECTOR(graph->os)[ (long int)vs->stdata[0]+1 ]) {
    vs->stdata[3] ++;
  }
}

bool_t igraph_vs_end_adj(igraph_t *graph, igraph_vs_t *vs) {
  return (vs->stdata[2] >= VECTOR(graph->os)[ (long int) vs->stdata[0]+1 ] &&
	  vs->stdata[3] >= VECTOR(graph->is)[ (long int) vs->stdata[0]+1 ]);
}

void igraph_vs_reset_adj(igraph_t *graph, igraph_vs_t *vs) {
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

integer_t igraph_vs_get_adj(igraph_t *graph, igraph_vs_t *vs) {
  if (vs->stdata[2] < VECTOR(graph->os)[ (long int)vs->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)vs->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)vs->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

int igraph_vs_unfold_adj(igraph_t *graph, igraph_vs_t *vs) {
  if (vs->v==0) {
    vs->v=Calloc(1, vector_t);
    if (vs->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vs->v); /* free it if cannot allocate vector */
    IGRAPH_CHECK(vector_init(vs->v, 0));
    IGRAPH_CHECK(igraph_neighbors(graph, vs->v, vs->stdata[0], vs->stdata[1]));
    IGRAPH_FINALLY_CLEAN(1);
    vs->vdestroy=1;
  }
  return 0;
}

void igraph_vs_adj_set(igraph_t *graph, igraph_vs_t *vs,
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
}

void igraph_vs_destroy_adj(igraph_vs_t *vs) {
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }
}  

/* -------------------------------------------------- */
/* Random walker                                      */
/* -------------------------------------------------- */

void igraph_vs_next_rw(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_rw(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_rw(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_rw(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_rw(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_rw(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_rw_table = {
  igraph_vs_next_rw, igraph_vs_end_rw, igraph_vs_reset_rw,
  igraph_vs_get_rw, igraph_vs_unfold_rw, igraph_vs_destroy_rw
};

igraph_vs_t igraph_vs_rw(igraph_t *graph,
			 integer_t vid, igraph_neimode_t mode) {
  igraph_vs_t vs;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs.type=IGRAPH_ITERATOR_VS_RW;
  
  vs.stdata[0]=vid;		/* actual vertex */
  vs.stdata[1]=vid;		/* start vertex  */
  vs.stdata[2]=mode;		/* mode */
  vs.stdata[3]=0;		/* number of steps so far */

  vs.table=&igraph_i_vs_rw_table;
  vs.v=0;
  vs.vdestroy=0;
  vs.view=0;
  return vs;
}

void igraph_vs_next_rw(igraph_t *graph, igraph_vs_t *vs) {
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

bool_t igraph_vs_end_rw(igraph_t *graph, igraph_vs_t *vs) {
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

void igraph_vs_reset_rw(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
  vs->stdata[3]=0;
}

integer_t igraph_vs_get_rw(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_rw(igraph_t *graph, igraph_vs_t *vs) {
  IGRAPH_FERROR("attempt to unfold random walker", IGRAPH_EUNFOLDINF);
  return 0;			/* return unneccesary */
}

long int igraph_vs_rw_length(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[3];
}

void igraph_vs_destroy_rw(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }
}

/* -------------------------------------------------- */
/* Random walker with one unit memory                 */
/* -------------------------------------------------- */

void igraph_vs_next_rw1(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_rw1(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_rw1(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_rw1(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_rw1(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_rw1(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_rw1_table = {
  igraph_vs_next_rw1, igraph_vs_end_rw1, igraph_vs_reset_rw1,
  igraph_vs_get_rw1, igraph_vs_unfold_rw1, igraph_vs_destroy_rw1
};

igraph_vs_t igraph_vs_rw1(igraph_t *graph,
			    integer_t vid, igraph_neimode_t mode) {
  igraph_vs_t vs;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  vs.type=IGRAPH_ITERATOR_VS_RW1;
  
  vs.stdata[0]=vid;			/* actual vertex */
  vs.stdata[1]=vid;			/* start vertex  */
  vs.stdata[2]=mode;			/* mode */
  vs.stdata[3]=0;			/* number of steps so far */
  vs.stdata[4]=-1;		        /* the previous vertex */

  vs.table=&igraph_i_vs_rw1_table;
  vs.v=0;
  vs.vdestroy=0;
  vs.view=0;
  return vs;
}

void igraph_vs_next_rw1(igraph_t *graph, igraph_vs_t *vs) {

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

bool_t igraph_vs_end_rw1(igraph_t *graph, igraph_vs_t *vs) {
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

integer_t igraph_vs_get_rw1(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0];
}

void igraph_vs_reset_rw1(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
  vs->stdata[3]=0;
  vs->stdata[5]=-1;
}

int igraph_vs_unfold_rw1(igraph_t *graph, igraph_vs_t *vs) {
  IGRAPH_FERROR("attempt to unfold random walker", IGRAPH_EUNFOLDINF);
  return 0;
}

long int igraph_vs_rw1_length(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[3];
}

void igraph_vs_destroy_rw1(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }  
}

/* -------------------------------------------------- */
/* Empty vertex set                                   */
/* -------------------------------------------------- */

void igraph_vs_next_none(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_none(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_none(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_none(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_none(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_none(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_none_table = {
  igraph_vs_next_none, igraph_vs_end_none, igraph_vs_reset_none,
  igraph_vs_get_none, igraph_vs_unfold_none, igraph_vs_destroy_none
};

igraph_vs_t igraph_vs_none(igraph_t *graph) {
  igraph_vs_t vs;
  vs.type=IGRAPH_ITERATOR_VS_NONE;
  vs.table=&igraph_i_vs_none_table;
  vs.v=0;
  vs.vdestroy=0;
  vs.view=0;
  return vs;
}

void igraph_vs_next_none(igraph_t *graph, igraph_vs_t *vs) {
  /* nothing to do */
}

bool_t igraph_vs_end_none(igraph_t *graph, igraph_vs_t *vs) {
  return 1;
}

void igraph_vs_reset_none(igraph_t *graph, igraph_vs_t *vs) {
  /* nothing to do */
}

integer_t igraph_vs_get_none(igraph_t *graph, igraph_vs_t *vs) {
  /* ooops this is an error, no way to signal it though... */
  return -1;
}

int igraph_vs_unfold_none(igraph_t *graph, igraph_vs_t *vs) {
  if (vs->v==0) {		/* otherwise already unfolded */
    vs->v=Calloc(1, vector_t);
    if (vs->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vs->v); /* free it if cannot allocate vector */
    IGRAPH_CHECK(vector_init(vs->v, 0));
    IGRAPH_FINALLY_CLEAN(1);
    vs->vdestroy=1;
  } 
  return 0;
}

void igraph_vs_destroy_none(igraph_vs_t *vs) {
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }
}

/* -------------------------------------------------- */
/* vertex set with single vertex                      */
/* -------------------------------------------------- */

void igraph_vs_next_1(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_1(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_1(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_1(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_1(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_1(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_1_table = {
  igraph_vs_next_1, igraph_vs_end_1, igraph_vs_reset_1,
  igraph_vs_get_1, igraph_vs_unfold_1, igraph_vs_destroy_1
};


igraph_vs_t igraph_vs_1(igraph_t *igraph, integer_t vid) {
  igraph_vs_t vs;
  vs.type=IGRAPH_ITERATOR_VS_1;
  vs.stdata[0]=vid;
  vs.stdata[1]=0;		/* write 1 here if end */
  vs.table=&igraph_i_vs_1_table;
  vs.v=0;
  vs.vdestroy=0;
  vs.view=0;
  return vs;
}

void igraph_vs_next_1(igraph_t *graph, igraph_vs_t *vs) {
  /* signal end */
  vs->stdata[1]=1;
}

bool_t igraph_vs_end_1(igraph_t *graph, igraph_vs_t *vs) {
  return (vs->stdata[1]==1);
}

void igraph_vs_reset_1(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[1]=0;
}

integer_t igraph_vs_get_1(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_1(igraph_t *graph, igraph_vs_t *vs) {
  if (vs->v==0) {
    vs->v=Calloc(1, vector_t);
    if (vs->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vs->v);
    IGRAPH_CHECK(vector_init(vs->v, 1));
    VECTOR(*vs->v)[0]=vs->stdata[0];
    IGRAPH_FINALLY_CLEAN(1);
    vs->vdestroy=1;
  }
  return 0;
}

void igraph_vs_destroy_1(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }
}

/* -------------------------------------------------- */
/* Vertex set with sequence of vertices               */
/* -------------------------------------------------- */

void igraph_vs_next_seq(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_seq(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_seq(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_seq(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_seq(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_seq(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_seq_table = {
  igraph_vs_next_seq, igraph_vs_end_seq, igraph_vs_reset_seq,
  igraph_vs_get_seq, igraph_vs_unfold_seq, igraph_vs_destroy_seq
};

igraph_vs_t igraph_vs_seq(igraph_t *igraph, integer_t from,
			   integer_t to) {
  igraph_vs_t vs;
  vs.type=IGRAPH_ITERATOR_VS_SEQ;
  vs.stdata[0]=from;
  vs.stdata[1]=from;
  vs.stdata[2]=to;
  vs.table=&igraph_i_vs_seq_table;
  vs.v=0;
  vs.vdestroy=0;
  vs.view=0;
  return vs;
}

void igraph_vs_next_seq(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_seq(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0] > vs->stdata[2];
}

void igraph_vs_reset_seq(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=vs->stdata[1];
}

integer_t igraph_vs_get_seq(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0];
}

int igraph_vs_unfold_seq(igraph_t *graph, igraph_vs_t *vs) {
  if (vs->v==0) {
    vs->v=Calloc(1, vector_t);
    if (vs->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, vs->v);
    IGRAPH_CHECK(vector_init_seq(vs->v, vs->stdata[1], vs->stdata[2]));
    IGRAPH_FINALLY_CLEAN(1);
    vs->vdestroy=1;
  }
  return 0;
}

void igraph_vs_destroy_seq(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }
}

/* -------------------------------------------------- */
/* Vertex ids in a vector, this is a view             */
/* -------------------------------------------------- */

void igraph_vs_next_vector(igraph_t *graph, igraph_vs_t *vs);
bool_t igraph_vs_end_vector(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_reset_vector(igraph_t *graph, igraph_vs_t *vs);
integer_t igraph_vs_get_vector(igraph_t *graph, igraph_vs_t *vs);
int igraph_vs_unfold_vector(igraph_t *graph, igraph_vs_t *vs);
void igraph_vs_destroy_vector(igraph_vs_t *vs);

igraph_i_vstable_t igraph_i_vs_vector_table = {
  igraph_vs_next_vector, igraph_vs_end_vector, igraph_vs_reset_vector,
  igraph_vs_get_vector, igraph_vs_unfold_vector, igraph_vs_destroy_vector
};

igraph_vs_t igraph_vs_vector(igraph_t *igraph, vector_t *vids) {
  igraph_vs_t vs;
  vs.type=IGRAPH_ITERATOR_VS_VECTOR;
  vs.stdata[0]=0;
  vs.table=&igraph_i_vs_vector_table;
  vs.v=vids;
  vs.vdestroy=0;
  vs.view=0;
  return vs;
}

void igraph_vs_next_vector(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0] ++;
}

bool_t igraph_vs_end_vector(igraph_t *graph, igraph_vs_t *vs) {
  return vs->stdata[0] >= vector_size(vs->v);
}

void igraph_vs_reset_vector(igraph_t *graph, igraph_vs_t *vs) {
  vs->stdata[0]=0;
}

integer_t igraph_vs_get_vector(igraph_t *graph, igraph_vs_t *vs) {
  return VECTOR(*vs->v)[ (long int) (vs->stdata[0]) ];
}

int igraph_vs_unfold_vector(igraph_t *graph, igraph_vs_t *vs) {
  /* nothing to do */
  return 0;
}

void igraph_vs_destroy_vector(igraph_vs_t *pvs) {
  igraph_vs_t *vs=(igraph_vs_t*)pvs;
  if (vs->vdestroy) {
    vector_destroy(vs->v);
    Free(vs->v);
  }
}

/* -------------------------------------------------- */
/* Edge iterator generics                             */
/* -------------------------------------------------- */

void igraph_es_next(igraph_t *graph, igraph_es_t *es) {
  es->table->next(graph, es);
}

bool_t igraph_es_end(igraph_t *graph, igraph_es_t *es) {
  return es->table->end(graph, es);
}

void igraph_es_reset(igraph_t *graph, igraph_es_t *es) {
  return es->table->reset(graph, es);
}

integer_t igraph_es_from(igraph_t *graph, igraph_es_t *es) {
  return es->table->from(graph, es);
}

integer_t igraph_es_to(igraph_t *graph, igraph_es_t *es) {
  return es->table->to(graph, es);
}

integer_t igraph_es_get(igraph_t *graph, igraph_es_t *es) {
  return es->table->get(graph, es);
}

int igraph_es_unfold(igraph_t *graph, igraph_es_t *es) {
  return es->table->unfold(graph, es);
}

void igraph_es_destroy(igraph_es_t *es) {
  es->table->destroy(es);
}

/* -------------------------------------------------- */
/* Simple edge iterator                               */
/* -------------------------------------------------- */

void igraph_es_next_all(igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_all(igraph_t *graph, igraph_es_t *es);
void igraph_es_reset_all(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_all(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_from_all(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_to_all(igraph_t *graph, igraph_es_t *es);
int igraph_es_unfold_all(igraph_t *graph, igraph_es_t *es);
void igraph_es_destroy_all(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_all_table = {
  igraph_es_next_all, igraph_es_end_all, igraph_es_reset_all,
  igraph_es_get_all, igraph_es_from_all, igraph_es_to_all,
  igraph_es_unfold_all, igraph_es_destroy_all
};

igraph_es_t igraph_es_all(igraph_t *graph) {
  igraph_es_t es;
  es.type=IGRAPH_ITERATOR_ES_ALL;
  es.stdata[0]=0;
  es.table=&igraph_i_es_all_table;
  es.v=0;
  es.vdestroy=0;
  es.view=0;
  return es;
}

void igraph_es_next_all(igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_all(igraph_t *graph, igraph_es_t *es) {
  return es->stdata[0] >= igraph_ecount(graph);
}

void igraph_es_reset_all(igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_all(igraph_t *graph, igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_all(igraph_t *graph, igraph_es_t *es) {
  return VECTOR(graph->from)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_to_all(igraph_t *graph, igraph_es_t *es) {
  return VECTOR(graph->to)[ (long int) es->stdata[0] ];
}

int igraph_es_unfold_all(igraph_t *graph, igraph_es_t *es) {
  long int n;
  if (es->v==0) {		/* otherwise already unfolded */
    n=igraph_ecount(graph);
    es->v=Calloc(1, vector_t);
    if (es->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, es->v); /* free it if cannot allocate vector */
    IGRAPH_CHECK(vector_init_seq(es->v, 0, n-1));
    IGRAPH_FINALLY_CLEAN(1);
    es->vdestroy=1;
  }
  return 0;
}

void igraph_es_destroy_all(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->vdestroy) {
    vector_destroy(es->v);
    Free(es->v);
  }
}

/* -------------------------------------------------- */
/* Ordered edge iterator                              */
/* -------------------------------------------------- */

void igraph_es_next_fromorder(igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_fromorder(igraph_t *graph, igraph_es_t *es);
void igraph_es_reset_fromorder(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_fromorder(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_from_fromorder(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_to_fromorder(igraph_t *graph, igraph_es_t *es);
int igraph_es_unfold_fromorder(igraph_t *graph, igraph_es_t *es);
void igraph_es_destroy_fromorder(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_fromorder_table = {
  igraph_es_next_fromorder, igraph_es_end_fromorder, igraph_es_reset_fromorder,
  igraph_es_get_fromorder, igraph_es_from_fromorder, igraph_es_to_fromorder,
  igraph_es_unfold_fromorder, igraph_es_destroy_fromorder
};
  
igraph_es_t igraph_es_fromorder(igraph_t *graph) {
  igraph_es_t es;
  es.type=IGRAPH_ITERATOR_ES_FROMORDER;
  es.stdata[0]=0;
  es.table=&igraph_i_es_fromorder_table;
  es.v=0;
  es.vdestroy=0;
  es.view=0;
  return es;
}

void igraph_es_next_fromorder(igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_fromorder(igraph_t *graph, igraph_es_t *es) {
  return es->stdata[0] >= igraph_ecount(graph);
}

void igraph_es_reset_fromorder(igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_fromorder(igraph_t *graph, igraph_es_t *es) {
  return VECTOR(graph->oi)[ (long int) es->stdata[0] ];
}

integer_t igraph_es_from_fromorder(igraph_t *graph, igraph_es_t *es) {
  long int idx=VECTOR(graph->oi)[ (long int) es->stdata[0] ];
  return VECTOR(graph->from)[ idx ];
}

integer_t igraph_es_to_fromorder(igraph_t *graph, igraph_es_t *es) {
  long int idx=VECTOR(graph->oi)[ (long int) es->stdata[0] ];
  return VECTOR(graph->to)[ idx ];
}

int igraph_es_unfold_fromorder(igraph_t *graph, igraph_es_t *es) {
  if (es->v==0) {
    long int i=0, n=igraph_ecount(graph);  
    es->v=Calloc(1, vector_t);
    if (es->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, es->v);
    IGRAPH_CHECK(vector_init(es->v, n));
    for (i=0; i<n; i++) {
      VECTOR(*es->v)[i]=VECTOR(graph->oi)[i];
    }
    IGRAPH_FINALLY_CLEAN(1);
    es->vdestroy=1;
  }
  return 0;
}

void igraph_es_destroy_fromorder(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->vdestroy) {
    vector_destroy(es->v);
    Free(es->v);
  }
}

/* -------------------------------------------------- */
/* Adjacent edges of a vertex                         */
/* -------------------------------------------------- */

void igraph_es_next_adj(igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_adj(igraph_t *graph, igraph_es_t *es);
void igraph_es_reset_adj(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_adj(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_from_adj(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_to_adj(igraph_t *graph, igraph_es_t *es);
int igraph_es_unfold_adj(igraph_t *graph, igraph_es_t *es);
void igraph_es_destroy_adj(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_adj_table = {
  igraph_es_next_adj, igraph_es_end_adj, igraph_es_reset_adj,
  igraph_es_get_adj, igraph_es_from_adj, igraph_es_to_adj,
  igraph_es_unfold_adj, igraph_es_destroy_adj
};
  
igraph_es_t igraph_es_adj(igraph_t *graph,
			    integer_t vid, igraph_neimode_t mode) {
  igraph_es_t es;
  if (!igraph_is_directed(graph)) {
    mode=IGRAPH_ALL;
  }
  
  es.type=IGRAPH_ITERATOR_ES_ADJ;

  es.stdata[0]=vid;
  es.stdata[1]=mode;
  if (mode & IGRAPH_OUT) {
    es.stdata[2]=VECTOR(graph->os)[(long int)vid];
  } else {
    es.stdata[2]=igraph_ecount(graph);
  }
  if (mode & IGRAPH_IN) {
    es.stdata[3]=VECTOR(graph->is)[(long int)vid];
  } else {
    es.stdata[3]=igraph_ecount(graph);
  }

  es.table=&igraph_i_es_adj_table;
  es.v=0;
  es.vdestroy=0;
  es.view=0;
  return es;
}

void igraph_es_next_adj(igraph_t *graph, igraph_es_t *es) {
  es->stdata[2] ++;
  if (es->stdata[2] > VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    es->stdata[3] ++;
  }
}

bool_t igraph_es_end_adj(igraph_t *graph, igraph_es_t *es) {
  return (es->stdata[2] >= VECTOR(graph->os)[ (long int) es->stdata[0]+1 ] &&
	  es->stdata[3] >= VECTOR(graph->is)[ (long int) es->stdata[0]+1 ]);
}

void igraph_es_reset_adj(igraph_t *graph, igraph_es_t *es) {
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

integer_t igraph_es_get_adj(igraph_t *graph, igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    return VECTOR(graph->oi)[(long int)es->stdata[2]];
  } else {
    return VECTOR(graph->ii)[(long int)es->stdata[3]];
  }
}

integer_t igraph_es_from_adj(igraph_t *graph, igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->from)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->from)[idx];
  }
}

integer_t igraph_es_to_adj(igraph_t *graph, igraph_es_t *es) {
  if (es->stdata[2] < VECTOR(graph->os)[ (long int)es->stdata[0]+1 ]) {
    long int idx=VECTOR(graph->oi)[(long int)es->stdata[2]];
    return VECTOR(graph->to)[idx];
  } else {
    long int idx=VECTOR(graph->ii)[(long int)es->stdata[3]];
    return VECTOR(graph->to)[idx];
  }
}

int igraph_es_unfold_adj(igraph_t *graph, igraph_es_t *es) {
  if (es->v==0) {
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
  
    es->v=Calloc(1, vector_t);
    if (es->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, es->v);
    IGRAPH_CHECK(vector_init(es->v, length));
    IGRAPH_FINALLY_CLEAN(1);
    es->vdestroy=1;
  
    if (mode & IGRAPH_OUT) {
      for (i=VECTOR(graph->os)[node]; i<VECTOR(graph->os)[node+1]; i++) {
	VECTOR(*es->v)[idx++] = VECTOR(graph->oi)[i];
      }
    }
    if (mode & IGRAPH_IN) {
      for (i=VECTOR(graph->is)[node]; i<VECTOR(graph->is)[node+1]; i++) {
	VECTOR(*es->v)[idx++] = VECTOR(graph->ii)[i];
      }
    }
  }

  return 0;
}

void igraph_es_destroy_adj(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->vdestroy) {
    vector_destroy(es->v);
    Free(es->v);
  }
}

void igraph_es_adj_set(igraph_t *graph, igraph_es_t *es,
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

integer_t igraph_es_adj_vertex(igraph_t *graph, igraph_es_t *es) {
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

void igraph_es_next_none(igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_none(igraph_t *graph, igraph_es_t *es);
void igraph_es_reset_none(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_none(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_from_none(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_to_none(igraph_t *graph, igraph_es_t *es);
int igraph_es_unfold_none(igraph_t *graph, igraph_es_t *es);
void igraph_es_destroy_none(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_none_table = {
  igraph_es_next_none, igraph_es_end_none, igraph_es_reset_none,
  igraph_es_get_none, igraph_es_from_none, igraph_es_to_none,
  igraph_es_unfold_none, igraph_es_destroy_none
};

igraph_es_t igraph_es_none(igraph_t *graph) {
  igraph_es_t es;
  es.type=IGRAPH_ITERATOR_ES_NONE;
  
  es.table=&igraph_i_es_none_table;
  es.v=0;
  es.vdestroy=0;
  es.view=0;
  return es;
}

void igraph_es_next_none(igraph_t *graph, igraph_es_t *es) {
  /* nothing to do */
}

bool_t igraph_es_end_none(igraph_t *graph, igraph_es_t *es) {
  return 1;
}

void igraph_es_reset_none(igraph_t *graph, igraph_es_t *es) {
  /* nothing to do */
}

integer_t igraph_es_get_none(igraph_t *graph, igraph_es_t *es) {
  /* ooops this is an error, no way to signal it though... */
  return -1;
}

integer_t igraph_es_from_none(igraph_t *graph, igraph_es_t *es) {
  /* error */
  return 0;
}

integer_t igraph_es_to_none(igraph_t *graph, igraph_es_t *es) {
  /* error */
  return 0;
}

int igraph_es_unfold_none(igraph_t *graph, igraph_es_t *es) {
  if (es->v==0) {
    es->v=Calloc(1, vector_t);
    if (es->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, es->v);
    IGRAPH_CHECK(vector_init(es->v, 0));
    IGRAPH_FINALLY_CLEAN(1);
    es->vdestroy=1;
  }
  return 0;
}

void igraph_es_destroy_none(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->vdestroy) {
    vector_destroy(es->v);
    Free(es->v);
  }
}

/* -------------------------------------------------- */
/* edge iterator, single edge                         */
/* -------------------------------------------------- */

void igraph_es_next_1(igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_1(igraph_t *graph, igraph_es_t *es);
void igraph_es_reset_1(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_1(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_from_1(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_to_1(igraph_t *graph, igraph_es_t *es);
int igraph_es_unfold_1(igraph_t *graph, igraph_es_t *es);
void igraph_es_destroy_1(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_1_table = {
  igraph_es_next_1, igraph_es_end_1, igraph_es_reset_1,
  igraph_es_get_1, igraph_es_from_1, igraph_es_to_1,
  igraph_es_unfold_1, igraph_es_destroy_1
};

igraph_es_t igraph_es_1(igraph_t *igraph, integer_t eid) {
  igraph_es_t es;
  es.type=IGRAPH_ITERATOR_ES_1;
  es.stdata[0]=eid;
  es.stdata[1]=0;		/* write 1 here if end */
  es.table=&igraph_i_es_1_table;
  es.v=0;
  es.vdestroy=0;
  es.view=0;
  return es;
}

void igraph_es_next_1(igraph_t *graph, igraph_es_t *es) {
  /* signal end */
  es->stdata[1]=1;
}

bool_t igraph_es_end_1(igraph_t *graph, igraph_es_t *es) {
  return (es->stdata[1]==1);
}

void igraph_es_reset_1(igraph_t *graph, igraph_es_t *es) {
  es->stdata[1]=0;
}

integer_t igraph_es_get_1(igraph_t *graph, igraph_es_t *es) {
  return es->stdata[0];
}

integer_t igraph_es_from_1(igraph_t *graph, igraph_es_t *es) {
  /* error */
  return 0;
}

integer_t igraph_es_to_1(igraph_t *graph, igraph_es_t *es) {
  /* error */
  return 0;
}

int igraph_es_unfold_1(igraph_t *graph, igraph_es_t *es) {
  if (es->v==0) {
    es->v=Calloc(1, vector_t);
    if (es->v==0) {
      IGRAPH_FERROR("Cannot unfold iterator", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, es->v);
    IGRAPH_CHECK(vector_init(es->v, 1));
    VECTOR(*es->v)[0]=es->stdata[0];
    IGRAPH_FINALLY_CLEAN(1);
    es->vdestroy=1;
  }
  return 0;
}

void igraph_es_destroy_1(igraph_es_t *pvs) {
  igraph_es_t *es=(igraph_es_t*)pvs;
  if (es->vdestroy) {
    vector_destroy(es->v);
    Free(es->v);
  }
}

/* -------------------------------------------------- */
/* Edge ids in a vector, this is a view               */
/* -------------------------------------------------- */

void igraph_es_next_vector(igraph_t *graph, igraph_es_t *es);
bool_t igraph_es_end_vector(igraph_t *graph, igraph_es_t *es);
void igraph_es_reset_vector(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_get_vector(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_from_vector(igraph_t *graph, igraph_es_t *es);
integer_t igraph_es_to_vector(igraph_t *graph, igraph_es_t *es);
int igraph_es_unfold_vector(igraph_t *graph, igraph_es_t *es);
void igraph_es_destroy_vector(igraph_es_t *es);

igraph_i_estable_t igraph_i_es_vector_table = {
  igraph_es_next_vector, igraph_es_end_vector, igraph_es_reset_vector,
  igraph_es_get_vector, igraph_es_from_vector, igraph_es_to_vector,
  igraph_es_unfold_vector, igraph_es_destroy_vector
};

igraph_es_t igraph_es_vector(igraph_t *igraph, vector_t *eids) {
  igraph_es_t es;
  es.type=IGRAPH_ITERATOR_ES_VECTOR;
  es.stdata[0]=0;
  es.table=&igraph_i_es_vector_table;
  es.v=eids;
  es.vdestroy=0;
  es.view=0;
  return es;
}

void igraph_es_next_vector(igraph_t *graph, igraph_es_t *es) {
  es->stdata[0] ++;
}

bool_t igraph_es_end_vector(igraph_t *graph, igraph_es_t *es) {
  return es->stdata[0] >= vector_size(es->v);
}

void igraph_es_reset_vector(igraph_t *graph, igraph_es_t *es) {
  es->stdata[0]=0;
}

integer_t igraph_es_get_vector(igraph_t *graph, igraph_es_t *es) {
  return VECTOR(*es->v)[ (long int) (es->stdata[0]) ];
}

integer_t igraph_es_from_vector(igraph_t *graph, igraph_es_t *es) {
  long int id=VECTOR(*es->v)[ (long int) (es->stdata[0]) ];
  return VECTOR(graph->from) [id];
}

integer_t igraph_es_to_vector(igraph_t *graph, igraph_es_t *es) {
  long int id=VECTOR(*es->v)[ (long int) (es->stdata[0]) ];
  return VECTOR(graph->to) [id];
}

int igraph_es_unfold_vector(igraph_t *graph, igraph_es_t *es) {
  /* nothing to do */
  return 0;
}

void igraph_es_destroy_vector(igraph_es_t *pes) {
  igraph_es_t *es=(igraph_es_t*)pes;
  if (es->vdestroy) {
    vector_destroy(es->v);
    Free(es->v);
  }
}

/* /\* -------------------------------------------------- *\/ */
/* /\* edge iterator, all edges between two vertex sets   *\/ */
/* /\* -------------------------------------------------- *\/   */

/* void igraph_es_next_fromto(igraph_t *graph, igraph_es_t *es); */
/* bool_t igraph_es_end_fromto(igraph_t *graph, igraph_es_t *es); */
/* void igraph_es_reset_fromto(igraph_t *graph, igraph_es_t *es); */
/* integer_t igraph_es_get_fromto(igraph_t *graph, igraph_es_t *es); */
/* int igraph_es_unfold_fromto(igraph_t *graph, igraph_es_t *es,  */
/* 			   vector_t *vids); */

/* igraph_es_t igraph_es_fromto(igraph_t *igraph, igraph_vs_t from, */
/* 			       igraph_vs_t to, igraph_neimode_t mode) { */
/*   igraph_es_t es; */
/*   es.type=IGRAPH_ITERATOR_FROMTO; */
/*   /\* todo: stdata *\/ */
/*   es.next=igraph_es_next_fromto; */
/*   es.end=igraph_es_end_fromto; */
/*   es.reset=igraph_es_reset_fromto; */
/*   es.get=igraph_es_get_fromto; */
/*   es.unfold=igraph_es_unfold_fromto; */
/*   return es;   */
/* } */

/* void igraph_es_next_fromto(igraph_t *graph, igraph_es_t *es) { */
/*   /\* TODO *\/ */
/* } */

/* bool_t igraph_es_end_fromto(igraph_t *graph, igraph_es_t *es) { */
/*   /\* TODO *\/ */
/*   return 1; */
/* } */


/* void igraph_es_reset_fromto(igraph_t *graph, igraph_es_t *es) { */
/*   /\* TODO *\/ */
/* } */


/* integer_t igraph_es_get_fromto(igraph_t *graph, igraph_es_t *es) { */
/*   /\* TODO *\/ */
/*   return 0; */
/* } */


/* int igraph_es_unfold_fromto(igraph_t *graph, igraph_es_t *es,  */
/* 			   vector_t *vids) { */
/*   /\* TODO *\/ */
/*   return 0; */
/* } */


/* -------------------------------------------------- */
/* FUNCTION TABLES                                    */
/* -------------------------------------------------- */

igraph_i_vstable_t *igraph_i_vstable[8] = {
  &igraph_i_vs_all_table,    &igraph_i_vs_none_table, &igraph_i_vs_adj_table,
  &igraph_i_vs_vector_table, &igraph_i_vs_1_table,    &igraph_i_vs_seq_table,
  &igraph_i_vs_rw_table,     &igraph_i_vs_rw1_table
};

igraph_i_estable_t *igraph_i_estable[6] = {
  &igraph_i_es_all_table,  &igraph_i_es_fromorder_table, 
  &igraph_i_es_adj_table,  &igraph_i_es_none_table, 
  &igraph_i_es_1_table,    &igraph_i_es_vector_table
};
