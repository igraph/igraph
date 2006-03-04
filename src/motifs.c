/* -*- mode: C -*-  */
/* 
   IGraph library.
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
#include "random.h"

#include <string.h>

int igraph_motifs_randesu(const igraph_t *graph, igraph_vector_t *hist, 
			  int size, igraph_vector_t *cut_prob) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_i_adjlist_t allneis;
  igraph_vector_t *neis;
  long int father;
  long int i, j, s;
  long int motifs=0;

  igraph_vector_t vids;		/* this is G */
  igraph_vector_t adjverts;	/* this is V_E */
  igraph_stack_t stack;		/* this is S */
  long int *added;
  char *subg;
  
  long int histlen;
  unsigned int *arr_idx, *arr_code;
  int class=0, code=0;
  unsigned char mul, idx;

  igraph_vector_t deg;

  IGRAPH_VECTOR_INIT_FINALLY(&deg, no_of_nodes);
  igraph_degree(graph, &deg, IGRAPH_VS_ALL(graph), IGRAPH_OUT, 1);
  
  if (size != 3 && size != 4) {
    IGRAPH_ERROR("Only 3 and 4 vertex motifs are implemented",
		 IGRAPH_EINVAL);
  }
  if (size==3) {
    mul=3;
    if (igraph_is_directed(graph)) {
      histlen=16;
      arr_idx=igraph_i_isoclass_3_idx;
      arr_code=igraph_i_isoclass2_3;
    } else {
      histlen=4;
      arr_idx=igraph_i_isoclass_3u_idx;
      arr_code=igraph_i_isoclass2_3u;
    }
  } else {
    mul=4;
    if (igraph_is_directed(graph)) {
      histlen=218;
      arr_idx=igraph_i_isoclass_4_idx;
      arr_code=igraph_i_isoclass2_4;
    } else {
      histlen=11;
      arr_idx=igraph_i_isoclass_4u_idx;
      arr_code=igraph_i_isoclass2_4u;
    }
  }

  IGRAPH_CHECK(igraph_vector_resize(hist, histlen));
  igraph_vector_null(hist);
  
  added=Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot find motifs", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);

  subg=Calloc(no_of_nodes, char);
  if (subg==0) {
    IGRAPH_ERROR("Cannot find motifs", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, subg);

  igraph_i_adjlist_init(graph, &allneis, IGRAPH_ALL);
  IGRAPH_FINALLY(igraph_i_adjlist_destroy, &allneis);  

  IGRAPH_VECTOR_INIT_FINALLY(&vids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjverts, 0);
  IGRAPH_CHECK(igraph_stack_init(&stack, 0));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);

  for (father=0; father<no_of_nodes; father++) {
    long int level;
    
    /* init G */
    igraph_vector_clear(&vids); level=0;
    IGRAPH_CHECK(igraph_vector_push_back(&vids, father));
    subg[father]=1; added[father] += 1; level += 1;
    
    /* init V_E */
    igraph_vector_clear(&adjverts);
    neis=igraph_i_adjlist_get(&allneis, father);
    s=igraph_vector_size(neis);
    for (i=0; i<s; i++) {
      long int nei=VECTOR(*neis)[i];
      if (!added[nei] && nei > father) {
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, father));
      }
      added[nei] += 1;
    }
    
    /* init S */
    igraph_stack_clear(&stack);

    while (level != 1 || !igraph_vector_empty(&adjverts)) {
      real_t cp=VECTOR(*cut_prob)[level];
      
      if (level==size-1) {
	s=igraph_vector_size(&adjverts)/2;
	for (i=0; i<s; i++) {
	  long int k, s2;
	  long int last;

	  if (cp!=0 && RNG_UNIF01() < cp) { continue; }
	  motifs+=1;
	  
	  last=VECTOR(adjverts)[2*i];
	  IGRAPH_CHECK(igraph_vector_push_back(&vids, last));
	  subg[last]=size;

	  code=0; idx=0;
	  for (k=0; k<size; k++) {
	    long int from=VECTOR(vids)[k];
 	    neis=igraph_i_adjlist_get(&allneis, from);
	    s2=VECTOR(deg)[from];
	    for (j=0; j<s2; j++) {
	      long int nei=VECTOR(*neis)[j];
	      if (subg[nei] && k != subg[nei]-1) {
		idx=mul*k+(subg[nei]-1);
		code |= arr_idx[idx];
	      }
	    }
	  }

	  igraph_vector_pop_back(&vids);
	  subg[last]=0;
	  VECTOR(*hist)[arr_code[code]] += 1;
	}
      }

      /* can we step down? */
      if (level < size-1 && 
	  !igraph_vector_empty(&adjverts) && 
	  (cp==0 || RNG_UNIF01() > cp)) {
	/* yes, step down */
	long int neifather=igraph_vector_pop_back(&adjverts);
	long int nei=igraph_vector_pop_back(&adjverts);
	IGRAPH_CHECK(igraph_vector_push_back(&vids, nei));
	subg[nei] = level+1; added[nei] += 1; level += 1;

	IGRAPH_CHECK(igraph_stack_push(&stack, neifather));
	IGRAPH_CHECK(igraph_stack_push(&stack, nei));
	IGRAPH_CHECK(igraph_stack_push(&stack, level));
	
	neis=igraph_i_adjlist_get(&allneis, nei);
	s=igraph_vector_size(neis);
	for (i=0; i<s; i++) {
	  long int nei2=VECTOR(*neis)[i];
	  if (!added[nei2] && nei2 > father) {
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei2));
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	  }
	  added[nei2] += 1;
	}
      } else {
	/* no, step back */
	long int nei, neifather;
	while (!igraph_stack_empty(&stack) &&
	       level==igraph_stack_top(&stack)-1) {
	  igraph_stack_pop(&stack);
	  nei=igraph_stack_pop(&stack);
	  neifather=igraph_stack_pop(&stack);
	  igraph_vector_push_back(&adjverts, nei);
	  igraph_vector_push_back(&adjverts, neifather);
	}

	nei=igraph_vector_pop_back(&vids);
	subg[nei]=0; added[nei] -= 1; level -= 1;
	neis=igraph_i_adjlist_get(&allneis, nei);
	s=igraph_vector_size(neis);
	for (i=0; i<s; i++) {
	  added[ (long int) VECTOR(*neis)[i] ] -= 1;
	}
	while (!igraph_vector_empty(&adjverts) && 
	       igraph_vector_tail(&adjverts)==nei) {
	  igraph_vector_pop_back(&adjverts);
	  igraph_vector_pop_back(&adjverts);
	}
      }
      
    } /* while */

    /* clear the added vector */
    added[father] -= 1;
    subg[father] = 0;
    neis=igraph_i_adjlist_get(&allneis, father);
    s=igraph_vector_size(neis);
    for (i=0; i<s; i++) {
      added[ (long int) VECTOR(*neis)[i] ] -= 1;
    }

  } /* for father */

  Free(added);
  Free(subg);
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&adjverts);
  igraph_i_adjlist_destroy(&allneis);
  igraph_stack_destroy(&stack);
  igraph_vector_destroy(&deg);
  IGRAPH_FINALLY_CLEAN(7);
  return 0;
}

int igraph_motifs_randesu_estimate(const igraph_t *graph, integer_t *est,
				   int size, igraph_vector_t *cut_prob, 
				   integer_t sample_size, 
				   igraph_vector_t *parsample) {

  long int no_of_nodes=igraph_vcount(graph);
  igraph_vector_t neis;
    
  igraph_vector_t vids;		/* this is G */
  igraph_vector_t adjverts;	/* this is V_E */
  igraph_stack_t stack;		/* this is S */
  long int *added;
  igraph_vector_t sample;
  long int sam;
  long int i;

  added=Calloc(no_of_nodes, long int);
  if (added==0) {
    IGRAPH_ERROR("Cannot find motifs", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, added);

  IGRAPH_VECTOR_INIT_FINALLY(&vids, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&adjverts, 0);
  IGRAPH_CHECK(igraph_stack_init(&stack, 0));
  IGRAPH_FINALLY(igraph_stack_destroy, &stack);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);

  if (igraph_vector_size(parsample)==0) {
    IGRAPH_VECTOR_INIT_FINALLY(&sample, 0);  
    IGRAPH_CHECK(igraph_random_sample(&sample, 0, no_of_nodes-1, sample_size));
  } else {
    sample=*parsample;
    sample_size=igraph_vector_size(&sample);
  }

  *est=0;

  for (sam=0; sam<sample_size; sam++) {
    long int father=VECTOR(sample)[sam];
    long int level, s;

    /* init G */
    igraph_vector_clear(&vids); level=0;
    IGRAPH_CHECK(igraph_vector_push_back(&vids, father));
    added[father] += 1; level += 1;
    
    /* init V_E */
    igraph_vector_clear(&adjverts);
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
    s=igraph_vector_size(&neis);
    for (i=0; i<s; i++) {
      long int nei=VECTOR(neis)[i];
      if (!added[nei]) {
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	IGRAPH_CHECK(igraph_vector_push_back(&adjverts, father));
      }
      added[nei] += 1;
    }
    
    /* init S */
    igraph_stack_clear(&stack);

    while (level != 1 || !igraph_vector_empty(&adjverts)) {
      real_t cp=VECTOR(*cut_prob)[level];
      
      if (level==size-1) {
	s=igraph_vector_size(&adjverts)/2;
	for (i=0; i<s; i++) {
	  if (cp!=0 && RNG_UNIF01() < cp) { continue; }
	  (*est) += 1;
	}
      }

      if (level < size-1 && 
	  !igraph_vector_empty(&adjverts) && 
	  (cp==0 || RNG_UNIF01() > cp)) {
	/* yes, step down */
	long int neifather=igraph_vector_pop_back(&adjverts);
	long int nei=igraph_vector_pop_back(&adjverts);
	IGRAPH_CHECK(igraph_vector_push_back(&vids, nei));
	added[nei] += 1; level += 1;

	IGRAPH_CHECK(igraph_stack_push(&stack, neifather));
	IGRAPH_CHECK(igraph_stack_push(&stack, nei));
	IGRAPH_CHECK(igraph_stack_push(&stack, level));

	IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
	s=igraph_vector_size(&neis);
	for (i=0; i<s; i++) {
	  long int nei2=VECTOR(neis)[i];
	  if (!added[nei2]) {
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei2));
	    IGRAPH_CHECK(igraph_vector_push_back(&adjverts, nei));
	  }
	  added[nei2] += 1;
	}
      } else {
	/* no, step back */
	long int nei, neifather;
	while (!igraph_stack_empty(&stack) &&
	       level==igraph_stack_top(&stack)-1) {
	  igraph_stack_pop(&stack);
	  nei=igraph_stack_pop(&stack);
	  neifather=igraph_stack_pop(&stack);
	  igraph_vector_push_back(&adjverts, nei);
	  igraph_vector_push_back(&adjverts, neifather);
	}

	nei=igraph_vector_pop_back(&vids);
	added[nei] -= 1; level -= 1;
	IGRAPH_CHECK(igraph_neighbors(graph, &neis, nei, IGRAPH_ALL));
	s=igraph_vector_size(&neis);
	for (i=0; i<s; i++) {
	  added[ (long int) VECTOR(neis)[i] ] -= 1;
	}
	while (!igraph_vector_empty(&adjverts) && 
	       igraph_vector_tail(&adjverts)==nei) {
	  igraph_vector_pop_back(&adjverts);
	  igraph_vector_pop_back(&adjverts);
	}
      }
      
    } /* while */

    /* clear the added vector */
    added[father] -= 1;
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, father, IGRAPH_ALL));
    s=igraph_vector_size(&neis);
    for (i=0; i<s; i++) {
      added[ (long int) VECTOR(neis)[i] ] -= 1;
    }

  } /* for father */

  (*est) *= ((double)no_of_nodes/sample_size);
  (*est) /= size;

  if (igraph_vector_size(parsample)==0) {
    igraph_vector_destroy(&sample);
    IGRAPH_FINALLY_CLEAN(1);
  }

  Free(added);
  igraph_vector_destroy(&vids);
  igraph_vector_destroy(&adjverts);
  igraph_stack_destroy(&stack);
  igraph_vector_destroy(&neis);
  IGRAPH_FINALLY_CLEAN(5);
  return 0;
}
