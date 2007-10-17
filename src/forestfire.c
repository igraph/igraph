/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2007  Gabor Csardi <csardi@rmki.kfki.hu>
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

typedef struct igraph_i_forest_fire_data_t {
  igraph_vector_t *inneis;
  igraph_vector_t *outneis;
  long int no_of_nodes;
} igraph_i_forest_fire_data_t;
  

void igraph_i_forest_fire_free(igraph_i_forest_fire_data_t *data) {
  long int i;
  for (i=0; i<data->no_of_nodes; i++) {
    igraph_vector_destroy(data->inneis+i);
    igraph_vector_destroy(data->outneis+i);
  }
}

int igraph_forest_fire_game(igraph_t *graph, igraph_integer_t nodes,
			    igraph_real_t fw_prob, igraph_real_t bw_factor,
			    igraph_integer_t pambs, igraph_bool_t directed) {
  
  igraph_vector_long_t visited;
  long int no_of_nodes=nodes, actnode, i;
  igraph_vector_t edges;
  igraph_vector_t *inneis, *outneis;
  igraph_i_forest_fire_data_t data;
  igraph_dqueue_t neiq;
  long int ambs=pambs;
  igraph_real_t param_geom_out=1-fw_prob;
  igraph_real_t param_geom_in=1-fw_prob*bw_factor;
  
  if (fw_prob < 0) {
    IGRAPH_ERROR("Forest fire model: 'fw_prob' should be between non-negative", 
		 IGRAPH_EINVAL);
  }
  if (bw_factor < 0) {
    IGRAPH_ERROR("Forest fire model: 'bw_factor' should be non-negative",
		 IGRAPH_EINVAL);
  }
  if (ambs < 0) {
    IGRAPH_ERROR("Number of ambassadors ('ambs') should be non-negative",
		 IGRAPH_EINVAL);
  }
  
  if (fw_prob == 0 || ambs == 0) {
    IGRAPH_WARNING("'fw_prob or ambs is zero, creating empty graph");
    IGRAPH_CHECK(igraph_empty(graph, nodes, directed));
    return 0;
  }
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  inneis=Calloc(no_of_nodes, igraph_vector_t);
  if (!inneis) {
    IGRAPH_ERROR("Cannot run forest fire model", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, inneis);
  outneis=Calloc(no_of_nodes, igraph_vector_t);
  if (!outneis) {
    IGRAPH_ERROR("Cannot run forest fire model", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, outneis);  
  data.inneis=inneis; 
  data.outneis=outneis;
  data.no_of_nodes=no_of_nodes;
  IGRAPH_FINALLY(igraph_i_forest_fire_free, &data);
  for (i=0; i<no_of_nodes; i++) {
    IGRAPH_CHECK(igraph_vector_init(inneis+i, 0));
    IGRAPH_CHECK(igraph_vector_init(outneis+i, 0));
  }  

  IGRAPH_CHECK(igraph_vector_long_init(&visited, no_of_nodes));
  IGRAPH_FINALLY(igraph_vector_long_destroy, &visited);
  IGRAPH_DQUEUE_INIT_FINALLY(&neiq, 10);

  RNG_BEGIN();

#define ADD_EDGE_TO(nei) \
      if (VECTOR(visited)[(nei)] != actnode+1) {                     \
	VECTOR(visited)[(nei)] = actnode+1;                          \
	IGRAPH_CHECK(igraph_dqueue_push(&neiq, nei));                \
	IGRAPH_CHECK(igraph_vector_push_back(&edges, actnode));      \
	IGRAPH_CHECK(igraph_vector_push_back(&edges, nei));          \
	IGRAPH_CHECK(igraph_vector_push_back(outneis+actnode, nei)); \
	IGRAPH_CHECK(igraph_vector_push_back(inneis+nei, actnode));  \
      }
  
  igraph_progress("Forest fire: ", 0.0, NULL);
  
  for (actnode=1; actnode < no_of_nodes; actnode++) {

    igraph_progress("Forest fire: ", 100.0*actnode/no_of_nodes, NULL);

    IGRAPH_ALLOW_INTERRUPTION();    
    
    /* We don't want to visit the current vertex */
    VECTOR(visited)[actnode] = actnode+1;

    /* Choose ambassador(s) */
    for (i=0; i<ambs; i++) {
      long int a=RNG_INTEGER(0, actnode-1);
      ADD_EDGE_TO(a);
    }
    
    while (!igraph_dqueue_empty(&neiq)) {
      long int actamb=igraph_dqueue_pop(&neiq);
      igraph_vector_t *outv=outneis+actamb;
      igraph_vector_t *inv=inneis+actamb;
      long int no_in=igraph_vector_size(inv);
      long int no_out=igraph_vector_size(outv);
      long int neis_out=RNG_GEOM(param_geom_out);
      long int neis_in=RNG_GEOM(param_geom_in);
      /* outgoing neighbors */
      if (neis_out >= no_out) {
	for (i=0; i<no_out; i++) {
	  long int nei=VECTOR(*outv)[i];
	  ADD_EDGE_TO(nei);
	}
      } else {
	long int oleft=no_out;
	for (i=0; i<neis_out && oleft > 0; ) {
	  long int which=RNG_INTEGER(0, oleft-1);
	  long int nei=VECTOR(*outv)[which];
	  VECTOR(*outv)[which] = VECTOR(*outv)[oleft-1];
	  VECTOR(*outv)[oleft-1] = nei;
	  if (VECTOR(visited)[nei] != actnode+1) {
	    ADD_EDGE_TO(nei);
	    i++;
	  }
	  oleft--;
	}
      }
      /* incoming neighbors */
      if (neis_in >= no_in) {
	for (i=0; i<no_in; i++) {
	  long int nei=VECTOR(*inv)[i];
	  ADD_EDGE_TO(nei);
	}
      } else {
	long int ileft=no_in;
	for (i=0; i<neis_in && ileft > 0; ) {
	  long int which=RNG_INTEGER(0, ileft-1);
	  long int nei=VECTOR(*inv)[which];
	  VECTOR(*inv)[which] = VECTOR(*inv)[ileft-1];
	  VECTOR(*inv)[ileft-1] = which;
	  if (VECTOR(visited)[nei] != actnode+1) {
	    ADD_EDGE_TO(nei);
	    i++;
	  }
	  ileft--;
	}
      }
      
    } /* while neiq not empty */

  } /* actnode < no_of_nodes */

#undef ADD_EDGE_TO  

  RNG_END();

  igraph_progress("Forest fire: ", 100.0, NULL);
  
  igraph_dqueue_destroy(&neiq);
  igraph_vector_long_destroy(&visited);
  igraph_i_forest_fire_free(&data);
  igraph_free(outneis);
  igraph_free(inneis);  
  IGRAPH_FINALLY_CLEAN(5);

  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}
