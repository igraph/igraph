/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <Gabor.Csardi@unil.ch>
   UNIL DGM, Rue de Bugnon 27, CH-1005 Lausanne, Switzerland
   
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

#include "igraph_mixing.h"
#include "igraph_interface.h"

int igraph_assortativity_nominal(const igraph_t *graph, 
				 const igraph_vector_long_t *types,
				 igraph_real_t *res,
				 igraph_bool_t directed) {

  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int no_of_types;
  igraph_vector_t ai, bi, eii;
  long int e, i;
  igraph_real_t sumaibi=0.0, sumeii=0.0;

  if (igraph_vector_long_size(types) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `types' vector length", IGRAPH_EINVAL);
  }

  directed = directed && igraph_is_directed(graph);

  no_of_types=igraph_vector_long_max(types)+1;
  IGRAPH_VECTOR_INIT_FINALLY(&ai, no_of_types);
  IGRAPH_VECTOR_INIT_FINALLY(&bi, no_of_types);
  IGRAPH_VECTOR_INIT_FINALLY(&eii, no_of_types);

  for (e=0; e<no_of_edges; e++) {
    long int from=IGRAPH_FROM(graph, e);
    long int to=IGRAPH_TO(graph, e);
    long int from_type = VECTOR(*types)[from];
    long int to_type = VECTOR(*types)[to];
    
    VECTOR(ai)[from_type] += 1;
    VECTOR(bi)[to_type] += 1;
    if (from_type == to_type) {
      VECTOR(eii)[from_type] += 1;
    }
    if (!directed) {
      if (from_type == to_type) {
	VECTOR(eii)[from_type] += 1;
      }  
      VECTOR(ai)[to_type] += 1;
      VECTOR(bi)[from_type] += 1;
    }
  }
  
  for (i=0; i<no_of_types; i++) {
    sumaibi += (VECTOR(ai)[i]/no_of_edges) * (VECTOR(bi)[i]/no_of_edges);
    sumeii  += (VECTOR(eii)[i]/no_of_edges);
  }
  
  if (!directed) { 
    sumaibi /= 4.0;
    sumeii  /= 2.0;
  }

  *res = (sumeii - sumaibi) / (1.0 - sumaibi);
  
  igraph_vector_destroy(&eii);
  igraph_vector_destroy(&bi);
  igraph_vector_destroy(&ai);
  IGRAPH_FINALLY_CLEAN(3);
  
  return 0;
}
    
int igraph_assortativity(const igraph_t *graph,
			 const igraph_vector_t *types1,
			 const igraph_vector_t *types2,
			 igraph_real_t *res,
			 igraph_bool_t directed) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  long int e;

  directed = directed && igraph_is_directed(graph);

  if (!directed && types2) { 
    IGRAPH_WARNING("Only `types1' is used for undirected case");
  }

  if (igraph_vector_size(types1) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `types1' vector length", IGRAPH_EINVAL);
  }

  if (types2 && igraph_vector_size(types2) != no_of_nodes) {
    IGRAPH_ERROR("Invalid `types2' vector length", IGRAPH_EINVAL);
  }

  if (!directed) {
    igraph_real_t num1=0.0, num2=0.0, den1=0.0;
  
    for (e=0; e<no_of_edges; e++) {
      long int from=IGRAPH_FROM(graph, e);
      long int to=IGRAPH_TO(graph, e);
      igraph_real_t from_type=VECTOR(*types1)[from];
      igraph_real_t to_type=VECTOR(*types1)[to];
      
      num1 += from_type * to_type;
      num2 += from_type + to_type;
      den1 += from_type * from_type + to_type * to_type;
    }
    
    num1 /= no_of_edges;
    den1 /= no_of_edges * 2;
    num2 /= no_of_edges * 2;
    num2 = num2 * num2;
    
    *res = (num1-num2) / (den1-num2);

  } else {
    igraph_real_t num1=0.0, num2=0.0, num3=0.0,
      den1=0.0, den2=0.0;
    igraph_real_t num, den;

    if (!types2) { types2=types1; }
    
    for (e=0; e<no_of_edges; e++) {
      long int from=IGRAPH_FROM(graph, e);
      long int to=IGRAPH_TO(graph, e);
      igraph_real_t from_type=VECTOR(*types1)[from];
      igraph_real_t to_type=VECTOR(*types2)[to];
      
      num1 += from_type * to_type;
      num2 += from_type;
      num3 += to_type;
      den1 += from_type * from_type;
      den2 += to_type * to_type;
    }
    
    num = num1 - num2*num3/no_of_edges;
    den = sqrt(den1 - num2*num2/no_of_edges) * 
      sqrt(den2 - num3*num3/no_of_edges);
    
    *res = num /den;
  }
      
  return 0;
}

int igraph_assortativity_degree(const igraph_t *graph,
				igraph_real_t *res, 
				igraph_bool_t directed) {

  directed = directed && igraph_is_directed(graph);

  if (directed) {
    igraph_vector_t indegree, outdegree;
    igraph_vector_init(&indegree, 0);
    igraph_vector_init(&outdegree, 0);
    igraph_degree(graph, &indegree, igraph_vss_all(), IGRAPH_IN, /*loops=*/ 1);
    igraph_degree(graph, &outdegree, igraph_vss_all(), IGRAPH_OUT, /*loops=*/ 1);
    igraph_vector_add_constant(&indegree, -1);
    igraph_vector_add_constant(&outdegree, -1);
    igraph_assortativity(graph, &outdegree, &indegree, res, /*directed=*/ 1);
    igraph_vector_destroy(&indegree);
    igraph_vector_destroy(&outdegree);
  } else {
    igraph_vector_t degree;
    igraph_vector_init(&degree, 0);
    igraph_degree(graph, &degree, igraph_vss_all(), IGRAPH_ALL, /*loops=*/ 1);
    igraph_vector_add_constant(&degree, -1);
    igraph_assortativity(graph, &degree, 0, res, /*directed=*/ 0);
    igraph_vector_destroy(&degree);
  }
  
  return 0;
}
