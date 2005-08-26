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

int igraph_diameter(igraph_t *graph, integer_t *res, 
		    bool_t directed, bool_t unconn) {

  long int no_of_nodes=igraph_vcount(graph);
  long int i, j;
  long int *already_added;
  long int nodes_reached;

  dqueue_t q;
  vector_t neis;
  integer_t dirmode;

  *res=0;  
  if (directed) { dirmode=1; } else { dirmode=3; }
  already_added=Calloc(no_of_nodes, long int);
  dqueue_init(&q, 100);
  vector_init(&neis, 0);
  
  for (i=0; i<no_of_nodes; i++) {
    nodes_reached=1;
    dqueue_push(&q, i);
    dqueue_push(&q, 0);
    already_added[i]=i+1;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      if (actdist>*res) { *res=actdist; }
      
      igraph_neighbors(graph, &neis, actnode, dirmode);
      for (j=0; j<vector_size(&neis); j++) {
	long int neighbor=VECTOR(neis)[j];
	if (already_added[neighbor] == i+1) { continue; }
	already_added[neighbor]=i+1;
	nodes_reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    } /* while !dqueue_empty */
    
    /* not connected, return largest possible */
    if (nodes_reached != no_of_nodes && !unconn) {
      *res=no_of_nodes;
      break;
    }
  } /* for i<no_of_nodes */
  
  /* clean */
  Free(already_added);
  dqueue_destroy(&q);
  vector_destroy(&neis);

  return 0;
}

int igraph_minimum_spanning_tree_unweighted(igraph_t *graph, igraph_t *mst) {
  /* TODO */
  return 0;
}

int igraph_minimum_spanning_tree_prim(igraph_t *graph, igraph_t *mst,
				      vector_t *weights) {
  /* TODO */
  return 0;
}

int igraph_closeness(igraph_t *graph, vector_t *res, vector_t *vids, 
		     integer_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  long int *already_counted;
  long int i, j;
  long int nodes_reached;

  dqueue_t q;
  
  long int nodes_to_calc=vector_size(vids);
  vector_t tmp;

  vector_resize(res, nodes_to_calc);
  vector_null(res);
  
  already_counted=Calloc(no_of_nodes, long int);
  vector_init(&tmp, 0);
  dqueue_init(&q, 100);

  for (i=0; i<nodes_to_calc; i++) {
    dqueue_push(&q, VECTOR(*vids)[i]);
    dqueue_push(&q, 0);
    nodes_reached=1;
    already_counted[(long int)VECTOR(*vids)[i]]=i+1;
    
    while (!dqueue_empty(&q)) {
      long int act=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      VECTOR(*res)[i] += actdist;

      igraph_neighbors(graph, &tmp, act, mode);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (already_counted[neighbor] == i+1) { continue; }
	already_counted[neighbor] = i+1;
	nodes_reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    }
    VECTOR(*res)[i] += (no_of_nodes * (no_of_nodes-nodes_reached));
    VECTOR(*res)[i] = (no_of_nodes-1) / VECTOR(*res)[i];
  }
  
  /* Clean */
  dqueue_destroy(&q);
  vector_destroy(&tmp);
  Free(already_counted);

  return 0;
}

int igraph_shortest_paths(igraph_t *graph, matrix_t *res, 
			  vector_t *from, integer_t mode) {

/* SEXP REST_shortest_paths(SEXP interface, SEXP graph, SEXP from, SEXP pmode) { */
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_from=vector_size(from);
  long int *already_counted;
  
  dqueue_t q;

  long int i, j;
  vector_t tmp;

  already_counted=Calloc(no_of_nodes, long int);
  vector_init(&tmp, 0);
  matrix_resize(res, no_of_from, no_of_nodes);
  matrix_null(res);

  dqueue_init(&q, 100);

  for (i=0; i<no_of_from; i++) {
    long int reached=1;
    dqueue_push(&q, VECTOR(*from)[i]);
    dqueue_push(&q, 0);
    already_counted[ (long int) VECTOR(*from)[i] ] = i+1;
    
    while (!dqueue_empty(&q)) {
      long int act=dqueue_pop(&q);
      long int actdist=dqueue_pop(&q);
      MATRIX(*res, i, act)=actdist;
      
      igraph_neighbors(graph, &tmp, act, mode);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (already_counted[neighbor] == i+1) { continue; }
	already_counted[neighbor] = i+1;
	reached++;
	dqueue_push(&q, neighbor);
	dqueue_push(&q, actdist+1);
      }
    }
    /* Plus the unreachable nodes */
    j=0;
    while (reached < no_of_nodes) {
      if (MATRIX(*res, i, j) == 0 && j != VECTOR(*from)[i]) {
	MATRIX(*res, i, j)=no_of_nodes;
	reached++;
      }
      j++;
    }
  }

  /* Clean */
  vector_destroy(&tmp);
  Free(already_counted);
  dqueue_destroy(&q);

  return 0;
}

int igraph_get_shortest_paths(igraph_t *graph, vector_t *res,
			      vector_t *from, integer_t mode) {
  /* TODO */
  return 0;
}

int igraph_subcomponent(igraph_t *graph, vector_t *res, real_t vertex, 
			integer_t mode) {

  long int no_of_nodes=igraph_vcount(graph);
  dqueue_t q;
  char *already_added;
  long int i,j;
  vector_t tmp;
  
  vector_init(&tmp, 0);
  already_added=Calloc(no_of_nodes, char);
  dqueue_init(&q, 100);
  
  dqueue_push(&q, vertex);
  vector_push_back(res, vertex);
  already_added[(long int)vertex]=1;
  
  while (!dqueue_empty(&q)) {
    long int actnode=dqueue_pop(&q);

    igraph_neighbors(graph, &tmp, actnode, mode);
    for (i=0; i<vector_size(&tmp); i++) {
      long int neighbor=VECTOR(tmp)[i];
      
      if (already_added[neighbor]) { continue; }
      already_added[neighbor]=1;
      vector_push_back(res, neighbor);
      dqueue_push(&q, neighbor);
    }
  }

  dqueue_destroy(&q);
  vector_destroy(&tmp);

  return 0;
}

int igraph_subgraph(igraph_t *graph, igraph_t *res, integer_t  mode) {
  /* TODO */
  return 0;
}

int igraph_simplify(igraph_t *graph, bool_t remove_loops, 
		    bool_t remove_multiple) {
  /* TODO */
  return 0;
}

int igraph_betweenness (igraph_t *graph, vector_t *res, vector_t *vids, 
			bool_t directed) {

  long int no_of_nodes=igraph_vcount(graph);
  dqueue_t q;
  long int *distance;
  long int *nrgeo;
  double *tmpscore;
  stack_t stack;
  long int source;
  long int j;
  vector_t tmp;
  integer_t modein, modeout;
  real_t *tmpres;

  if (directed) { modeout=1; modein=2; } else { modeout=modein=3; }

  distance=Calloc(no_of_nodes, long int);
  nrgeo=Calloc(no_of_nodes, long int);
  tmpscore=Calloc(no_of_nodes, double);
  tmpres=Calloc(no_of_nodes, real_t);

  vector_init(&tmp, 0);
  vector_resize(res, vector_size(vids));
  vector_null(res);

  dqueue_init(&q, 100);
  stack_init(&stack, no_of_nodes);

  /* here we go */
  
  for (source=0; source<no_of_nodes; source++) {
    
    memset(distance, 0, no_of_nodes*sizeof(long int));
    memset(nrgeo, 0, no_of_nodes*sizeof(long int));
    memset(tmpscore, 0, no_of_nodes*sizeof(double));
    stack_clear(&stack); /* it should be empty anyway... */
    
    dqueue_push(&q, source);
    nrgeo[source]=1;
    distance[source]=0;
    
    while (!dqueue_empty(&q)) {
      long int actnode=dqueue_pop(&q);

      igraph_neighbors(graph, &tmp, actnode, modeout);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (nrgeo[neighbor] != 0) {
	  /* we've already seen this node, another shortest path? */
	  if (distance[neighbor]==distance[actnode]+1) {
	    nrgeo[neighbor]++;
	  }
	} else {
	  /* we haven't seen this node yet */
	  nrgeo[neighbor]++;
	  distance[neighbor]=distance[actnode]+1;
	  dqueue_push(&q, neighbor);
	  stack_push(&stack, neighbor);
	}
      }
    } /* while !dqueue_empty */
    
    /* Ok, we've the distance of each node and also the number of
       shortest paths to them. Now we do an inverse search, starting
       with the farthest nodes. */
    while (!stack_empty(&stack)) {
      long int actnode=stack_pop(&stack);
      long int friends=0;
      if (distance[actnode]<=1) { continue; } /* skip source node */
      
      /* search for the neighbors on the geodesics */
      igraph_neighbors(graph, &tmp, actnode, modein);
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (distance[neighbor]==distance[actnode]-1) { friends++; }
      }

      /* set the temporary score of the friends */
      for (j=0; j<vector_size(&tmp); j++) {
	long int neighbor=VECTOR(tmp)[j];
	if (distance[neighbor]==distance[actnode]-1) {
	  tmpscore[neighbor] += (tmpscore[actnode]+1)/friends;
	}
      }
    }
    
    /* Ok, we've the scores for this source */
    for (j=0; j<vector_size(vids); j++) {
      long int node=VECTOR(*vids)[j];
      VECTOR(*res)[node] += tmpscore[node];
      tmpscore[node] = 0.0; /* in case a node is in vids multiple times */
    }

  } /* for source < no_of_nodes */

  /* divide by 2 for undirected graph */
  if (!directed || !igraph_is_directed(graph)) {
    for (j=0; j<vector_size(res); j++) {
      VECTOR(*res)[j] /= 2.0;
    }
  }
  
  /* clean  */
  Free(distance);
  Free(nrgeo);
  Free(tmpscore);
  dqueue_destroy(&q);
  stack_destroy(&stack);
  vector_destroy(&tmp);

  return 0;
}

int igraph_edge_betweenness (igraph_t *graph, vector_t *result, 
			     bool_t directed) {
  /* TODO */
  return 0;
}


/* SEXP REST_edge_betweenness (SEXP interface, SEXP graph, SEXP pdirected) { */

/*   REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph); */

/*   SEXP result; */
/*   SEXP dim; */
  
/*   long int no_of_nodes; */
/*   long int no_of_edges; */
/*   dqueue_t q; */
/*   long int *distance; */
/*   long int *nrgeo; */
/*   double *tmpscore; */
/*   stack_t stack; */
/*   long int source; */
/*   long int j; */

/*   SEXP tmp, modein, modeout; */

/*   no_of_nodes=R(ptrtable.vcount(interface, graph)); */
/*   no_of_edges=R(ptrtable.ecount(interface, graph)); */

/*   if (LOGICAL(pdirected)[0]) { */
/*     PROTECT(modeout=ScalarString(CREATE_STRING_VECTOR("out"))); */
/*     PROTECT(modein =ScalarString(CREATE_STRING_VECTOR("in"))); */
/*   } else { */
/*     PROTECT(modeout=ScalarString(CREATE_STRING_VECTOR("all"))); */
/*     PROTECT(modein =ScalarString(CREATE_STRING_VECTOR("all"))); */
/*   }         */
  
/*   distance=(long int*) R_alloc(no_of_nodes, sizeof(long int)); */
/*   nrgeo=(long int*) R_alloc(no_of_nodes, sizeof(long int)); */
/*   tmpscore=(double*) R_alloc(no_of_nodes, sizeof(double)); */
/*   /\* others memsetting later *\/ */

/*   PROTECT(result=NEW_NUMERIC(no_of_nodes*no_of_nodes)); */
/*   memset(REAL(result), 0, no_of_nodes*no_of_nodes*sizeof(double)); */
/*   PROTECT(dim=NEW_INTEGER(2)); */
/*   INTEGER(dim)[0]=no_of_nodes; */
/*   INTEGER(dim)[1]=no_of_nodes; */
/*   SET_DIM(result, dim); */

/*   dqueue_init(&q, 100); */
/*   stack_init(&stack, no_of_nodes); */

/*   /\* here we go *\/ */
  
/*   for (source=1; source<=no_of_nodes; source++) { */
    
/*     memset(distance, 0, no_of_nodes*sizeof(long int)); */
/*     memset(nrgeo, 0, no_of_nodes*sizeof(long int)); */
/*     memset(tmpscore, 0, no_of_nodes*sizeof(double)); */
/*     stack_clear(&stack); /\* it should be empty anyway... *\/ */
    
/*     dqueue_push(&q, source); */
/*     nrgeo[source-1]=1; */
/*     distance[source-1]=0; */
    
/*     while (!dqueue_empty(&q)) { */
/*       long int actnode=dqueue_pop(&q); */
      
/*       tmp=ptrtable.neighbors(interface, graph, actnode, modeout); */
/*       for (j=0; j<GET_LENGTH(tmp); j++) { */
/* 	long int neighbor=REAL(tmp)[j]; */
/* 	if (nrgeo[neighbor-1] != 0) { */
/* 	  /\* we've already seen this node, another shortest path? *\/ */
/* 	  if (distance[neighbor-1]==distance[actnode-1]+1) { */
/* 	    nrgeo[neighbor-1]++; */
/* 	  } */
/* 	} else { */
/* 	  /\* we haven't seen this node yet *\/ */
/* 	  nrgeo[neighbor-1]++; */
/* 	  distance[neighbor-1]=distance[actnode-1]+1; */
/* 	  dqueue_push(&q, neighbor); */
/* 	  stack_push(&stack, neighbor); */
/* 	} */
/*       } */
/*     } /\* while !dqueue_empty *\/ */
    
/*     /\* Ok, we've the distance of each node and also the number of  */
/*        shortest paths to them. Now we do an inverse search, starting */
/*        with the farthest nodes. *\/ */
/*     while (!stack_empty(&stack)) { */
/*       long int actnode=stack_pop(&stack); */
/*       long int friends=0; */
/*       if (distance[actnode-1]<1) { continue; } /\* skip source node *\/ */
      
/*       /\* search for the neighbors on the geodesics *\/ */
/*       tmp=ptrtable.neighbors(interface, graph, actnode, modein); */
/*       for (j=0; j<GET_LENGTH(tmp); j++) { */
/* 	long int neighbor=REAL(tmp)[j]; */
/* 	if (distance[neighbor-1]==distance[actnode-1]-1) { friends++; } */
/*       } */

/*       /\* set the temporary score of the friends *\/ */
/*       for (j=0; j<GET_LENGTH(tmp); j++) { */
/* 	long int neighbor=REAL(tmp)[j]; */
/* 	if (distance[neighbor-1]==distance[actnode-1]-1) { 	   */
/* 	  tmpscore[neighbor-1] += (tmpscore[actnode-1]+1)/friends; */
/* 	  RMATRIX(result,neighbor, actnode) += (tmpscore[actnode-1]+1)/friends; */
/* 	} */
/*       } */
/*     } */
    
/*     /\* Ok, we've the scores for this source *\/ */
/*     R_CheckUserInterrupt(); */

/*   } /\* for source <= no_of_nodes *\/ */
  
  
/*   /\* clean and return *\/ */
/*   dqueue_destroy(&q); */
/*   stack_destroy(&stack); */
/*   UNPROTECT(4); */
/*   return result; */

/* } */

/* SEXP REST_get_shortest_paths(SEXP interface, SEXP graph, SEXP pfrom, SEXP pmode) { */
  
/*   REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph); */

/*   SEXP result; */
  
/*   long int no_of_nodes; */
/*   long int *father; */
/*   long int from; */
  
/*   dqueue_t q; */

/*   long int i, j; */
/*   SEXP tmp; */

/*   no_of_nodes=R(ptrtable.vcount(interface, graph)); */
/*   from=R(pfrom); */
  
/*   father=(long int*) R_alloc(no_of_nodes, sizeof(long int)); */
/*   memset(father, 0, no_of_nodes*sizeof(long int)); */

/*   dqueue_init(&q, 100); */

/*   dqueue_push(&q, from); */
/*   father[ from-1 ] = from; */
  
/*   while (!dqueue_empty(&q)) { */
/* 	  long int act=dqueue_pop(&q); */
	  
/* 	  tmp=ptrtable.neighbors(interface,graph, act, pmode); */
/* 	  for (j=0; j<GET_LENGTH(tmp); j++) { */
/* 		  long int neighbor=REAL(tmp)[j]; */
/* 		  if (father[neighbor-1] != 0) { continue; } */
/* 		  father[neighbor-1] = act; */
/* 		  dqueue_push(&q, neighbor); */
/* 	  } */
/*   } */
    
/*   /\* Create result object *\/ */
/*   PROTECT(result=NEW_LIST(no_of_nodes)); */
/*   for (j=0; j<no_of_nodes; j++) { */
/* 	  if (father[j]==0) { */
/* 		  SET_VECTOR_ELT(result, j, NEW_NUMERIC(0)); */
/* 	  } else { */
/* 		  long int act=j+1; */
/* 		  while (father[act-1] != act) { */
/* 			  dqueue_push(&q, father[act-1]); */
/* 			  act=father[act-1]; */
/* 		  } */
/* 		  act=dqueue_size(&q); */
/* 		  SET_VECTOR_ELT(result, j, NEW_NUMERIC(act+1)); */
/* 		  REAL(VECTOR_ELT(result, j))[0] = j+1; */
/* 		  for (i=1; i<=act; i++) { */
/* 			  REAL(VECTOR_ELT(result, j))[i] = dqueue_pop(&q); */
/* 		  } */
/* 	  } */
/*   } */
  
/*   /\* Clean *\/ */
/*   dqueue_destroy(&q); */
/*   UNPROTECT(1); */
/*   return result; */
/* } */



/* SEXP REST_minimum_spanning_tree_unweighted(SEXP interface,  */
/* 					   SEXP graph) { */
  
/*   REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph); */
  
/*   SEXP result; */
/*   long int no_of_nodes; */
/*   char *already_added; */
  
/*   dqueue_t q; */
/*   dqueue_t edges; */
/*   SEXP mode, tmp; */
/*   long int i, j; */

/*   no_of_nodes=R(ptrtable.vcount(interface, graph)); */
/*   already_added=(char*) R_alloc(no_of_nodes, sizeof(char)); */
/*   memset(already_added, 0, no_of_nodes*sizeof(char)); */

/*   PROTECT(mode=NEW_CHARACTER(1)); */
/*   SET_STRING_ELT(mode, 0, CREATE_STRING_VECTOR("all")); */

/*   dqueue_init(&q, 100); */
/*   dqueue_init(&edges, (no_of_nodes-1)*2); */
  
/*   for (i=1; i<=no_of_nodes; i++) { */
/*     if (already_added[i-1]>0) { continue; } */

/*     already_added[i-1]=1; */
/*     dqueue_push(&q, i); */
/*     while (! dqueue_empty(&q)) { */
/*       long int act_node=dqueue_pop(&q); */
/*       tmp=ptrtable.neighbors(interface, graph, act_node, mode); */
/*       for (j=0; j<GET_LENGTH(tmp); j++) { */
/* 	long int neighbor=REAL(tmp)[j]; */
/* 	if (already_added[neighbor-1]==0) { */
/* 	  already_added[neighbor-1]=1; */
/* 	  dqueue_push(&q, neighbor); */
/* 	  dqueue_push(&edges, act_node); */
/* 	  dqueue_push(&edges, neighbor); */
/* 	} */
/*       } */
/*     } */
/*     R_CheckUserInterrupt(); */
/*   } */
  
/*   /\* Save the result *\/ */
  
/*   dqueue_destroy(&q); */
/*   j=dqueue_size(&edges); */
/*   PROTECT(result=NEW_NUMERIC(j)); */
/*   for (i=0; i<j; i++) { */
/*     REAL(result)[i] = dqueue_pop(&edges); */
/*   } */
    
/*   dqueue_destroy(&edges); */
/*   UNPROTECT(2); */
/*   return result; */
/* } */

/* SEXP REST_minimum_spanning_tree_prim(SEXP interface,  */
/* 				     SEXP graph) { */

/*   REST_i_ptrtable_t ptrtable = REST_i_getptrtable(graph); */

/*   SEXP result; */
/*   long int no_of_nodes; */
/*   char *already_added; */

/*   d_indheap_t heap; */
/*   dqueue_t edges; */

/*   SEXP mode, weight, tmp, tmp2; */
/*   long int i, j; */
  
/*   no_of_nodes=R(ptrtable.vcount(interface, graph)); */
/*   already_added=(char*) R_alloc(no_of_nodes, sizeof(char)); */
/*   memset(already_added, 0, no_of_nodes*sizeof(char)); */
  
/*   PROTECT(mode=ScalarString(CREATE_STRING_VECTOR("all"))); */
/*   PROTECT(weight=ScalarString(CREATE_STRING_VECTOR("weight"))); */
  
/*   d_indheap_init(&heap, 0); */
/*   dqueue_init(&edges, (no_of_nodes-1)*2); */
  
/*   for (i=1; i<=no_of_nodes; i++) { */
/*     if (already_added[i-1]>0) { continue; } */

/*     already_added[i-1]=1; */
/*     tmp=ptrtable.neighbors(interface, graph, i, mode); */
/*     tmp2=ptrtable.get_edge_attribute(interface, graph, weight,  */
/* 				     ScalarReal((double)i), */
/* 				     NULL_USER_OBJECT); */
/*     /\* add all edges of the first vertex *\/ */
/*     for (j=0; j<GET_LENGTH(tmp); j++) { */
/*       long int neighbor=REAL(tmp)[j]; */
/*       if (already_added[neighbor-1] == 0) { */
/* 	d_indheap_push(&heap, -REAL(tmp2)[j], i, neighbor); */
/*       } */
/*     } */

/*     while(! d_indheap_empty(&heap)) { */
/*       /\* Get minimal edge *\/ */
/*       long int from, to; */
/*       d_indheap_max_index(&heap, &from, &to); */

/*       /\* Erase it *\/ */
/*       d_indheap_delete_max(&heap); */
      
/*       /\* Does it point to a visited node? *\/ */
/*       if (already_added[to-1]==0) { */
/* 	already_added[to-1]=1; */
/* 	dqueue_push(&edges, from); */
/* 	dqueue_push(&edges, to); */
/* 	tmp=ptrtable.neighbors(interface, graph, to, mode); */
/* 	tmp2=ptrtable.get_edge_attribute(interface, graph, weight,  */
/* 					 ScalarReal((double)to),  */
/* 					 NULL_USER_OBJECT); */
/* 	/\* add all outgoing edges *\/ */
/* 	for (j=0; j<GET_LENGTH(tmp); j++) { */
/* 	  long int neighbor=REAL(tmp)[j]; */
/* 	  if (already_added[neighbor-1] == 0) { */
/* 	    d_indheap_push(&heap, -REAL(tmp2)[j], to, neighbor); */
/* 	  } */
/* 	} /\* for *\/ */
/*       } /\* if !already_added *\/ */
/*     } /\* while in the same component *\/ */
/*   } /\* for all nodes *\/ */

/*   d_indheap_destroy(&heap); */
/*   j=dqueue_size(&edges); */
/*   PROTECT(result=NEW_NUMERIC(j)); */
/*   for (i=0; i<j; i++) { */
/*     REAL(result)[i] = dqueue_pop(&edges); */
/*   }   */

/*   dqueue_destroy(&edges); */
/*   UNPROTECT(3); */
/*   return result; */
/* } */
