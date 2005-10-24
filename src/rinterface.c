/* -*- mode: C -*-  */
/* 
   IGraph library R interface.
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

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP R_matrix_to_SEXP(matrix_t *m) {

  SEXP result, dim; 
  
  PROTECT(result=NEW_NUMERIC(matrix_size(m)));
  matrix_copy_to(m, REAL(result));
  PROTECT(dim=NEW_INTEGER(2));
  INTEGER(dim)[0]=matrix_nrow(m);
  INTEGER(dim)[1]=matrix_ncol(m);
  SET_DIM(result, dim);

  UNPROTECT(2);
  return result;
}

SEXP R_strarray_to_SEXP(igraph_strarray_t *sa) {
  
  SEXP result;
  Rbyte *ptr;
  
  PROTECT(result=allocVector(RAWSXP, sa->sa_end-sa->sa_begin));
  memcpy(RAW(result), sa->sa_begin, sizeof(char)*(sa->sa_end-sa->sa_begin));
  for (ptr = RAW(result); ptr < RAW(result)+(sa->sa_end-sa->sa_begin);
       ptr++) {
    if (*ptr == '\0') {
      *ptr = '\n';
    }
  }  
  
  UNPROTECT(1);
  return result;
}

SEXP R_attributes_to_SEXP(igraph_attribute_list_t *al) {

  SEXP result;

  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, allocVector(RAWSXP, al->sa_end-al->sa_begin));
  memcpy(RAW(VECTOR_ELT(result, 0)), al->sa_begin, 
	 sizeof(char)*(al->sa_end-al->sa_begin));
  SET_VECTOR_ELT(result, 1, R_matrix_to_SEXP(&al->numattrs));

  UNPROTECT(1);
  return result;
}

int R_SEXP_to_matrix(SEXP pakl, matrix_t *akl) {
  akl->data=vector_as_vector(REAL(pakl), GET_LENGTH(pakl));
  akl->nrow=INTEGER(GET_DIM(pakl))[0];
  akl->ncol=INTEGER(GET_DIM(pakl))[1];

  return 0;
}

int R_SEXP_to_matrix_copy(SEXP pakl, matrix_t *akl) {
  vector_init_copy(&akl->data, REAL(pakl), GET_LENGTH(pakl));
  akl->nrow=INTEGER(GET_DIM(pakl))[0];
  akl->ncol=INTEGER(GET_DIM(pakl))[1];

  return 0;
}
 
int R_SEXP_to_attributes(SEXP attr, igraph_attribute_list_t *al) {
  
  al->len=INTEGER(GET_DIM(VECTOR_ELT(attr, 1)))[0];
  al->nstr=INTEGER(GET_DIM(VECTOR_ELT(attr, 1)))[1];
  al->sa_begin=RAW(VECTOR_ELT(attr, 0));
  al->sa_end=al->sa_begin+GET_LENGTH(VECTOR_ELT(attr, 0));
  R_SEXP_to_matrix(VECTOR_ELT(attr,1), &al->numattrs);
  
  return 0;
}

int R_SEXP_to_attributes_copy(SEXP attr, igraph_attribute_list_t *al) {

  al->len=INTEGER(GET_DIM(VECTOR_ELT(attr, 1)))[0];
  al->nstr=INTEGER(GET_DIM(VECTOR_ELT(attr, 1)))[1];
  al->sa_begin=Calloc(GET_LENGTH(VECTOR_ELT(attr, 0)), char);
  al->sa_end=al->sa_begin+GET_LENGTH(VECTOR_ELT(attr, 0));
  memcpy(al->sa_begin, RAW(VECTOR_ELT(attr, 0)), 
	 sizeof(char)*GET_LENGTH(VECTOR_ELT(attr, 0)));
  R_SEXP_to_matrix_copy(VECTOR_ELT(attr,1), &al->numattrs);
  
  return 0;
}

SEXP R_igraph_to_SEXP(igraph_t *graph) {
  
  SEXP result;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  
  PROTECT(result=NEW_LIST(10));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(1));
  SET_VECTOR_ELT(result, 1, NEW_LOGICAL(1));
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 3, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 4, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 5, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 6, NEW_NUMERIC(no_of_nodes+1));
  SET_VECTOR_ELT(result, 7, NEW_NUMERIC(no_of_nodes+1));

  SET_VECTOR_ELT(result, 8, R_attributes_to_SEXP(&graph->gal));
  SET_VECTOR_ELT(result, 9, R_attributes_to_SEXP(&graph->val));
  
  REAL(VECTOR_ELT(result, 0))[0]=no_of_nodes;
  LOGICAL(VECTOR_ELT(result, 1))[0]=graph->directed;
  memcpy(REAL(VECTOR_ELT(result, 2)), graph->from.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 3)), graph->to.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 4)), graph->oi.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 5)), graph->ii.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 6)), graph->os.stor_begin, 
	 sizeof(real_t)*(no_of_nodes+1));
  memcpy(REAL(VECTOR_ELT(result, 7)), graph->is.stor_begin, 
	 sizeof(real_t)*(no_of_nodes+1));
  
  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("igraph")));
  
  UNPROTECT(1);
  return result;
}

int R_SEXP_to_igraph(SEXP graph, igraph_t *res) {
  
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  res->from=vector_as_vector(REAL(VECTOR_ELT(graph, 2)), 
			    GET_LENGTH(VECTOR_ELT(graph, 2)));
  res->to=vector_as_vector(REAL(VECTOR_ELT(graph, 3)), 
			   GET_LENGTH(VECTOR_ELT(graph, 3)));
  res->oi=vector_as_vector(REAL(VECTOR_ELT(graph, 4)), 
			   GET_LENGTH(VECTOR_ELT(graph, 4)));
  res->ii=vector_as_vector(REAL(VECTOR_ELT(graph, 5)), 
			   GET_LENGTH(VECTOR_ELT(graph, 5)));
  res->os=vector_as_vector(REAL(VECTOR_ELT(graph, 6)), 
			   GET_LENGTH(VECTOR_ELT(graph, 6)));
  res->is=vector_as_vector(REAL(VECTOR_ELT(graph, 7)), 
			   GET_LENGTH(VECTOR_ELT(graph, 7)));
  
  R_SEXP_to_attributes(VECTOR_ELT(graph, 8), &res->gal);
  R_SEXP_to_attributes(VECTOR_ELT(graph, 9), &res->val);
  
  return 0;
}

int R_SEXP_to_igraph_copy(SEXP graph, igraph_t *res) {
  
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  vector_init_copy(&res->from, REAL(VECTOR_ELT(graph, 2)), 
		   GET_LENGTH(VECTOR_ELT(graph, 2)));
  vector_init_copy(&res->to, REAL(VECTOR_ELT(graph, 3)), 
		   GET_LENGTH(VECTOR_ELT(graph, 3)));
  vector_init_copy(&res->oi, REAL(VECTOR_ELT(graph, 4)), 
		   GET_LENGTH(VECTOR_ELT(graph, 4)));
  vector_init_copy(&res->ii, REAL(VECTOR_ELT(graph, 5)), 
		   GET_LENGTH(VECTOR_ELT(graph, 5)));
  vector_init_copy(&res->os, REAL(VECTOR_ELT(graph, 6)), 
		   GET_LENGTH(VECTOR_ELT(graph, 6)));
  vector_init_copy(&res->is, REAL(VECTOR_ELT(graph, 7)),
		   GET_LENGTH(VECTOR_ELT(graph, 7)));

  R_SEXP_to_attributes_copy(VECTOR_ELT(graph, 8), &res->gal);
  R_SEXP_to_attributes_copy(VECTOR_ELT(graph, 9), &res->val);
  
  return 0;
}

SEXP R_igraph_iterator_to_SEXP(igraph_t *g, igraph_iterator_t *it) {
  
  SEXP result;
  long int datalen=1;
  
  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(1));
  REAL(VECTOR_ELT(result, 0))[0] = it->type;
  switch ((long int)it->type) {
  case IGRAPH_ITERATOR_VID:
  case IGRAPH_ITERATOR_EID:
  case IGRAPH_ITERATOR_EFROMORDER:
    datalen=1;
    break;
  case IGRAPH_ITERATOR_ENEIS: 
  case IGRAPH_ITERATOR_VNEIS:
    datalen=4;
    break;
  default:
    break;
  }
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(datalen));
  memcpy(REAL(VECTOR_ELT(result, 1)), (real_t *)it->data, 
	 datalen * sizeof (real_t));

  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("igraph.iterator")));  

  UNPROTECT(1);
  return result;
}

int R_SEXP_to_igraph_iterator(SEXP rit, igraph_iterator_t *it) {
  
  it->type=REAL(VECTOR_ELT(rit, 0))[0];
  it->data=REAL(VECTOR_ELT(rit, 1));
  switch ( (long int) it->type ) {
  case IGRAPH_ITERATOR_VID:
    it->next=igraph_next_vid;
    it->prev=igraph_prev_vid;
    it->end=igraph_end_vid;
    it->reset=igraph_reset_vid;
    it->getvertex=igraph_get_vertex_vid;
    it->getvertexfrom=it->getvertexto=it->getedge=0;
    it->getvertexnei=0;    
    break;
  case IGRAPH_ITERATOR_EID:
    it->next=igraph_next_eid;
    it->prev=igraph_prev_eid;
    it->end=igraph_end_eid;
    it->reset=igraph_reset_eid;
    it->getvertex=0;
    it->getvertexfrom=igraph_get_vertex_from_eid;
    it->getvertexto=igraph_get_vertex_to_eid;
    it->getedge=igraph_get_edge_eid;
    it->getvertexnei=0;
    break;
  case IGRAPH_ITERATOR_EFROMORDER:
    it->next=igraph_next_efromorder;
    it->prev=igraph_prev_efromorder;
    it->end=igraph_end_efromorder;
    it->reset=igraph_reset_efromorder;
    it->getvertex=0;
    it->getvertexfrom=igraph_get_vertex_from_efromorder;
    it->getvertexto=igraph_get_vertex_to_efromorder;
    it->getedge=igraph_get_edge_efromorder;
    it->getvertexnei=0;
    break;
  case IGRAPH_ITERATOR_ENEIS:
    it->next=igraph_next_eneis;
    it->prev=0;
    it->end=igraph_end_eneis;
    it->reset=igraph_reset_eneis;
    it->getvertex=0;
    it->getvertexfrom=igraph_get_vertex_from_eneis;
    it->getvertexto=igraph_get_vertex_to_eneis;
    it->getedge=igraph_get_edge_eneis;
    it->getvertexnei=igraph_get_vertex_nei_eneis;
    break;
  case IGRAPH_ITERATOR_VNEIS:
    it->type=IGRAPH_ITERATOR_VNEIS;
    it->next=igraph_next_vneis;
    it->prev=0;
    it->end=igraph_end_vneis;
    it->reset=igraph_reset_vneis;
    it->getvertex=igraph_get_vertex_vneis;
    it->getvertexfrom=0;
    it->getvertexto=0;
    it->getedge=0;
    it->getvertexnei=0;
    break;
  }
  
  return 0;
}

int R_SEXP_to_igraph_iterator_copy(SEXP rit, igraph_iterator_t *it) {
  
  it->type=REAL(VECTOR_ELT(rit, 0))[0];
  long int datalen=1;
  switch ( (long int) it->type ) {
  case IGRAPH_ITERATOR_VID:
    datalen=1;
    it->next=igraph_next_vid;
    it->prev=igraph_prev_vid;
    it->end=igraph_end_vid;
    it->reset=igraph_reset_vid;
    it->getvertex=igraph_get_vertex_vid;
    it->getvertexfrom=it->getvertexto=it->getedge=0;
    it->getvertexnei=0;    
    break;
  case IGRAPH_ITERATOR_EID:
    datalen=1;
    it->next=igraph_next_eid;
    it->prev=igraph_prev_eid;
    it->end=igraph_end_eid;
    it->reset=igraph_reset_eid;
    it->getvertex=0;
    it->getvertexfrom=igraph_get_vertex_from_eid;
    it->getvertexto=igraph_get_vertex_to_eid;
    it->getedge=igraph_get_edge_eid;
    it->getvertexnei=0;
    break;
  case IGRAPH_ITERATOR_EFROMORDER:
    datalen=1;
    it->next=igraph_next_efromorder;
    it->prev=igraph_prev_efromorder;
    it->end=igraph_end_efromorder;
    it->reset=igraph_reset_efromorder;
    it->getvertex=0;
    it->getvertexfrom=igraph_get_vertex_from_efromorder;
    it->getvertexto=igraph_get_vertex_to_efromorder;
    it->getedge=igraph_get_edge_efromorder;
    it->getvertexnei=0;
    break;
  case IGRAPH_ITERATOR_ENEIS:
    datalen=4;
    it->next=igraph_next_eneis;
    it->prev=0;
    it->end=igraph_end_eneis;
    it->reset=igraph_reset_eneis;
    it->getvertex=0;
    it->getvertexfrom=igraph_get_vertex_from_eneis;
    it->getvertexto=igraph_get_vertex_to_eneis;
    it->getedge=igraph_get_edge_eneis;
    it->getvertexnei=igraph_get_vertex_nei_eneis;
    break;
  case IGRAPH_ITERATOR_VNEIS:
    datalen=4;
    it->type=IGRAPH_ITERATOR_VNEIS;
    it->next=igraph_next_vneis;
    it->prev=0;
    it->end=igraph_end_vneis;
    it->reset=igraph_reset_vneis;
    it->getvertex=igraph_get_vertex_vneis;
    it->getvertexfrom=0;
    it->getvertexto=0;
    it->getedge=0;
    it->getvertexnei=0;
    break;
  default: 
    break;
  }
  it->data=Calloc(datalen, real_t);
  memcpy(it->data, REAL(VECTOR_ELT(rit, 1)), datalen * sizeof(real_t));
  
  return 0;
}

/*******************************************************************/

SEXP R_igraph_empty(SEXP pn, SEXP pdirected) {
  
  SEXP result;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  igraph_t g;

  igraph_empty(&g, n, directed);  
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_edges(SEXP graph, SEXP edges) {
  
  vector_t v;			/* do NOT destroy! */
  igraph_t g;
  SEXP result;

  v=vector_as_vector(REAL(edges), GET_LENGTH(edges));
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_edges(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_vertices(SEXP graph, SEXP pnv) {
  
  integer_t nv;
  igraph_t g;
  SEXP result;
  
  nv=REAL(pnv)[0];
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_vertices(&g, nv);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_vcount(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=igraph_vcount(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_ecount(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=igraph_ecount(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_neighbors(SEXP graph, SEXP pvid, SEXP pmode) {
  
  igraph_t g;
  vector_t neis;
  SEXP result;
  real_t vid;
  integer_t mode;
  
  vector_init(&neis, 0);
  vid=REAL(pvid)[0];
  mode=REAL(pmode)[0];
  R_SEXP_to_igraph(graph, &g);
  igraph_neighbors(&g, &neis, vid, mode);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&neis)));
  vector_copy_to(&neis, REAL(result));
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_delete_edges(SEXP graph, SEXP edges) {
  
  vector_t v;
  igraph_t g;
  SEXP result;
  
  v=vector_as_vector(REAL(edges), GET_LENGTH(edges));
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_delete_edges(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_delete_vertices(SEXP graph, SEXP vertices) {
  
  vector_t v;
  igraph_t g;
  SEXP result;
  
  v=vector_as_vector(REAL(vertices), GET_LENGTH(vertices));
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_delete_vertices(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_is_directed(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=igraph_is_directed(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_create(SEXP edges, SEXP pn, SEXP pdirected) {
  
  igraph_t g;
  vector_t v;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  v=vector_as_vector(REAL(edges), GET_LENGTH(edges));
  igraph_create(&g, &v, n, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_degree(SEXP graph, SEXP vids, SEXP pmode, SEXP ploops) {
  
  igraph_t g;
  vector_t v;
  vector_t res;
  integer_t mode=REAL(pmode)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;

  v=vector_as_vector(REAL(vids), GET_LENGTH(vids));
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, vector_size(&v));
  igraph_degree(&g, &res, &v, mode, loops);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  
  vector_destroy(&res);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_clusters(SEXP graph, SEXP pmode) {

  igraph_t g;
  integer_t mode=REAL(pmode)[0];
  vector_t membership;
  vector_t csize;
  SEXP result, names;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&membership, 0);
  vector_init(&csize, 0);
  igraph_clusters(&g, &membership, &csize, mode);
  
  PROTECT(result=NEW_LIST(2));
  PROTECT(names=NEW_CHARACTER(2));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(vector_size(&membership)));
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(vector_size(&csize)));
  SET_STRING_ELT(names, 0, CREATE_STRING_VECTOR("membership"));
  SET_STRING_ELT(names, 1, CREATE_STRING_VECTOR("csize"));
  SET_NAMES(result, names);
  vector_copy_to(&membership, REAL(VECTOR_ELT(result, 0)));
  vector_copy_to(&csize, REAL(VECTOR_ELT(result, 1)));
  
  vector_destroy(&membership);
  vector_destroy(&csize);
  UNPROTECT(2);
  return result;
}

SEXP R_igraph_is_connected(SEXP graph, SEXP pmode) {

  igraph_t g;
  integer_t mode=REAL(pmode)[0];
  bool_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  igraph_is_connected(&g, &res, mode);
  
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_diameter(SEXP graph, SEXP pdirected, SEXP punconnected) {
  
  igraph_t g;
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t unconnected=LOGICAL(punconnected)[0];
  integer_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  igraph_diameter(&g, &res, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_closeness(SEXP graph, SEXP pvids, SEXP pmode) {
  
  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  integer_t mode=REAL(pmode)[0];
  vector_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_closeness(&g, &res, &vids, mode);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_subcomponent(SEXP graph, SEXP pvertex, SEXP pmode) {
  
  igraph_t g;
  real_t vertex=REAL(pvertex)[0];
  integer_t mode=REAL(pmode)[0];
  vector_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_subcomponent(&g, &res, vertex, mode);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_betweenness(SEXP graph, SEXP pvids, SEXP pdirected) {
  
  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  bool_t directed=LOGICAL(pdirected)[0];
  vector_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_betweenness(&g, &res, &vids, directed);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_running_mean(SEXP pdata, SEXP pbinwidth) {
  
  vector_t data=vector_as_vector(REAL(pdata), GET_LENGTH(pdata));
  integer_t binwidth=REAL(pbinwidth)[0];
  vector_t res;
  SEXP result;
  
  vector_init(&res, 0);
  igraph_running_mean(&data, &res, binwidth);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_cocitation(SEXP graph, SEXP pvids) {

  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  matrix_t m;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&m, 0, 0);
  igraph_cocitation(&g, &m, &vids);
  
  PROTECT(result=R_matrix_to_SEXP(&m));
  matrix_destroy(&m);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_bibcoupling(SEXP graph, SEXP pvids) {

  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  matrix_t m;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&m, 0, 0);
  igraph_bibcoupling(&g, &m, &vids);
  
  PROTECT(result=R_matrix_to_SEXP(&m));
  matrix_destroy(&m);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_growing_random_game(SEXP pn, SEXP pm, SEXP pdirected,
				  SEXP pcitation) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t m=REAL(pm)[0];
  bool_t citation=LOGICAL(pcitation)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  igraph_growing_random_game(&g, n, m, directed, citation);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_shortest_paths(SEXP graph, SEXP pvids, SEXP pmode) {
  
  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  integer_t mode=REAL(pmode)[0];
  matrix_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_shortest_paths(&g, &res, &vids, mode);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_lattice(SEXP pdimvector, SEXP pnei, SEXP pdirected,
		      SEXP pmutual, SEXP pcircular) {
  
  igraph_t g;
  vector_t dimvector=vector_as_vector(REAL(pdimvector), 
				      GET_LENGTH(pdimvector));
  integer_t nei=REAL(pnei)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t mutual=LOGICAL(pmutual)[0];
  bool_t circular=LOGICAL(pcircular)[0];  
  SEXP result;
  
  igraph_lattice(&g, &dimvector, nei, directed, mutual, circular);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_barabasi_game(SEXP pn, SEXP pm, SEXP poutseq,
			    SEXP poutpref, SEXP pdirected) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t m=REAL(pm)[0]; 
  vector_t outseq=vector_as_vector(REAL(poutseq), GET_LENGTH(poutseq));
  bool_t outpref=LOGICAL(poutpref)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  igraph_barabasi_game(&g, n, m, &outseq, outpref, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_layout_kamada_kawai(SEXP graph, SEXP pniter, SEXP pinitemp, 
				  SEXP pcoolexp, SEXP pkkconst, SEXP psigma) {
  
  igraph_t g;
  integer_t niter=REAL(pniter)[0];
  real_t initemp=REAL(pinitemp)[0];
  real_t coolexp=REAL(pcoolexp)[0];
  real_t kkconst=REAL(pkkconst)[0];
  real_t sigma=REAL(psigma)[0];
  matrix_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_layout_kamada_kawai(&g, &res, niter, sigma, 
			     initemp, coolexp, kkconst);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_minimum_spanning_tree_unweighted(SEXP graph) {
  
  igraph_t g;
  igraph_t mst;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_minimum_spanning_tree_unweighted(&g, &mst);
  PROTECT(result=R_igraph_to_SEXP(&mst));
  igraph_destroy(&mst);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_minimum_spanning_tree_prim(SEXP graph, SEXP pweights) {
  
  igraph_t g;
  igraph_t mst;
  vector_t weights=vector_as_vector(REAL(pweights), GET_LENGTH(pweights));
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_minimum_spanning_tree_prim(&g, &mst, &weights);
  PROTECT(result=R_igraph_to_SEXP(&mst));
  igraph_destroy(&mst);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_edge_betweenness(SEXP graph, SEXP pdirected) {
  
  igraph_t g;
  vector_t res;
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_edge_betweenness(&g, &res, directed);
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  
  vector_destroy(&res);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage(SEXP graph, SEXP pst, SEXP pagebins,
				     SEXP pmaxind, SEXP plsd) {
  
  igraph_t g;
  matrix_t akl, sd;
  vector_t st=vector_as_vector(REAL(pst), GET_LENGTH(pst));
  integer_t agebins=REAL(pagebins)[0];
  integer_t maxind=REAL(pmaxind)[0];
  bool_t lsd=LOGICAL(plsd)[0];
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&akl, 0, 0);
  matrix_init(&sd, 0, 0);
  igraph_measure_dynamics_idage(&g, &akl, &sd, &st, agebins, maxind, lsd);
  
  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, R_matrix_to_SEXP(&akl));
  matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_matrix_to_SEXP(&sd));
  matrix_destroy(&sd);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage_debug(SEXP graph, SEXP pst, 
					   SEXP pagebins, SEXP pmaxind, 
					   SEXP plsd, SEXP pest_ind, 
					   SEXP pest_age) {
  igraph_t g;
  matrix_t akl, sd;
  vector_t st=vector_as_vector(REAL(pst), GET_LENGTH(pst));
  integer_t agebins=REAL(pagebins)[0];
  integer_t maxind=REAL(pmaxind)[0];
  bool_t lsd=LOGICAL(plsd)[0];
  integer_t est_ind=REAL(pest_ind)[0];
  integer_t est_age=REAL(pest_age)[0];
  vector_t estimates;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&akl, 0, 0);
  matrix_init(&sd, 0, 0);
  vector_init(&estimates, 0);
  igraph_measure_dynamics_idage_debug(&g, &akl, &sd, &st, agebins, maxind, lsd,
				      &estimates, est_ind, est_age);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_matrix_to_SEXP(&akl));
  matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_matrix_to_SEXP(&sd));
  matrix_destroy(&sd);
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(vector_size(&estimates)));
  vector_copy_to(&estimates, REAL(VECTOR_ELT(result, 2)));
  vector_destroy(&estimates);

  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_measure_dynamics_idage_st(SEXP graph, SEXP pakl) {
  
  igraph_t g;
  matrix_t akl;
  vector_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_matrix(pakl, &akl);
  vector_init(&res, 0);
  igraph_measure_dynamics_idage_st(&g, &res, &akl);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  
  vector_destroy(&res);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_shortest_paths(SEXP graph, SEXP pfrom, SEXP pmode) {

  igraph_t g;
  integer_t from=REAL(pfrom)[0];
  integer_t mode=REAL(pmode)[0];
  vector_t *vects;
  long int no_of_nodes, i;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  no_of_nodes=igraph_vcount(&g);
  vects=Calloc(no_of_nodes, vector_t);
  for (i=0; i<no_of_nodes; i++) {
    vector_init(&vects[i], 0);
  }
  igraph_get_shortest_paths(&g, vects, from, mode);
  PROTECT(result=NEW_LIST(no_of_nodes));
  for (i=0; i<no_of_nodes; i++) {
    SET_VECTOR_ELT(result, i, NEW_NUMERIC(vector_size(&vects[i])));
    vector_copy_to(&vects[i], REAL(VECTOR_ELT(result, i)));
    vector_destroy(&vects[i]);
  }
  
  Free(vects);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_are_connected(SEXP graph, SEXP pv1, SEXP pv2) {
  
  igraph_t g;
  integer_t v1=REAL(pv1)[0];
  integer_t v2=REAL(pv2)[0];
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=igraph_are_connected(&g, v1, v2);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_graph_adjacency(SEXP adjmatrix, SEXP pmode) {
  
  igraph_t g;
  matrix_t adjm;
  integer_t mode=REAL(pmode)[0];
  SEXP result;
  
  R_SEXP_to_matrix(adjmatrix, &adjm);
  igraph_adjacency(&g, &adjm, mode);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_average_path_length(SEXP graph, SEXP pdirected, 
				  SEXP punconnected) {
  
  igraph_t g;
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t unconnected=LOGICAL(punconnected)[0];
  real_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  igraph_average_path_length(&g, &res, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_star(SEXP pn, SEXP pmode, SEXP pcenter) {

  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t mode=REAL(pmode)[0];
  integer_t center=REAL(pcenter)[0];
  SEXP result;
  
  igraph_star(&g, n, mode, center);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_ring(SEXP pn, SEXP pdirected, SEXP pmutual, SEXP pcircular) {

  igraph_t g;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t mutual=LOGICAL(pmutual)[0];
  bool_t circular=LOGICAL(pcircular)[0];
  SEXP result;
  
  igraph_ring(&g, n, directed, mutual, circular);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_tree(SEXP pn, SEXP pchildren, SEXP pmode) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t children=REAL(pchildren)[0];
  integer_t mode=REAL(pmode)[0];
  SEXP result;

  igraph_tree(&g, n, children, mode);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_subgraph(SEXP graph, SEXP pvids) {
  
  igraph_t g;
  igraph_t sub;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_subgraph(&g, &sub, &vids);
  PROTECT(result=R_igraph_to_SEXP(&sub));
  igraph_destroy(&sub);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_random(SEXP graph) {
  
  igraph_t g;
  matrix_t res;
  SEXP result=R_NilValue;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_layout_random(&g, &res);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_circle(SEXP graph) {
  
  igraph_t g;
  matrix_t res;
  SEXP result=R_NilValue;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_layout_circle(&g, &res);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_erdos_renyi_game(SEXP pn, SEXP ptype,
			       SEXP pporm, SEXP pdirected, SEXP ploops) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t type=REAL(ptype)[0];
  real_t porm=REAL(pporm)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;
  
  igraph_erdos_renyi_game(&g, type, n, porm, directed, loops);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_full(SEXP pn, SEXP pdirected, SEXP ploops) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;
  
  igraph_full(&g, n, directed, loops);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_random_sample(SEXP plow, SEXP phigh, SEXP plength) {
  
  vector_t res;
  integer_t low=REAL(plow)[0];
  integer_t high=REAL(phigh)[0];
  integer_t length=REAL(plength)[0];
  SEXP result;

  vector_init(&res, 0);
  igraph_random_sample(&res, low, high, length);
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_edgelist(SEXP graph, SEXP pbycol) {
  
  igraph_t g;
  vector_t res;
  bool_t bycol=LOGICAL(pbycol)[0];
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_get_edgelist(&g, &res, bycol);
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_get_adjacency(SEXP graph, SEXP ptype) {
  
  igraph_t g;
  matrix_t res;
  integer_t type=REAL(ptype)[0];
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_get_adjacency(&g, &res, type);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_simplify(SEXP graph, SEXP pmultiple, SEXP ploops) {
  
  igraph_t g;
  bool_t multiple=LOGICAL(pmultiple)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_simplify(&g, multiple, loops);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_fruchterman_reingold(SEXP graph, SEXP pniter, 
					  SEXP pmaxdelta, SEXP parea,
					  SEXP pcoolexp, SEXP prepulserad) {
  igraph_t g;
  integer_t niter=REAL(pniter)[0];
  real_t maxdelta=REAL(pmaxdelta)[0];
  real_t area=REAL(parea)[0];
  real_t coolexp=REAL(pcoolexp)[0];
  real_t repulserad=REAL(prepulserad)[0];
  matrix_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_layout_fruchterman_reingold(&g, &res, niter, maxdelta, area, 
				     coolexp, repulserad, 0);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_degree_sequence_game(SEXP pout_seq, SEXP pin_seq, 
				   SEXP pmethod) {
  igraph_t g;
  vector_t outseq=vector_as_vector(REAL(pout_seq), GET_LENGTH(pout_seq));
  vector_t inseq=vector_as_vector(REAL(pin_seq), GET_LENGTH(pin_seq));
  integer_t method=REAL(pmethod)[0];
  SEXP result;

  igraph_degree_sequence_game(&g, &outseq, &inseq, method);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_transitivity(SEXP graph, SEXP ptype) {
  
  igraph_t g;
  vector_t res;
  integer_t type=REAL(ptype)[0];
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_transitivity(&g, &res, type);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_add_graph_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  SEXP result;

  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_graph_attribute(&g, name);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_remove_graph_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  SEXP result;

  R_SEXP_to_igraph_copy(graph, &g);
  igraph_remove_graph_attribute(&g, name);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_graph_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  real_t value;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  igraph_get_graph_attribute(&g, name, &value);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=value;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_graph_attribute(SEXP graph, SEXP pname, SEXP pvalue) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  real_t value=REAL(pvalue)[0];
  SEXP result;

  R_SEXP_to_igraph_copy(graph, &g);
  igraph_set_graph_attribute(&g, name, value);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_vertex_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  SEXP result;

  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_vertex_attribute(&g, name);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_remove_vertex_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  SEXP result;

  R_SEXP_to_igraph_copy(graph, &g);
  igraph_remove_vertex_attribute(&g, name);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_vertex_attribute(SEXP graph, SEXP pname, SEXP pv) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  long int v=REAL(pv)[0];
  real_t value;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  igraph_get_vertex_attribute(&g, name, v, &value);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=value;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_vertex_attribute(SEXP graph, SEXP pname, SEXP pv, 
				   SEXP pvalue) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  long int v=REAL(pv)[0];
  real_t value=REAL(pvalue)[0];
  SEXP result;

  R_SEXP_to_igraph_copy(graph, &g);
  igraph_set_vertex_attribute(&g, name, v, value);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_vertex_attributes(SEXP graph, SEXP pname, SEXP pv) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  vector_t v=vector_as_vector(REAL(pv), GET_LENGTH(pv));
  vector_t value;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  vector_init(&value, 0);
  igraph_get_vertex_attributes(&g, name, &v, &value);
  PROTECT(result=NEW_NUMERIC(vector_size(&value)));
  vector_copy_to(&value, REAL(result));
  vector_destroy(&value);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_vertex_attributes(SEXP graph, SEXP pname, SEXP pv, 
				    SEXP pvalue) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  vector_t v=vector_as_vector(REAL(pv), GET_LENGTH(pv));
  vector_t value=vector_as_vector(REAL(pvalue), GET_LENGTH(pvalue));
  SEXP result;

  R_SEXP_to_igraph_copy(graph, &g);
  igraph_set_vertex_attributes(&g, name, &v, &value);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_list_graph_attributes(SEXP graph) {
  igraph_t g;
  igraph_strarray_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_strarray_init(&res);
  igraph_list_graph_attributes(&g, &res);
  PROTECT(result=R_strarray_to_SEXP(&res));
  igraph_strarray_destroy(&res);
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_list_vertex_attributes(SEXP graph) {
  igraph_t g;
  igraph_strarray_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_strarray_init(&res);
  igraph_list_vertex_attributes(&g, &res);
  PROTECT(result=R_strarray_to_SEXP(&res));
  igraph_strarray_destroy(&res);
  
  UNPROTECT(1);
  return result;  

}

SEXP R_igraph_iterator_vid(SEXP graph) {
  
  igraph_t g;
  igraph_iterator_t it;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_iterator_vid(&g, &it);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_eid(SEXP graph) {

  igraph_t g;
  igraph_iterator_t it;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_iterator_eid(&g, &it);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_efromorder(SEXP graph) {

  igraph_t g;
  igraph_iterator_t it;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_iterator_efromorder(&g, &it);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_eneis(SEXP graph, SEXP pvid, SEXP pmode) {

  igraph_t g;
  igraph_iterator_t it;
  integer_t vid=REAL(pvid)[0];
  integer_t mode=REAL(pmode)[0];
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_iterator_eneis(&g, &it, vid, mode);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_vneis(SEXP graph, SEXP pvid, SEXP pmode) {

  igraph_t g;
  igraph_iterator_t it;
  integer_t vid=REAL(pvid)[0];
  integer_t mode=REAL(pmode)[0];
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_iterator_vneis(&g, &it, vid, mode);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_next(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator_copy(pit, &it);
  igraph_next(&g, &it);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_prev(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator_copy(pit, &it);
  igraph_prev(&g, &it);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_reset(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator_copy(pit, &it);
  igraph_reset(&g, &it);
  PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it));
  igraph_iterator_destroy(&g, &it);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_end(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  bool_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator(pit, &it);
  res=igraph_end(&g, &it);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_get_vertex_nei(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  integer_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator(pit, &it);
  res=igraph_get_vertex_nei(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_get_vertex(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  integer_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator(pit, &it);
  res=igraph_get_vertex(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_from(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  integer_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator(pit, &it);
  res=igraph_get_vertex_from(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_to(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  integer_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator(pit, &it);
  res=igraph_get_vertex_to(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_iterator_edge(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_iterator_t it;
  integer_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_iterator(pit, &it);
  res=igraph_get_edge(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}
