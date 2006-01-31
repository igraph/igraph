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
#include "error.h"

#include "config.h"

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <stdio.h>

SEXP R_igraph_matrix_to_SEXP(igraph_matrix_t *m);
SEXP R_igraph_strvector_to_SEXP(igraph_strvector_t *m);
SEXP R_attributes_to_SEXP(igraph_attribute_list_t *al);
SEXP R_igraph_vs_to_SEXP(igraph_t *g, igraph_vs_t *vs);
SEXP R_igraph_es_to_SEXP(igraph_t *g, igraph_es_t *vs);
SEXP R_igraph_to_SEXP(igraph_t *graph);

int R_igraph_SEXP_to_strvector(SEXP rval, igraph_strvector_t *sv);
int R_igraph_SEXP_to_strvector_copy(SEXP rval, igraph_strvector_t *sv);
int R_SEXP_to_vector(SEXP sv, igraph_vector_t *v);
int R_SEXP_to_vector_copy(SEXP sv, igraph_vector_t *v);
int R_SEXP_to_matrix(SEXP pakl, igraph_matrix_t *akl);
int R_SEXP_to_igraph_matrix_copy(SEXP pakl, igraph_matrix_t *akl);
int R_SEXP_to_attributes(SEXP attr, igraph_attribute_list_t *al);
int R_SEXP_to_attributes_copy(SEXP attr, igraph_attribute_list_t *al);
int R_SEXP_to_igraph(SEXP graph, igraph_t *res);
int R_SEXP_to_igraph_copy(SEXP graph, igraph_t *res);
int R_SEXP_to_igraph_attr(SEXP graph, igraph_t *res);
int R_SEXP_to_igraph_vs_copy(SEXP rit, igraph_t *graph, igraph_vs_t *it);
int R_SEXP_to_igraph_es_copy(SEXP rit, igraph_t *graph, igraph_es_t *it);

/******************************************************
 * things to do before and after                      *
 * calling an interface function                      *
 *****************************************************/

igraph_error_handler_t *R_igraph_oldhandler;

void R_igraph_myhandler (const char *reason, const char *file,
			 int line, int igraph_errno) {
  IGRAPH_FINALLY_FREE();
  error("At %s:%i : %s, %s", file, line, reason, 
	igraph_strerror(igraph_errno));
}

R_INLINE void R_igraph_before() {
  R_igraph_oldhandler=igraph_set_error_handler(R_igraph_myhandler);
}

R_INLINE void R_igraph_after() {
  igraph_set_error_handler(R_igraph_oldhandler);
}

/******************************************************
 * functions to convert igraph objects to SEXP
 *****************************************************/

SEXP R_igraph_matrix_to_SEXP(igraph_matrix_t *m) {

  SEXP result, dim; 
  
  PROTECT(result=NEW_NUMERIC(igraph_matrix_size(m)));
  igraph_matrix_copy_to(m, REAL(result));
  PROTECT(dim=NEW_INTEGER(2));
  INTEGER(dim)[0]=igraph_matrix_nrow(m);
  INTEGER(dim)[1]=igraph_matrix_ncol(m);
  SET_DIM(result, dim);

  UNPROTECT(2);
  return result;
}

SEXP R_igraph_strvector_to_SEXP(igraph_strvector_t *m) {
  SEXP result;
  long int i;
  char *str;
  long int len;
  
  len=igraph_strvector_size(m);
  PROTECT(result=NEW_CHARACTER(len));
  for (i=0; i<len; i++) {
    igraph_strvector_get(m, i, &str);
    SET_STRING_ELT(result, i, CREATE_STRING_VECTOR(str));
  }
  
  UNPROTECT(1);
  return result;
}

SEXP R_attributes_to_SEXP(igraph_attribute_list_t *al) {
  long int nattr, i;
  SEXP result;

  nattr=igraph_attribute_list_size(al);
  PROTECT(result=NEW_LIST(nattr+3));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(1));
  REAL(VECTOR_ELT(result, 0))[0]=al->len;
  SET_VECTOR_ELT(result, 1, R_igraph_strvector_to_SEXP(&al->names));
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(nattr));
  igraph_vector_copy_to(&al->types, REAL(VECTOR_ELT(result, 2)));
  for (i=0; i<nattr; i++) {
    if (VECTOR(al->types)[i]==IGRAPH_ATTRIBUTE_NUM) {
      igraph_vector_t *data=VECTOR(al->data)[i];
      SET_VECTOR_ELT(result, i+3, NEW_NUMERIC(igraph_vector_size(data)));
      igraph_vector_copy_to(data, REAL(VECTOR_ELT(result, i+3)));
    } else if (VECTOR(al->types)[i]==IGRAPH_ATTRIBUTE_STR) {
      igraph_strvector_t *data=VECTOR(al->data)[i];
      SET_VECTOR_ELT(result, i+3, R_igraph_strvector_to_SEXP(data));
    }
  }

  UNPROTECT(1);
  return result;
}

/*
 * All types of iterators are converted to a vector or a sequence
 * Vector right now. :)
 */
SEXP R_igraph_vs_to_SEXP(igraph_t *g, igraph_vs_t *it) {

  SEXP result;
  igraph_vs_t vectorvs;
  const igraph_vector_t *v;
  long int vlen;

  igraph_vs_vectorview_it(g, it, &vectorvs);
  IGRAPH_FINALLY(igraph_vs_destroy, &vectorvs);
  v=igraph_vs_vector_getvector(g, &vectorvs);
  vlen=igraph_vector_size(v);
  
  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(3));
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(vlen));
  REAL(VECTOR_ELT(result, 0))[0]=1;
  REAL(VECTOR_ELT(result, 0))[1]=1;
  REAL(VECTOR_ELT(result, 0))[2]=vlen;
  igraph_vector_copy_to(v, REAL(VECTOR_ELT(result, 1)));
  
  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("igraphvsvector")));
  
  igraph_vs_destroy(&vectorvs);
  IGRAPH_FINALLY_CLEAN(1);
  UNPROTECT(1);
  return result;
}

/* 
 * All types of iterators are converted to a vector or a sequence
 * Vector right now. :)
 */
SEXP R_igraph_es_to_SEXP(igraph_t *g, igraph_es_t *it) {

  SEXP result;
  igraph_es_t vectores;
  const igraph_vector_t *v;
  long int vlen;

  igraph_es_vectorview_it(g, it, &vectores);
  IGRAPH_FINALLY(igraph_es_destroy, &vectores);
  v=igraph_es_vector_getvector(g, &vectores);
  vlen=igraph_vector_size(v);
  
  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(3));
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(vlen));
  REAL(VECTOR_ELT(result, 0))[0]=1;
  REAL(VECTOR_ELT(result, 0))[1]=1;
  REAL(VECTOR_ELT(result, 0))[2]=vlen;
  igraph_vector_copy_to(v, REAL(VECTOR_ELT(result, 1)));
  
  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("igraphesvector")));
  
  igraph_es_destroy(&vectores);
  IGRAPH_FINALLY_CLEAN(1);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_to_SEXP(igraph_t *graph) {
  
  SEXP result;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  
  PROTECT(result=NEW_LIST(11));
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
  SET_VECTOR_ELT(result,10, R_attributes_to_SEXP(&graph->eal));
  
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

int R_igraph_SEXP_to_strvector(SEXP rval, igraph_strvector_t *sv) {
  long int i;
  sv->len=GET_LENGTH(rval);
  sv->data=(char**) R_alloc(sv->len, sizeof(char*));
  for (i=0; i<sv->len; i++) {
    sv->data[i]=CHAR(STRING_ELT(rval, i));
  }

  return 0;
}

int R_igraph_SEXP_to_strvector_copy(SEXP rval, igraph_strvector_t *sv) {
  long int i;
  sv->len=GET_LENGTH(rval);
  sv->data=Calloc(sv->len, char*);
  for (i=0; i<sv->len; i++) {
    igraph_strvector_set(sv, i, CHAR(STRING_ELT(rval, i)));
  }
  
  return 0;
}

int R_SEXP_to_vector(SEXP sv, igraph_vector_t *v) {
  v->stor_begin=REAL(sv);
  v->stor_end=v->stor_begin+GET_LENGTH(sv);
  v->end=v->stor_end;
  return 0;
}

int R_SEXP_to_vector_copy(SEXP sv, igraph_vector_t *v) {
  return igraph_vector_init_copy(v, REAL(sv), GET_LENGTH(sv));  
}

int R_SEXP_to_matrix(SEXP pakl, igraph_matrix_t *akl) {
  R_SEXP_to_vector(pakl, &akl->data);
  akl->nrow=INTEGER(GET_DIM(pakl))[0];
  akl->ncol=INTEGER(GET_DIM(pakl))[1];

  return 0;
}

int R_SEXP_to_igraph_matrix_copy(SEXP pakl, igraph_matrix_t *akl) {
  igraph_vector_init_copy(&akl->data, REAL(pakl), GET_LENGTH(pakl));
  akl->nrow=INTEGER(GET_DIM(pakl))[0];
  akl->ncol=INTEGER(GET_DIM(pakl))[1];

  return 0;
}

int R_SEXP_to_attributes(SEXP attr, igraph_attribute_list_t *al) {
  long int i, nattr;
  nattr=GET_LENGTH(VECTOR_ELT(attr, 1));
  al->len=REAL(VECTOR_ELT(attr, 0))[0];
  R_igraph_SEXP_to_strvector(VECTOR_ELT(attr, 1), &al->names);
  R_SEXP_to_vector(VECTOR_ELT(attr, 2), &al->types);
  igraph_vector_ptr_view(&al->data, (void*) R_alloc(nattr, sizeof(void*)),
		  nattr);
  for (i=0; i<nattr; i++) {
    if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
      igraph_vector_t *data=(igraph_vector_t *)R_alloc(1, sizeof(igraph_vector_t));
      VECTOR(al->data)[i] = data;
      R_SEXP_to_vector(VECTOR_ELT(attr, i+3), data);
    } else if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_STR) {
      igraph_strvector_t *data=
	(igraph_strvector_t*) R_alloc(1, sizeof(igraph_strvector_t));
      VECTOR(al->data)[i] = data;
      R_igraph_SEXP_to_strvector(VECTOR_ELT(attr, i+3), data);
    }
  }
  
  return 0;
}

int R_SEXP_to_attributes_copy(SEXP attr, igraph_attribute_list_t *al) {
  long int i, nattr;
  nattr=GET_LENGTH(VECTOR_ELT(attr, 1));
  al->len=REAL(VECTOR_ELT(attr, 0))[0];
  R_igraph_SEXP_to_strvector_copy(VECTOR_ELT(attr, 1), &al->names);
  igraph_vector_init_copy(&al->types, REAL(VECTOR_ELT(attr, 2)), 
		   GET_LENGTH(VECTOR_ELT(attr, 2)));
  igraph_vector_ptr_init(&al->data, nattr);  

  for (i=0; i<nattr; i++) {
    if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_NUM) {      
      igraph_vector_t *data=Calloc(1, igraph_vector_t);
      VECTOR(al->data)[i]=data;
      igraph_vector_init_copy(data, REAL(VECTOR_ELT(attr, i+3)),
		       GET_LENGTH(VECTOR_ELT(attr, i+3)));
    } else if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_STR) {
      igraph_strvector_t *data=Calloc(1, igraph_strvector_t);
      VECTOR(al->data)[i]=data;
      R_igraph_SEXP_to_strvector_copy(VECTOR_ELT(attr, i+3), data);
    }
  }

  return 0;
}

int R_SEXP_to_igraph(SEXP graph, igraph_t *res) {
  
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  R_SEXP_to_vector(VECTOR_ELT(graph, 2), &res->from);
  R_SEXP_to_vector(VECTOR_ELT(graph, 3), &res->to);
  R_SEXP_to_vector(VECTOR_ELT(graph, 4), &res->oi);
  R_SEXP_to_vector(VECTOR_ELT(graph, 5), &res->ii);
  R_SEXP_to_vector(VECTOR_ELT(graph, 6), &res->os);
  R_SEXP_to_vector(VECTOR_ELT(graph, 7), &res->is);
  
  return 0;
}

int R_SEXP_to_igraph_copy(SEXP graph, igraph_t *res) {
  
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  igraph_vector_init_copy(&res->from, REAL(VECTOR_ELT(graph, 2)), 
		   GET_LENGTH(VECTOR_ELT(graph, 2)));
  igraph_vector_init_copy(&res->to, REAL(VECTOR_ELT(graph, 3)), 
		   GET_LENGTH(VECTOR_ELT(graph, 3)));
  igraph_vector_init_copy(&res->oi, REAL(VECTOR_ELT(graph, 4)), 
		   GET_LENGTH(VECTOR_ELT(graph, 4)));
  igraph_vector_init_copy(&res->ii, REAL(VECTOR_ELT(graph, 5)), 
		   GET_LENGTH(VECTOR_ELT(graph, 5)));
  igraph_vector_init_copy(&res->os, REAL(VECTOR_ELT(graph, 6)), 
		   GET_LENGTH(VECTOR_ELT(graph, 6)));
  igraph_vector_init_copy(&res->is, REAL(VECTOR_ELT(graph, 7)),
		   GET_LENGTH(VECTOR_ELT(graph, 7)));

  R_SEXP_to_attributes_copy(VECTOR_ELT(graph, 8), &res->gal);
  R_SEXP_to_attributes_copy(VECTOR_ELT(graph, 9), &res->val);
  R_SEXP_to_attributes_copy(VECTOR_ELT(graph,10), &res->eal);
  
  return 0;
}

int R_SEXP_to_igraph_attr(SEXP graph, igraph_t *res) {

  R_SEXP_to_igraph(graph, res);
  
  R_SEXP_to_attributes(VECTOR_ELT(graph, 8), &res->gal);
  R_SEXP_to_attributes(VECTOR_ELT(graph, 9), &res->val);
  R_SEXP_to_attributes(VECTOR_ELT(graph,10), &res->eal);
  
  return 0;
}

/* 
 * We have only seq and vector types
 */

int R_SEXP_to_igraph_vs_copy(SEXP rit, igraph_t *graph, igraph_vs_t *it) {

  if (!strcmp(CHAR(STRING_ELT(GET_CLASS(rit),0)), "igraphvsseq")) {
    long int from=REAL(VECTOR_ELT(rit, 0))[1];
    long int to=REAL(VECTOR_ELT(rit, 0))[2];    
    igraph_vs_seq(graph, it, from, to);
  } else if (!strcmp(CHAR(STRING_ELT(GET_CLASS(rit),0)), "igraphvsvector")) {
    igraph_vector_t *tmpv=(igraph_vector_t*)R_alloc(1, sizeof(igraph_vector_t));
    igraph_vs_vectorview(graph, it, 
			 igraph_vector_view(tmpv, REAL(VECTOR_ELT(rit,1)), 
				     GET_LENGTH(VECTOR_ELT(rit,1))));
  } else {
    error("unknown vertex set type");
  }
  return 0;
}

/* 
 * We have only seq and vector types
 */

int R_SEXP_to_igraph_es_copy(SEXP rit, igraph_t *graph, igraph_es_t *it) {
  
  if (!strcmp(CHAR(STRING_ELT(GET_CLASS(rit),0)), "igraphesseq")) {
    long int from=REAL(VECTOR_ELT(rit, 0))[1];
    long int to=REAL(VECTOR_ELT(rit, 0))[2];    
    igraph_es_seq(graph, it, from, to);
  } else if (!strcmp(CHAR(STRING_ELT(GET_CLASS(rit),0)), "igraphesvector")) {
    igraph_vector_t *tmpv=(igraph_vector_t*)R_alloc(1, sizeof(igraph_vector_t));
    igraph_es_vectorview(graph, it, 
			 igraph_vector_view(tmpv, REAL(VECTOR_ELT(rit,1)), 
				     GET_LENGTH(VECTOR_ELT(rit,1))));
  } else {
    error("unknown edge set type");
  }    
  return 0;
}

/*******************************************************************/

SEXP R_igraph_empty(SEXP pn, SEXP pdirected) {
  
  SEXP result;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  igraph_t g;

  R_igraph_before();
  
  igraph_empty(&g, n, directed);  
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_edges(SEXP graph, SEXP edges) {
  
  igraph_vector_t v;			/* do NOT destroy! */
  igraph_t g;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_vector(edges, &v);
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_edges(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_vertices(SEXP graph, SEXP pnv) {
  
  integer_t nv;
  igraph_t g;
  SEXP result;
  
  R_igraph_before();
  
  nv=REAL(pnv)[0];
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_vertices(&g, nv);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_vcount(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=igraph_vcount(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_ecount(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=igraph_ecount(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_neighbors(SEXP graph, SEXP pvid, SEXP pmode) {
  
  igraph_t g;
  igraph_vector_t neis;
  SEXP result;
  real_t vid;
  integer_t mode;
  
  R_igraph_before();
  
  igraph_vector_init(&neis, 0);
  vid=REAL(pvid)[0];
  mode=REAL(pmode)[0];
  R_SEXP_to_igraph(graph, &g);
  igraph_neighbors(&g, &neis, vid, mode);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&neis)));
  igraph_vector_copy_to(&neis, REAL(result));

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_delete_edges(SEXP graph, SEXP edges) {
  
  igraph_vector_t v;
  igraph_t g;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_vector(edges, &v);
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_delete_edges(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_delete_vertices(SEXP graph, SEXP vertices) {
  
  igraph_vs_t vs;
  igraph_t g;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  R_SEXP_to_igraph_vs_copy(vertices, &g, &vs);
  igraph_delete_vertices(&g, &vs);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  igraph_vs_destroy(&vs);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_is_directed(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=igraph_is_directed(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_create(SEXP edges, SEXP pn, SEXP pdirected) {
  
  igraph_t g;
  igraph_vector_t v;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_vector(edges, &v);
  igraph_create(&g, &v, n, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_degree(SEXP graph, SEXP vids, SEXP pmode, SEXP ploops) {
  
  igraph_t g;
  igraph_vs_t vs;
  igraph_vector_t res;
  integer_t mode=REAL(pmode)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(vids, &g, &vs);
  igraph_vector_init(&res, 0);
  igraph_degree(&g, &res, &vs, mode, loops);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  igraph_vs_destroy(&vs);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_clusters(SEXP graph, SEXP pmode) {

  igraph_t g;
  integer_t mode=REAL(pmode)[0];
  igraph_vector_t membership;
  igraph_vector_t csize;
  SEXP result, names;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vector_init(&membership, 0);
  igraph_vector_init(&csize, 0);
  igraph_clusters(&g, &membership, &csize, mode);
  
  PROTECT(result=NEW_LIST(2));
  PROTECT(names=NEW_CHARACTER(2));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(igraph_vector_size(&membership)));
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(igraph_vector_size(&csize)));
  SET_STRING_ELT(names, 0, CREATE_STRING_VECTOR("membership"));
  SET_STRING_ELT(names, 1, CREATE_STRING_VECTOR("csize"));
  SET_NAMES(result, names);
  igraph_vector_copy_to(&membership, REAL(VECTOR_ELT(result, 0)));
  igraph_vector_copy_to(&csize, REAL(VECTOR_ELT(result, 1)));
  
  igraph_vector_destroy(&membership);
  igraph_vector_destroy(&csize);

  R_igraph_after();
  
  UNPROTECT(2);
  return result;
}

SEXP R_igraph_is_connected(SEXP graph, SEXP pmode) {

  igraph_t g;
  integer_t mode=REAL(pmode)[0];
  bool_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_is_connected(&g, &res, mode);
  
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=res;
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_diameter(SEXP graph, SEXP pdirected, SEXP punconnected) {
  
  igraph_t g;
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t unconnected=LOGICAL(punconnected)[0];
  integer_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_diameter(&g, &res, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_closeness(SEXP graph, SEXP pvids, SEXP pmode) {
  
  igraph_t g;
  igraph_vs_t vs;
  integer_t mode=REAL(pmode)[0];
  igraph_vector_t res;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_vector_init(&res, 0);
  igraph_closeness(&g, &res, &vs, mode);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  igraph_vs_destroy(&vs);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_subcomponent(SEXP graph, SEXP pvertex, SEXP pmode) {
  
  igraph_t g;
  real_t vertex=REAL(pvertex)[0];
  integer_t mode=REAL(pmode)[0];
  igraph_vector_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vector_init(&res, 0);
  igraph_subcomponent(&g, &res, vertex, mode);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_betweenness(SEXP graph, SEXP pvids, SEXP pdirected) {
  
  igraph_t g;
  igraph_vs_t vs;
  bool_t directed=LOGICAL(pdirected)[0];
  igraph_vector_t res;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_vector_init(&res, 0);
  igraph_betweenness(&g, &res, &vs, directed);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  igraph_vs_destroy(&vs);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_running_mean(SEXP pdata, SEXP pbinwidth) {
  
  igraph_vector_t data;
  integer_t binwidth=REAL(pbinwidth)[0];
  igraph_vector_t res;
  SEXP result;
  
  R_igraph_before();
  R_SEXP_to_vector(pdata, &data);
  
  igraph_vector_init(&res, 0);
  igraph_running_mean(&data, &res, binwidth);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_cocitation(SEXP graph, SEXP pvids) {

  igraph_t g;
  igraph_vs_t vs;
  igraph_matrix_t m;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_matrix_init(&m, 0, 0);
  igraph_cocitation(&g, &m, &vs);
  
  PROTECT(result=R_igraph_matrix_to_SEXP(&m));
  igraph_matrix_destroy(&m);
  igraph_vs_destroy(&vs);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_bibcoupling(SEXP graph, SEXP pvids) {

  igraph_t g;
  igraph_vs_t vs;
  igraph_matrix_t m;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_matrix_init(&m, 0, 0);
  igraph_bibcoupling(&g, &m, &vs);
  
  PROTECT(result=R_igraph_matrix_to_SEXP(&m));
  igraph_matrix_destroy(&m);
  igraph_vs_destroy(&vs);
  
  R_igraph_after();
  
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
  
  R_igraph_before();
  
  igraph_growing_random_game(&g, n, m, directed, citation);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_shortest_paths(SEXP graph, SEXP pvids, SEXP pmode) {
  
  igraph_t g;
  igraph_vs_t vs;
  integer_t mode=REAL(pmode)[0];
  igraph_matrix_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_matrix_init(&res, 0, 0);
  igraph_shortest_paths(&g, &res, &vs, mode);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  igraph_vs_destroy(&vs);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_lattice(SEXP pdimvector, SEXP pnei, SEXP pdirected,
		      SEXP pmutual, SEXP pcircular) {
  
  igraph_t g;
  igraph_vector_t dimvector;
  integer_t nei=REAL(pnei)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t mutual=LOGICAL(pmutual)[0];
  bool_t circular=LOGICAL(pcircular)[0];  
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_vector(pdimvector, &dimvector);
  
  igraph_lattice(&g, &dimvector, nei, directed, mutual, circular);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_barabasi_game(SEXP pn, SEXP pm, SEXP poutseq,
			    SEXP poutpref, SEXP pdirected) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t m=REAL(pm)[0]; 
  igraph_vector_t outseq;
  bool_t outpref=LOGICAL(poutpref)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_vector(poutseq, &outseq);
  
  igraph_barabasi_game(&g, n, m, &outseq, outpref, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
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
  igraph_matrix_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_kamada_kawai(&g, &res, niter, sigma, 
			     initemp, coolexp, kkconst);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_lgl(SEXP graph, SEXP pmaxiter, SEXP pmaxdelta,
			 SEXP parea, SEXP pcoolexp, SEXP prepulserad,
			 SEXP pcellsize, SEXP proot) {
  
  igraph_t g;
  igraph_matrix_t res;
  integer_t maxiter=REAL(pmaxiter)[0];
  real_t maxdelta=REAL(pmaxdelta)[0];
  real_t area=REAL(parea)[0];
  real_t coolexp=REAL(pcoolexp)[0];
  real_t repulserad=REAL(prepulserad)[0];
  real_t cellsize=REAL(pcellsize)[0];
  integer_t root=REAL(proot)[0];
  SEXP result;

  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_lgl(&g, &res, maxiter, maxdelta, area, coolexp, repulserad,
		    cellsize, root);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}  

SEXP R_igraph_layout_fruchterman_reingold_grid(SEXP graph, SEXP mat,
					       SEXP pniter,
					       SEXP pmaxdelta, SEXP parea,
					       SEXP pcoolexp, SEXP prepulserad,
					       SEXP pcellsize, SEXP puseseed) {

  igraph_t g;
  igraph_matrix_t res;
  integer_t niter=REAL(pniter)[0];
  real_t maxdelta=REAL(pmaxdelta)[0];
  real_t area=REAL(parea)[0];
  real_t coolexp=REAL(pcoolexp)[0];
  real_t repulserad=REAL(prepulserad)[0];
  real_t cellsize=REAL(pcellsize)[0];
  bool_t use_seed=LOGICAL(puseseed)[0];
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  if (use_seed) {
    R_SEXP_to_igraph_matrix_copy(mat, &res);
  } else {
    igraph_matrix_init(&res, 0, 0);
  }
  igraph_layout_grid_fruchterman_reingold(&g, &res, niter, maxdelta, area,
					  coolexp, repulserad, cellsize, 
					  use_seed);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_minimum_spanning_tree_unweighted(SEXP graph) {
  
  igraph_t g;
  igraph_t mst;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_minimum_spanning_tree_unweighted(&g, &mst);
  PROTECT(result=R_igraph_to_SEXP(&mst));
  igraph_destroy(&mst);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_minimum_spanning_tree_prim(SEXP graph, SEXP pweights) {
  
  igraph_t g;
  igraph_t mst;
  igraph_vector_t weights;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_vector(pweights, &weights);
  
  R_SEXP_to_igraph(graph, &g);
  igraph_minimum_spanning_tree_prim(&g, &mst, &weights);
  PROTECT(result=R_igraph_to_SEXP(&mst));
  igraph_destroy(&mst);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_edge_betweenness(SEXP graph, SEXP pdirected) {
  
  igraph_t g;
  igraph_vector_t res;
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vector_init(&res, 0);
  igraph_edge_betweenness(&g, &res, directed);
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage(SEXP graph, SEXP pst, SEXP pagebins,
				     SEXP pmaxind, SEXP plsd) {
  
  igraph_t g;
  igraph_matrix_t akl, sd;
  igraph_vector_t st;
  integer_t agebins=REAL(pagebins)[0];
  integer_t maxind=REAL(pmaxind)[0];
  bool_t lsd=LOGICAL(plsd)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_vector(pst, &st);

  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&akl, 0, 0);
  igraph_matrix_init(&sd, 0, 0);
  igraph_measure_dynamics_idage(&g, &akl, &sd, &st, agebins, maxind, lsd);
  
  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, R_igraph_matrix_to_SEXP(&akl));
  igraph_matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_igraph_matrix_to_SEXP(&sd));
  igraph_matrix_destroy(&sd);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage_debug(SEXP graph, SEXP pst, 
					   SEXP pagebins, SEXP pmaxind, 
					   SEXP plsd, SEXP pest_ind, 
					   SEXP pest_age) {
  igraph_t g;
  igraph_matrix_t akl, sd;
  igraph_vector_t st;
  integer_t agebins=REAL(pagebins)[0];
  integer_t maxind=REAL(pmaxind)[0];
  bool_t lsd=LOGICAL(plsd)[0];
  integer_t est_ind=REAL(pest_ind)[0];
  integer_t est_age=REAL(pest_age)[0];
  igraph_vector_t estimates;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_vector(pst, &st);
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&akl, 0, 0);
  igraph_matrix_init(&sd, 0, 0);
  igraph_vector_init(&estimates, 0);
  igraph_measure_dynamics_idage_debug(&g, &akl, &sd, &st, agebins, maxind, lsd,
				      &estimates, est_ind, est_age);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_igraph_matrix_to_SEXP(&akl));
  igraph_matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_igraph_matrix_to_SEXP(&sd));
  igraph_matrix_destroy(&sd);
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(igraph_vector_size(&estimates)));
  igraph_vector_copy_to(&estimates, REAL(VECTOR_ELT(result, 2)));
  igraph_vector_destroy(&estimates);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_measure_dynamics_idage_st(SEXP graph, SEXP pakl) {
  
  igraph_t g;
  igraph_matrix_t akl;
  igraph_vector_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_matrix(pakl, &akl);
  igraph_vector_init(&res, 0);
  igraph_measure_dynamics_idage_st(&g, &res, &akl);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  
  igraph_vector_destroy(&res);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_shortest_paths(SEXP graph, SEXP pfrom, SEXP pmode) {

  igraph_t g;
  integer_t from=REAL(pfrom)[0];
  integer_t mode=REAL(pmode)[0];
  igraph_vector_t *vects;
  long int no_of_nodes, i;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  no_of_nodes=igraph_vcount(&g);
  vects=Calloc(no_of_nodes, igraph_vector_t);
  for (i=0; i<no_of_nodes; i++) {
    igraph_vector_init(&vects[i], 0);
  }
  igraph_get_shortest_paths(&g, vects, from, mode);
  PROTECT(result=NEW_LIST(no_of_nodes));
  for (i=0; i<no_of_nodes; i++) {
    SET_VECTOR_ELT(result, i, NEW_NUMERIC(igraph_vector_size(&vects[i])));
    igraph_vector_copy_to(&vects[i], REAL(VECTOR_ELT(result, i)));
    igraph_vector_destroy(&vects[i]);
  }
  
  Free(vects);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_are_connected(SEXP graph, SEXP pv1, SEXP pv2) {
  
  igraph_t g;
  integer_t v1=REAL(pv1)[0];
  integer_t v2=REAL(pv2)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=igraph_are_connected(&g, v1, v2);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_graph_adjacency(SEXP adjmatrix, SEXP pmode) {
  
  igraph_t g;
  igraph_matrix_t adjm;
  integer_t mode=REAL(pmode)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_matrix(adjmatrix, &adjm);
  igraph_adjacency(&g, &adjm, mode);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
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

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_average_path_length(&g, &res, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_star(SEXP pn, SEXP pmode, SEXP pcenter) {

  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t mode=REAL(pmode)[0];
  integer_t center=REAL(pcenter)[0];
  SEXP result;
  
  R_igraph_before();
  
  igraph_star(&g, n, mode, center);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
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
  
  R_igraph_before();
  
  igraph_ring(&g, n, directed, mutual, circular);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_tree(SEXP pn, SEXP pchildren, SEXP pmode) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t children=REAL(pchildren)[0];
  integer_t mode=REAL(pmode)[0];
  SEXP result;

  R_igraph_before();
  
  igraph_tree(&g, n, children, mode);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_subgraph(SEXP graph, SEXP pvids) {
  
  igraph_t g;
  igraph_t sub;
  igraph_vs_t vs;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_subgraph(&g, &sub, &vs);
  PROTECT(result=R_igraph_to_SEXP(&sub));
  igraph_destroy(&sub);
  igraph_vs_destroy(&vs);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_random(SEXP graph) {
  
  igraph_t g;
  igraph_matrix_t res;
  SEXP result=R_NilValue;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_random(&g, &res);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_circle(SEXP graph) {
  
  igraph_t g;
  igraph_matrix_t res;
  SEXP result=R_NilValue;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_circle(&g, &res);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  R_igraph_after();
  
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
  
  R_igraph_before();
  
  igraph_erdos_renyi_game(&g, type, n, porm, directed, loops);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_full(SEXP pn, SEXP pdirected, SEXP ploops) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;
  
  R_igraph_before();
  
  igraph_full(&g, n, directed, loops);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_random_sample(SEXP plow, SEXP phigh, SEXP plength) {
  
  igraph_vector_t res;
  integer_t low=REAL(plow)[0];
  integer_t high=REAL(phigh)[0];
  integer_t length=REAL(plength)[0];
  SEXP result;

  R_igraph_before();
  
  igraph_vector_init(&res, 0);
  igraph_random_sample(&res, low, high, length);
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_edgelist(SEXP graph, SEXP pbycol) {
  
  igraph_t g;
  igraph_vector_t res;
  bool_t bycol=LOGICAL(pbycol)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vector_init(&res, 0);
  igraph_get_edgelist(&g, &res, bycol);
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_get_adjacency(SEXP graph, SEXP ptype) {
  
  igraph_t g;
  igraph_matrix_t res;
  integer_t type=REAL(ptype)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_get_adjacency(&g, &res, type);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_simplify(SEXP graph, SEXP pmultiple, SEXP ploops) {
  
  igraph_t g;
  bool_t multiple=LOGICAL(pmultiple)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_simplify(&g, multiple, loops);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
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
  igraph_matrix_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_fruchterman_reingold(&g, &res, niter, maxdelta, area, 
				     coolexp, repulserad, 0);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_degree_sequence_game(SEXP pout_seq, SEXP pin_seq, 
				   SEXP pmethod) {
  igraph_t g;
  igraph_vector_t outseq;
  igraph_vector_t inseq;
  integer_t method=REAL(pmethod)[0];
  SEXP result;

  R_igraph_before();

  R_SEXP_to_vector(pout_seq, &outseq);
  R_SEXP_to_vector(pin_seq, &inseq);
  igraph_degree_sequence_game(&g, &outseq, &inseq, method);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_transitivity(SEXP graph, SEXP ptype) {
  
  igraph_t g;
  igraph_vector_t res;
  integer_t type=REAL(ptype)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vector_init(&res, 0);
  igraph_transitivity(&g, &res, type);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_add_graph_attribute(SEXP graph, SEXP pname, SEXP ptype) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  igraph_attribute_type_t type=REAL(ptype)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_graph_attribute(&g, name, type);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_remove_graph_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_remove_graph_attribute(&g, name);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_graph_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  void *value;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  igraph_get_graph_attribute(&g, name, &value, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    PROTECT(result=NEW_NUMERIC(1));
    REAL(result)[0]=*(real_t*)value;
  } else {
    PROTECT(result=NEW_CHARACTER(1));
    SET_STRING_ELT(result, 0, CREATE_STRING_VECTOR((char*)value));
  }
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_graph_attribute(SEXP graph, SEXP pname, SEXP pvalue) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  void *value;
  igraph_attribute_type_t type;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_get_graph_attribute_type(&g, name, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    value=REAL(AS_NUMERIC(pvalue));
  } else {
    value=CHAR(STRING_ELT(AS_CHARACTER(pvalue), 0));
  }
  igraph_set_graph_attribute(&g, name, value);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_vertex_attribute(SEXP graph, SEXP pname, SEXP ptype) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  igraph_attribute_type_t type=REAL(ptype)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_vertex_attribute(&g, name, type);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_remove_vertex_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_remove_vertex_attribute(&g, name);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_vertex_attribute(SEXP graph, SEXP pname, SEXP pv) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  long int v=REAL(pv)[0];
  void *value;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  igraph_get_vertex_attribute(&g, name, v, &value, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    PROTECT(result=NEW_NUMERIC(1));
    REAL(result)[0]=*(real_t*)value;
  } else {
    PROTECT(result=NEW_CHARACTER(1));
    SET_STRING_ELT(result, 0, CREATE_STRING_VECTOR((char*)value));
  }    
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_vertex_attribute(SEXP graph, SEXP pname, SEXP pv, 
				   SEXP pvalue) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  long int v=REAL(pv)[0];
  void *value;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_get_vertex_attribute_type(&g, name, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    value=REAL(AS_NUMERIC(pvalue));
  } else {
    value=CHAR(STRING_ELT(AS_CHARACTER(pvalue), 0));
  }    
  igraph_set_vertex_attribute(&g, name, v, value);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_vertex_attributes(SEXP graph, SEXP pname, SEXP pv) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  igraph_vs_t vs;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();

  R_SEXP_to_igraph_attr(graph, &g);
  R_SEXP_to_igraph_vs_copy(pv, &g, &vs);

  igraph_get_vertex_attribute_type(&g, name, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t value;
    void *valueptr=&value;
    igraph_vector_init(&value, 0);
    igraph_get_vertex_attributes(&g, name, &vs, &valueptr);
    PROTECT(result=NEW_NUMERIC(igraph_vector_size(&value)));
    igraph_vector_copy_to(&value, REAL(result));
    igraph_vector_destroy(&value);
  } else {
    igraph_strvector_t value;
    void *valueptr=&value;
    igraph_strvector_init(&value, 0);
    igraph_get_vertex_attributes(&g, name, &vs, &valueptr);
    PROTECT(result=R_igraph_strvector_to_SEXP(&value));
    igraph_strvector_destroy(&value);
  }
  igraph_vs_destroy(&vs);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_vertex_attributes(SEXP graph, SEXP pname, SEXP pv, 
				    SEXP pvalue) {
  
  igraph_t g;
  igraph_vs_t vs;
  const char *name=CHAR(STRING_ELT(pname, 0));
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  R_SEXP_to_igraph_vs_copy(pv, &g, &vs);
  igraph_get_vertex_attribute_type(&g, name, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) { 
    igraph_vector_t value;
    R_SEXP_to_vector(AS_NUMERIC(pvalue), &value);
    igraph_set_vertex_attributes(&g, name, &vs, &value);
  } else {
    igraph_strvector_t value;
    R_igraph_SEXP_to_strvector_copy(AS_CHARACTER(pvalue), &value);
    igraph_set_vertex_attributes(&g, name, &vs, &value);
  }
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  igraph_vs_destroy(&vs);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_list_graph_attributes(SEXP graph) {
  igraph_t g;
  igraph_strvector_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  igraph_strvector_init(&res, 0);
  igraph_list_graph_attributes(&g, &res, 0);
  PROTECT(result=R_igraph_strvector_to_SEXP(&res));
  igraph_strvector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_list_vertex_attributes(SEXP graph) {
  igraph_t g;
  igraph_strvector_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  igraph_strvector_init(&res, 0);
  igraph_list_vertex_attributes(&g, &res, 0);
  PROTECT(result=R_igraph_strvector_to_SEXP(&res));
  igraph_strvector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  

}

SEXP R_igraph_vs_all(SEXP graph) {
  
  igraph_t g;
  igraph_vs_t it;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vs_all(&g, &it);
  PROTECT(result=R_igraph_vs_to_SEXP(&g, &it));
  igraph_vs_destroy(&it);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_all(SEXP graph) {

  igraph_t g;
  igraph_es_t it;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_es_all(&g, &it);
  PROTECT(result=R_igraph_es_to_SEXP(&g, &it));
  igraph_es_destroy(&it);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_fromorder(SEXP graph) {

  igraph_t g;
  igraph_es_t it;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_es_fromorder(&g, &it);
  PROTECT(result=R_igraph_es_to_SEXP(&g, &it));
  igraph_es_destroy(&it);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_adj(SEXP graph, SEXP pvid, SEXP pmode) {

  igraph_t g;
  igraph_es_t it;
  integer_t vid=REAL(pvid)[0];
  integer_t mode=REAL(pmode)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_es_adj(&g, &it, vid, mode);
  PROTECT(result=R_igraph_es_to_SEXP(&g, &it));
  igraph_es_destroy(&it);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_vs_adj(SEXP graph, SEXP pvid, SEXP pmode) {

  igraph_t g;
  igraph_vs_t it;
  integer_t vid=REAL(pvid)[0];
  integer_t mode=REAL(pmode)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vs_adj(&g, &it, vid, mode);
  PROTECT(result=R_igraph_vs_to_SEXP(&g, &it));
  igraph_vs_destroy(&it);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_vs_vector(SEXP graph, SEXP vec) {
  
  igraph_t g;
  igraph_vs_t it;
  igraph_vector_t v;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_vector(vec, &v);
  igraph_vs_vectorview(&g, &it, &v);
  PROTECT(result=R_igraph_vs_to_SEXP(&g, &it));
  igraph_vs_destroy(&it);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_vector(SEXP graph, SEXP vec) {
  
  igraph_t g;
  igraph_es_t it;
  igraph_vector_t v;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_vector(vec, &v);
  igraph_es_vectorview(&g, &it, &v);
  PROTECT(result=R_igraph_es_to_SEXP(&g, &it));
  igraph_es_destroy(&it);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_vs_next(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_vs_t it;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pit, &g, &it);
  igraph_vs_next(&g, &it);
  PROTECT(result=R_igraph_vs_to_SEXP(&g, &it));
  igraph_vs_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_next(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_es_t it;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_es_copy(pit, &g, &it);
  R_SEXP_to_igraph(graph, &g);
  igraph_es_next(&g, &it);
  PROTECT(result=R_igraph_es_to_SEXP(&g, &it));
  igraph_es_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_vs_reset(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_vs_t it;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pit, &g, &it);
  igraph_vs_reset(&g, &it);
  PROTECT(result=R_igraph_vs_to_SEXP(&g, &it));
  igraph_vs_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_reset(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_es_t it;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_es_copy(pit, &g, &it);
  igraph_es_reset(&g, &it);
  PROTECT(result=R_igraph_es_to_SEXP(&g, &it));
  igraph_es_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_vs_end(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_vs_t it;
  bool_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pit, &g, &it);
  res=igraph_vs_end(&g, &it);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=res;
  igraph_vs_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_end(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_es_t it;
  bool_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_es_copy(pit, &g, &it);
  res=igraph_es_end(&g, &it);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=res;
  igraph_es_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_vs_get(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_vs_t it;
  integer_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pit, &g, &it);
  res=igraph_vs_get(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  igraph_vs_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_get(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_es_t it;
  integer_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_es_copy(pit, &g, &it);
  res=igraph_es_get(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  igraph_es_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_from(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_es_t it;
  integer_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_es_copy(pit, &g, &it);
  res=igraph_es_from(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  igraph_es_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_es_to(SEXP graph, SEXP pit) {

  igraph_t g;
  igraph_es_t it;
  integer_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_es_copy(pit, &g, &it);
  res=igraph_es_to(&g, &it);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  igraph_es_destroy(&it);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_edge_attribute(SEXP graph, SEXP pname, SEXP ptype) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  igraph_attribute_type_t type=REAL(ptype)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_edge_attribute(&g, name, type);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_remove_edge_attribute(SEXP graph, SEXP pname) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_remove_edge_attribute(&g, name);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_edge_attribute(SEXP graph, SEXP pname, SEXP pv) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  long int v=REAL(pv)[0];
  void *value;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  igraph_get_edge_attribute(&g, name, v, &value, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    PROTECT(result=NEW_NUMERIC(1));
    REAL(result)[0]=*(real_t*)value;
  } else {
    PROTECT(result=NEW_CHARACTER(1));
    SET_STRING_ELT(result, 0, CREATE_STRING_VECTOR((char*)value));
  }    
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_edge_attribute(SEXP graph, SEXP pname, SEXP pv, 
				 SEXP pvalue) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  long int v=REAL(pv)[0];
  void *value;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_get_edge_attribute_type(&g, name, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    value=REAL(AS_NUMERIC(pvalue));
  } else {
    value=CHAR(STRING_ELT(AS_CHARACTER(pvalue), 0));
  }  
  igraph_set_edge_attribute(&g, name, v, value);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_edge_attributes(SEXP graph, SEXP pname, SEXP pv) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  igraph_es_t es;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();

  R_SEXP_to_igraph_attr(graph, &g);
  R_SEXP_to_igraph_es_copy(pv, &g, &es);

  igraph_get_edge_attribute_type(&g, name, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t value;
    void *valueptr=&value;
    igraph_vector_init(&value, 0);
    igraph_get_edge_attributes(&g, name, &es, &valueptr);
    PROTECT(result=NEW_NUMERIC(igraph_vector_size(&value)));
    igraph_vector_copy_to(&value, REAL(result));
    igraph_vector_destroy(&value);
  } else {
    igraph_strvector_t value;
    void *valueptr=&value;
    igraph_strvector_init(&value, 0);
    igraph_get_edge_attributes(&g, name, &es, &valueptr);
    PROTECT(result=R_igraph_strvector_to_SEXP(&value));
    igraph_strvector_destroy(&value);
  }    
  igraph_es_destroy(&es);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_set_edge_attributes(SEXP graph, SEXP pname, SEXP pv, 
				  SEXP pvalue) {
  
  igraph_t g;
  const char *name=CHAR(STRING_ELT(pname, 0));
  igraph_es_t es;
  igraph_attribute_type_t type;
  SEXP result;

  R_igraph_before();

  R_SEXP_to_igraph_copy(graph, &g);
  R_SEXP_to_igraph_es_copy(pv, &g, &es);

  igraph_get_edge_attribute_type(&g, name, &type);
  if (type==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t value;
    R_SEXP_to_vector(AS_NUMERIC(pvalue), &value);
    igraph_set_edge_attributes(&g, name, &es, &value);
  } else {
    igraph_strvector_t value;
    R_igraph_SEXP_to_strvector(AS_CHARACTER(pvalue), &value);
    igraph_set_edge_attributes(&g, name, &es, &value);
  }    
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  igraph_es_destroy(&es);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_list_edge_attributes(SEXP graph) {
  igraph_t g;
  igraph_strvector_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  igraph_strvector_init(&res, 0);
  igraph_list_edge_attributes(&g, &res, 0);
  PROTECT(result=R_igraph_strvector_to_SEXP(&res));
  igraph_strvector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  

}

SEXP R_igraph_read_graph_edgelist(SEXP pvfile, SEXP pn, SEXP pdirected) {
  igraph_t g;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  FILE *file;
  SEXP result;
  
  R_igraph_before();
  
#if HAVE_FMEMOPEN == 1
  file=fmemopen(RAW(pvfile), GET_LENGTH(pvfile), "r");
#else 
  file=fopen(CHAR(STRING_ELT(pvfile, 0)), "r");
#endif
  if (file==0) { igraph_error("Cannot read edgelist", __FILE__, __LINE__,
			      IGRAPH_EFILE); }
  igraph_read_graph_edgelist(&g, file, n, directed);
  fclose(file);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_write_graph_edgelist(SEXP graph, SEXP file) {
  igraph_t g;
  FILE *stream;
  char *bp;
  size_t size;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
#if HAVE_OPEN_MEMSTREAM == 1
  stream=open_memstream(&bp, &size);
#else
  stream=fopen(CHAR(STRING_ELT(file, 0)), "w");
#endif
  if (stream==0) { igraph_error("Cannot write edgelist", __FILE__, __LINE__,
				IGRAPH_EFILE); }
  igraph_write_graph_edgelist(&g, stream);
  fclose(stream);
#if HAVE_OPEN_MEMSTREAM == 1
  PROTECT(result=allocVector(RAWSXP, size));
  memcpy(RAW(result), bp, sizeof(char)*size);
  free(bp);
#else 
  PROTECT(result=NEW_NUMERIC(0));
#endif
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_read_graph_ncol(SEXP pvfile, SEXP pnames, SEXP pweights) {
  igraph_t g;
  bool_t names=LOGICAL(pnames)[0];
  bool_t weights=LOGICAL(pweights)[0];
  FILE *file;
  SEXP result;
  
  R_igraph_before();
  
#if HAVE_FMEMOPEN == 1
  file=fmemopen(RAW(pvfile), GET_LENGTH(pvfile), "r");
#else 
  file=fopen(CHAR(STRING_ELT(pvfile, 0)), "r");
#endif
  if (file==0) { igraph_error("Cannot read edgelist", __FILE__, __LINE__,
			      IGRAPH_EFILE); }
  igraph_read_graph_ncol(&g, file, names, weights);
  fclose(file);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_write_graph_ncol(SEXP graph, SEXP file, SEXP pnames, 
			       SEXP pweights) {
  igraph_t g;
  FILE *stream;
  char *bp;
  size_t size;
  const char *names, *weights;
  SEXP result;

  R_igraph_before();
  
  if (isNull(pnames)) {
    names=0; 
  } else {
    names=CHAR(STRING_ELT(pnames, 0));
  } 
  if (isNull(pweights)) {
    weights=0; 
  } else {
    weights=CHAR(STRING_ELT(pweights, 0));
  }   

  R_SEXP_to_igraph_attr(graph, &g);
#if HAVE_OPEN_MEMSTREAM == 1
  stream=open_memstream(&bp, &size);
#else 
  stream=fopen(CHAR(STRING_ELT(file,0)), "w");
#endif
  if (stream==0) { igraph_error("Cannot write .ncol file", __FILE__, __LINE__,
				IGRAPH_EFILE); }
  igraph_write_graph_ncol(&g, stream, names, weights);
  fclose(stream);
#if HAVE_OPEN_MEMSTREAM == 1
  PROTECT(result=allocVector(RAWSXP, size));
  memcpy(RAW(result), bp, sizeof(char)*size);
  free(bp);
#else
  PROTECT(result=NEW_NUMERIC(0));
#endif
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_read_graph_lgl(SEXP pvfile, SEXP pnames, SEXP pweights) {
  igraph_t g;
  bool_t names=LOGICAL(pnames)[0];
  bool_t weights=LOGICAL(pweights)[0];
  FILE *file;
  SEXP result;
  
  R_igraph_before();
  
#if HAVE_FMEMOPEN == 1
  file=fmemopen(RAW(pvfile), GET_LENGTH(pvfile), "r");
#else 
  file=fopen(CHAR(STRING_ELT(pvfile, 0)), "r");
#endif
  if (file==0) { igraph_error("Cannot read edgelist", __FILE__, __LINE__,
			      IGRAPH_EFILE); }
  igraph_read_graph_lgl(&g, file, names, weights);
  fclose(file);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_write_graph_lgl(SEXP graph, SEXP file, SEXP pnames, 
			      SEXP pweights, SEXP pisolates) {
  igraph_t g;
  FILE *stream;
  char *bp;
  size_t size;
  const char *names, *weights;
  bool_t isolates=LOGICAL(pisolates)[0];
  SEXP result;

  R_igraph_before();
  
  if (isNull(pnames)) {
    names=0; 
  } else {
    names=CHAR(STRING_ELT(pnames, 0));
  } 
  if (isNull(pweights)) {
    weights=0; 
  } else {
    weights=CHAR(STRING_ELT(pweights, 0));
  }   

  R_SEXP_to_igraph_attr(graph, &g);
#if HAVE_OPEN_MEMSTREAM == 1
  stream=open_memstream(&bp, &size);
#else
  stream=fopen(CHAR(STRING_ELT(file, 0)), "w");
#endif
  igraph_write_graph_lgl(&g, stream, names, weights, isolates);
  fclose(stream);
#if HAVE_OPEN_MEMSTREAM == 1
  PROTECT(result=allocVector(RAWSXP, size));
  memcpy(RAW(result), bp, sizeof(char)*size);
  free(bp);
#else
  PROTECT(result=NEW_NUMERIC(0));
#endif
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

/* SEXP R_igraph_iterator_randomwalk(SEXP graph, SEXP pvid, SEXP pmode) { */
  
/*   igraph_t g; */
/*   igraph_iterator_t it; */
/*   integer_t vid=REAL(pvid)[0]; */
/*   integer_t mode=REAL(pmode)[0]; */
/*   SEXP result;  */
  
/*   R_igraph_before(); */
  
/*   R_SEXP_to_igraph(graph, &g); */
/*   igraph_iterator_randomwalk(&g, &it, vid, mode); */
/*   PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it)); */
/*   igraph_iterator_destroy(&g, &it); */
  
/*   R_igraph_after(); */
  
/*   UNPROTECT(1); */
/*   return result; */
/* } */
  
/* SEXP R_igraph_iterator_randomwalk1(SEXP graph, SEXP pvid, SEXP pmode) { */
  
/*   igraph_t g; */
/*   igraph_iterator_t it; */
/*   integer_t vid=REAL(pvid)[0]; */
/*   integer_t mode=REAL(pmode)[0]; */
/*   SEXP result;  */
  
/*   R_igraph_before(); */
  
/*   R_SEXP_to_igraph(graph, &g); */
/*   igraph_iterator_randomwalk1(&g, &it, vid, mode); */
/*   PROTECT(result=R_igraph_iterator_to_SEXP(&g, &it)); */
/*   igraph_iterator_destroy(&g, &it); */
  
/*   R_igraph_after(); */
  
/*   UNPROTECT(1); */
/*   return result; */
/* } */
  
SEXP R_igraph_decompose(SEXP graph, SEXP pmode, SEXP pmaxcompno, 
			SEXP pminelements) {

  igraph_t g;
  integer_t mode=REAL(pmode)[0];
  integer_t maxcompno=REAL(pmaxcompno)[0];
  integer_t minelements=REAL(pminelements)[0];
  igraph_vector_ptr_t comps;
  SEXP result;
  long int i;
  
  R_igraph_before();
  
  R_SEXP_to_igraph_attr(graph, &g);
  igraph_vector_ptr_init(&comps, 0);
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &comps);
  igraph_decompose(&g, &comps, mode, maxcompno, minelements);
  PROTECT(result=NEW_LIST(igraph_vector_ptr_size(&comps)));
  for (i=0; i<igraph_vector_ptr_size(&comps); i++) {
    SET_VECTOR_ELT(result, i, R_igraph_to_SEXP(VECTOR(comps)[i]));
    igraph_destroy(VECTOR(comps)[i]);
    igraph_free(VECTOR(comps)[i]);
  }
  igraph_vector_ptr_destroy(&comps);
  IGRAPH_FINALLY_CLEAN(1);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_atlas(SEXP pno) {
  
  int no=REAL(pno)[0];
  igraph_t g;
  SEXP result;
  
  R_igraph_before();

  igraph_atlas(&g, no);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
