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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA 

*/

#include "igraph.h"
#include "error.h"

#include "config.h"

#define USE_RINTERNALS
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <stdio.h>

SEXP R_igraph_matrix_to_SEXP(igraph_matrix_t *m);
SEXP R_igraph_strvector_to_SEXP(igraph_strvector_t *m);
SEXP R_igraph_to_SEXP(igraph_t *graph);

int R_igraph_SEXP_to_strvector(SEXP rval, igraph_strvector_t *sv);
int R_igraph_SEXP_to_strvector_copy(SEXP rval, igraph_strvector_t *sv);
int R_SEXP_to_vector(SEXP sv, igraph_vector_t *v);
int R_SEXP_to_vector_copy(SEXP sv, igraph_vector_t *v);
int R_SEXP_to_matrix(SEXP pakl, igraph_matrix_t *akl);
int R_SEXP_to_igraph_matrix_copy(SEXP pakl, igraph_matrix_t *akl);
int R_SEXP_to_igraph(SEXP graph, igraph_t *res);
int R_SEXP_to_igraph_copy(SEXP graph, igraph_t *res);
int R_SEXP_to_igraph_vs_copy(SEXP rit, igraph_t *graph, igraph_vs_t *it);
int R_SEXP_to_igraph_es_copy(SEXP rit, igraph_t *graph, igraph_es_t *it);

/* get the list element named str, or return NULL */
/* from the R Manual */

SEXP R_igraph_getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

/******************************************************
 * Attributes                                         *
 *****************************************************/
int R_igraph_attribute_init(igraph_t *graph) {
  SEXP result;
  long int i;
  PROTECT(result=NEW_LIST(4));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(2));
  REAL(VECTOR_ELT(result, 0))[0]=0; /* R objects */
  REAL(VECTOR_ELT(result, 0))[1]=1; /* igraph_t objects */
  for (i=0; i<3; i++) {
    SET_VECTOR_ELT(result, i+1, NEW_LIST(0)); /* gal, val, eal */
  }
  graph->attr=result;
  return 0;
}

void R_igraph_attribute_destroy(igraph_t *graph) {
  SEXP attr=graph->attr;
  REAL(VECTOR_ELT(attr, 0))[1] -= 1; /* refcount for igraph_t */
  if (REAL(VECTOR_ELT(attr, 0))[1]==0) {
    UNPROTECT_PTR(attr);
  }
  graph->attr=0;
}

int R_igraph_attribute_copy(igraph_t *to, const igraph_t *from) {
  SEXP fromattr=from->attr;
  to->attr=from->attr;
  REAL(VECTOR_ELT(fromattr, 0))[1] += 1; /* refcount only */
  if (REAL(VECTOR_ELT(fromattr, 0))[1] == 1) {
    PROTECT(to->attr);
  }
  return 0;
}

int R_igraph_attribute_add_vertices(igraph_t *graph, long int nv, 
				    igraph_vector_ptr_t *nattr) {
  SEXP attr=graph->attr;
  SEXP val, rep=0, names, newnames;
  igraph_vector_t news;
  long int valno, i, origlen, nattrno, newattrs;
  if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
    SEXP newattr;
    PROTECT(newattr=duplicate(attr));
    REAL(VECTOR_ELT(attr, 0))[1] -= 1;
    if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
      UNPROTECT_PTR(attr);
    }
    REAL(VECTOR_ELT(newattr, 0))[0] = 0;
    REAL(VECTOR_ELT(newattr, 0))[1] = 1;
    attr=graph->attr=newattr;
  }

  val=VECTOR_ELT(attr, 2);
  valno=GET_LENGTH(val);  
  names=GET_NAMES(val);
  if (nattr==NULL) { 
    nattrno=0;
  } else {
    nattrno=igraph_vector_ptr_size(nattr); 
  }
  origlen=igraph_vcount(graph);

  /* First add the new attributes, if any */
  newattrs=0;
  IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
  for (i=0; i<nattrno; i++) {
    igraph_i_attribute_record_t *nattr_entry=VECTOR(*nattr)[i];
    const char *nname=nattr_entry->name;
    long int j; 
    igraph_bool_t l=0;
    for (j=0; !l && j<valno; j++) {
      l=!strcmp(nname, CHAR(STRING_ELT(names, j)));
    }
    if (!l) {
      newattrs++;
      IGRAPH_CHECK(igraph_vector_push_back(&news, i));
    }
  }
  if (newattrs != 0) {
    SEXP app, newval;
    PROTECT(app=NEW_LIST(newattrs));
    PROTECT(newnames=NEW_CHARACTER(newattrs));
    PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL), 
			   ScalarInteger(origlen))));
    for (i=0; i<newattrs; i++) {
      igraph_i_attribute_record_t *tmp=
	VECTOR(*nattr)[(long int)VECTOR(news)[i]];
      SET_VECTOR_ELT(app, i, rep);
      SET_STRING_ELT(newnames, i, CREATE_STRING_VECTOR(tmp->name));
    }
    UNPROTECT(1); 		/* rep */
    PROTECT(newval=EVAL(lang3(install("c"), val, app)));
    PROTECT(newnames=EVAL(lang3(install("c"), names, newnames)));
    SET_NAMES(newval, newnames);
    SET_VECTOR_ELT(attr, 2, newval);
    val=VECTOR_ELT(attr, 2);    
    valno=GET_LENGTH(val);  
    names=GET_NAMES(val);
    UNPROTECT(4);
    rep=0;
  }
  igraph_vector_destroy(&news);
  IGRAPH_FINALLY_CLEAN(1);	/* news */

  /* Now append the new values */
  for (i=0; i<valno; i++) {
    SEXP oldva=VECTOR_ELT(val, i), newva;
    const char *sexpname=CHAR(STRING_ELT(names,i));
    igraph_bool_t l=0;
    long int j;
    for (j=0; !l && j<nattrno; j++) {
      igraph_i_attribute_record_t *tmp=VECTOR(*nattr)[j];
      l=!strcmp(sexpname, tmp->name);
    }
    if (l) {
      /* This attribute is present in nattr */
      SEXP app;
      igraph_i_attribute_record_t *tmprec=VECTOR(*nattr)[j-1];
      switch (tmprec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	if (nv != igraph_vector_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=NEW_NUMERIC(nv));
	igraph_vector_copy_to(tmprec->value, REAL(app));
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	if (nv != igraph_strvector_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=R_igraph_strvector_to_SEXP(tmprec->value));
	break;
      case IGRAPH_ATTRIBUTE_R_OBJECT:
	/* TODO */
	IGRAPH_ERROR("R_objects not implemented yet", IGRAPH_UNIMPLEMENTED);
	break;
      default:
	warning("Ignoring unknown attribute type");
	break;
      }
      PROTECT(newva=EVAL(lang3(install("c"), oldva, app)));
      SET_VECTOR_ELT(val, i, newva);
      UNPROTECT(2);		/* app & newva */
    } else {
      /* No such attribute, append NA's */
      if (rep==0) {
	PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL), 
			       ScalarInteger(nv))));
      }
      PROTECT(newva=EVAL(lang3(install("c"), oldva, rep)));
      SET_VECTOR_ELT(val, i, newva);
      UNPROTECT(1); 		/* newva */
    }
  }
  if (rep != 0) {
    UNPROTECT(1);
  } 
  
  return 0;
}

void R_igraph_attribute_delete_vertices(igraph_t *graph, 
					const igraph_vector_t *eidx,
					const igraph_vector_t *vidx) {
  SEXP attr=graph->attr;
  SEXP eal, val;
  long int valno, ealno, i;
  if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
    SEXP newattr;
    PROTECT(newattr=duplicate(attr));
    REAL(VECTOR_ELT(attr, 0))[1] -= 1;
    if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
      UNPROTECT_PTR(attr);
    }
    REAL(VECTOR_ELT(newattr, 0))[0] = 0;
    REAL(VECTOR_ELT(newattr, 0))[1] = 1;
    attr=graph->attr=newattr;
  }

  /* Vertices */
  val=VECTOR_ELT(attr, 2);
  valno=GET_LENGTH(val);
  for (i=0; i<valno; i++) {
    SEXP oldva=VECTOR_ELT(val, i), newva, ss;
    long int origlen=GET_LENGTH(oldva);
    long int newlen=0, j;
    for (j=0; j<igraph_vector_size(vidx); j++) {
      if (VECTOR(*vidx)[j] > 0) {
	newlen++;
      }
    }
    PROTECT(ss=NEW_NUMERIC(newlen));
    for (j=0; j<origlen; j++) {
      if (VECTOR(*vidx)[j]>0) {
	REAL(ss)[(long int)VECTOR(*vidx)[j]-1]=j+1;
      }
    }
    PROTECT(newva=EVAL(lang3(install("["), oldva, ss)));
    SET_VECTOR_ELT(val, i, newva);
    UNPROTECT(2);
  }    

  /* Edges */
  eal=VECTOR_ELT(attr, 3);
  ealno=GET_LENGTH(eal);
  for (i=0; i<ealno; i++) {
    SEXP oldea=VECTOR_ELT(eal, i), newea, ss;
    long int origlen=GET_LENGTH(oldea);
    long int newlen=0, j;
    /* calculate new length */
    for (j=0; j<origlen; j++) {
      if (VECTOR(*eidx)[j] > 0) {
	newlen++;
      }
    }    
    PROTECT(ss=NEW_NUMERIC(newlen));
    for (j=0; j<origlen; j++) {
      if (VECTOR(*eidx)[j]>0) {
	REAL(ss)[(long int)VECTOR(*eidx)[j]-1]=j+1;
      }
    }
    PROTECT(newea=EVAL(lang3(install("["), oldea, ss)));
    SET_VECTOR_ELT(eal, i, newea);
    UNPROTECT(2);
  }
}


int R_igraph_attribute_add_edges(igraph_t *graph, 
				 const igraph_vector_t *edges,
				 igraph_vector_ptr_t *nattr) {
  SEXP attr=graph->attr;
  SEXP eal, rep=0, names, newnames;
  igraph_vector_t news;
  long int ealno, i, origlen, nattrno, newattrs;  
  long int ne=igraph_vector_size(edges)/2;
  if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
    SEXP newattr;
    PROTECT(newattr=duplicate(attr));
    REAL(VECTOR_ELT(attr, 0))[1] -= 1;
    if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
      UNPROTECT_PTR(attr);
    }
    REAL(VECTOR_ELT(newattr, 0))[0] = 0;
    REAL(VECTOR_ELT(newattr, 0))[1] = 1;
    attr=graph->attr=newattr;
  }

  eal=VECTOR_ELT(attr, 3);
  ealno=GET_LENGTH(eal);
  names=GET_NAMES(eal);
  if (nattr==NULL) {
    nattrno=0; 
  } else {
    nattrno=igraph_vector_ptr_size(nattr);
  }
  origlen=igraph_ecount(graph)-ne;

  /* First add the new attributes, if any */
  newattrs=0;
  IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
  for (i=0; i<nattrno; i++) {
    igraph_i_attribute_record_t *nattr_entry=VECTOR(*nattr)[i];
    const char *nname=nattr_entry->name;
    long int j;
    igraph_bool_t l=0;
    for (j=0; !l && j<ealno; j++) {
      l=!strcmp(nname, CHAR(STRING_ELT(names, j)));
    }
    if (!l) {
      newattrs++;
      IGRAPH_CHECK(igraph_vector_push_back(&news, i));
    }
  }
  if (newattrs != 0) {
    SEXP app, neweal;
    PROTECT(app=NEW_LIST(newattrs));
    PROTECT(newnames=NEW_CHARACTER(newattrs));
    PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL),
			   ScalarInteger(origlen))));
    for (i=0; i<newattrs; i++) {
      igraph_i_attribute_record_t *tmp=
	VECTOR(*nattr)[ (long int) VECTOR(news)[i]];
      SET_VECTOR_ELT(app, i, rep);
      SET_STRING_ELT(newnames, i, CREATE_STRING_VECTOR(tmp->name));
    }
    UNPROTECT(1);		/* rep */
    PROTECT(neweal=EVAL(lang3(install("c"), eal, app)));
    PROTECT(newnames=EVAL(lang3(install("c"), names, newnames)));
    SET_NAMES(neweal, newnames);
    SET_VECTOR_ELT(attr, 3, neweal);
    eal=VECTOR_ELT(attr, 3);
    ealno=GET_LENGTH(eal);
    names=GET_NAMES(eal);
    UNPROTECT(4);
    rep=0;
  }
  igraph_vector_destroy(&news);
  IGRAPH_FINALLY_CLEAN(1);

  /* Now append the new values */
  for (i=0; i<ealno; i++) {
    SEXP oldea=VECTOR_ELT(eal, i), newea;
    const char *sexpname=CHAR(STRING_ELT(names, i));
    igraph_bool_t l=0;
    long int j;
    for (j=0; !l && j<nattrno; j++) {
      igraph_i_attribute_record_t *tmp=VECTOR(*nattr)[j];
      l=!strcmp(sexpname, tmp->name);
    }
    if (l) {
      /* This attribute is present in nattr */
      SEXP app;
      igraph_i_attribute_record_t *tmprec=VECTOR(*nattr)[j-1];
      switch (tmprec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	if (ne != igraph_vector_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=NEW_NUMERIC(ne));
	igraph_vector_copy_to(tmprec->value, REAL(app));
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	if (ne != igraph_vector_ptr_size(tmprec->value)) {
	  IGRAPH_ERROR("Invalid attribute length", IGRAPH_EINVAL);
	}
	PROTECT(app=R_igraph_strvector_to_SEXP(tmprec->value));
	break;
      case IGRAPH_ATTRIBUTE_R_OBJECT:
	/* TODO */
	IGRAPH_ERROR("R objects not implemented yet", IGRAPH_UNIMPLEMENTED);
	break;
      default:
	warning("Ignoring unknown attribute type");
	break;
      }
      PROTECT(newea=EVAL(lang3(install("c"), oldea, app)));
      SET_VECTOR_ELT(eal, i, newea);
      UNPROTECT(2);		/* app & newea */
    } else {
      /* No such attribute, append NA's */
      if (rep==0) {
	PROTECT(rep=EVAL(lang3(install("rep"), ScalarLogical(NA_LOGICAL),
			       ScalarInteger(ne))));
      }
      PROTECT(newea=EVAL(lang3(install("c"), oldea, rep)));
      SET_VECTOR_ELT(eal, i, newea);
      UNPROTECT(1);		/* newea */
    }
  }
  if (rep != 0) {
    UNPROTECT(1);
  }

  return 0;
}

void R_igraph_attribute_delete_edges(igraph_t *graph, 
				     const igraph_vector_t *idx) {
  SEXP attr=graph->attr;
  SEXP eal;
  long int ealno, i;
  if (REAL(VECTOR_ELT(attr, 0))[0]+REAL(VECTOR_ELT(attr, 0))[1] > 1) {
    SEXP newattr;
    PROTECT(newattr=duplicate(attr));
    REAL(VECTOR_ELT(attr, 0))[1] -= 1;
    if (REAL(VECTOR_ELT(attr, 0))[1] == 0) {
      UNPROTECT_PTR(attr);
    }
    REAL(VECTOR_ELT(newattr, 0))[0] = 0;
    REAL(VECTOR_ELT(newattr, 0))[1] = 1;
    attr=graph->attr=newattr;
  }

  eal=VECTOR_ELT(attr, 3);
  ealno=GET_LENGTH(eal);
  for (i=0; i<ealno; i++) {
    SEXP oldea=VECTOR_ELT(eal, i), newea, ss;
    long int origlen=GET_LENGTH(oldea);
    long int newlen=0, j, k;
    /* create subscript vector */
    for (j=0; j<origlen; j++) {
      if (VECTOR(*idx)[j] > 0) {
	newlen++;
      }
    }
    PROTECT(ss=NEW_NUMERIC(newlen));
    for (j=0; j<origlen; j++) {
      if (VECTOR(*idx)[j] > 0) {
	REAL(ss)[(long int)VECTOR(*idx)[j]-1] = j+1;
      }
    }
    PROTECT(newea=EVAL(lang3(install("["), oldea, ss)));
    SET_VECTOR_ELT(eal, i, newea);
    UNPROTECT(2);
  }
}

int R_igraph_attribute_get_info(const igraph_t *graph,
				igraph_strvector_t *gnames,
				igraph_vector_t *gtypes,
				igraph_strvector_t *vnames,
				igraph_vector_t *vtypes,
				igraph_strvector_t *enames,
				igraph_vector_t *etypes) {
  igraph_strvector_t *names[3] = { gnames, vnames, enames };
  igraph_vector_t *types[3] = { gtypes, vtypes, etypes };
  long int i, j;

  SEXP attr=graph->attr;

  for (i=0; i<3; i++) {
    igraph_strvector_t *n=names[i];
    igraph_vector_t *t=types[i];
    SEXP al=VECTOR_ELT(attr, i+1);
    
    if (n) {			/* return names */
      SEXP names=GET_NAMES(al);
      R_igraph_SEXP_to_strvector_copy(names, n);
    }

    if (t) {			/* return types */
      for (j=0; j<GET_LENGTH(al); j++) {
	SEXP a=VECTOR_ELT(al, j);
	if (IS_NUMERIC(a)) {
	  VECTOR(*t)[j]=IGRAPH_ATTRIBUTE_NUMERIC;
	} else if (IS_CHARACTER(a)) {
	  VECTOR(*t)[j]=IGRAPH_ATTRIBUTE_STRING;
	} else {
	  VECTOR(*t)[j]=IGRAPH_ATTRIBUTE_R_OBJECT;
	}
      }
    }
  }

  return 0;
}

igraph_bool_t R_igraph_attribute_has_attr(const igraph_t *graph,
					  igraph_attribute_elemtype_t type,
					  const char *name) {
  long int attrnum;
  SEXP res;
  
  switch (type) {
  case IGRAPH_ATTRIBUTE_GRAPH:
    attrnum=1;
    break;
  case IGRAPH_ATTRIBUTE_VERTEX:
    attrnum=2;
    break;
  case IGRAPH_ATTRIBUTE_EDGE:
    attrnum=3;
    break;
  default:
    IGRAPH_ERROR("Unkwown attribute element type", IGRAPH_EINVAL);
    break;
  }
  
  res=R_igraph_getListElement(VECTOR_ELT(graph->attr, attrnum), name);
  return res != NULL;
}

int R_igraph_attribute_get_numeric_graph_attr(const igraph_t *graph,
					      const char *name, 
					      igraph_vector_t *value) {
  SEXP gal=VECTOR_ELT(graph->attr, 1);
  SEXP ga=R_igraph_getListElement(gal, name);
  
  if (ga == NULL) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_vector_resize(value, 1));
  VECTOR(*value)[0]=REAL(ga)[0];

  return 0;
}

int R_igraph_attribute_get_string_graph_attr(const igraph_t *graph,
					     const char *name,
					     igraph_strvector_t *value) {
  /* TODO: serialization */
  SEXP gal=VECTOR_ELT(graph->attr, 1);
  SEXP ga=R_igraph_getListElement(gal, name);
  
  if (ga == NULL) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_strvector_resize(value, 1));
  IGRAPH_CHECK(igraph_strvector_set(value, 0, CHAR(STRING_ELT(ga, 0))));
  
  return 0;
}

int R_igraph_attribute_get_numeric_vertex_attr(const igraph_t *graph, 
					       const char *name,
					       igraph_vs_t vs,
					       igraph_vector_t *value) {
  /* TODO: handle vs, serialization */
  SEXP val=VECTOR_ELT(graph->attr, 1);
  SEXP va=R_igraph_getListElement(val, name);
  igraph_vector_t newvalue;
  
  if (va == NULL) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  
  R_SEXP_to_vector_copy(va, &newvalue);
  igraph_vector_destroy(value);
  *value=newvalue;

  return 0;
}

int R_igraph_attribute_get_string_vertex_attr(const igraph_t *graph, 
					      const char *name,
					      igraph_vs_t vs,
					      igraph_strvector_t *value) {
  /* TODO: handle vs, serialization */
  SEXP val=VECTOR_ELT(graph->attr, 1);
  SEXP va=R_igraph_getListElement(val, name);
  igraph_vector_t newvalue;
  
  if (va == NULL) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  
  R_igraph_SEXP_to_strvector_copy(va, value);
  
  return 0;
}

int R_igraph_attribute_get_numeric_edge_attr(const igraph_t *graph,
					     const char *name,
					     igraph_es_t es,
					     igraph_vector_t *value) {
  /* TODO: handle es, serialization */
  SEXP eal=VECTOR_ELT(graph->attr, 1);
  SEXP ea=R_igraph_getListElement(eal, name);
  igraph_vector_t newvalue;
  
  if (ea == NULL) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  
  R_SEXP_to_vector_copy(ea, &newvalue);
  igraph_vector_destroy(value);
  *value=newvalue;

  return 0;
}

int R_igraph_attribute_get_string_edge_attr(const igraph_t *graph,
					    const char *name,
					    igraph_es_t es,
					    igraph_strvector_t *value) {
  /* TODO: handle es, serialization */
  SEXP eal=VECTOR_ELT(graph->attr, 1);
  SEXP ea=R_igraph_getListElement(eal, name);
  igraph_vector_t newvalue;
  
  if (ea == NULL) {
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  }
  
  R_igraph_SEXP_to_strvector_copy(ea, value);
  
  return 0;
}

igraph_attribute_table_t R_igraph_attribute_table={
  &R_igraph_attribute_init, &R_igraph_attribute_destroy,
  &R_igraph_attribute_copy, &R_igraph_attribute_add_vertices,
  &R_igraph_attribute_delete_vertices, &R_igraph_attribute_add_edges,
  &R_igraph_attribute_delete_edges, &R_igraph_attribute_get_info,
  &R_igraph_attribute_has_attr, &R_igraph_attribute_get_numeric_graph_attr,
  &R_igraph_attribute_get_string_graph_attr,
  &R_igraph_attribute_get_numeric_vertex_attr,
  &R_igraph_attribute_get_string_vertex_attr,
  &R_igraph_attribute_get_numeric_edge_attr,
  &R_igraph_attribute_get_string_edge_attr
};

igraph_attribute_table_t *R_igraph_attribute_oldtable;

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

igraph_interruption_handler_t *R_igraph_oldinterrupt;

extern int R_interrupts_pending;

int R_igraph_interrupt_handler(void *data) {
#if  ( defined(HAVE_AQUA) || defined(Win32) )
  /* TODO!!!: proper interrupt handling for these systems, 
     right now memory is not deallocated!!! */
  R_CheckUserInterrupt();
#else
  if (R_interrupts_pending) {
    IGRAPH_FINALLY_FREE();
    R_CheckUserInterrupt();
  }
#endif
  return 0;
}

R_INLINE void R_igraph_before() {
  R_igraph_oldhandler=igraph_set_error_handler(R_igraph_myhandler);
  R_igraph_oldinterrupt=
    igraph_set_interruption_handler(R_igraph_interrupt_handler);
  R_igraph_attribute_oldtable=
    igraph_i_set_attribute_table(&R_igraph_attribute_table);
}

R_INLINE void R_igraph_after() {
  igraph_set_error_handler(R_igraph_oldhandler);
  igraph_set_interruption_handler(R_igraph_oldinterrupt);
  igraph_i_set_attribute_table(R_igraph_attribute_oldtable);
}

int R_igraph_progress_handler(const char *message, igraph_real_t percent,
			      void* data) {
  static igraph_real_t last=0;
  if (percent == 0) {
    last=0;
  } 
  while (last < percent) {
    fputc('#', stderr);
    last+=2;
  }
  return 0;
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

SEXP R_igraph_to_SEXP(igraph_t *graph) {
  
  SEXP result;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  SEXP attr;
  
  PROTECT(result=NEW_LIST(9));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(1));
  SET_VECTOR_ELT(result, 1, NEW_LOGICAL(1));
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 3, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 4, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 5, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 6, NEW_NUMERIC(no_of_nodes+1));
  SET_VECTOR_ELT(result, 7, NEW_NUMERIC(no_of_nodes+1));

  REAL(VECTOR_ELT(result, 0))[0]=no_of_nodes;
  LOGICAL(VECTOR_ELT(result, 1))[0]=graph->directed;
  memcpy(REAL(VECTOR_ELT(result, 2)), graph->from.stor_begin, 
	 sizeof(igraph_real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 3)), graph->to.stor_begin, 
	 sizeof(igraph_real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 4)), graph->oi.stor_begin, 
	 sizeof(igraph_real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 5)), graph->ii.stor_begin, 
	 sizeof(igraph_real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 6)), graph->os.stor_begin, 
	 sizeof(igraph_real_t)*(no_of_nodes+1));
  memcpy(REAL(VECTOR_ELT(result, 7)), graph->is.stor_begin, 
	 sizeof(igraph_real_t)*(no_of_nodes+1));
  
  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("igraph")));

  /* Attributes */
  SET_VECTOR_ELT(result, 8, graph->attr);
  REAL(VECTOR_ELT(graph->attr, 0))[0] += 1;
  
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

int R_SEXP_to_igraph(SEXP graph, igraph_t *res) {
  
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  R_SEXP_to_vector(VECTOR_ELT(graph, 2), &res->from);
  R_SEXP_to_vector(VECTOR_ELT(graph, 3), &res->to);
  R_SEXP_to_vector(VECTOR_ELT(graph, 4), &res->oi);
  R_SEXP_to_vector(VECTOR_ELT(graph, 5), &res->ii);
  R_SEXP_to_vector(VECTOR_ELT(graph, 6), &res->os);
  R_SEXP_to_vector(VECTOR_ELT(graph, 7), &res->is);
  
  /* attributes */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[0] = 1; /* R objects refcount */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[1] = 0; /* igraph_t objects */
  res->attr=VECTOR_ELT(graph, 8);
  
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

  /* attributes */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[0] = 1; /* R objects */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[1] = 1; /* igraph_t objects */
  PROTECT(res->attr=VECTOR_ELT(graph, 8));  

  return 0;
}

/* 
 * We have only seq and vector types
 */

int R_SEXP_to_igraph_vs_copy(SEXP rit, igraph_t *graph, igraph_vs_t *it) {

  if (!strcmp(CHAR(STRING_ELT(GET_CLASS(rit),0)), "igraphvsseq")) { 
    
    long int from=REAL(VECTOR_ELT(rit, 0))[1]; 
    long int to=REAL(VECTOR_ELT(rit, 0))[2];
    igraph_vs_seq(it, from, to); 
  } else if (!strcmp(CHAR(STRING_ELT(GET_CLASS(rit),0)), "igraphvsvector")) { 
    igraph_vector_t *tmpv=(igraph_vector_t*)R_alloc(1,sizeof(igraph_vector_t));
    igraph_vs_vector(it, igraph_vector_view(tmpv, REAL(VECTOR_ELT(rit,1)), 
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
    igraph_es_seq(it, from, to);
  } else if (!strcmp(CHAR(STRING_ELT(GET_CLASS(rit),0)), "igraphesvector")) {
    igraph_vector_t *tmpv=(igraph_vector_t*)R_alloc(1,sizeof(igraph_vector_t));
    igraph_es_vector(it, igraph_vector_view(tmpv, REAL(VECTOR_ELT(rit,1)),
					    GET_LENGTH(VECTOR_ELT(rit,1))));
  } else {
    error("unknown edge set type");
  }
  return 0;
}

/*******************************************************************/

SEXP R_igraph_empty(SEXP pn, SEXP pdirected) {
  
  SEXP result;
  igraph_integer_t n=REAL(pn)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
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
  igraph_add_edges(&g, &v, 0);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_vertices(SEXP graph, SEXP pnv) {
  
  igraph_integer_t nv;
  igraph_t g;
  SEXP result;
  
  R_igraph_before();
  
  nv=REAL(pnv)[0];
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_vertices(&g, nv, 0);
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
  igraph_real_t vid;
  igraph_integer_t mode;
  
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
  igraph_delete_vertices(&g, vs);
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
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
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_bool_t loops=LOGICAL(ploops)[0];
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(vids, &g, &vs);
  igraph_vector_init(&res, 0);
  igraph_degree(&g, &res, vs, mode, loops);
  
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
  igraph_integer_t mode=REAL(pmode)[0];
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
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_bool_t res;
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
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t unconnected=LOGICAL(punconnected)[0];
  igraph_integer_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_diameter(&g, &res, 0, 0, 0, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_diameter(SEXP graph, SEXP pdirected, SEXP punconnected) {
  
  igraph_t g;
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t unconnected=LOGICAL(punconnected)[0];
  igraph_vector_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_vector_init(&res, 0);
  igraph_diameter(&g, 0, 0, 0, &res, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_farthest_points(SEXP graph, SEXP pdirected, SEXP punconnected) {
  
  igraph_t g;
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t unconnected=LOGICAL(punconnected)[0];
  igraph_integer_t from, to, len;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_diameter(&g, &len, &from, &to, 0, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(3));
  if (from < 0) {
    REAL(result)[0]=REAL(result)[1]=REAL(result)[2]=NA_REAL;
  } else {
    REAL(result)[0]=from;
    REAL(result)[1]=to;
    REAL(result)[2]=len;
  }
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_closeness(SEXP graph, SEXP pvids, SEXP pmode) {
  
  igraph_t g;
  igraph_vs_t vs;
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_vector_t res;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_vector_init(&res, 0);
  igraph_closeness(&g, &res, vs, mode);
  
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
  igraph_real_t vertex=REAL(pvertex)[0];
  igraph_integer_t mode=REAL(pmode)[0];
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
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_vector_t res;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_vector_init(&res, 0);
  igraph_betweenness(&g, &res, vs, directed);
  
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
  igraph_integer_t binwidth=REAL(pbinwidth)[0];
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
  igraph_cocitation(&g, &m, vs);
  
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
  igraph_bibcoupling(&g, &m, vs);
  
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_integer_t m=REAL(pm)[0];
  igraph_bool_t citation=LOGICAL(pcitation)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
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
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_matrix_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_matrix_init(&res, 0, 0);
  igraph_shortest_paths(&g, &res, vs, mode);
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
  igraph_integer_t nei=REAL(pnei)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t mutual=LOGICAL(pmutual)[0];
  igraph_bool_t circular=LOGICAL(pcircular)[0];  
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_integer_t m=REAL(pm)[0]; 
  igraph_vector_t outseq;
  igraph_bool_t outpref=LOGICAL(poutpref)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
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
				  SEXP pcoolexp, SEXP pkkconst, SEXP psigma, 
				  SEXP pverbose) {
  
  igraph_t g;
  igraph_integer_t niter=REAL(pniter)[0];
  igraph_real_t initemp=REAL(pinitemp)[0];
  igraph_real_t coolexp=REAL(pcoolexp)[0];
  igraph_real_t kkconst=REAL(pkkconst)[0];
  igraph_real_t sigma=REAL(psigma)[0];
  igraph_matrix_t res;
  igraph_bool_t verbose=LOGICAL(pverbose)[0];
  igraph_progress_handler_t *oldprogress;
  SEXP result;
  
  R_igraph_before();
  if (verbose) {
    fprintf(stderr, "KK layout: ");
    oldprogress=igraph_set_progress_handler(R_igraph_progress_handler);
  }
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_kamada_kawai(&g, &res, niter, sigma, 
			     initemp, coolexp, kkconst);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  if (verbose) {
    igraph_set_progress_handler(oldprogress);
    fputc('\n', stderr);
  }
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_kamada_kawai_3d(SEXP graph, SEXP pniter, SEXP pinitemp,
				     SEXP pcoolexp, SEXP pkkconst, 
				     SEXP psigma, SEXP pverbose) {
  igraph_t g;
  igraph_integer_t niter=REAL(pniter)[0];
  igraph_real_t initemp=REAL(pinitemp)[0];
  igraph_real_t coolexp=REAL(pcoolexp)[0];
  igraph_real_t kkconst=REAL(pkkconst)[0];
  igraph_real_t sigma=REAL(psigma)[0];
  igraph_matrix_t res;
  igraph_bool_t verbose=LOGICAL(pverbose)[0];
  igraph_progress_handler_t *oldprogress;
  SEXP result;
  
  R_igraph_before();
  if (verbose) {
    fprintf(stderr, "3D KK layout: ");
    oldprogress=igraph_set_progress_handler(R_igraph_progress_handler);
  }
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_kamada_kawai_3d(&g, &res, niter, sigma, 
				initemp, coolexp, kkconst);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  if (verbose) {
    igraph_set_progress_handler(oldprogress);
    fputc('\n', stderr);
  }
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}


SEXP R_igraph_layout_lgl(SEXP graph, SEXP pmaxiter, SEXP pmaxdelta,
			 SEXP parea, SEXP pcoolexp, SEXP prepulserad,
			 SEXP pcellsize, SEXP proot) {
  
  igraph_t g;
  igraph_matrix_t res;
  igraph_integer_t maxiter=REAL(pmaxiter)[0];
  igraph_real_t maxdelta=REAL(pmaxdelta)[0];
  igraph_real_t area=REAL(parea)[0];
  igraph_real_t coolexp=REAL(pcoolexp)[0];
  igraph_real_t repulserad=REAL(prepulserad)[0];
  igraph_real_t cellsize=REAL(pcellsize)[0];
  igraph_integer_t root=REAL(proot)[0];
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
					       SEXP pcellsize, SEXP puseseed,
					       SEXP pverbose) {

  igraph_t g;
  igraph_matrix_t res;
  igraph_integer_t niter=REAL(pniter)[0];
  igraph_real_t maxdelta=REAL(pmaxdelta)[0];
  igraph_real_t area=REAL(parea)[0];
  igraph_real_t coolexp=REAL(pcoolexp)[0];
  igraph_real_t repulserad=REAL(prepulserad)[0];
  igraph_real_t cellsize=REAL(pcellsize)[0];
  igraph_bool_t use_seed=LOGICAL(puseseed)[0];
  igraph_bool_t verbose=LOGICAL(pverbose)[0];
  igraph_progress_handler_t *oldprogress;
  SEXP result;
  
  R_igraph_before();
  if (verbose) {
    fprintf(stderr, "Grid-FR layout: ");
    oldprogress=igraph_set_progress_handler(R_igraph_progress_handler);
  }

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
  
  if (verbose) {
    igraph_set_progress_handler(oldprogress);
    fputc('\n', stderr);
  }
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
  igraph_bool_t directed=LOGICAL(pdirected)[0];
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

SEXP R_igraph_measure_dynamics_id(SEXP graph, SEXP pstartvertex,
				  SEXP pst, SEXP pmaxind, 
				  SEXP psign, SEXP pno) {
  
  igraph_t g;
  igraph_matrix_t ak, sd, confint, no;
  igraph_vector_t st;
  igraph_integer_t startvertex=REAL(pstartvertex)[0];
  igraph_integer_t maxind=REAL(pmaxind)[0];
  igraph_real_t sign=REAL(psign)[0];
  igraph_bool_t lno=LOGICAL(pno)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_vector(pst, &st);
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&ak, 0, 0);
  igraph_matrix_init(&sd, 0, 0);
  igraph_matrix_init(&confint, 0, 0);
  igraph_matrix_init(&no, 0, 0);
  igraph_measure_dynamics_id(&g, startvertex, &ak, &sd, &confint, &no, &st,
			     maxind, sign, lno);
  PROTECT(result=NEW_LIST(4));
  SET_VECTOR_ELT(result, 0, R_igraph_matrix_to_SEXP(&ak));
  igraph_matrix_destroy(&ak);
  SET_VECTOR_ELT(result, 1, R_igraph_matrix_to_SEXP(&sd));
  igraph_matrix_destroy(&sd);
  SET_VECTOR_ELT(result, 2, R_igraph_matrix_to_SEXP(&confint));
  igraph_matrix_destroy(&confint);
  SET_VECTOR_ELT(result, 3, R_igraph_matrix_to_SEXP(&no));
  igraph_matrix_destroy(&no);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_id_st(SEXP graph, SEXP pak) {
  
  igraph_t g;
  igraph_matrix_t ak;
  igraph_vector_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_matrix(pak, &ak);
  igraph_vector_init(&res, 0);
  igraph_measure_dynamics_id_st(&g, &res, &ak);
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  
  igraph_vector_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage(SEXP graph, SEXP pstartvertex,
				     SEXP pst, SEXP pagebins,
				     SEXP pmaxind, SEXP psign, SEXP pno) {

  igraph_t g;
  igraph_matrix_t akl, sd, confint, no;
  igraph_vector_t st;
  igraph_integer_t startvertex=REAL(pstartvertex)[0];
  igraph_integer_t agebins=REAL(pagebins)[0];
  igraph_integer_t maxind=REAL(pmaxind)[0];
  igraph_real_t sign=REAL(psign)[0];
  igraph_bool_t lno=LOGICAL(pno)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_vector(pst, &st);

  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&akl, 0, 0);
  igraph_matrix_init(&sd, 0, 0);
  igraph_matrix_init(&confint, 0, 0);
  igraph_matrix_init(&no, 0, 0);
  igraph_measure_dynamics_idage(&g, startvertex, &akl, &sd, &confint, &no, &st,
				agebins, maxind, sign, lno);
    
  PROTECT(result=NEW_LIST(4));  
  SET_VECTOR_ELT(result, 0, R_igraph_matrix_to_SEXP(&akl));
  igraph_matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_igraph_matrix_to_SEXP(&sd));
  igraph_matrix_destroy(&sd);
  SET_VECTOR_ELT(result, 2, R_igraph_matrix_to_SEXP(&confint));
  igraph_matrix_destroy(&confint);
  SET_VECTOR_ELT(result, 3, R_igraph_matrix_to_SEXP(&no));
  igraph_matrix_destroy(&no);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage_debug(SEXP graph, SEXP pst, 
					   SEXP pagebins, SEXP pmaxind, 
					   SEXP psign, SEXP pest_ind, 
					   SEXP pest_age, SEXP pno) {
  igraph_t g;
  igraph_matrix_t akl, sd, confint, no;
  igraph_vector_t st;
  igraph_integer_t agebins=REAL(pagebins)[0];
  igraph_integer_t maxind=REAL(pmaxind)[0];
  igraph_real_t sign=REAL(psign)[0];
  igraph_bool_t lno=LOGICAL(pno)[0];
  igraph_integer_t est_ind=REAL(pest_ind)[0];
  igraph_integer_t est_age=REAL(pest_age)[0];
  igraph_vector_t estimates;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_vector(pst, &st);
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&akl, 0, 0);
  igraph_matrix_init(&sd, 0, 0);
  igraph_matrix_init(&confint, 0, 0);
  igraph_matrix_init(&no, 0, 0);
  igraph_vector_init(&estimates, 0);
  igraph_measure_dynamics_idage_debug(&g, &akl, &sd, &no, &confint, 
				      &st, agebins, 
				      maxind, sign, &estimates, est_ind, 
				      est_age, lno);
    
  PROTECT(result=NEW_LIST(5));
  SET_VECTOR_ELT(result, 0, R_igraph_matrix_to_SEXP(&akl));
  igraph_matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_igraph_matrix_to_SEXP(&sd));
  igraph_matrix_destroy(&sd);
  SET_VECTOR_ELT(result, 2, R_igraph_matrix_to_SEXP(&confint));
  igraph_matrix_destroy(&confint);
  SET_VECTOR_ELT(result, 3, R_igraph_matrix_to_SEXP(&no));
  igraph_matrix_destroy(&no);
  SET_VECTOR_ELT(result, 4, NEW_NUMERIC(igraph_vector_size(&estimates)));
  igraph_vector_copy_to(&estimates, REAL(VECTOR_ELT(result, 4)));
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

SEXP R_igraph_get_shortest_paths(SEXP graph, SEXP pfrom, SEXP pto, 
				 SEXP pmode, SEXP pno) {

  igraph_t g;
  igraph_integer_t from=REAL(pfrom)[0];
  igraph_vs_t to;
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_vector_t *vects;
  long int i;
  igraph_vector_ptr_t ptrvec;
  SEXP result;
  
  long int no=REAL(pno)[0];
  
  R_igraph_before();  

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pto, &g, &to);

  igraph_vector_ptr_init(&ptrvec, no);
  vects=(igraph_vector_t*) R_alloc(GET_LENGTH(pto), sizeof(igraph_vector_t));
  for (i=0; i<no; i++) {
    igraph_vector_init(&vects[i], 0);
    VECTOR(ptrvec)[i]=&vects[i];
  }
  igraph_get_shortest_paths(&g, &ptrvec, from, to, mode);
  PROTECT(result=NEW_LIST(no));
  for (i=0; i<no; i++) {
    SET_VECTOR_ELT(result, i, NEW_NUMERIC(igraph_vector_size(&vects[i])));
    igraph_vector_copy_to(&vects[i], REAL(VECTOR_ELT(result, i)));
    igraph_vector_destroy(&vects[i]);
  }
  igraph_vector_ptr_destroy(&ptrvec);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_are_connected(SEXP graph, SEXP pv1, SEXP pv2) {
  
  igraph_t g;
  igraph_integer_t v1=REAL(pv1)[0];
  igraph_integer_t v2=REAL(pv2)[0];
  igraph_bool_t res;
  SEXP result;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_LOGICAL(1));
  igraph_are_connected(&g, v1, v2, &res);
  LOGICAL(result)[0]=res;
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_graph_adjacency(SEXP adjmatrix, SEXP pmode) {
  
  igraph_t g;
  igraph_matrix_t adjm;
  igraph_integer_t mode=REAL(pmode)[0];
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
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t unconnected=LOGICAL(punconnected)[0];
  igraph_real_t res;
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_integer_t center=REAL(pcenter)[0];
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t mutual=LOGICAL(pmutual)[0];
  igraph_bool_t circular=LOGICAL(pcircular)[0];
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_integer_t children=REAL(pchildren)[0];
  igraph_integer_t mode=REAL(pmode)[0];
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
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pvids, &g, &vs);
  igraph_subgraph(&g, &sub, vs);
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_integer_t type=REAL(ptype)[0];
  igraph_real_t porm=REAL(pporm)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t loops=LOGICAL(ploops)[0];
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
  igraph_integer_t n=REAL(pn)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  igraph_bool_t loops=LOGICAL(ploops)[0];
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
  igraph_integer_t low=REAL(plow)[0];
  igraph_integer_t high=REAL(phigh)[0];
  igraph_integer_t length=REAL(plength)[0];
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
  igraph_bool_t bycol=LOGICAL(pbycol)[0];
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
  igraph_integer_t type=REAL(ptype)[0];
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
  igraph_bool_t multiple=LOGICAL(pmultiple)[0];
  igraph_bool_t loops=LOGICAL(ploops)[0];
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
					  SEXP pcoolexp, SEXP prepulserad, 
					  SEXP pverbose) {
  igraph_t g;
  igraph_integer_t niter=REAL(pniter)[0];
  igraph_real_t maxdelta=REAL(pmaxdelta)[0];
  igraph_real_t area=REAL(parea)[0];
  igraph_real_t coolexp=REAL(pcoolexp)[0];
  igraph_real_t repulserad=REAL(prepulserad)[0];
  igraph_matrix_t res;
  igraph_bool_t verbose=LOGICAL(pverbose)[0];
  igraph_progress_handler_t *oldprogress;
  SEXP result;
  
  R_igraph_before();
  if (verbose) {
    fprintf(stderr, "FR layout: ");
    oldprogress=igraph_set_progress_handler(R_igraph_progress_handler);
  }
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_fruchterman_reingold(&g, &res, niter, maxdelta, area, 
				     coolexp, repulserad, 0);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  if (verbose) {
    igraph_set_progress_handler(oldprogress);
    fputc('\n', stderr);
  }
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_fruchterman_reingold_3d(SEXP graph, SEXP pniter, 
					     SEXP pmaxdelta, SEXP parea,
					     SEXP pcoolexp, SEXP prepulserad,
					     SEXP pverbose) {
  igraph_t g;
  igraph_integer_t niter=REAL(pniter)[0];
  igraph_real_t maxdelta=REAL(pmaxdelta)[0];
  igraph_real_t area=REAL(parea)[0];
  igraph_real_t coolexp=REAL(pcoolexp)[0];
  igraph_real_t repulserad=REAL(prepulserad)[0];
  igraph_matrix_t res;
  igraph_bool_t verbose=LOGICAL(pverbose)[0];
  igraph_progress_handler_t *oldprogress;
  SEXP result;
  
  R_igraph_before();
  if (verbose) {
    fprintf(stderr, "3D FR layout: ");
    oldprogress=igraph_set_progress_handler(R_igraph_progress_handler);
  }
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_fruchterman_reingold_3d(&g, &res, niter, maxdelta, area, 
				     coolexp, repulserad, 0);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  if (verbose) {
    igraph_set_progress_handler(oldprogress);
    fputc('\n', stderr);
  }
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_degree_sequence_game(SEXP pout_seq, SEXP pin_seq, 
				   SEXP pmethod) {
  igraph_t g;
  igraph_vector_t outseq;
  igraph_vector_t inseq;
  igraph_integer_t method=REAL(pmethod)[0];
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
  igraph_integer_t type=REAL(ptype)[0];
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
  
SEXP R_igraph_read_graph_edgelist(SEXP pvfile, SEXP pn, SEXP pdirected) {
  igraph_t g;
  igraph_integer_t n=REAL(pn)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
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

SEXP R_igraph_read_graph_ncol(SEXP pvfile, SEXP ppredef,
			      SEXP pnames, SEXP pweights,
			      SEXP pdirected) {
  igraph_t g;
  igraph_bool_t names=LOGICAL(pnames)[0];
  igraph_bool_t weights=LOGICAL(pweights)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  FILE *file;
  igraph_strvector_t predef, *predefptr=0;  
  SEXP result;
  
  R_igraph_before();
  
#if HAVE_FMEMOPEN == 1
  file=fmemopen(RAW(pvfile), GET_LENGTH(pvfile), "r");
#else 
  file=fopen(CHAR(STRING_ELT(pvfile, 0)), "r");
#endif
  if (file==0) { igraph_error("Cannot read edgelist", __FILE__, __LINE__,
			      IGRAPH_EFILE); }
  if (GET_LENGTH(ppredef)>0) {
    R_igraph_SEXP_to_strvector(ppredef, &predef);
    predefptr=&predef;
  } 
  igraph_read_graph_ncol(&g, file, predefptr, names, weights, directed);
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

  R_SEXP_to_igraph(graph, &g);
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
  igraph_bool_t names=LOGICAL(pnames)[0];
  igraph_bool_t weights=LOGICAL(pweights)[0];
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
  igraph_bool_t isolates=LOGICAL(pisolates)[0];
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

  R_SEXP_to_igraph(graph, &g);
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

SEXP R_igraph_decompose(SEXP graph, SEXP pmode, SEXP pmaxcompno, 
			SEXP pminelements) {

  igraph_t g;
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_integer_t maxcompno=REAL(pmaxcompno)[0];
  igraph_integer_t minelements=REAL(pminelements)[0];
  igraph_vector_ptr_t comps;
  SEXP result;
  long int i;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
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

SEXP R_igraph_layout_random_3d(SEXP graph) {
  
  igraph_t g;
  igraph_matrix_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_random_3d(&g, &res);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);  
  
  R_igraph_after();

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_sphere(SEXP graph) {

  igraph_t g;
  igraph_matrix_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_sphere(&g, &res);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_isoclass_34(SEXP graph) {
  
  igraph_t g;
  int class;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  igraph_isoclass(&g, &class);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=class;
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_isomorphic_34(SEXP graph1, SEXP graph2) {
  
  igraph_t g1, g2;
  igraph_bool_t res;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_igraph(graph1, &g1);
  R_SEXP_to_igraph(graph2, &g2);
  igraph_isomorphic(&g1, &g2, &res);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=(res>0);
    
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_callaway_traits_game(SEXP pnodes, SEXP ptypes, 
				  SEXP pepers, SEXP ptype_dist,
				  SEXP pmatrix, SEXP pdirected) {

  igraph_t g;
  igraph_integer_t nodes=REAL(pnodes)[0];
  igraph_integer_t types=REAL(ptypes)[0];
  igraph_integer_t epers=REAL(pepers)[0];
  igraph_vector_t type_dist;
  igraph_matrix_t matrix;
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  SEXP result; 

  R_igraph_before();
  
  R_SEXP_to_vector(ptype_dist, &type_dist);
  R_SEXP_to_matrix(pmatrix, &matrix);
  igraph_callaway_traits_game(&g, nodes, types, epers, &type_dist,
			      &matrix, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_establishment_game(SEXP pnodes, SEXP ptypes, SEXP pk,
				 SEXP ptype_dist, SEXP pmatrix,
				 SEXP pdirected) {
  igraph_t g;
  igraph_integer_t nodes=REAL(pnodes)[0];
  igraph_integer_t types=REAL(ptypes)[0];
  igraph_integer_t k=REAL(pk)[0];
  igraph_vector_t type_dist;
  igraph_matrix_t matrix;
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_vector(ptype_dist, &type_dist);
  R_SEXP_to_matrix(pmatrix, &matrix);
  igraph_establishment_game(&g, nodes, types, k, &type_dist, &matrix,
			    directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
			    
SEXP R_igraph_motifs_randesu(SEXP graph, SEXP psize, SEXP pcutprob) {
  igraph_t g;
  igraph_integer_t size=REAL(psize)[0];
  igraph_vector_t cutprob;
  igraph_vector_t res;
  SEXP result;

  R_igraph_before();

  R_SEXP_to_vector(pcutprob, &cutprob);
  R_SEXP_to_igraph(graph, &g);
  igraph_vector_init(&res, 0);
  igraph_motifs_randesu(&g, &res, size, &cutprob);
  PROTECT(result=NEW_NUMERIC(igraph_vector_size(&res)));
  igraph_vector_copy_to(&res, REAL(result));
  igraph_vector_destroy(&res);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_motifs_randesu_no(SEXP graph, SEXP psize, SEXP pcutprob) {
  igraph_t g;
  igraph_integer_t size=REAL(psize)[0];
  igraph_vector_t cutprob;
  igraph_integer_t res;
  SEXP result;

  R_igraph_before();

  R_SEXP_to_vector(pcutprob, &cutprob);
  R_SEXP_to_igraph(graph, &g);
  igraph_motifs_randesu_no(&g, &res, size, &cutprob);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_motifs_randesu_estimate(SEXP graph, SEXP psize, SEXP pcutprob,
				      SEXP psamplesize, SEXP psample) {
  igraph_t g;
  igraph_integer_t size=REAL(psize)[0];
  igraph_vector_t cutprob;
  igraph_integer_t samplesize=REAL(psamplesize)[0];
  igraph_vector_t sample;
  igraph_integer_t res;
  SEXP result;
  
  R_igraph_before();

  R_SEXP_to_vector(pcutprob, &cutprob);
  R_SEXP_to_vector(psample, &sample);
  R_SEXP_to_igraph(graph, &g);
  igraph_motifs_randesu_estimate(&g, &res, size, &cutprob, samplesize,
				 &sample);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_all_shortest_paths(SEXP graph, SEXP pfrom, SEXP pto,
				     SEXP pmode) {
  
  igraph_t g;
  igraph_integer_t from=REAL(pfrom)[0];
  igraph_vs_t to;
  igraph_integer_t mode=REAL(pmode)[0];
  igraph_vector_ptr_t res;
  SEXP result;
  long int i;

  R_igraph_before();
  
  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_igraph_vs_copy(pto, &g, &to);
  igraph_vector_ptr_init(&res, 0);
  igraph_get_all_shortest_paths(&g, &res, 0, from, to, mode);
  PROTECT(result=NEW_LIST(igraph_vector_ptr_size(&res)));
  for (i=0; i<igraph_vector_ptr_size(&res); i++) {
    long int len=igraph_vector_size(VECTOR(res)[i]);
    SET_VECTOR_ELT(result, i, NEW_NUMERIC(len));
    igraph_vector_copy_to(VECTOR(res)[i], REAL(VECTOR_ELT(result, i)));
    igraph_vector_destroy(VECTOR(res)[i]);
  }
  igraph_vector_ptr_destroy(&res);
  igraph_vs_destroy(&to);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_isoclass_create(SEXP psize, SEXP pnumber, SEXP pdirected) {

  igraph_t g;
  igraph_integer_t size=REAL(psize)[0];
  igraph_integer_t number=REAL(pnumber)[0];
  igraph_bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  R_igraph_before();

  igraph_isoclass_create(&g, size, number, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_layout_merge_dla(SEXP graphs, SEXP layouts, SEXP pverbose) {

  igraph_vector_ptr_t graphvec;
  igraph_vector_ptr_t ptrvec;
  igraph_t *gras;
  igraph_matrix_t *mats;
  igraph_matrix_t res;
  long int i;
  igraph_bool_t verbose=LOGICAL(pverbose)[0];
  igraph_progress_handler_t *oldprogress;
  SEXP result;
  
  R_igraph_before();
  if (verbose) {
    fprintf(stderr, "Merge layouts: ");
    oldprogress=igraph_set_progress_handler(R_igraph_progress_handler);
  }
  
  igraph_vector_ptr_init(&graphvec, GET_LENGTH(graphs));
  igraph_vector_ptr_init(&ptrvec, GET_LENGTH(layouts));
  gras=(igraph_t*)R_alloc(GET_LENGTH(graphs), sizeof(igraph_t));
  mats=(igraph_matrix_t*)R_alloc(GET_LENGTH(layouts),
				 sizeof(igraph_matrix_t));
  for (i=0; i<GET_LENGTH(graphs); i++) {
    R_SEXP_to_igraph(VECTOR_ELT(graphs, i), &gras[i]);
    VECTOR(graphvec)[i]=&gras[i];
  }
  for (i=0; i<GET_LENGTH(layouts); i++) {
    R_SEXP_to_matrix(VECTOR_ELT(layouts, i), &mats[i]);
    VECTOR(ptrvec)[i]=&mats[i];
  }
  igraph_matrix_init(&res, 0, 0);
  igraph_layout_merge_dla(&graphvec, &ptrvec, &res);
  igraph_vector_ptr_destroy(&graphvec);
  igraph_vector_ptr_destroy(&ptrvec);
  PROTECT(result=R_igraph_matrix_to_SEXP(&res));
  igraph_matrix_destroy(&res);
  
  if (verbose) {
    igraph_set_progress_handler(oldprogress);
    fputc('\n', stderr);
  }
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_disjoint_union(SEXP pgraphs) {
  
  igraph_vector_ptr_t ptrvec;
  igraph_t *graphs;
  igraph_t res;
  long int i;
  SEXP result;

  R_igraph_before();
  
  igraph_vector_ptr_init(&ptrvec, GET_LENGTH(pgraphs));
  graphs=(igraph_t *)R_alloc(GET_LENGTH(pgraphs),
			     sizeof(igraph_t));
  for (i=0; i<GET_LENGTH(pgraphs); i++) {
    R_SEXP_to_igraph(VECTOR_ELT(pgraphs, i), &graphs[i]);
    VECTOR(ptrvec)[i]=&graphs[i];
  }
  igraph_disjoint_union_many(&res, &ptrvec);
  igraph_vector_ptr_destroy(&ptrvec);
  PROTECT(result=R_igraph_to_SEXP(&res));
  igraph_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_union(SEXP pgraphs) {
  
  igraph_vector_ptr_t ptrvec;
  igraph_t *graphs;
  igraph_t res;
  long int i;
  SEXP result;

  R_igraph_before();
  
  igraph_vector_ptr_init(&ptrvec, GET_LENGTH(pgraphs));
  graphs=(igraph_t *)R_alloc(GET_LENGTH(pgraphs),
			     sizeof(igraph_t));
  for (i=0; i<GET_LENGTH(pgraphs); i++) {
    R_SEXP_to_igraph(VECTOR_ELT(pgraphs, i), &graphs[i]);
    VECTOR(ptrvec)[i]=&graphs[i];
  }
  igraph_union_many(&res, &ptrvec);
  igraph_vector_ptr_destroy(&ptrvec);
  PROTECT(result=R_igraph_to_SEXP(&res));
  igraph_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_intersection(SEXP pgraphs) {
  
  igraph_vector_ptr_t ptrvec;
  igraph_t *graphs;
  igraph_t res;
  long int i;
  SEXP result;
  
  R_igraph_before();
  
  igraph_vector_ptr_init(&ptrvec, GET_LENGTH(pgraphs));
  graphs=(igraph_t *)R_alloc(GET_LENGTH(pgraphs),
			     sizeof(igraph_t));
  for (i=0; i<GET_LENGTH(pgraphs); i++) {
    R_SEXP_to_igraph(VECTOR_ELT(pgraphs, i), &graphs[i]);
    VECTOR(ptrvec)[i]=&graphs[i];
  }
  igraph_intersection_many(&res, &ptrvec);
  igraph_vector_ptr_destroy(&ptrvec);
  PROTECT(result=R_igraph_to_SEXP(&res));
  igraph_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_difference(SEXP pleft, SEXP pright) {
  
  igraph_t left, right;
  igraph_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(pleft, &left);
  R_SEXP_to_igraph(pright, &right);
  igraph_difference(&res, &left, &right);
  PROTECT(result=R_igraph_to_SEXP(&res));
  igraph_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_complementer(SEXP pgraph, SEXP ploops) {
  
  igraph_t g;
  igraph_t res;
  igraph_bool_t loops=LOGICAL(ploops)[0];
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(pgraph, &g);
  igraph_complementer(&res, &g, loops);
  PROTECT(result=R_igraph_to_SEXP(&res));
  igraph_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_compose(SEXP pleft, SEXP pright) {
  
  igraph_t left, right;
  igraph_t res;
  SEXP result;
  
  R_igraph_before();
  
  R_SEXP_to_igraph(pleft, &left);
  R_SEXP_to_igraph(pright, &right);
  igraph_compose(&res, &left, &right);
  PROTECT(result=R_igraph_to_SEXP(&res));
  igraph_destroy(&res);
  
  R_igraph_after();
  
  UNPROTECT(1);
  return result;
}
