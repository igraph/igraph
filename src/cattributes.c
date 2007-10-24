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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "attributes.h"
#include "memory.h"
#include "igraph.h"

#include <string.h>

/* An attribute is either a numeric vector (vector_t) or a string
   vector (strvector_t). The attribute itself is stored in a 
   struct igraph_i_attribute_record_t, there is one such object for each
   attribute. The igraph_t has a pointer to an array of three
   vector_ptr_t's which contains pointers to
   igraph_i_cattribute_t's. Graph attributes are first, then vertex 
   and edge attributes. */

igraph_bool_t igraph_i_cattribute_find(const igraph_vector_ptr_t *ptrvec, 
				       const char *name, long int *idx) {
  long int i, n=igraph_vector_ptr_size(ptrvec);
  igraph_bool_t l=0;
  for (i=0; !l && i<n; i++) {
    igraph_i_attribute_record_t *rec=VECTOR(*ptrvec)[i];
    l= !strcmp(rec->name, name);
  }
  if (idx) { *idx=i-1; }
  return l;
}

typedef struct igraph_i_cattributes_t {
  igraph_vector_ptr_t gal;
  igraph_vector_ptr_t val;
  igraph_vector_ptr_t eal;
} igraph_i_cattributes_t;

int igraph_i_cattribute_init(igraph_t *graph, igraph_vector_ptr_t *attr) {
  igraph_i_cattributes_t *nattr=igraph_Calloc(1, igraph_i_cattributes_t);
  if (!nattr) {
    IGRAPH_ERROR("Can't init attributes", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, nattr);
  IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->gal, 0));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &nattr->gal);
  IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->val, 0));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &nattr->gal);
  IGRAPH_CHECK(igraph_vector_ptr_init(&nattr->eal, 0));
  
  IGRAPH_FINALLY_CLEAN(3);
  graph->attr=nattr;
  return 0;
}

void igraph_i_cattribute_destroy(igraph_t *graph) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *als[3]= { &attr->gal, &attr->val, &attr->eal };
  long int i, n, a;
  igraph_vector_t *num;
  igraph_strvector_t *str;
  igraph_i_attribute_record_t *rec;
  for (a=0; a<3; a++) {
    n=igraph_vector_ptr_size(als[a]);
    for (i=0; i<n; i++) {
      rec=VECTOR(*als[a])[i];
      if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
	num=(igraph_vector_t*)rec->value;
	igraph_vector_destroy(num);
	igraph_free(num);
      } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
	str=(igraph_strvector_t*)rec->value;
	igraph_strvector_destroy(str);
	igraph_free(str);
      }
      igraph_free((char*)rec->name);
      igraph_free(rec);
    }
  }
  graph->attr=0;
}

/* Almost the same as destroy, but we might have null pointers */

void igraph_i_cattribute_copy_free(igraph_i_cattributes_t *attr) {
  igraph_vector_ptr_t *als[3] = { &attr->gal, &attr->val, &attr->eal };
  long int i, n, a;
  igraph_vector_t *num;
  igraph_strvector_t *str;
  igraph_i_attribute_record_t *rec;
  for (a=0; a<3; a++) {
    n=igraph_vector_ptr_size(als[a]);
    for (i=0; i<n; i++) {
      rec=VECTOR(*als[a])[i];
      if (!rec) { continue; }
      if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
	num=(igraph_vector_t*)rec->value;
	igraph_vector_destroy(num);
	igraph_free(num);
      } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
	str=(igraph_strvector_t*)rec->value;
	igraph_strvector_destroy(str);
	igraph_free(str);
      }
      igraph_free((char*)rec->name);
      igraph_free(rec);
    }
  }  
}

int igraph_i_cattributes_copy_attribute_record(igraph_i_attribute_record_t **newrec, 
					       const igraph_i_attribute_record_t *rec) {
  igraph_vector_t *num, *newnum;
  igraph_strvector_t *str, *newstr;
  
  *newrec=igraph_Calloc(1, igraph_i_attribute_record_t);
  if (!(*newrec)) { IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); }
  IGRAPH_FINALLY(igraph_free, *newrec);
  (*newrec)->type=rec->type;
  (*newrec)->name=strdup(rec->name);
  if (!(*newrec)->name) { IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); }
  IGRAPH_FINALLY(igraph_free, (void*)(*newrec)->name);
  if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
    num=(igraph_vector_t *)rec->value;
    newnum=igraph_Calloc(1, igraph_vector_t);
    if (!newnum) { 
      IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); 
    }
    IGRAPH_FINALLY(igraph_free, newnum);
    IGRAPH_CHECK(igraph_vector_copy(newnum, num));
    IGRAPH_FINALLY(igraph_vector_destroy, newnum);
    (*newrec)->value=newnum;
  } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
    str=(igraph_strvector_t*)rec->value;
    newstr=igraph_Calloc(1, igraph_strvector_t);
    if (!newstr) { 
      IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM); 
    }
    IGRAPH_FINALLY(igraph_free, newstr);
    IGRAPH_CHECK(igraph_strvector_copy(newstr, str));
    IGRAPH_FINALLY(igraph_strvector_destroy, newstr);
    (*newrec)->value=newstr;
  }

  IGRAPH_FINALLY_CLEAN(4);
  return 0;
}


/* No reference counting here. If you use attributes in C you should
   know what you're doing. */

int igraph_i_cattribute_copy(igraph_t *to, const igraph_t *from, 
			     igraph_bool_t ga, igraph_bool_t va, igraph_bool_t ea) {
  igraph_i_cattributes_t *attrfrom=from->attr, *attrto;
  igraph_vector_ptr_t *alto[3], *alfrom[3]={ &attrfrom->gal, &attrfrom->val, 
					     &attrfrom->eal };
  long int i, n, a;
  igraph_bool_t copy[3] = { ga, va, ea };
  to->attr=attrto=igraph_Calloc(1, igraph_i_cattributes_t);
  if (!attrto) { 
    IGRAPH_ERROR("Cannot copy attributes", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, attrto);
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->gal, 0);
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->val, 0);
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&attrto->eal, 0);
  IGRAPH_FINALLY_CLEAN(3);
  IGRAPH_FINALLY(igraph_i_cattribute_copy_free, attrto);

  alto[0]=&attrto->gal; alto[1]=&attrto->val; alto[2]=&attrto->eal;
  for (a=0; a<3; a++) {
    if (copy[a]) {
      n=igraph_vector_ptr_size(alfrom[a]);
      IGRAPH_CHECK(igraph_vector_ptr_resize(alto[a], n));
      igraph_vector_ptr_null(alto[a]);
      for (i=0; i<n; i++) {
	igraph_i_attribute_record_t *newrec;
	IGRAPH_CHECK(igraph_i_cattributes_copy_attribute_record(&newrec, 
								VECTOR(*alfrom[a])[i]));
	VECTOR(*alto[a])[i]=newrec;
      }
    }
  }
  
  IGRAPH_FINALLY_CLEAN(2);
  return 0;
}

int igraph_i_cattribute_add_vertices(igraph_t *graph, long int nv,
				     igraph_vector_ptr_t *nattr) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int length=igraph_vector_ptr_size(val);
  long int nattrno=nattr==NULL ? 0 : igraph_vector_ptr_size(nattr);
  long int origlen=igraph_vcount(graph)-nv;
  long int newattrs=0, i;
  igraph_vector_t news;
  
  /* First add the new attributes if any */
  newattrs=0;
  IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
  for (i=0; i<nattrno; i++) {
    igraph_i_attribute_record_t *nattr_entry=VECTOR(*nattr)[i];
    const char *nname=nattr_entry->name;
    long int j;
    igraph_bool_t l=igraph_i_cattribute_find(val, nname, &j);
    if (!l) {
      newattrs++;
      IGRAPH_CHECK(igraph_vector_push_back(&news, i));
    } else {
      /* check types */
      if (nattr_entry->type != 
	  ((igraph_i_attribute_record_t*)VECTOR(*val)[j])->type) {
	IGRAPH_ERROR("You cannot mix attribute types", IGRAPH_EINVAL);
      }
    }
  }

  /* Add NA/empty string vectors for the existing vertices */
  if (newattrs != 0) {
    for (i=0; i<newattrs; i++) {
      igraph_i_attribute_record_t *tmp=VECTOR(*nattr)[(long int)VECTOR(news)[i]];
      igraph_i_attribute_record_t *newrec=igraph_Calloc(1, igraph_i_attribute_record_t);
      igraph_attribute_type_t type=tmp->type;
      if (!newrec) { 
	IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, newrec);
      newrec->type=type;
      newrec->name=strdup(tmp->name);
      if (!newrec->name) { 
	IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, (char*)newrec->name);
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_vector_t *newnum=igraph_Calloc(1, igraph_vector_t);
	if (!newnum) {
	  IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
	}
	IGRAPH_FINALLY(igraph_free, newnum);
	IGRAPH_VECTOR_INIT_FINALLY(newnum, origlen);
	newrec->value=newnum;
	igraph_vector_fill(newnum, IGRAPH_NAN);
      } else if (type==IGRAPH_ATTRIBUTE_STRING) {
	igraph_strvector_t *newstr=igraph_Calloc(1, igraph_strvector_t);
	if (!newstr) {
	  IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
	}
	IGRAPH_FINALLY(igraph_free, newstr);
	IGRAPH_STRVECTOR_INIT_FINALLY(newstr, origlen);
	newrec->value=newstr;
      }
      IGRAPH_CHECK(igraph_vector_ptr_push_back(val, newrec));
      IGRAPH_FINALLY_CLEAN(4);
    }
    length=igraph_vector_ptr_size(val);
  }
  
  /* Now append the new values */
  for (i=0; i<length; i++) {
    igraph_i_attribute_record_t *oldrec=VECTOR(*val)[i];
    igraph_i_attribute_record_t *newrec=0;
    const char *name=oldrec->name;
    long int j;
    igraph_bool_t l=0;
    if (nattr) { l=igraph_i_cattribute_find(nattr, name, &j); }
    if (l) {
      /* This attribute is present in nattr */
      igraph_vector_t *oldnum, *newnum;
      igraph_strvector_t *oldstr, *newstr;
      newrec=VECTOR(*nattr)[j];
      oldnum=(igraph_vector_t*)oldrec->value;
      newnum=(igraph_vector_t*)newrec->value;
      oldstr=(igraph_strvector_t*)oldrec->value;
      newstr=(igraph_strvector_t*)newrec->value;
      if (oldrec->type != newrec->type) {
	IGRAPH_ERROR("Attribute types do not match", IGRAPH_EINVAL);
      }
      switch (oldrec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	if (nv != igraph_vector_size(newnum)) {
	  IGRAPH_ERROR("Invalid numeric attribute length", IGRAPH_EINVAL);
	}
	IGRAPH_CHECK(igraph_vector_append(oldnum, newnum));
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	if (nv != igraph_strvector_size(newstr)) {
	  IGRAPH_ERROR("Invalid string attribute length", IGRAPH_EINVAL);
	}
	IGRAPH_CHECK(igraph_strvector_append(oldstr, newstr));
	break;
      default:
	IGRAPH_WARNING("Invalid attribute type");	
	break;
      }
    } else {
      /* No such attribute, append NA's */
      igraph_vector_t *oldnum=(igraph_vector_t *)oldrec->value;
      igraph_strvector_t *oldstr=(igraph_strvector_t*)oldrec->value;
      switch (oldrec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	IGRAPH_CHECK(igraph_vector_resize(oldnum, origlen+nv));
	for (j=origlen; j<origlen+nv; j++) {
	  VECTOR(*oldnum)[j]=IGRAPH_NAN;
	}
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	IGRAPH_CHECK(igraph_strvector_resize(oldstr, origlen+nv));
	break;
      default:
	IGRAPH_WARNING("Invalid attribute type");
	break;
      }
    }
  }  
  
  return 0;
}

void igraph_i_cattribute_delete_vertices(igraph_t *graph,
				       const igraph_vector_t *eidx,
				       const igraph_vector_t *vidx) {
  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int valno=igraph_vector_ptr_size(val);
  long int ealno=igraph_vector_ptr_size(eal);
  long int i;
  long int origlen, newlen;
  
  /* Vertices */
  origlen=igraph_vector_size(vidx);
  newlen=0;
  for (i=0; i<origlen; i++) {
    if (VECTOR(*vidx)[i]>0) {
      newlen++;
    }
  }
  for (i=0; i<valno; i++) {
    igraph_i_attribute_record_t *oldrec=VECTOR(*val)[i];
    igraph_attribute_type_t type=oldrec->type;
    igraph_vector_t *num=(igraph_vector_t*)oldrec->value;
    igraph_strvector_t *str=(igraph_strvector_t*)oldrec->value;
    switch (type) {
    case IGRAPH_ATTRIBUTE_NUMERIC:
      igraph_vector_permdelete(num, vidx, origlen-newlen);
      break;
    case IGRAPH_ATTRIBUTE_STRING:
      igraph_strvector_permdelete(str, vidx, origlen-newlen);
      break;
    default:
      IGRAPH_WARNING("Unknown vertex attribute ignored");
    }
  }

  /* Edges */
  origlen=igraph_vector_size(eidx);
  newlen=0;
  for (i=0; i<origlen; i++) {
    if (VECTOR(*eidx)[i]>0) {
      newlen++;
    }
  }
  for (i=0; i<ealno; i++) {
    igraph_i_attribute_record_t *oldrec=VECTOR(*eal)[i];
    igraph_attribute_type_t type=oldrec->type;
    igraph_vector_t *num=(igraph_vector_t*)oldrec->value;
    igraph_strvector_t *str=(igraph_strvector_t*)oldrec->value;
    switch (type) {
    case IGRAPH_ATTRIBUTE_NUMERIC:
      igraph_vector_permdelete(num, eidx, origlen-newlen);
      break;
    case IGRAPH_ATTRIBUTE_STRING:
      igraph_strvector_permdelete(str, eidx, origlen-newlen);
      break;
    default:
      IGRAPH_WARNING("Unknown edge attribute ignored");
    }
  }
}

int igraph_i_cattribute_add_edges(igraph_t *graph, const igraph_vector_t *edges,
				 igraph_vector_ptr_t *nattr) {
  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int ealno=igraph_vector_ptr_size(eal);
  long int ne=igraph_vector_size(edges)/2;
  long int origlen=igraph_ecount(graph)-ne;
  long int nattrno= nattr == 0 ? 0 : igraph_vector_ptr_size(nattr);
  igraph_vector_t news;
  long int newattrs, i;
  
  /* First add the new attributes if any */
  newattrs=0;
  IGRAPH_VECTOR_INIT_FINALLY(&news, 0);
  for (i=0; i<nattrno; i++) {
    igraph_i_attribute_record_t *nattr_entry=VECTOR(*nattr)[i];
    const char *nname=nattr_entry->name;
    long int j;
    igraph_bool_t l=igraph_i_cattribute_find(eal, nname, &j);
    if (!l) {
      newattrs++;
      IGRAPH_CHECK(igraph_vector_push_back(&news, i));
    } else {
      /* check types */
      if (nattr_entry->type != 
	  ((igraph_i_attribute_record_t*)VECTOR(*eal)[j])->type) {
	IGRAPH_ERROR("You cannot mix attribute types", IGRAPH_EINVAL);
      }
    }
  }

  /* Add NA/empty string vectors for the existing vertices */
  if (newattrs != 0) {
    for (i=0; i<newattrs; i++) {
      igraph_i_attribute_record_t *tmp=VECTOR(*nattr)[(long int)VECTOR(news)[i]];
      igraph_i_attribute_record_t *newrec=igraph_Calloc(1, igraph_i_attribute_record_t);
      igraph_attribute_type_t type=tmp->type;
      if (!newrec) { 
	IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, newrec);
      newrec->type=type;
      newrec->name=strdup(tmp->name);
      if (!newrec->name) { 
	IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
      }
      IGRAPH_FINALLY(igraph_free, (char*)newrec->name);
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_vector_t *newnum=igraph_Calloc(1, igraph_vector_t);
	if (!newnum) {
	  IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
	}
	IGRAPH_FINALLY(igraph_free, newnum);
	IGRAPH_VECTOR_INIT_FINALLY(newnum, origlen);
	newrec->value=newnum;
	igraph_vector_fill(newnum, IGRAPH_NAN);
      } else if (type==IGRAPH_ATTRIBUTE_STRING) {
	igraph_strvector_t *newstr=igraph_Calloc(1, igraph_strvector_t);
	if (!newstr) {
	  IGRAPH_ERROR("Cannot add attributes", IGRAPH_ENOMEM);
	}
	IGRAPH_FINALLY(igraph_free, newstr);
	IGRAPH_STRVECTOR_INIT_FINALLY(newstr, origlen);
	newrec->value=newstr;	
      }
      IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, newrec));
      IGRAPH_FINALLY_CLEAN(4);
    }
    ealno=igraph_vector_ptr_size(eal);
  }
  
  /* Now append the new values */
  for (i=0; i<ealno; i++) {
    igraph_i_attribute_record_t *oldrec=VECTOR(*eal)[i];
    igraph_i_attribute_record_t *newrec=0;
    const char *name=oldrec->name;
    long int j;
    igraph_bool_t l=0;
    if (nattr) { l=igraph_i_cattribute_find(nattr, name, &j); }
    if (l) {
      /* This attribute is present in nattr */
      igraph_vector_t *oldnum, *newnum;
      igraph_strvector_t *oldstr, *newstr;
      newrec=VECTOR(*nattr)[j];
      oldnum=(igraph_vector_t*)oldrec->value;
      newnum=(igraph_vector_t*)newrec->value;
      oldstr=(igraph_strvector_t*)oldrec->value;
      newstr=(igraph_strvector_t*)newrec->value;
      if (oldrec->type != newrec->type) {
	IGRAPH_ERROR("Attribute types do not match", IGRAPH_EINVAL);
      }
      switch (oldrec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	if (ne != igraph_vector_size(newnum)) {
	  IGRAPH_ERROR("Invalid numeric attribute length", IGRAPH_EINVAL);
	}
	IGRAPH_CHECK(igraph_vector_append(oldnum, newnum));
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	if (ne != igraph_strvector_size(newstr)) {
	  IGRAPH_ERROR("Invalid string attribute length", IGRAPH_EINVAL);
	}
	IGRAPH_CHECK(igraph_strvector_append(oldstr, newstr));
	break;
      default:
	IGRAPH_WARNING("Invalid attribute type");	
	break;
      }
    } else {
      /* No such attribute, append NA's */
      igraph_vector_t *oldnum=(igraph_vector_t *)oldrec->value;
      igraph_strvector_t *oldstr=(igraph_strvector_t*)oldrec->value;
      switch (oldrec->type) {
      case IGRAPH_ATTRIBUTE_NUMERIC:
	IGRAPH_CHECK(igraph_vector_resize(oldnum, origlen+ne));
	for (j=origlen; j<origlen+ne; j++) {
	  VECTOR(*oldnum)[j]=IGRAPH_NAN;
	}
	break;
      case IGRAPH_ATTRIBUTE_STRING:
	IGRAPH_CHECK(igraph_strvector_resize(oldstr, origlen+ne));
	break;
      default:
	IGRAPH_WARNING("Invalid attribute type");
	break;
      }
    }
  }  
   
  return 0;
}

void igraph_i_cattribute_delete_edges(igraph_t *graph, const igraph_vector_t *idx) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int ealno=igraph_vector_ptr_size(eal);
  long int i;
  long int origlen=igraph_vector_size(idx), newlen;

  newlen=0;
  for (i=0; i<origlen; i++) {
    if (VECTOR(*idx)[i]>0) {
      newlen++;
    }
  }
  for (i=0; i<ealno; i++) {
    igraph_i_attribute_record_t *oldrec=VECTOR(*eal)[i];
    igraph_attribute_type_t type=oldrec->type;
    igraph_vector_t *num=(igraph_vector_t*)oldrec->value;
    igraph_strvector_t *str=(igraph_strvector_t*)oldrec->value;
    switch (type) {
    case IGRAPH_ATTRIBUTE_NUMERIC:
      igraph_vector_permdelete(num, idx, origlen-newlen);
      break;
    case IGRAPH_ATTRIBUTE_STRING:
      igraph_strvector_permdelete(str, idx, origlen-newlen);
      break;
    default:
      IGRAPH_WARNING("Unknown edge attribute ignored");
    }
  }
  
}

/* Is this just another name for this? */

int igraph_i_cattribute_permute_edges(igraph_t *graph,
				    const igraph_vector_t *idx) {
  
  igraph_i_cattribute_delete_edges(graph, idx);
  return 0;
}



int igraph_i_cattribute_get_info(const igraph_t *graph,
				 igraph_strvector_t *gnames,
				 igraph_vector_t *gtypes,
				 igraph_strvector_t *vnames,
				 igraph_vector_t *vtypes,
				 igraph_strvector_t *enames,
				 igraph_vector_t *etypes) {
  
  igraph_strvector_t *names[3] = { gnames, vnames, enames };
  igraph_vector_t *types[3] = { gtypes, vtypes, etypes };
  igraph_i_cattributes_t *at=graph->attr;
  igraph_vector_ptr_t *attr[3]={ &at->gal, &at->val, &at->eal };
  long int i,j;

  for (i=0; i<3; i++) {
    igraph_strvector_t *n=names[i];
    igraph_vector_t *t=types[i];
    igraph_vector_ptr_t *al=attr[i];
    long int len=igraph_vector_ptr_size(al);
    
    if (n) {
      IGRAPH_CHECK(igraph_strvector_resize(n, len));
    }
    if (t) {
      IGRAPH_CHECK(igraph_vector_resize(t, len));
    }

    for (j=0; j<len; j++) {
      igraph_i_attribute_record_t *rec=VECTOR(*al)[j];
      const char *name=rec->name;
      igraph_attribute_type_t type=rec->type;      
      if (n) {
	IGRAPH_CHECK(igraph_strvector_set(n, j, name));
      }
      if (t) {
	VECTOR(*t)[j]=type;
      }
    }
  }
  
  return 0;
}

igraph_bool_t igraph_i_cattribute_has_attr(const igraph_t *graph,
					 igraph_attribute_elemtype_t type,
					 const char *name) {
  igraph_i_cattributes_t *at=graph->attr;
  igraph_vector_ptr_t *attr[3]={ &at->gal, &at->val, &at->eal };
  long int attrnum;

  switch (type) {
  case IGRAPH_ATTRIBUTE_GRAPH:
    attrnum=0;
    break;
  case IGRAPH_ATTRIBUTE_VERTEX:
    attrnum=1;
    break;
  case IGRAPH_ATTRIBUTE_EDGE:
    attrnum=2;
    break;
  default:
    IGRAPH_ERROR("Unknown attribute element type", IGRAPH_EINVAL);
    break;
  }

  return igraph_i_cattribute_find(attr[attrnum], name, 0);
}

int igraph_i_cattribute_gettype(const igraph_t *graph,
			      igraph_attribute_type_t *type,
			      igraph_attribute_elemtype_t elemtype,
			      const char *name) {
  long int attrnum;
  igraph_i_attribute_record_t *rec;  
  igraph_i_cattributes_t *at=graph->attr;
  igraph_vector_ptr_t *attr[3]={ &at->gal, &at->val, &at->eal };
  igraph_vector_ptr_t *al;
  long int j;
  igraph_bool_t l=0;
  
  switch (elemtype) {
  case IGRAPH_ATTRIBUTE_GRAPH:
    attrnum=0;
    break;
  case IGRAPH_ATTRIBUTE_VERTEX:
    attrnum=1;
    break;
  case IGRAPH_ATTRIBUTE_EDGE:
    attrnum=2;
    break;
  default:
    IGRAPH_ERROR("Unknown attribute element type", IGRAPH_EINVAL);
    break;
  }

  al=attr[attrnum];
  l=igraph_i_cattribute_find(al, name, &j);  
  if (!l) {
    IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
  }  
  rec=VECTOR(*al)[j];
  *type=rec->type;

  return 0;
}

int igraph_i_cattribute_get_numeric_graph_attr(const igraph_t *graph,
					      const char *name,
					      igraph_vector_t *value) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *gal=&attr->gal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_vector_t *num;
  igraph_bool_t l=igraph_i_cattribute_find(gal, name, &j);
  
  if (!l) {
    IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
  }

  rec=VECTOR(*gal)[j];
  num=(igraph_vector_t*)rec->value;
  IGRAPH_CHECK(igraph_vector_resize(value, 1));
  VECTOR(*value)[0]=VECTOR(*num)[0];
  
  return 0;
}

int igraph_i_cattribute_get_string_graph_attr(const igraph_t *graph,
					     const char *name,
					     igraph_strvector_t *value) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *gal=&attr->gal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_strvector_t *str;
  igraph_bool_t l=igraph_i_cattribute_find(gal, name, &j);
  
  if (!l) {
    IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
  }

  rec=VECTOR(*gal)[j];
  str=(igraph_strvector_t*)rec->value;
  IGRAPH_CHECK(igraph_strvector_resize(value, 1));
  IGRAPH_CHECK(igraph_strvector_set(value, 0, STR(*str,0)));

  return 0;
}

int igraph_i_cattribute_get_numeric_vertex_attr(const igraph_t *graph,
					      const char *name,
					      igraph_vs_t vs,
					      igraph_vector_t *value) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_vector_t *num;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  if (!l) {
    IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
  }

  rec=VECTOR(*val)[j];
  num=(igraph_vector_t*)rec->value;
  if (igraph_vs_is_all(&vs)) {
    igraph_vector_clear(value);
    IGRAPH_CHECK(igraph_vector_append(value, num));
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_VIT_SIZE(it)));
    for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
      long int v=IGRAPH_VIT_GET(it);
      VECTOR(*value)[i]=VECTOR(*num)[v];
    }
    igraph_vit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

int igraph_i_cattribute_get_string_vertex_attr(const igraph_t *graph,
					     const char *name,
					     igraph_vs_t vs,
					     igraph_strvector_t *value) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_strvector_t *str;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  if (!l) {
    IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
  }

  rec=VECTOR(*val)[j];
  str=(igraph_strvector_t*)rec->value;
  if (igraph_vs_is_all(&vs)) {
    igraph_strvector_resize(value, 0);
    IGRAPH_CHECK(igraph_strvector_append(value, str));
  } else {
    igraph_vit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_vit_create(graph, vs, &it));
    IGRAPH_FINALLY(igraph_vit_destroy, &it);
    IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_VIT_SIZE(it)));
    for (; !IGRAPH_VIT_END(it); IGRAPH_VIT_NEXT(it), i++) {
      long int v=IGRAPH_VIT_GET(it);
      char *s;
      igraph_strvector_get(str, v, &s);
      IGRAPH_CHECK(igraph_strvector_set(str, i, s));
    }
    igraph_vit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

int igraph_i_cattribute_get_numeric_edge_attr(const igraph_t *graph,
					    const char *name,
					    igraph_es_t es,
					    igraph_vector_t *value) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_vector_t *num;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  if (!l) {
    IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
  }

  rec=VECTOR(*eal)[j];
  num=(igraph_vector_t*)rec->value;
  if (igraph_es_is_all(&es)) {
    igraph_vector_clear(value);
    IGRAPH_CHECK(igraph_vector_append(value, num));
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_vector_resize(value, IGRAPH_EIT_SIZE(it)));
    for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
      long int e=IGRAPH_EIT_GET(it);
      VECTOR(*value)[i]=VECTOR(*num)[e];
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }

  return 0;
}

int igraph_i_cattribute_get_string_edge_attr(const igraph_t *graph,
					   const char *name,
					   igraph_es_t es,
					   igraph_strvector_t *value) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_strvector_t *str;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  if (!l) {
    IGRAPH_ERROR("Unknown attribute", IGRAPH_EINVAL);
  }

  rec=VECTOR(*eal)[j];
  str=(igraph_strvector_t*)rec->value;
  if (igraph_es_is_all(&es)) {
    igraph_strvector_resize(value, 0);
    IGRAPH_CHECK(igraph_strvector_append(value, str));
  } else {
    igraph_eit_t it;
    long int i=0;
    IGRAPH_CHECK(igraph_eit_create(graph, es, &it));
    IGRAPH_FINALLY(igraph_eit_destroy, &it);
    IGRAPH_CHECK(igraph_strvector_resize(value, IGRAPH_EIT_SIZE(it)));
    for (; !IGRAPH_EIT_END(it); IGRAPH_EIT_NEXT(it), i++) {
      long int e=IGRAPH_EIT_GET(it);
      char *s;
      igraph_strvector_get(str, e, &s);
      IGRAPH_CHECK(igraph_strvector_set(str, i, s));
    }
    igraph_eit_destroy(&it);
    IGRAPH_FINALLY_CLEAN(1);
  }
  
  return 0;
}

/* -------------------------------------- */

igraph_attribute_table_t igraph_cattribute_table={
  &igraph_i_cattribute_init, &igraph_i_cattribute_destroy,
  &igraph_i_cattribute_copy, &igraph_i_cattribute_add_vertices,
  &igraph_i_cattribute_delete_vertices, &igraph_i_cattribute_add_edges,
  &igraph_i_cattribute_delete_edges, &igraph_i_cattribute_permute_edges,
  &igraph_i_cattribute_get_info,
  &igraph_i_cattribute_has_attr, &igraph_i_cattribute_gettype,
  &igraph_i_cattribute_get_numeric_graph_attr,
  &igraph_i_cattribute_get_string_graph_attr,
  &igraph_i_cattribute_get_numeric_vertex_attr,
  &igraph_i_cattribute_get_string_vertex_attr,
  &igraph_i_cattribute_get_numeric_edge_attr,
  &igraph_i_cattribute_get_string_edge_attr
};

/* -------------------------------------- */

igraph_real_t igraph_cattribute_GAN(const igraph_t *graph, const char *name) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *gal=&attr->gal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_vector_t *num;
  igraph_bool_t l=igraph_i_cattribute_find(gal, name, &j);
  
  if (!l) {
    igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
    return 0;
  }

  rec=VECTOR(*gal)[j];
  num=(igraph_vector_t*)rec->value;
  return VECTOR(*num)[0];
}

const char* igraph_cattribute_GAS(const igraph_t *graph, const char *name) {
  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *gal=&attr->gal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_strvector_t *str;
  igraph_bool_t l=igraph_i_cattribute_find(gal, name, &j);
  
  if (!l) {
    igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
    return 0;
  }

  rec=VECTOR(*gal)[j];
  str=(igraph_strvector_t*)rec->value;
  return STR(*str, 0);  
}

igraph_real_t igraph_cattribute_VAN(const igraph_t *graph, const char *name,
				      igraph_integer_t vid) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_vector_t *num;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  if (!l) {
    igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
    return 0;
  }

  rec=VECTOR(*val)[j];
  num=(igraph_vector_t*)rec->value;
  return VECTOR(*num)[(long int)vid];
}

const char* igraph_cattribute_VAS(const igraph_t *graph, const char *name,
				    igraph_integer_t vid) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_strvector_t *str;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  if (!l) {
    igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
    return 0;
  }

  rec=VECTOR(*val)[j];
  str=(igraph_strvector_t*)rec->value;
  return STR(*str, (long int)vid);  
}

igraph_real_t igraph_cattribute_EAN(const igraph_t *graph, const char *name,
				      igraph_integer_t eid) {  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_vector_t *num;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  if (!l) {
    igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
    return 0;
  }

  rec=VECTOR(*eal)[j];
  num=(igraph_vector_t*)rec->value;
  return VECTOR(*num)[(long int)eid];
}

const char* igraph_cattribute_EAS(const igraph_t *graph, const char *name,
				    igraph_integer_t eid) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_i_attribute_record_t *rec;
  igraph_strvector_t *str;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  if (!l) {
    igraph_error("Unknown attribute", __FILE__, __LINE__, IGRAPH_EINVAL);
    return 0;
  }
  
  rec=VECTOR(*eal)[j];
  str=(igraph_strvector_t*)rec->value;
  return STR(*str, (long int)eid);  
}

int igraph_cattribute_list(const igraph_t *graph,
			   igraph_strvector_t *gnames, igraph_vector_t *gtypes,
			   igraph_strvector_t *vnames, igraph_vector_t *vtypes,
			   igraph_strvector_t *enames, igraph_vector_t *etypes) {
  return igraph_i_cattribute_get_info(graph, gnames, gtypes, vnames, vtypes,
				      enames, etypes);
}

int igraph_cattribute_GAN_set(igraph_t *graph, const char *name, 
			      igraph_real_t value) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *gal=&attr->gal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(gal, name, &j);
  
  if (l) {
    igraph_i_attribute_record_t *rec=VECTOR(*gal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
      IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
    } else {
      igraph_vector_t *num=(igraph_vector_t *)rec->value;
      VECTOR(*num)[0]=value;
    }
  } else {
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_vector_t *num;
    if (!rec) {
      IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    rec->type=IGRAPH_ATTRIBUTE_NUMERIC;
    num=igraph_Calloc(1, igraph_vector_t);
    if (!num) { 
      IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, num);
    IGRAPH_VECTOR_INIT_FINALLY(num, 1);
    VECTOR(*num)[0]=value;
    rec->value=num;
    IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }

  return 0;
}

int igraph_cattribute_GAS_set(igraph_t *graph, const char *name, 
			      const char *value) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *gal=&attr->gal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(gal, name, &j);
  
  if (l) {
    igraph_i_attribute_record_t *rec=VECTOR(*gal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
      IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
    } else {
      igraph_strvector_t *str=(igraph_strvector_t*)rec->value;
      IGRAPH_CHECK(igraph_strvector_set(str, 0, value));
    }
  } else {
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_strvector_t *str;
    if (!rec) {
      IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    rec->type=IGRAPH_ATTRIBUTE_STRING;
    str=igraph_Calloc(1, igraph_strvector_t);
    if (!str) {
      IGRAPH_ERROR("Cannot add graph attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, str);
    IGRAPH_STRVECTOR_INIT_FINALLY(str, 1);
    IGRAPH_CHECK(igraph_strvector_set(str, 0, value));
    rec->value=str;
    IGRAPH_CHECK(igraph_vector_ptr_push_back(gal, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }

  return 0;
}

int igraph_cattribute_VAN_set(igraph_t *graph, const char *name, 
			      igraph_integer_t vid, igraph_real_t value) {
  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  if (l) {
    igraph_i_attribute_record_t *rec=VECTOR(*val)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
      IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
    } else {
      igraph_vector_t *num=(igraph_vector_t*)rec->value;
      VECTOR(*num)[(long int)vid]=value;
    }
  } else {
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_vector_t *num;
    if (!rec) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    rec->type=IGRAPH_ATTRIBUTE_NUMERIC;
    num=igraph_Calloc(1, igraph_vector_t);
    if (!num) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, num);
    IGRAPH_VECTOR_INIT_FINALLY(num, igraph_vcount(graph));
    igraph_vector_fill(num, IGRAPH_NAN);
    VECTOR(*num)[(long int)vid]=value;
    rec->value=num;
    IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }
  
  return 0;
}

int igraph_cattribute_VAS_set(igraph_t *graph, const char *name, 
			      igraph_integer_t vid, const char *value) {
  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  if (l) {
    igraph_i_attribute_record_t *rec=VECTOR(*val)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
      IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
    } else {
      igraph_strvector_t *str=(igraph_strvector_t*)rec->value;
      IGRAPH_CHECK(igraph_strvector_set(str, vid, value));
    }
  } else {
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_strvector_t *str;
    if (!rec) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    rec->type=IGRAPH_ATTRIBUTE_STRING;
    str=igraph_Calloc(1, igraph_strvector_t);
    if (!str) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, str);
    IGRAPH_STRVECTOR_INIT_FINALLY(str, igraph_vcount(graph));
    IGRAPH_CHECK(igraph_strvector_set(str, vid, value));
    rec->value=str;
    IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }
  
  return 0;
}

int igraph_cattribute_EAN_set(igraph_t *graph, const char *name, 
			      igraph_integer_t eid, igraph_real_t value) {
  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  if (l) {
    igraph_i_attribute_record_t *rec=VECTOR(*eal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
      IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
    } else {
      igraph_vector_t *num=(igraph_vector_t*)rec->value;
      VECTOR(*num)[(long int)eid]=value;
    }
  } else {
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_vector_t *num;
    if (!rec) {
      IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    rec->type=IGRAPH_ATTRIBUTE_NUMERIC;
    num=igraph_Calloc(1, igraph_vector_t);
    if (!num) {
      IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, num);
    IGRAPH_VECTOR_INIT_FINALLY(num, igraph_ecount(graph));
    igraph_vector_fill(num, IGRAPH_NAN);
    VECTOR(*num)[(long int)eid]=value;
    rec->value=num;
    IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }
  
  return 0;
}

int igraph_cattribute_EAS_set(igraph_t *graph, const char *name, 
			      igraph_integer_t eid, const char *value) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  if (l) {
    igraph_i_attribute_record_t *rec=VECTOR(*eal)[j];
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
      IGRAPH_ERROR("Invalid attribute type", IGRAPH_EINVAL);
    } else {
      igraph_strvector_t *str=(igraph_strvector_t*)rec->value;
      IGRAPH_CHECK(igraph_strvector_set(str, eid, value));
    }
  } else {
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_strvector_t *str;
    if (!rec) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    rec->type=IGRAPH_ATTRIBUTE_STRING;
    str=igraph_Calloc(1, igraph_strvector_t);
    if (!str) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, str);
    IGRAPH_STRVECTOR_INIT_FINALLY(str, igraph_vcount(graph));
    IGRAPH_CHECK(igraph_strvector_set(str, eid, value));
    rec->value=str;
    IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }

  return 0;
}

int igraph_cattribute_VAN_setv(igraph_t *graph, const char *name, 
			       const igraph_vector_t *v) {
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);

  /* Check length first */
  if (igraph_vector_size(v) != igraph_vcount(graph)) {
    IGRAPH_ERROR("Invalid vertex attribute vector length", IGRAPH_EINVAL);
  }
  
  if (l) {
    /* Already present, check type */
    igraph_i_attribute_record_t *rec=VECTOR(*val)[j];
    igraph_vector_t *num=(igraph_vector_t *)rec->value;
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
      IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
    }
    igraph_vector_clear(num);
    IGRAPH_CHECK(igraph_vector_append(num, v));
  } else {
    /* Add it */
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_vector_t *num;
    if (!rec) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->type=IGRAPH_ATTRIBUTE_NUMERIC;
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    num=igraph_Calloc(1, igraph_vector_t);
    if (!num) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, num);
    rec->value=num;
    IGRAPH_CHECK(igraph_vector_copy(num, v));
    IGRAPH_FINALLY(igraph_vector_destroy, num);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }

  return 0;
}

int igraph_cattribute_VAS_setv(igraph_t *graph, const char *name,
			       const igraph_strvector_t *sv) {
  
  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  /* Check length first */
  if (igraph_strvector_size(sv) != igraph_vcount(graph)) {
    IGRAPH_ERROR("Invalid vertex attribute vector length", IGRAPH_EINVAL);
  }

  if (l) {
    /* Already present, check type */
    igraph_i_attribute_record_t *rec=VECTOR(*val)[j];
    igraph_strvector_t *str=(igraph_strvector_t *)rec->value;
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
      IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
    }
    igraph_strvector_clear(str);
    IGRAPH_CHECK(igraph_strvector_append(str, sv));
  } else { 
    /* Add it */
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_strvector_t *str;
    if (!rec) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->type=IGRAPH_ATTRIBUTE_STRING;
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    str=igraph_Calloc(1, igraph_strvector_t);
    if (!str) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, str);
    rec->value=str;
    IGRAPH_CHECK(igraph_strvector_copy(str, sv));
    IGRAPH_FINALLY(igraph_strvector_destroy, str);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(val, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }

  return 0;
}

int igraph_cattribute_EAN_setv(igraph_t *graph, const char *name, 
			       const igraph_vector_t *v) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);

  /* Check length first */
  if (igraph_vector_size(v) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid edge attribute vector length", IGRAPH_EINVAL);
  }
  
  if (l) {
    /* Already present, check type */
    igraph_i_attribute_record_t *rec=VECTOR(*eal)[j];
    igraph_vector_t *num=(igraph_vector_t *)rec->value;
    if (rec->type != IGRAPH_ATTRIBUTE_NUMERIC) {
      IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
    }
    igraph_vector_clear(num);
    IGRAPH_CHECK(igraph_vector_append(num, v));
  } else {
    /* Add it */
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_vector_t *num;
    if (!rec) {
      IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->type=IGRAPH_ATTRIBUTE_NUMERIC;
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    num=igraph_Calloc(1, igraph_vector_t);
    if (!num) {
      IGRAPH_ERROR("Cannot add edge attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, num);
    rec->value=num;
    IGRAPH_CHECK(igraph_vector_copy(num, v));
    IGRAPH_FINALLY(igraph_vector_destroy, num);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }

  return 0;
}

int igraph_cattribute_EAS_setv(igraph_t *graph, const char *name,
			       const igraph_strvector_t *sv) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  /* Check length first */
  if (igraph_strvector_size(sv) != igraph_ecount(graph)) {
    IGRAPH_ERROR("Invalid edge attribute vector length", IGRAPH_EINVAL);
  }

  if (l) {
    /* Already present, check type */
    igraph_i_attribute_record_t *rec=VECTOR(*eal)[j];
    igraph_strvector_t *str=(igraph_strvector_t *)rec->value;
    if (rec->type != IGRAPH_ATTRIBUTE_STRING) {
      IGRAPH_ERROR("Attribute type mismatch", IGRAPH_EINVAL);
    }
    igraph_strvector_clear(str);
    IGRAPH_CHECK(igraph_strvector_append(str, sv));
  } else { 
    /* Add it */
    igraph_i_attribute_record_t *rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    igraph_strvector_t *str;
    if (!rec) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, rec);
    rec->type=IGRAPH_ATTRIBUTE_STRING;
    rec->name=strdup(name);
    if (!rec->name) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, (char*)rec->name);
    str=igraph_Calloc(1, igraph_strvector_t);
    if (!str) {
      IGRAPH_ERROR("Cannot add vertex attribute", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, str);
    rec->value=str;
    IGRAPH_CHECK(igraph_strvector_copy(str, sv));
    IGRAPH_FINALLY(igraph_strvector_destroy, str);
    IGRAPH_CHECK(igraph_vector_ptr_push_back(eal, rec));
    IGRAPH_FINALLY_CLEAN(4);
  }

  return 0;
}

void igraph_i_cattribute_free_rec(igraph_i_attribute_record_t *rec) {
  
  if (rec->type==IGRAPH_ATTRIBUTE_NUMERIC) {
    igraph_vector_t *num=(igraph_vector_t*)rec->value;
    igraph_vector_destroy(num);
  } else if (rec->type==IGRAPH_ATTRIBUTE_STRING) {
    igraph_strvector_t *str=(igraph_strvector_t*)rec->value;
    igraph_strvector_destroy(str);
  }
  igraph_Free(rec->name);
  igraph_Free(rec->value);
  igraph_Free(rec);
}

void igraph_cattribute_remove_g(igraph_t *graph, const char *name) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *gal=&attr->gal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(gal, name, &j);
  
  if (l) {
    igraph_i_cattribute_free_rec(VECTOR(*gal)[j]);
    igraph_vector_ptr_remove(gal, j);
  } else {
    IGRAPH_WARNING("Cannot remove non-existant graph attribute");
  }  
}

void igraph_cattribute_remove_v(igraph_t *graph, const char *name) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *val=&attr->val;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(val, name, &j);
  
  if (l) {
    igraph_i_cattribute_free_rec(VECTOR(*val)[j]);
    igraph_vector_ptr_remove(val, j);
  } else {
    IGRAPH_WARNING("Cannot remove non-existant graph attribute");
  }  
}

void igraph_cattribute_remove_e(igraph_t *graph, const char *name) {

  igraph_i_cattributes_t *attr=graph->attr;
  igraph_vector_ptr_t *eal=&attr->eal;
  long int j;
  igraph_bool_t l=igraph_i_cattribute_find(eal, name, &j);
  
  if (l) {
    igraph_i_cattribute_free_rec(VECTOR(*eal)[j]);
    igraph_vector_ptr_remove(eal, j);
  } else {
    IGRAPH_WARNING("Cannot remove non-existant graph attribute");
  }  
}

void igraph_cattribute_remove_all(igraph_t *graph, igraph_bool_t g,
				  igraph_bool_t v, igraph_bool_t e) {

  igraph_i_cattributes_t *attr=graph->attr;

  if (g) {
    igraph_vector_ptr_t *gal=&attr->gal;
    long int i, n=igraph_vector_ptr_size(gal);
    for (i=0;i<n;i++) {
      igraph_i_cattribute_free_rec(VECTOR(*gal)[i]);
    }
    igraph_vector_ptr_clear(gal);
  }
  if (v) {
    igraph_vector_ptr_t *val=&attr->val;
    long int i, n=igraph_vector_ptr_size(val);
    for (i=0;i<n;i++) {
      igraph_i_cattribute_free_rec(VECTOR(*val)[i]);
    }
    igraph_vector_ptr_clear(val);    
  }
  if (e) {
    igraph_vector_ptr_t *eal=&attr->eal;
    long int i, n=igraph_vector_ptr_size(eal);
    for (i=0;i<n;i++) {
      igraph_i_cattribute_free_rec(VECTOR(*eal)[i]);
    }
    igraph_vector_ptr_clear(eal);
  }
}
