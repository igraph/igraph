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

#include "attributes.h"
#include "memory.h"
#include "igraph.h"

#include <string.h>

/**
 * \ingroup internal
 */

long int igraph_i_attribute_list_get_pos(const igraph_attribute_list_t *al, 
					 const char *name) {
  long int pos=-1;
  long int n=igraph_strvector_size(&al->names);
  bool_t l=0;
  char *str;
  
  while(!l && pos < n-1) {
    igraph_strvector_get(&al->names, pos+1, &str);
    l=!strcmp(name, str);
    pos++;
  }
  
  if (!l) {
    pos = -1;
  }

  return pos;
}

/**
 * \ingroup internal
 */

void igraph_i_attribute_list_free(igraph_attribute_list_t *al, long int i) {
  if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t *numv=VECTOR(al->data)[i];
    igraph_vector_destroy(numv);
  } else if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_STR) {
    igraph_strvector_t *strv=VECTOR(al->data)[i];
    igraph_strvector_destroy(strv);
  }
}

/**
 * \ingroup internal
 * \brief Initializes an attribute list
 */

int igraph_attribute_list_init(igraph_attribute_list_t *al, long int len) {
  al->len=len;
  IGRAPH_STRVECTOR_INIT_FINALLY(&al->names, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&al->types, 0);
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&al->data, 0);
  IGRAPH_FINALLY_CLEAN(3);
  return 0;
}

/**
 * \ingroup internal
 * \brief Frees the memory allocated for an attribute list
 */

void igraph_attribute_list_destroy(igraph_attribute_list_t *al) {
  long int i;

  for (i=0; i<igraph_vector_size(&al->types); i++) {
    igraph_i_attribute_list_free(al, i);
  }

  igraph_strvector_destroy(&al->names);
  igraph_vector_destroy(&al->types);
  igraph_vector_ptr_destroy_all(&al->data);
}

/**
 * \ingroup internal
 * \brief Adds a new attribute to an attribute list.
 */

int igraph_attribute_list_add(igraph_attribute_list_t *al,
			      const char *name, igraph_attribute_type_t type){
  long int pos;
  void *data=NULL;

  /* Checks */
  if (strlen(name)==0) {
    IGRAPH_ERROR("invalid attribute name", IGRAPH_EINVAL);
  }
  pos=igraph_i_attribute_list_get_pos(al, name);
  if (pos >= 0) {
    IGRAPH_ERROR("attribute already exists", IGRAPH_EXISTS);
  }

  if (type==IGRAPH_ATTRIBUTE_NUM) {
    data=(void*) Calloc(1, igraph_vector_t);
    if (data != 0) { 
      IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)data, al->len); 
    } else {
      IGRAPH_ERROR("cannot add attribute", IGRAPH_ENOMEM);
    }
  } else /* if (type==IGRAPH_ATTRIBUTE_STR) */ {
    data=(void*)Calloc(1, igraph_strvector_t);
    if (data != 0) { 
      IGRAPH_STRVECTOR_INIT_FINALLY((igraph_strvector_t*)data, al->len);
    } else {
      IGRAPH_ERROR("cannot add attribute", IGRAPH_ENOMEM);
    }
  }

  IGRAPH_CHECK(igraph_vector_ptr_reserve(&al->data, igraph_vector_ptr_size(&al->data)+1));
  IGRAPH_CHECK(igraph_vector_reserve(&al->types, igraph_vector_size(&al->types)+1));
  IGRAPH_CHECK(igraph_strvector_add(&al->names, name));
    
  igraph_vector_ptr_push_back(&al->data, data);
  igraph_vector_push_back(&al->types, type);
  /* Space is allocated already, no need to check errors... */

  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup internal
 * \brief Removes an attribute from an attribute list.
 */

int igraph_attribute_list_remove(igraph_attribute_list_t *al, 
				 const char *name) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);
  void* ptr;
  
  if (pos < 0) {
    IGRAPH_ERROR("no such attribute", IGRAPH_EINVAL);
  }

  igraph_i_attribute_list_free(al, pos);
  ptr=igraph_vector_ptr_e(&al->data, pos);
  Free(ptr);
  igraph_vector_ptr_remove(&al->data, pos);
  igraph_vector_remove(&al->types, pos);
  igraph_strvector_remove(&al->names, pos);
  return 0;
}

/**
 * \ingroup internal
 * \brief Returns an attribute for a single element 
 */

int igraph_attribute_list_get(const igraph_attribute_list_t *al, 
			      const char *name, 
			      long int idx, void **value, 
			      igraph_attribute_type_t *type) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);
  igraph_attribute_type_t atype;

  if (pos < 0) {
    IGRAPH_ERROR("no such attribute", IGRAPH_EINVAL);
  }

  atype=VECTOR(al->types)[pos];
  if (type != 0) {
    *type = atype;
  }
  if (atype==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t *data=VECTOR(al->data)[pos];
    *value=(void*) igraph_vector_e_ptr(data, idx);
  } else if (atype==IGRAPH_ATTRIBUTE_STR) {
    igraph_strvector_t *data=VECTOR(al->data)[pos];
    igraph_strvector_get(data, idx, (char**)value);
  }
  return 0;
}

/**
 * \ingroup internal
 * \brief Sets an attribute for a single element
 */

int igraph_attribute_list_set(igraph_attribute_list_t *al, const char *name,
			      long int idx, const void *value) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);
  igraph_attribute_type_t atype;

  if (pos < 0) {
    IGRAPH_ERROR("no such attribute", IGRAPH_EINVAL);
  }

  atype=VECTOR(al->types)[pos];
  if (atype==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t *data=VECTOR(al->data)[pos];
    igraph_vector_set(data, idx, *(real_t*)(value));
  } else if (atype==IGRAPH_ATTRIBUTE_STR) {
    igraph_strvector_t *data=VECTOR(al->data)[pos];
    IGRAPH_CHECK(igraph_strvector_set(data, idx, (char*)value));
  }
  return 0;
}

/**
 * \ingroup internal
 * \brief Returns an attribute for many elements
 */

int igraph_attribute_list_get_many(const igraph_attribute_list_t *al, 
				   const char *name,
				   const igraph_vector_t *idx, void **value) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);
  igraph_attribute_type_t atype;
  long int i;

  if (pos < 0) {
    IGRAPH_ERROR("no such attribute", IGRAPH_EINVAL);
  }

  atype=VECTOR(al->types)[pos];
  if (atype==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t *data=VECTOR(al->data)[pos];
    igraph_vector_t *nvalue=*value;
    IGRAPH_CHECK(igraph_vector_resize(nvalue, igraph_vector_size(idx)));
    for (i=0; i<igraph_vector_size(idx); i++) {
      VECTOR(*nvalue)[i] = VECTOR(*data)[ (long int)VECTOR(*idx)[i] ];
    }
  } else if (atype==IGRAPH_ATTRIBUTE_STR) {
    igraph_strvector_t *data=VECTOR(al->data)[pos];
    igraph_strvector_t *svalue=*value;
    IGRAPH_CHECK(igraph_strvector_resize(svalue, igraph_vector_size(idx)));
    for (i=0; i<igraph_vector_size(idx); i++) {
      char *str;
      igraph_strvector_get(data, VECTOR(*idx)[i], &str);
      IGRAPH_CHECK(igraph_strvector_set(svalue, i, str));
    }    
  }
  return 0;
}

/**
 * \ingroup internal
 * \brief Sets an attribute for many elements
 */

int igraph_attribute_list_set_many(igraph_attribute_list_t *al, 
				   const char *name,
				   const igraph_vector_t *idx, const void *value) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);
  igraph_attribute_type_t atype;
  long int i;

  if (pos < 0) {
    IGRAPH_ERROR("no such attribute", IGRAPH_EINVAL);
  }

  atype=VECTOR(al->types)[pos];
  if (atype==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t *data=VECTOR(al->data)[pos];
    const igraph_vector_t *nvalue=value;
    long int idxlen=igraph_vector_size(nvalue);
    long int j=0;
    for (i=0; i<igraph_vector_size(idx); i++) {
      VECTOR(*data)[ (long int)VECTOR(*idx)[i] ] = VECTOR(*nvalue)[j++];
      if (j>=idxlen) { j=0; }
    }
  } else if (atype==IGRAPH_ATTRIBUTE_STR) {
    igraph_strvector_t bak;
    igraph_strvector_t *data=VECTOR(al->data)[pos];
    const igraph_strvector_t *svalue=value;
    long int idxlen=igraph_strvector_size(svalue);
    long int j=0;
    IGRAPH_CHECK(igraph_strvector_copy(&bak, data));
    IGRAPH_FINALLY(igraph_strvector_destroy, &bak);
    for (i=0; i<igraph_vector_size(idx); i++) {
      char *str;
      igraph_strvector_get(svalue, j++, &str);
      IGRAPH_CHECK(igraph_strvector_set(&bak, VECTOR(*idx)[i], str));
      if (j>=idxlen) { j=0; }
    }
    IGRAPH_FINALLY_CLEAN(1);
    igraph_strvector_destroy(data);
    *data=bak;
  }

  return 0;
}

/**
 * \ingroup internal
 * \brief Returns an attribute for all elements (untested!)
 */

int igraph_attribute_list_get_all(const igraph_attribute_list_t *al, 
				  const char *name, void **value, 
				  igraph_attribute_type_t *type) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);
  igraph_attribute_type_t atype=IGRAPH_ATTRIBUTE_NUM;

  if (pos < 0) {
    IGRAPH_ERROR("no such attribute", IGRAPH_EINVAL);
  }
  
  if (type != 0) {
    *type=atype;
  }
  
  atype=VECTOR(al->types)[pos];
  if (atype==IGRAPH_ATTRIBUTE_NUM) {
    igraph_vector_t *data=VECTOR(al->data)[pos];
    igraph_vector_t *nvalue=*value;
    igraph_vector_t tmp=*nvalue;
    IGRAPH_CHECK(igraph_vector_copy(nvalue, data));
    igraph_vector_destroy(&tmp);
  } else if (atype==IGRAPH_ATTRIBUTE_STR) {
    igraph_strvector_t *data=VECTOR(al->data)[pos];
    igraph_strvector_t *svalue=*value;
    igraph_strvector_t tmp=*svalue;
    IGRAPH_CHECK(igraph_strvector_copy(svalue, data));
    igraph_strvector_destroy(&tmp);
  }

  return 0;
}

/**
 * \ingroup internal
 * \brief Returns the number of attributes in an attribute list
 */

long int igraph_attribute_list_size(const igraph_attribute_list_t *al) {
  return igraph_vector_size(&al->types);
}

/**
 * \ingroup internal
 * \brief Adds new elements to an attribute list (not attributes, elements!)
 */

int igraph_attribute_list_add_elem(igraph_attribute_list_t *al, long int ne) {
  long int i;
  int ret;
  bool_t error=0;
  igraph_error_handler_t *oldhandler;
  
  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);
  for (i=0; i<igraph_vector_size(&al->types); i++) {
    if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
      igraph_vector_t *data=VECTOR(al->data)[i];
      ret=igraph_vector_resize(data, al->len+ne);
      if (ret != 0) {
	error=1; break;
      }
    } else if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_STR) {
      igraph_strvector_t *data=VECTOR(al->data)[i];
      ret=igraph_strvector_resize(data, al->len+ne);
      if (ret != 0) {
	error=1; break;
      }
    }
  }
  al->len += ne;

  if (error) {
    al->len -= ne;
    for (i=0; i<igraph_vector_size(&al->types); i++) {
      if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
	igraph_vector_t *data=VECTOR(al->data)[i];
	igraph_vector_resize(data, al->len);
      } else if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_STR) {
	igraph_strvector_t *data=VECTOR(al->data)[i];
	igraph_strvector_resize(data, al->len);
      }
    }
  }

  igraph_set_error_handler(oldhandler);
  return 0;
}

/**
 * \ingroup internal
 * \brief Returns names of the attributes in an attribute list
 */

int igraph_attribute_list_names(const igraph_attribute_list_t *al,
				igraph_strvector_t *names, igraph_vector_t *types) {
  if (names != 0) {
    igraph_strvector_t tmp=*names;
    IGRAPH_CHECK(igraph_strvector_copy(names, &al->names));
    igraph_strvector_destroy(&tmp);
  }
  if (types != 0) {
    igraph_vector_t tmp=*types;
    IGRAPH_CHECK(igraph_vector_copy(types, &al->types));
    igraph_vector_destroy(&tmp);
  }
  
  return 0;
}

/**
 * \ingroup internal
 * \brief Creates a (deep) copy of an attribute list
 */

int igraph_attribute_list_copy(igraph_attribute_list_t *to,
			       const igraph_attribute_list_t *from) {
  long int i;
  igraph_error_handler_t *oldhandler;
  bool_t error=0;

  to->len=from->len;
  IGRAPH_CHECK(igraph_strvector_copy(&to->names, &from->names));
  IGRAPH_FINALLY(&to->names, igraph_strvector_destroy);  
  IGRAPH_CHECK(igraph_vector_copy(&to->types, &from->types));
  IGRAPH_FINALLY(&to->types, igraph_vector_destroy);
  IGRAPH_CHECK(igraph_vector_ptr_copy(&to->data, &from->data));
  igraph_vector_ptr_null(&to->data);

  oldhandler=igraph_set_error_handler(igraph_error_handler_ignore);
  for (i=0; i<igraph_vector_size(&from->types); i++) {
    int ret;
    if (VECTOR(from->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
      igraph_vector_t *data=VECTOR(from->data)[i];
      igraph_vector_t *ndata=Calloc(1, igraph_vector_t);
      if (ndata==0) {
	error=1; break;
      }
      ret=igraph_vector_copy(ndata, data);
      if (ret != 0) {
	error=1; break;
      }
      VECTOR(to->data)[i]=ndata;
    } else if (VECTOR(from->types)[i] == IGRAPH_ATTRIBUTE_STR) {
      igraph_strvector_t *data=VECTOR(from->data)[i];
      igraph_strvector_t *ndata=Calloc(1, igraph_strvector_t);
      if (ndata==0) {
	error=1; break;
      }
      ret=igraph_strvector_copy(ndata, data);
      if (ret != 0) {
	error=1; break;
      }
      VECTOR(to->data)[i]=ndata;
    }
  }
  
  if (error != 0) {
    /* names & types are deleted already */
    for (i=0; i<igraph_vector_size(&from->types); i++) {
      if (VECTOR(from->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
	igraph_vector_t *data=VECTOR(to->data)[i];
	if (data != 0) {
	  igraph_vector_destroy(data);
	  Free(data);
	}
      } else if (VECTOR(from->types)[i] == IGRAPH_ATTRIBUTE_STR) {
	igraph_strvector_t *data=VECTOR(to->data)[i];
	if (data != 0) {
	  igraph_strvector_destroy(data);
	  Free(data);
	}
      }
    }
    igraph_vector_ptr_destroy(&to->data);
  }
  
  igraph_set_error_handler(oldhandler);

  if (!error) { IGRAPH_FINALLY_CLEAN(2); }
  return 0;
}

/**
 * \ingroup internal
 * \brief Returns 0 if the given attribute exists and returns
 * its type in the third argument. Returns an error code otherwise.
 */

int igraph_attribute_list_get_type(const igraph_attribute_list_t *al, 
				   const char *name,
				   igraph_attribute_type_t *type) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);

  if (pos < 0) {
    IGRAPH_ERROR("no such attribute", IGRAPH_EINVAL);
  }

  if (type) *type=VECTOR(al->types)[pos];
  return 0;
}

/**
 * \ingroup internal
 * \brief Removes elements (not attributes!) from an attribute list
 * (in a weird way)
 */

void igraph_attribute_list_remove_elem_idx(igraph_attribute_list_t *al, 
					  long int *index, long int nremove) {
  long int i;
  al->len -= nremove;
  
  for (i=0; i<igraph_vector_size(&al->types); i++) {
    if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
      igraph_vector_t *data=VECTOR(al->data)[i];
      igraph_vector_permdelete(data, index, nremove);
    } else if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_STR) {
      igraph_strvector_t *data=VECTOR(al->data)[i];
      igraph_strvector_permdelete(data, index, nremove);
    }
  }
}

/**
 * \ingroup internal
 * \brief Removes elements (not attributes!) from an attribute list
 * (in another weird way)
 */

void igraph_attribute_list_remove_elem_neg(igraph_attribute_list_t *al,
					   const igraph_vector_t *neg, 
					   long int nremove) {
  long int i;
  al->len -= nremove;
  
  for (i=0; i<igraph_vector_size(&al->types); i++) {
    if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_NUM) {
      igraph_vector_t *data=VECTOR(al->data)[i];
      igraph_vector_remove_negidx(data, neg, nremove);
    } else if (VECTOR(al->types)[i] == IGRAPH_ATTRIBUTE_STR) {
      igraph_strvector_t *data=VECTOR(al->data)[i];
      igraph_strvector_remove_negidx(data, neg, nremove);
    }
  }
}

/**
 * \ingroup internal
 * \brief Checks whether the list contains the named attribute
 */

bool_t igraph_attribute_list_has(const igraph_attribute_list_t *al, 
				 const char *name) {
  long int pos=igraph_i_attribute_list_get_pos(al, name);
  return (pos != -1);
}

/* ------------------------------------------------------------------------- */

/**
 * \section about_attributes
 * 
 * <para>Attributes are associated values of graph, vertices or edges. A
 * graph attribute can contain the date of its creations, vertex
 * attributes describe the color to be used for the vertices when the
 * graph is plotted, edge attributes can simply be the weights of the
 * edges in a weighted graph.</para>
 * 
 * <para>The name space for graph, vertex and edge attributes are different,
 * thus the <quote>color</quote> vertex attribute has nothing to do with the
 * <quote>color</quote> edge attribute.</para>
 *
 * <para>In order to use an attribute, first it has to be added by the
 * \ref igraph_add_graph_attribute(), \ref igraph_add_vertex_attribute() or
 * \ref igraph_add_edge_attribute() function. The value(s) of the added
 * attribute is/are undefined.</para>
 * 
 * <para>The value of the attribute can be set after it was added, and also
 * it can be requested, see the documentation of the specific
 * functions.</para>
 * 
 * <para>The attribute can be removed if it is not needed any more.</para>
 * 
 * <para>The attributes are handled properly if vertices and/or edges are
 * added or removed to or from a graph, the attributes are also copied
 * if the graph is copied by \ref igraph_copy().</para>
 * 
 * <para>There are two types of attributes, numeric
 * (<constant>IGRAPH_ATTRIBUTE_NUM</constant>) and  string
 * (<constant>IGRAPH_ATTRIBUTE_STR</constant>). These types cannot be
 * mixed for an 
 * attribute, ie. if the <quote>id</quote> vertex attribute is a
 * string attribute then it is a string attribute for all vertices.</para>
 *
 * <para>Attribute handling does not change the time complexity of any
 * functions in the igraph library.</para>
 */

/**
 * \ingroup attributes
 * \function igraph_add_graph_attribute
 * \brief Adds a graph attribute.
 *
 * Attributes have to be added by calling this function before setting
 * or getting them.
 * \param graph A graph object.
 * \param name The name of the attribute to install.
 * \param type Numeric constant giving the type of the attribute,
 *        either <constant>IGRAPH_ATTRIBUTE_NUM</constant> (numeric) or
 *        <constant>IGRAPH_ATTRIBUTE_STR</constant> (string).
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_EINVAL
 *           invalid (empty) attribute name
 *        \cli IGRAPH_EXISTS
 *           the attribute already exists
 *        \endclist
 *
 * Time complexity: O(1). (Assuming
 * the number of graph attributes of
 * graph is
 * O(1).) 
 */

int igraph_add_graph_attribute(igraph_t *graph, const char *name,
			       igraph_attribute_type_t type) {
  return igraph_attribute_list_add(&graph->gal, name, type);
}

/**
 * \ingroup attributes
 * \function igraph_remove_graph_attribute
 * \brief Removes a graph attribute.
 *
 * \param graph A graph object.
 * \param name The name of the attribute to remove.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1). (Assuming
 * the number of graph attributes of
 * graph is
 * O(1).) 
 */

int igraph_remove_graph_attribute(igraph_t *graph, const char *name) {
  return igraph_attribute_list_remove(&graph->gal, name);
}

/**
 * \ingroup attributes
 * \function igraph_get_graph_attribute
 * \brief Queries the value of a graph attribute.
 *
 * \param graph A graph object.
 * \param name The name of the attribute to query.
 * \param value Pointer to a typeless pointer. The address of the
 *        result will be stored here, a <type>real_t</type> pointer
 *        for numeric attributes or a <type>const char</type> pointer
 *        to string attributes.
 * \param type Pointer to the attribute type, it will be stored here
 *        if not <constant>NULL</constant>.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1). (Assuming
 * the number of graph attributes of
 * graph is
 * O(1).) 
 */

int igraph_get_graph_attribute(const igraph_t *graph, const char *name,
			       void **value, igraph_attribute_type_t *type) {
  return igraph_attribute_list_get(&graph->gal, name, 0, value, type);
}

/**
 * \ingroup attributes
 * \function igraph_set_graph_attribute
 * \brief Sets the value of a graph attribute.
 *
 * \param graph A graph object.
 * \param name The name of the attribute to set.
 * \param value Pointer to the new value of the attribute, either a
 *        <type>real_t</type> or a <type>const char</type> pointer.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1). (Assuming
 * the number of graph attributes of
 * graph is
 * O(1).) 
 */

int igraph_set_graph_attribute(igraph_t *graph, const char *name,
			       const void *value) {
  return igraph_attribute_list_set(&graph->gal, name, 0, value);
}

/**
 * \ingroup attributes
 * \function igraph_list_graph_attributes
 * \brief Queries the list of installed graph attributes.
 *
 * \param graph A graph object.
 * \param names This string vector will contain the names of the
 *        attributes. It should be initialized and will be resized.
 *        This parameter can be <constant>NULL</constant>, which means that it
 *        is ignored.
 * \param types Pointer to a vector, this will be set to the types of
 *        the attributes if the pointer is not <constant>NULL</constant>. The
 *        vector will be resized if neccessary.
 * \return Error code.
 *
 * Time complexity: O(1). (Assuming
 * the number of graph attributes of
 * graph is
 * O(1).) 
 */

int igraph_list_graph_attributes(const igraph_t *graph, 
				 igraph_strvector_t *names,
				 igraph_vector_t *types) {
  return igraph_attribute_list_names(&graph->gal, names, types);
}

/**
 * \ingroup attributes
 * \function igraph_add_vertex_attribute
 * \brief Adds a vertex attribute.
 *
 * \param graph The graph object.
 * \param name The name of the attribute to install.
 * \param type Numeric constant giving the type of the attribute,
 *        either <constant>IGRAPH_ATTRIBUTE_NUM</constant> (numeric) or
 *        <constant>IGRAPH_ATTRIBUTE_STR</constant> (string).
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_EINVAL 
 *          invalid (empty) attribute name
 *        \cli IGRAPH_EXISTS
 *          the attribute already exists
 *        \endclist
 *
 * Time complexity: O(|V|), the
 * number of vertices in the graph.
 */

int igraph_add_vertex_attribute(igraph_t *graph, const char *name,
				igraph_attribute_type_t type) {
  return igraph_attribute_list_add(&graph->val, name, type);
}

/**
 * \ingroup attributes
 * \function igraph_remove_vertex_attribute
 * \brief Removes a vertex attribute.
 *
 * \param graph A graph object.
 * \param name The name of the attribute to remove.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(|V|), assuming
 * that the graph has O(1) vertex
 * attributes. |V| is the number 
 * of vertices.
 */

int igraph_remove_vertex_attribute(igraph_t *graph, const char *name) {
  return igraph_attribute_list_remove(&graph->val, name);
}

/**
 * \ingroup attributes
 * \function igraph_get_vertex_attribute
 * \brief Queries the value of a vertex attribute for a single vertex
 *
 * \param graph The graph object.
 * \param name The name of the vertex attribute.
 * \param v The id of the vertex of which the attribute is requested.
 * \param value Pointer to a typeless pointer. The address of the
 *        result will be stored here, a <type>real_t</type> pointer
 *        for numeric attributes or a <type>const char</type> pointer
 *        to string attributes.
 * \param type If not <constant>NULL</constant> the type of the attribute will
 *        be stored here.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1), assuming
 * that the graph has  O(1) vertex
 * attributes installed. 
 */

int igraph_get_vertex_attribute(const igraph_t *graph, const char *name,
				long int v, void **value,
				igraph_attribute_type_t *type) {
  return igraph_attribute_list_get(&graph->val, name, v, value, type);
}

/**
 * \ingroup attributes
 * \function igraph_set_vertex_attribute
 * \brief Set the value of a vertex attribute for a single vertex.
 *
 * \param graph The graph object.
 * \param name Name of the vertex attribute.
 * \param v The id of the vertex of which the attribute is set.
 * \param value The new value of the attribute.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1), assuming
 * that the graph has O(1) vertex
 * attributes installed. 
 */

int igraph_set_vertex_attribute(igraph_t *graph, const char *name,
				long int v, const void* value) {
  return igraph_attribute_list_set(&graph->val, name, v, value);
}

/**
 * \ingroup attributes
 * \function igraph_get_vertex_attributes
 * \brief Query the value of a vertex attribute for many vertices.
 *
 * \param graph The graph object.
 * \param name The name of the attribute to get.
 * \param v Vector with the vertex ids of the vertices of which the
 *        attribute will be returned.
 * \param value Pointer to a typeless pointer, which contains the
 *        address of either a <type>igraph_vector_t</type> or a
 *        <type>igraph_strvector_t</type> depending on the type of the
 *        attribute.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(|v|), the
 * number of queried vertices, assuming the graph has
 * O(1) vertex attributes.
 */

int igraph_get_vertex_attributes(const igraph_t *graph, const char *name,
				 const igraph_vs_t *v, void **value) {
  int ret;
  igraph_vs_t myv;
  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, v, &myv));
  IGRAPH_FINALLY(igraph_vs_destroy, &myv);
  ret=igraph_attribute_list_get_many(&graph->val, name, 
				     igraph_vs_vector_getvector(graph, &myv), 
				     value);
  igraph_vs_destroy(&myv);
  IGRAPH_FINALLY_CLEAN(1);
  if (v->shorthand) { igraph_vs_destroy((igraph_vs_t*)v); }
  return ret;
}

/**
 * \ingroup attributes
 * \function igraph_set_vertex_attributes
 * \brief Set the value of a vertex attribute for many vertices.
 *
 * \param graph The graph object.
 * \param name The name of the attribute to set.
 * \param v Vector with the vertex ids of the vertices of which the
 *        attribute will be set.
 * \param value The new value(s) of the attribute. This may be of
 *        different length than <parameter>v</parameter>, if it is shorter it
 *        will be recycled (ie. after the last element the first one
 *        is used again), if it is longer the unneeded values are
 *        ignored. Thus it is easy to set an attribute to a single
 *        constant value for many vertices, just give a vector of
 *        length 1 here. This parameter is either a pointer to a
 *        <type>igraph_vector_t</type> or an <type>igraph_strvector_t</type>, 
 *        depending on the type of the attribute.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(|v|), the
 * number of affected vertices, assuming the graph has
 * O(1) vertex attributes.
 */

int igraph_set_vertex_attributes(igraph_t *graph, const char *name,
				 const igraph_vs_t *v, const void *value) {
  int ret;
  igraph_vs_t myv;
  IGRAPH_CHECK(igraph_vs_vectorview_it(graph, v, &myv));
  IGRAPH_FINALLY(igraph_vs_destroy, &myv);
  ret=igraph_attribute_list_set_many(&graph->val, name, 
				     igraph_vs_vector_getvector(graph, &myv), 
				     value);
  igraph_vs_destroy(&myv);
  IGRAPH_FINALLY_CLEAN(1);
  if (v->shorthand) { igraph_vs_destroy((igraph_vs_t*)v); }
  return ret;
}

/**
 * \ingroup attributes
 * \function igraph_list_vertex_attributes
 * \brief Queries the list of installed vertex attributes.
 *
 * \param graph A graph object.
 * \param l This string array will contain the names of the
 *        attributes. It should be initialized and will be resized.
 *        If <constant>NULL</constant> then this parameter is ignored.
 * \param types Pointer to a vector, this will be set to the types of
 *        the attributes if the pointer is not <constant>NULL</constant>. The
 *        vector will be resized if neccessary.
 * \return Error code.
 *
 * Time complexity: O(1). (Assuming
 * the number of vertex attributes of
 * graph is
 * O(1).) 
 */

int igraph_list_vertex_attributes(const igraph_t *graph, igraph_strvector_t *l,
				  igraph_vector_t *types) {
  return igraph_attribute_list_names(&graph->val, l, types);
}

/**
 * \ingroup attributes
 * \function igraph_add_edge_attribute
 * \brief Adds an edge attribute.
 *
 * \param graph The graph object.
 * \param name The name of the attribute to install.
 * \param type Numeric constant giving the type of the attribute,
 *        either <constant>IGRAPH_ATTRIBUTE_NUM</constant> (numeric) or
 *        <constant>IGRAPH_ATTRIBUTE_STR</constant> (string).
 * \return Error code:
 *        \clist
 *        \cli IGRAPH_EINVAL
 *          invalid (empty) attribute name
 *        \cli IGRAPH_EXISTS 
 *          the attribute already exists
 *        \endclist
 *
 * Time complexity: O(|E|), the
 * number of edges in the graph.
 */

int igraph_add_edge_attribute(igraph_t *graph, const char *name,
			      igraph_attribute_type_t type) {
  return igraph_attribute_list_add(&graph->eal, name, type);
}

/**
 * \ingroup attributes
 * \function igraph_remove_edge_attribute
 * \brief Removes an edge attribute.
 *
 * \param graph A graph object.
 * \param name The name of the attribute to remove.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(|E|), assuming
 * that the graph has O(1) edge
 * attributes. |E| is the number 
 * of edges.
 */

int igraph_remove_edge_attribute(igraph_t *graph, const char *name) {
  return igraph_attribute_list_remove(&graph->eal, name);
}

/**
 * \ingroup attributes
 * \function igraph_get_edge_attribute
 * \brief Queries the value of an edge attribute for a single edge
 *
 * It is easy to get the id of an edge by using an edge iterator.
 * \param graph The graph object.
 * \param name The name of the edge attribute.
 * \param e The id of the edge of which the attribute is requested.
 * \param value Pointer to a typeless pointer. The address of the
 *        result will be placed here, a <type>real_t</type> pointer
 *        for numeric attributes and a <type>const char</type> pointer
 *        for string attributes. 
 * \param type If not <constant>NULL</constant> then the type of the attribute
 *        will be stored here.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1), assuming
 * that the graph has  O(1) edge
 * attributes installed. 
 */

int igraph_get_edge_attribute(const igraph_t *graph, const char *name,
			      long int e, void **value,
			      igraph_attribute_type_t *type) {
  return igraph_attribute_list_get(&graph->eal, name, e, value, type);
}

/**
 * \ingroup attributes
 * \function igraph_set_edge_attribute
 * \brief Set the value of an edge attribute for a single edge.
 *
 * \param graph The graph object.
 * \param name Name of the edge attribute.
 * \param e The id of the edge of which the attribute is set.
 * \param value The new value of the attribute. Pointer to a
 *        <type>real_t</type> or a <type>const char</type>.
 * \return Error code:
 *        <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1), assuming
 * that the graph has O(1) edge
 * attributes installed. 
 */

int igraph_set_edge_attribute(igraph_t *graph, const char *name,
			      long int e, const void *value) {
  return igraph_attribute_list_set(&graph->eal, name, e, value);
}

/**
 * \ingroup attributes
 * \function igraph_get_edge_attributes
 * \brief Query the value of an edge attribute for many edges.
 *
 * \param graph The graph object.
 * \param name The name of the attribute to get.
 * \param e Vector with the edge ids of the edges of which the
 *        attribute will be returned.
 * \param value Pointer to a typeless pointer, which contains the
 *        address of either a <type>igraph_vector_t</type> or a
 *        <type>igraph_strvector_t</type> depending on the type of the
 *        attribute.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(|e|), the
 * number of queried edges, assuming the graph has
 * O(1) edge attributes.
 */

int igraph_get_edge_attributes(const igraph_t *graph, const char *name,
			       const igraph_vector_t *e, void **value) {
  return igraph_attribute_list_get_many(&graph->eal, name, e, value);
}

/**
 * \ingroup attributes
 * \function igraph_set_edge_attributes
 * \brief Set the value of an edge attribute for many edges.
 *
 * \param graph The graph object.
 * \param name The name of the attribute to set.
 * \param e Vector with the edge ids of the edges of which the
 *        attribute will be set.
 * \param value The new value(s) of the attribute. This may be of
 *        different length than <parameter>v</parameter>, if it is shorter it
 *        will be recycled (ie. after the last element the first one
 *        is used again), if it is longer the unneeded values are
 *        ignored. Thus it is easy to set an attribute to a single
 *        constant value for many vertices, just give a vector of
 *        length 1 here. This parameter is either a pointer to a
 *        <type>igraph_vector_t</type> or an <type>igraph_strvector_t</type>, 
 *        depending on the type of the attribute.
 * \return Error code:
 *        <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(|v|), the
 * number of affected edges, assuming the graph has
 * O(1) edge attributes.
 */

int igraph_set_edge_attributes(igraph_t *graph, const char *name,
			       const igraph_vector_t *e, const void *value) {
  return igraph_attribute_list_set_many(&graph->eal, name, e, value);
}

/**
 * \ingroup attributes
 * \function igraph_list_edge_attributes
 * \brief Queries the list of installed edge attributes.
 *
 * \param graph A graph object.
 * \param l This string array will contain the names of the
 *        attributes (if not <constant>NULL</constant>). It should be
 *        initialized and will be resized.
 * \param types Pointer to a vector, this will be set to the types of
 *        the attributes if the pointer is not <constant>NULL</constant>. The
 *        vector will be resized if neccessary.
 * \return Error code.
 *
 * Time complexity: O(1). (Assuming
 * the number of edge attributes of
 * graph is
 * O(1).) 
 */

int igraph_list_edge_attributes(const igraph_t *graph, igraph_strvector_t *l,
				igraph_vector_t *types) {
  return igraph_attribute_list_names(&graph->eal, l, types);
}

/**
 * \ingroup attributes
 * \function igraph_get_graph_attribute_type
 * \brief Queries the type of a graph attribute.
 * 
 * \param graph The graph object.
 * \param name The name of the attribute.
 * \param type Pointer to the attribute type.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1) assuming there are
 * O(1) graph attributes.
 */

int igraph_get_graph_attribute_type(const igraph_t *graph, const char *name,
				    igraph_attribute_type_t *type) {
  return igraph_attribute_list_get_type(&graph->gal, name, type);
}

/**
 * \ingroup attributes
 * \function igraph_get_vertex_attribute_type
 * \brief Queries the type of a vertex attribute.
 * 
 * \param graph The graph object.
 * \param name The name of the attribute.
 * \param type Pointer to the attribute type.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1) assuming there are
 * O(1) vertex attributes.
 */

int igraph_get_vertex_attribute_type(const igraph_t *graph, const char *name,
				     igraph_attribute_type_t *type) {
  return igraph_attribute_list_get_type(&graph->val, name, type);
}

/**
 * \ingroup attributes
 * \function igraph_get_edge_attribute_type
 * \brief Queries the type of an edge attribute.
 * 
 * \param graph The graph object.
 * \param name The name of the attribute.
 * \param type Pointer to the attribute type.
 * \return Error code:
 *         <constant>IGRAPH_EINVAL</constant>: the attribute does not exist.
 *
 * Time complexity: O(1) assuming there are
 * O(1) edge attributes.
 */

int igraph_get_edge_attribute_type(const igraph_t *graph, const char *name,
				   igraph_attribute_type_t *type) {
  return igraph_attribute_list_get_type(&graph->eal, name, type);
}

/**
 * \ingroup attributes
 * \function igraph_has_graph_attribute
 * \brief Checks whether a graph has the named graph attribute.
 * 
 * \param graph The graph object.
 * \param name The name of the (potential) attribute.
 * \return Non-zero (TRUE) value if the attribute is installed, zero
 *         (FALSE) otherwise.
 * 
 * Time complexity: O(1) assuming there are
 * O(1) edge attributes.
 */

bool_t igraph_has_graph_attribute(const igraph_t *graph, const char *name) {
  return igraph_attribute_list_has(&graph->gal, name); 
}

/**
 * \ingroup attributes
 * \function igraph_has_vertex_attribute
 * \brief Checks whether a graph has the named vertex attribute.
 * 
 * \param graph The graph object.
 * \param name The name of the (potential) attribute.
 * \return Non-zero (TRUE) value if the attribute is installed, zero
 *         (FALSE) otherwise.
 * 
 * Time complexity: O(1) assuming there are
 * O(1) edge attributes.
 */

bool_t igraph_has_vertex_attribute(const igraph_t *graph, const char *name) {
  return igraph_attribute_list_has(&graph->val, name); 
}

/**
 * \ingroup attributes
 * \function igraph_has_edge_attribute
 * \brief Checks whether a graph has the named edge attribute.
 * 
 * \param graph The graph object.
 * \param name The name of the (potential) attribute.
 * \return Non-zero (TRUE) value if the attribute is installed, zero
 *         (FALSE) otherwise.
 * 
 * Time complexity: O(1) assuming there are
 * O(1) edge attributes.
 */

bool_t igraph_has_edge_attribute(const igraph_t *graph, const char *name) {
  return igraph_attribute_list_has(&graph->eal, name); 
}
