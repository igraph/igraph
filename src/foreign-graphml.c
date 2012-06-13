/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
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

#include "igraph_foreign.h"
#include "config.h"
#include <math.h>               /* isnan */
#include "igraph_math.h"
#include "igraph_attributes.h"
#include "igraph_interface.h"
#include "igraph_types_internal.h"

#include <ctype.h>		/* isspace */
#include <string.h>
#include "igraph_memory.h"
#include <stdarg.h> 		/* va_start & co */

#if HAVE_LIBXML == 1
#include <libxml/encoding.h>
#include <libxml/parser.h>


xmlEntity blankEntityStruct = {
#ifndef XML_WITHOUT_CORBA
  0,
#endif
  XML_ENTITY_DECL,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  XML_EXTERNAL_GENERAL_PARSED_ENTITY,
  0,
  0,
  0,
  0
};

xmlEntityPtr blankEntity = &blankEntityStruct;

/* TODO: proper error handling */

typedef struct igraph_i_graphml_attribute_record_t {
  const char *id;         	/* GraphML id */
  enum { I_GRAPHML_BOOLEAN, I_GRAPHML_INTEGER, I_GRAPHML_LONG,
	 I_GRAPHML_FLOAT, I_GRAPHML_DOUBLE, I_GRAPHML_STRING,
	 I_GRAPHML_UNKNOWN_TYPE } type;	/* GraphML type */
  igraph_attribute_record_t record;
} igraph_i_graphml_attribute_record_t;

struct igraph_i_graphml_parser_state {
  enum { START, INSIDE_GRAPHML, INSIDE_GRAPH, INSIDE_NODE, INSIDE_EDGE,
      INSIDE_KEY, INSIDE_DEFAULT, INSIDE_DATA, FINISH, UNKNOWN, ERROR } st;
  igraph_t *g;
  igraph_trie_t node_trie;
  igraph_strvector_t edgeids;
  igraph_vector_t edgelist;
  unsigned int prev_state;
  unsigned int unknown_depth;
  int index;
  igraph_bool_t successful, edges_directed, destroyed;
  igraph_trie_t v_names;
  igraph_vector_ptr_t v_attrs;
  igraph_trie_t e_names;
  igraph_vector_ptr_t e_attrs;
  igraph_trie_t g_names;
  igraph_vector_ptr_t g_attrs;
  xmlChar *data_key;
  igraph_attribute_elemtype_t data_type;
  char *error_message;
  char *data_char;
};

void igraph_i_graphml_destroy_state(struct igraph_i_graphml_parser_state* state) {
  long int i;

  if (state->destroyed) return;
  state->destroyed=1;

  /* this is the easy part */
  igraph_trie_destroy(&state->node_trie);
  igraph_strvector_destroy(&state->edgeids);
  igraph_trie_destroy(&state->v_names);
  igraph_trie_destroy(&state->e_names);
  igraph_trie_destroy(&state->g_names);
  igraph_vector_destroy(&state->edgelist);
   
  if (state->error_message) { free(state->error_message); }
  if (state->data_key) { free(state->data_key); }
  if (state->data_char) { free(state->data_char); }
  
  for (i=0; i<igraph_vector_ptr_size(&state->v_attrs); i++) {
    igraph_i_graphml_attribute_record_t *rec=VECTOR(state->v_attrs)[i];
    if (rec->record.type==IGRAPH_ATTRIBUTE_NUMERIC) {
      if (rec->record.value != 0) {
	igraph_vector_destroy((igraph_vector_t*)rec->record.value);
	igraph_Free(rec->record.value);
      }
    } else if (rec->record.type==IGRAPH_ATTRIBUTE_STRING) {
      if (rec->record.value != 0) {
	igraph_strvector_destroy((igraph_strvector_t*)rec->record.value);
	igraph_Free(rec->record.value);
      }
    }
    if (rec->id != 0) igraph_Free(rec->id);
    if (rec->record.name != 0) igraph_Free(rec->record.name);
    igraph_Free(rec);
  }	 

  for (i=0; i<igraph_vector_ptr_size(&state->e_attrs); i++) {
    igraph_i_graphml_attribute_record_t *rec=VECTOR(state->e_attrs)[i];
    if (rec->record.type==IGRAPH_ATTRIBUTE_NUMERIC) {
      if (rec->record.value != 0) {
	igraph_vector_destroy((igraph_vector_t*)rec->record.value);
	igraph_Free(rec->record.value);
      }
    } else if (rec->record.type==IGRAPH_ATTRIBUTE_STRING) {
      if (rec->record.value != 0) {
	igraph_strvector_destroy((igraph_strvector_t*)rec->record.value);
	igraph_Free(rec->record.value);
      }
    }
    if (rec->id != 0) igraph_Free(rec->id);
    if (rec->record.name != 0) igraph_Free(rec->record.name);
    igraph_Free(rec);
  }

  for (i=0; i<igraph_vector_ptr_size(&state->g_attrs); i++) {
    igraph_i_graphml_attribute_record_t *rec=VECTOR(state->g_attrs)[i];
    if (rec->record.type==IGRAPH_ATTRIBUTE_NUMERIC) {
      if (rec->record.value != 0) {
	igraph_vector_destroy((igraph_vector_t*)rec->record.value);
	igraph_Free(rec->record.value);
      }
    } else if (rec->record.type==IGRAPH_ATTRIBUTE_STRING) {
      if (rec->record.value != 0) {
	igraph_strvector_destroy((igraph_strvector_t*)rec->record.value);
	igraph_Free(rec->record.value);
      }
    }
    if (rec->id != 0) igraph_Free(rec->id);
    if (rec->record.name != 0) igraph_Free(rec->record.name);
    igraph_Free(rec);
  }

  igraph_vector_ptr_destroy(&state->v_attrs);
  igraph_vector_ptr_destroy(&state->e_attrs);
  igraph_vector_ptr_destroy(&state->g_attrs);
  
  IGRAPH_FINALLY_CLEAN(1);
}

void igraph_i_graphml_sax_handler_error(void *state0, const char* msg, ...) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  va_list ap;
  
  va_start(ap, msg);
  
  if (state->error_message == 0)
    state->error_message=igraph_Calloc(4096, char);
   
  state->successful=0;
  state->st=ERROR;
  vsnprintf(state->error_message, 4096, msg, ap);
   
  va_end(ap);
}

xmlEntityPtr igraph_i_graphml_sax_handler_get_entity(void *state0,
						     const xmlChar* name) {
  xmlEntityPtr predef = xmlGetPredefinedEntity(name);
  if (predef != NULL) return predef;
  IGRAPH_WARNING("unknown XML entity found\n");
  return blankEntity;
}

void igraph_i_graphml_handle_unknown_start_tag(struct igraph_i_graphml_parser_state *state) {
  if (state->st != UNKNOWN) {
    state->prev_state=state->st;
    state->st=UNKNOWN;
    state->unknown_depth=1;
  } else state->unknown_depth++;
}

void igraph_i_graphml_sax_handler_start_document(void *state0) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  int ret;
  
  state->st=START;
  state->successful=1;
  state->edges_directed=0;
  state->destroyed=0;
  state->data_key=0;
  state->error_message=0;
  state->data_char=0;
  
  ret=igraph_vector_ptr_init(&state->v_attrs, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &state->v_attrs);
  ret=igraph_vector_ptr_init(&state->e_attrs, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &state->e_attrs);
  ret=igraph_vector_ptr_init(&state->g_attrs, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &state->g_attrs);
  ret=igraph_vector_init(&state->edgelist, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_vector_destroy, &state->edgelist);
  ret=igraph_trie_init(&state->node_trie, 1);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_trie_destroy, &state->node_trie);
  ret=igraph_strvector_init(&state->edgeids, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_strvector_destroy, &state->edgeids);
  ret=igraph_trie_init(&state->v_names, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_trie_destroy, &state->v_names);
  ret=igraph_trie_init(&state->e_names, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_trie_destroy, &state->e_names);
  ret=igraph_trie_init(&state->g_names, 0);
  if (ret) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_trie_destroy, &state->g_names);
  
  IGRAPH_FINALLY_CLEAN(9);
  IGRAPH_FINALLY(igraph_i_graphml_destroy_state, state);
}

void igraph_i_graphml_sax_handler_end_document(void *state0) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  long i, l;
  int r;
  igraph_attribute_record_t idrec, eidrec;
  const char *idstr="id";
  igraph_bool_t already_has_vertex_id=0, already_has_edge_id=0;

  if (!state->successful) return;

  if (state->index<0) {

    igraph_vector_ptr_t vattr, eattr, gattr;
    long int esize=igraph_vector_ptr_size(&state->e_attrs);
    const void **tmp;
    r=igraph_vector_ptr_init(&vattr, 
			     igraph_vector_ptr_size(&state->v_attrs)+1);
    if (r) {
      igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, r);
      igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
      return;
    }
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &vattr);
    if (igraph_strvector_size(&state->edgeids) != 0) {
      esize++;      
    }
    r=igraph_vector_ptr_init(&eattr, esize);
    if (r) {
      igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, r);
      igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
      return;
    }
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &eattr);
    r=igraph_vector_ptr_init(&gattr, igraph_vector_ptr_size(&state->g_attrs));
    if (r) {
      igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, r);
      igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
      return;
    }
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &gattr);

    for (i=0; i<igraph_vector_ptr_size(&state->v_attrs); i++) {
      igraph_i_graphml_attribute_record_t *graphmlrec=
	VECTOR(state->v_attrs)[i];
      igraph_attribute_record_t *rec=&graphmlrec->record;

      /* Check that the name of the vertex attribute is not 'id'.
	 If it is then we cannot the complimentary 'id' attribute. */
      if (! strcmp(rec->name, idstr)) {
	already_has_vertex_id=1;
      }

      if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_vector_t *vec=(igraph_vector_t*)rec->value;
	long int origsize=igraph_vector_size(vec);
	long int nodes=igraph_trie_size(&state->node_trie);
	igraph_vector_resize(vec, nodes);
	for (l=origsize; l<nodes; l++) {
	  VECTOR(*vec)[l]=IGRAPH_NAN;
	}
      } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
	igraph_strvector_t *strvec=(igraph_strvector_t*)rec->value;
	long int origsize=igraph_strvector_size(strvec);
	long int nodes=igraph_trie_size(&state->node_trie);
	igraph_strvector_resize(strvec, nodes);
	for (l=origsize; l<nodes; l++) {
	  igraph_strvector_set(strvec, l, "");
	}
      }
      VECTOR(vattr)[i]=rec;
    }
    if (!already_has_vertex_id) {
      idrec.name=idstr;
      idrec.type=IGRAPH_ATTRIBUTE_STRING;
      tmp=&idrec.value;
      igraph_trie_getkeys(&state->node_trie, (const igraph_strvector_t **)tmp);
      VECTOR(vattr)[i]=&idrec;
    } else {
      igraph_vector_ptr_pop_back(&vattr);
      IGRAPH_WARNING("Could not add vertex ids, "
		     "there is already an 'id' vertex attribute");
    }

    for (i=0; i<igraph_vector_ptr_size(&state->e_attrs); i++) {
      igraph_i_graphml_attribute_record_t *graphmlrec=
	VECTOR(state->e_attrs)[i];
      igraph_attribute_record_t *rec=&graphmlrec->record;

      if (! strcmp(rec->name, idstr)) {
	already_has_edge_id=1;
      }

      if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_vector_t *vec=(igraph_vector_t*)rec->value;
	long int origsize=igraph_vector_size(vec);
	long int edges=igraph_vector_size(&state->edgelist)/2;
	igraph_vector_resize(vec, edges);
	for (l=origsize; l<edges; l++) {
	  VECTOR(*vec)[l]=IGRAPH_NAN;
	}
      } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
	igraph_strvector_t *strvec=(igraph_strvector_t*)rec->value;
	long int origsize=igraph_strvector_size(strvec);
	long int edges=igraph_vector_size(&state->edgelist)/2;
	igraph_strvector_resize(strvec, edges);
	for (l=origsize; l<edges; l++) {
	  igraph_strvector_set(strvec, l, "");
	}
      }
      VECTOR(eattr)[i]=rec;
    }
    if (igraph_strvector_size(&state->edgeids) != 0) {
      if (!already_has_edge_id) {
	long int origsize=igraph_strvector_size(&state->edgeids);
	eidrec.name=idstr;
	eidrec.type=IGRAPH_ATTRIBUTE_STRING;
	igraph_strvector_resize(&state->edgeids, 
				igraph_vector_size(&state->edgelist)/2);
	for (; origsize < igraph_strvector_size(&state->edgeids); origsize++) {
	  igraph_strvector_set(&state->edgeids, origsize, "");
	}
	eidrec.value=&state->edgeids;
	VECTOR(eattr)[(long int)igraph_vector_ptr_size(&eattr)-1]=&eidrec;
      } else {
	igraph_vector_ptr_pop_back(&eattr);
	IGRAPH_WARNING("Could not add edge ids, "
		       "there is already an 'id' edge attribute");
      }
    }

    for (i=0; i<igraph_vector_ptr_size(&state->g_attrs); i++) {
      igraph_i_graphml_attribute_record_t *graphmlrec=
	VECTOR(state->g_attrs)[i];
      igraph_attribute_record_t *rec=&graphmlrec->record;
      if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_vector_t *vec=(igraph_vector_t*)rec->value;
	long int origsize=igraph_vector_size(vec);
	igraph_vector_resize(vec, 1);
	for (l=origsize; l<1; l++) {
	  VECTOR(*vec)[l]=IGRAPH_NAN;
	}
      } else if (rec->type == IGRAPH_ATTRIBUTE_STRING) {
	igraph_strvector_t *strvec=(igraph_strvector_t*)rec->value;
	long int origsize=igraph_strvector_size(strvec);
	igraph_strvector_resize(strvec, 1);
	for (l=origsize; l<1; l++) {
	  igraph_strvector_set(strvec, l, "");
	}
      }
      VECTOR(gattr)[i]=rec;
    }
    
    igraph_empty_attrs(state->g, 0, state->edges_directed, &gattr);
    igraph_add_vertices(state->g, igraph_trie_size(&state->node_trie),
			&vattr);
    igraph_add_edges(state->g, &state->edgelist, &eattr);

    igraph_vector_ptr_destroy(&vattr);
    igraph_vector_ptr_destroy(&eattr);
    igraph_vector_ptr_destroy(&gattr);
    IGRAPH_FINALLY_CLEAN(3);     
  }

  igraph_i_graphml_destroy_state(state);
}

#define toXmlChar(a)   (BAD_CAST(a))
#define fromXmlChar(a) ((char *)(a)) /* not the most elegant way... */

void igraph_i_graphml_add_attribute_key(const xmlChar** attrs, 
					struct igraph_i_graphml_parser_state *state) {
  xmlChar **it;
  igraph_trie_t *trie=0;
  igraph_vector_ptr_t *ptrvector=0;
  long int id;
  int ret;
  igraph_i_graphml_attribute_record_t *rec=
    igraph_Calloc(1, igraph_i_graphml_attribute_record_t);

  if (!state->successful) return;
   
   if (rec==0) { 
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, 
		 IGRAPH_ENOMEM);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  IGRAPH_FINALLY(igraph_free, rec);
  for (it=(xmlChar**)attrs; *it; it+=2) {
    if (xmlStrEqual(*it, toXmlChar("id"))) {
      const char *id=(const char*)(*(it+1));
      rec->id=strdup(id);
    } else if (xmlStrEqual(*it, toXmlChar("attr.name"))) {
      const char *name=fromXmlChar(*(it+1));
      rec->record.name=strdup(name);
    } else if (xmlStrEqual(*it, toXmlChar("attr.type"))) {
      if (xmlStrEqual(*(it+1), (xmlChar*)"boolean")) { 
	rec->type=I_GRAPHML_BOOLEAN;
	rec->record.type=IGRAPH_ATTRIBUTE_NUMERIC;	    
      } else if (xmlStrEqual(*(it+1), toXmlChar("string"))) {
	rec->type=I_GRAPHML_STRING;
	rec->record.type=IGRAPH_ATTRIBUTE_STRING;
      } else if (xmlStrEqual(*(it+1), toXmlChar("float"))) { 
	rec->type=I_GRAPHML_FLOAT;
	rec->record.type=IGRAPH_ATTRIBUTE_NUMERIC;
      } else if (xmlStrEqual(*(it+1), toXmlChar("double"))) { 
	rec->type=I_GRAPHML_DOUBLE;
	rec->record.type=IGRAPH_ATTRIBUTE_NUMERIC;
      } else if (xmlStrEqual(*(it+1), toXmlChar("int"))) {
	rec->type=I_GRAPHML_INTEGER;
	rec->record.type=IGRAPH_ATTRIBUTE_NUMERIC;
      } else if (xmlStrEqual(*(it+1), toXmlChar("long"))) {
	rec->type=I_GRAPHML_LONG;
	rec->record.type=IGRAPH_ATTRIBUTE_NUMERIC;
      } else {
	igraph_error("Cannot parse GraphML file, unknown attribute type", 
		     __FILE__, __LINE__, IGRAPH_PARSEERROR);
        igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file, unknown attribute type");
        return;
      }
    } else if (xmlStrEqual(*it, toXmlChar("for"))) {
      /* graph, vertex or edge attribute? */
      if (xmlStrEqual(*(it+1), toXmlChar("graph"))) { 
	trie=&state->g_names;
	ptrvector=&state->g_attrs;
      } else if (xmlStrEqual(*(it+1), toXmlChar("node"))) {
	trie=&state->v_names;
	ptrvector=&state->v_attrs;
      } else if (xmlStrEqual(*(it+1), toXmlChar("edge"))) {
	trie=&state->e_names;
	ptrvector=&state->e_attrs;
      } else {
	igraph_error("Cannot parse GraphML file, unknown attribute type",
		     __FILE__, __LINE__, IGRAPH_PARSEERROR);
        igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file, unknown attribute type");
        return;
      }
    }
  }

  /* in case of a missing attr.name attribute, use the id as the attribute name */
  if (rec->record.name == 0) {
    rec->record.name=strdup(rec->id);
  }

  if (trie == 0 && state->successful) {
    igraph_error("Cannot parse GraphML file, missing 'for' attribute", __FILE__, __LINE__, IGRAPH_PARSEERROR);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file, missing 'for' attribute");
    return;
  }
	
  /* add to trie, attribues */
  igraph_trie_get(trie, rec->id, &id);
  if (id != igraph_trie_size(trie)-1) {
    igraph_error("Cannot parse GraphML file, duplicate attribute", 
		 __FILE__, __LINE__, IGRAPH_PARSEERROR);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file, duplicate attribute");
    return;
  }
  ret=igraph_vector_ptr_push_back(ptrvector, rec);
  if (ret) {
    igraph_error("Cannot read GraphML file", __FILE__, __LINE__, ret);
    igraph_i_graphml_sax_handler_error(state, "Cannot read GraphML file");
    return;
  }

  /* create the attribute values */
  switch (rec->record.type) {
    igraph_vector_t *vec;
    igraph_strvector_t *strvec;
  case IGRAPH_ATTRIBUTE_NUMERIC:
    vec=igraph_Calloc(1, igraph_vector_t);
    if (vec==0) {
      igraph_error("Cannot parse GraphML file", __FILE__, __LINE__,
		   IGRAPH_ENOMEM);
      igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
      return;
    }
    rec->record.value=vec;
    igraph_vector_init(vec, 0);    
    break;
  case IGRAPH_ATTRIBUTE_STRING:
    strvec=igraph_Calloc(1, igraph_strvector_t);
    if (strvec==0) {
      igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, 
		   IGRAPH_ENOMEM);
      igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
      return;
    }
    rec->record.value=strvec;
    igraph_strvector_init(strvec, 0);
    break;
  default: break;
  }

  IGRAPH_FINALLY_CLEAN(1);	/* rec */
}

void igraph_i_graphml_attribute_data_setup(struct igraph_i_graphml_parser_state *state,
					   const xmlChar **attrs,
					   igraph_attribute_elemtype_t type) {
  xmlChar **it;
  for (it=(xmlChar**)attrs; *it; it+=2) {
    if (xmlStrEqual(*it, toXmlChar("key"))) {
      if (state->data_key) {
	free(state->data_key);
      }
      state->data_key=xmlStrdup(*(it+1));
      if (state->data_char) {
        free(state->data_char);
      }
      state->data_char=0;
      state->data_type=type;
    } else {
      /* ignore */
    }
  }
}

void igraph_i_graphml_attribute_data_add(struct igraph_i_graphml_parser_state *state,
					 const xmlChar *data, int len) {
  long int data_char_new_start=0;

  if (!state->successful) return;

  if (state->data_char) {
    data_char_new_start=strlen(state->data_char);
    state->data_char=igraph_Realloc(state->data_char, data_char_new_start+len+1, char);
  } else {
    state->data_char=igraph_Calloc(len+1, char);
  }
  if (state->data_char==0) {
    igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, 
		 IGRAPH_ENOMEM);
    igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
    return;
  }
  memcpy(state->data_char+data_char_new_start, data, len*sizeof(xmlChar));
  state->data_char[data_char_new_start+len]='\0';
}

void igraph_i_graphml_attribute_data_finish(struct igraph_i_graphml_parser_state *state) {
  const char *key=fromXmlChar(state->data_key);
  igraph_attribute_elemtype_t type=state->data_type;
  igraph_trie_t *trie=0;
  igraph_vector_ptr_t *ptrvector=0;
  igraph_i_graphml_attribute_record_t *graphmlrec;
  igraph_attribute_record_t *rec;
  long int recid, id=0;
  int ret;
  
  switch (type) {
  case IGRAPH_ATTRIBUTE_GRAPH:
    trie=&state->g_names;
    ptrvector=&state->g_attrs;
    id=0;
    break;
  case IGRAPH_ATTRIBUTE_VERTEX:
    trie=&state->v_names;
    ptrvector=&state->v_attrs;
    id=igraph_trie_size(&state->node_trie)-1; /* hack */
    break;
  case IGRAPH_ATTRIBUTE_EDGE:
    trie=&state->e_names;
    ptrvector=&state->e_attrs;
    id=igraph_vector_size(&state->edgelist)/2-1; /* hack */
    break;
  default:
    /* impossible */
    break;
  }
  
  igraph_trie_check(trie, key, &recid);
  if (recid < 0) {
    /* no such attribute key, issue a warning */
    IGRAPH_WARNING("unknown attribute key in GraphML file, ignoring attribute");
    igraph_Free(state->data_char);
    return;
  }
   
  graphmlrec=VECTOR(*ptrvector)[recid];
  rec=&graphmlrec->record;

  switch (rec->type) {
    igraph_vector_t *vec;
    igraph_strvector_t *strvec;
    igraph_real_t num;
    long int s, i;
  case IGRAPH_ATTRIBUTE_NUMERIC:
    vec=(igraph_vector_t *)rec->value;
    s=igraph_vector_size(vec);
    if (id >= s) {
      ret=igraph_vector_resize(vec, id+1);
      if (ret) {
	igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
        igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
        return;
      }
      for (i=s; i<id; i++) {
	VECTOR(*vec)[i]=IGRAPH_NAN;
      }
    }
    if (state->data_char)
      sscanf(state->data_char, "%lf", &num);
    else
      num=0;
    VECTOR(*vec)[id]=num;
    break;
  case IGRAPH_ATTRIBUTE_STRING:
    strvec=(igraph_strvector_t *)rec->value;
    s=igraph_strvector_size(strvec);
    if (id >= s) {
      ret=igraph_strvector_resize(strvec, id+1);
      if (ret) {
	igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
        igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
        return;
      }
      for (i=s;i<id;i++) {
	igraph_strvector_set(strvec, i, "");
      }
    }
    if (state->data_char)
      ret=igraph_strvector_set(strvec, id, (char*)state->data_char);
    else
      ret=igraph_strvector_set(strvec, id, "");
    if (ret) {
      igraph_error("Cannot parse GraphML file", __FILE__, __LINE__, ret);
      igraph_i_graphml_sax_handler_error(state, "Cannot parse GraphML file");
      return;
    }
    break;
  default:
    break;
  }

  if (state->data_char) igraph_Free(state->data_char);
}

void igraph_i_graphml_sax_handler_start_element(void *state0,
						const xmlChar* name,
						const xmlChar** attrs) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  xmlChar** it;
  long int id1, id2;

  switch (state->st) {
  case START:
    /* If we are in the START state and received a graphml tag,
     * change to INSIDE_GRAPHML state. Otherwise, change to UNKNOWN. */
    if (xmlStrEqual(name, toXmlChar("graphml")))
      state->st=INSIDE_GRAPHML;
    else
      igraph_i_graphml_handle_unknown_start_tag(state);
    break;
    
  case INSIDE_GRAPHML:
    /* If we are in the INSIDE_GRAPHML state and received a graph tag,
     * change to INSIDE_GRAPH state if the state->index counter reached
     * zero (this is to handle multiple graphs in the same file).
     * Otherwise, change to UNKNOWN. */
    if (xmlStrEqual(name, toXmlChar("graph"))) {
      if (state->index==0) {
	state->st=INSIDE_GRAPH;
	for (it=(xmlChar**)attrs; *it; it+=2) {
	  if (xmlStrEqual(*it, toXmlChar("edgedefault"))) {
	    if (xmlStrEqual(*(it+1), toXmlChar("directed"))) state->edges_directed=1;
	    else if (xmlStrEqual(*(it+1), toXmlChar("undirected"))) state->edges_directed=0;
	  }
	}
      }
      state->index--;
    } else if (xmlStrEqual(name, toXmlChar("key"))) {
      igraph_i_graphml_add_attribute_key(attrs, state);
      state->st=INSIDE_KEY;
    } else
      igraph_i_graphml_handle_unknown_start_tag(state);
    break;

  case INSIDE_KEY:
    /* If we are in the INSIDE_KEY state, check for default tag */
    if (xmlStrEqual(name, toXmlChar("default"))) state->st=INSIDE_DEFAULT;
    else igraph_i_graphml_handle_unknown_start_tag(state);
    break;

  case INSIDE_DEFAULT:
    /* If we are in the INSIDE_DEFAULT state, every further tag will be unknown */
    igraph_i_graphml_handle_unknown_start_tag(state);
    break;
    
  case INSIDE_GRAPH:
    /* If we are in the INSIDE_GRAPH state, check for node and edge tags */
    if (xmlStrEqual(name, toXmlChar("edge"))) {
      id1=-1; id2=-1; 
      for (it=(xmlChar**)attrs; *it; it+=2) {
	if (xmlStrEqual(*it, toXmlChar("source"))) {
	  igraph_trie_get(&state->node_trie, fromXmlChar(*(it+1)), &id1);
	}
	if (xmlStrEqual(*it, toXmlChar("target"))) {
	  igraph_trie_get(&state->node_trie, fromXmlChar(*(it+1)), &id2);
	}
	if (xmlStrEqual(*it, toXmlChar("id"))) {
	  long int edges=igraph_vector_size(&state->edgelist)/2+1;
	  long int origsize=igraph_strvector_size(&state->edgeids);
	  igraph_strvector_resize(&state->edgeids, edges);
	  for (;origsize < edges-1; origsize++) {
	    igraph_strvector_set(&state->edgeids, origsize, "");
	  }
	  igraph_strvector_set(&state->edgeids, edges-1, fromXmlChar(*(it+1)));
	}
      }
      if (id1>=0 && id2>=0) {
	igraph_vector_push_back(&state->edgelist, id1);
	igraph_vector_push_back(&state->edgelist, id2);
      } else {
	igraph_i_graphml_sax_handler_error(state, "Edge with missing source or target encountered");
	return;
      }
      state->st=INSIDE_EDGE;
    } else if (xmlStrEqual(name, toXmlChar("node"))) {
      for (it=(xmlChar**)attrs; *it; it+=2) {
	if (xmlStrEqual(*it, toXmlChar("id"))) {
	  it++;
	  igraph_trie_get(&state->node_trie, fromXmlChar(*it), &id1);
	  break;
	}
      }
      state->st=INSIDE_NODE;
    } else if (xmlStrEqual(name, toXmlChar("data"))) {
      igraph_i_graphml_attribute_data_setup(state, attrs, IGRAPH_ATTRIBUTE_GRAPH);
      state->prev_state=state->st;
      state->st=INSIDE_DATA;
    } else
      igraph_i_graphml_handle_unknown_start_tag(state);
    break;
    
  case INSIDE_NODE:
    if (xmlStrEqual(name, toXmlChar("data"))) {
      igraph_i_graphml_attribute_data_setup(state, attrs,
					    IGRAPH_ATTRIBUTE_VERTEX);
      state->prev_state=state->st;
      state->st=INSIDE_DATA;
    }
    break;
    
  case INSIDE_EDGE:
    if (xmlStrEqual(name, toXmlChar("data"))) {
      igraph_i_graphml_attribute_data_setup(state, attrs, 
					    IGRAPH_ATTRIBUTE_EDGE);
      state->prev_state=state->st;
      state->st=INSIDE_DATA;
    }
    break;
    
  default:
    break;
  }
}

void igraph_i_graphml_sax_handler_end_element(void *state0,
						const xmlChar* name) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  
  switch (state->st) {
  case INSIDE_GRAPHML:
    state->st=FINISH;
    break;
    
  case INSIDE_GRAPH:
    state->st=INSIDE_GRAPHML;
    break;
    
  case INSIDE_KEY:
    state->st=INSIDE_GRAPHML;
    break;

  case INSIDE_DEFAULT:
    state->st=INSIDE_KEY;
    break;
    
  case INSIDE_NODE:
    state->st=INSIDE_GRAPH;
    break;
    
  case INSIDE_EDGE:
    state->st=INSIDE_GRAPH;
    break;

  case INSIDE_DATA:
    igraph_i_graphml_attribute_data_finish(state);
    state->st=state->prev_state;
    break;
    
  case UNKNOWN:
    state->unknown_depth--;
    if (!state->unknown_depth) state->st=state->prev_state;
    break;
    
  default:
    break;
  }
}

void igraph_i_graphml_sax_handler_chars(void* state0, const xmlChar* ch, int len) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  
  switch (state->st) {
  case INSIDE_KEY:
  case INSIDE_DEFAULT:
    break;
    
  case INSIDE_DATA:
    igraph_i_graphml_attribute_data_add(state, ch, len);
    break;
    
  default:
    /* just ignore it */
    break;
  }
}

static xmlSAXHandler igraph_i_graphml_sax_handler={
  NULL, NULL, NULL, NULL, NULL,
    igraph_i_graphml_sax_handler_get_entity,
    NULL, NULL, NULL, NULL, NULL, NULL,
    igraph_i_graphml_sax_handler_start_document,
    igraph_i_graphml_sax_handler_end_document,
    igraph_i_graphml_sax_handler_start_element,
    igraph_i_graphml_sax_handler_end_element,
    NULL,
    igraph_i_graphml_sax_handler_chars,
    NULL, NULL, NULL,
    igraph_i_graphml_sax_handler_error,
    igraph_i_graphml_sax_handler_error,
    igraph_i_graphml_sax_handler_error,
};

#endif

int igraph_i_xml_escape(char* src, char** dest) {
  long int destlen=0;
  char *s, *d;
  for (s=src; *s; s++, destlen++) {
    if (*s == '&') destlen += 4;
    else if (*s == '<') destlen += 3;
    else if (*s == '>') destlen += 3;
    else if (*s == '"') destlen += 5;
    else if (*s == '\'') destlen += 5;
  }
  *dest=igraph_Calloc(destlen+1, char);
  if (!*dest) IGRAPH_ERROR("Not enough memory", IGRAPH_ENOMEM);
  for (s=src, d=*dest; *s; s++, d++) {
    switch (*s) {
    case '&':
      strcpy(d, "&amp;"); d+=4; break;
    case '<':
      strcpy(d, "&lt;"); d+=3; break;
    case '>':
      strcpy(d, "&gt;"); d+=3; break;
    case '"':
      strcpy(d, "&quot;"); d+=5; break;
    case '\'':
      strcpy(d, "&apos;"); d+=5; break;
    default:
      *d = *s;
    }
  }
  *d=0;
  return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_read_graph_graphml
 * \brief Reads a graph from a GraphML file.
 * 
 * </para><para>
 * GraphML is an XML-based file format for representing various types of
 * graphs. Currently only the most basic import functionality is implemented
 * in igraph: it can read GraphML files without nested graphs and hyperedges.
 * Attributes of the graph are loaded only if an attribute interface
 * is attached, ie. if you use igraph from R or Python.
 *
 * </para><para>
 * Graph attribute names are taken from the \c attr.name attributes of the
 * \c key tags in the GraphML file. Since \c attr.name is not mandatory,
 * igraph will fall back to the \c id attribute of the \c key tag if
 * \c attr.name is missing.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param instream A stream, it should be readable.
 * \param index If the GraphML file contains more than one graph, the one
 *              specified by this index will be loaded. Indices start from
 *              zero, so supply zero here if your GraphML file contains only
 *              a single graph.
 * 
 * \return Error code:
 *         \c IGRAPH_PARSEERROR: if there is a
 *         problem reading the file, or the file is syntactically
 *         incorrect.
 *         \c IGRAPH_UNIMPLEMENTED: the GraphML functionality was disabled
 *         at compile-time
 * 
 * \example examples/simple/graphml.c
 */
int igraph_read_graph_graphml(igraph_t *graph, FILE *instream,
			      int index) {

#if HAVE_LIBXML == 1
  xmlParserCtxtPtr ctxt;
  struct igraph_i_graphml_parser_state state;
  int res;
  char buffer[4096];

  if (index<0)
    IGRAPH_ERROR("Graph index must be non-negative", IGRAPH_EINVAL);

  xmlInitParser();

  /* Create a progressive parser context */
  state.g=graph;
  state.index=index<0?0:index;
  res=fread(buffer, 1, 4096, instream);
  ctxt=xmlCreatePushParserCtxt(&igraph_i_graphml_sax_handler,
			       &state,
			       buffer,
			       res,
			       NULL);
/*   ctxt=xmlCreateIOParserCtxt(&igraph_i_graphml_sax_handler, &state, */
/* 			     igraph_i_libxml2_read_callback, */
/* 			     igraph_i_libxml2_close_callback, */
/* 			     instream, XML_CHAR_ENCODING_NONE); */
  if (ctxt==NULL)
    IGRAPH_ERROR("Can't create progressive parser context", IGRAPH_PARSEERROR);

  /* Parse the file */
  while ((res=fread(buffer, 1, 4096, instream))>0) {
    xmlParseChunk(ctxt, buffer, res, 0);
    if (!state.successful) break;
  }
  xmlParseChunk(ctxt, buffer, res, 1);
  
  /* Free the context */
  xmlFreeParserCtxt(ctxt);
  if (!state.successful) {
    if (state.error_message != 0)
      IGRAPH_ERROR(state.error_message, IGRAPH_PARSEERROR);
    else
      IGRAPH_ERROR("Malformed GraphML file", IGRAPH_PARSEERROR);
  }
  if (state.index>=0)
    IGRAPH_ERROR("Graph index was too large", IGRAPH_EINVAL);
  
  return 0;
#else
  IGRAPH_ERROR("GraphML support is disabled", IGRAPH_UNIMPLEMENTED);
#endif
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_graphml
 * \brief Writes the graph to a file in GraphML format
 *
 * </para><para>
 * GraphML is an XML-based file format for representing various types of
 * graphs. See the GraphML Primer (http://graphml.graphdrawing.org/primer/graphml-primer.html)
 * for detailed format description.
 * 
 * \param graph The graph to write. 
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error
 *         writing the file. 
 *
 * Time complexity: O(|V|+|E|) otherwise. All
 * file operations are expected to have time complexity 
 * O(1). 
 * 
 * \example examples/simple/graphml.c
 */
int igraph_write_graph_graphml(const igraph_t *graph, FILE *outstream) {
  int ret;
  igraph_integer_t l, vc;
  igraph_eit_t it;
  igraph_strvector_t gnames, vnames, enames;
  igraph_vector_t gtypes, vtypes, etypes;
  long int i;
  igraph_vector_t numv;
  igraph_strvector_t strv;
  
  ret=fprintf(outstream, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  ret=fprintf(outstream, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  ret=fprintf(outstream, "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  ret=fprintf(outstream, "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  ret=fprintf(outstream, "         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  ret=fprintf(outstream, "<!-- Created by igraph -->\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);

  /* dump the <key> elements if any */

  IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
  IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);
  
  IGRAPH_STRVECTOR_INIT_FINALLY(&gnames, 0);
  IGRAPH_STRVECTOR_INIT_FINALLY(&vnames, 0);
  IGRAPH_STRVECTOR_INIT_FINALLY(&enames, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&gtypes, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vtypes, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&etypes, 0);
  igraph_i_attribute_get_info(graph,
			      &gnames, &gtypes,
			      &vnames, &vtypes,
			      &enames, &etypes);
  
  /* graph attributes */
  for (i=0; i<igraph_vector_size(&gtypes); i++) {
    char *name, *name_escaped;
    igraph_strvector_get(&gnames, i, &name);
    IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
    if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
      ret=fprintf(outstream, "  <key id=\"g_%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"string\"/>\n", name_escaped, name_escaped);
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);      
    } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
      ret=fprintf(outstream, "  <key id=\"g_%s\" for=\"graph\" attr.name=\"%s\" attr.type=\"double\"/>\n", name_escaped, name_escaped);
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);            
    }
    igraph_Free(name_escaped);
  }

  /* vertex attributes */
  for (i=0; i<igraph_vector_size(&vtypes); i++) {
    char *name, *name_escaped;
    igraph_strvector_get(&vnames, i, &name);
    IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
    if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
      ret=fprintf(outstream, "  <key id=\"v_%s\" for=\"node\" attr.name=\"%s\" attr.type=\"string\"/>\n", name_escaped, name_escaped);
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);      
    } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
      ret=fprintf(outstream, "  <key id=\"v_%s\" for=\"node\" attr.name=\"%s\" attr.type=\"double\"/>\n", name_escaped, name_escaped);
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);            
    }
    igraph_Free(name_escaped);
  }

  /* edge attributes */
  for (i=0; i<igraph_vector_size(&etypes); i++) {
    char *name, *name_escaped;
    igraph_strvector_get(&enames, i, &name);
    IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
    if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
      ret=fprintf(outstream, "  <key id=\"e_%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"string\"/>\n", name_escaped, name_escaped);
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
      ret=fprintf(outstream, "  <key id=\"e_%s\" for=\"edge\" attr.name=\"%s\" attr.type=\"double\"/>\n", name_escaped, name_escaped);
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
    igraph_Free(name_escaped);
  }

  ret=fprintf(outstream, "  <graph id=\"G\" edgedefault=\"%s\">\n", (igraph_is_directed(graph)?"directed":"undirected"));
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);

  /* Write the graph atributes before anything else */
  
  for (i=0; i<igraph_vector_size(&gtypes); i++) {
    char *name, *name_escaped;
    if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_strvector_get(&gnames, i, &name);
      IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
      if (!isnan(VECTOR(numv)[0])) {
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        ret=fprintf(outstream, "    <data key=\"g_%s\">%g</data>\n",
                    name_escaped, VECTOR(numv)[0]);
        igraph_Free(name_escaped);
        if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
    } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
      char *s, *s_escaped;
      igraph_strvector_get(&gnames, i, &name);
      IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
      ret=fprintf(outstream, "    <data key=\"g_%s\">", name_escaped);
      igraph_Free(name_escaped);
      IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
      igraph_strvector_get(&strv, 0, &s);
      IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
      ret=fprintf(outstream, "%s", s_escaped);
      igraph_Free(s_escaped);
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      ret=fprintf(outstream, "</data>\n");
      if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    }
  }
    
  /* Let's dump the nodes first */
  vc=igraph_vcount(graph);
  for (l=0; l<vc; l++) {
    char *name, *name_escaped;
    ret=fprintf(outstream, "    <node id=\"n%ld\">\n", (long)l);

    if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    
    for (i=0; i<igraph_vector_size(&vtypes); i++) {
      if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
        igraph_strvector_get(&vnames, i, &name);
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name,
                     igraph_vss_1(l), &numv));
        if (!isnan(VECTOR(numv)[0])) {
          IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
          ret=fprintf(outstream, "      <data key=\"v_%s\">%g</data>\n",
                      name_escaped, VECTOR(numv)[0]);
          igraph_Free(name_escaped);
          if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        }
      } else if (VECTOR(vtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
        char *s, *s_escaped;
        igraph_strvector_get(&vnames, i, &name);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        ret=fprintf(outstream, "      <data key=\"v_%s\">", name_escaped);
        igraph_Free(name_escaped);
        IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name,
                     igraph_vss_1(l), &strv));
        igraph_strvector_get(&strv, 0, &s);
        IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
        ret=fprintf(outstream, "%s", s_escaped);
        igraph_Free(s_escaped);
        if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        ret=fprintf(outstream, "</data>\n");
        if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
    }

    ret=fprintf(outstream, "    </node>\n");
    if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);    
  }
  
  /* Now the edges */
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &it));
  IGRAPH_FINALLY(igraph_eit_destroy, &it);
  while (!IGRAPH_EIT_END(it)) {
    igraph_integer_t from, to;
    char *name, *name_escaped;
    long int edge=IGRAPH_EIT_GET(it);
    igraph_edge(graph, edge, &from, &to);
    ret=fprintf(outstream, "    <edge source=\"n%ld\" target=\"n%ld\">\n", 
		(long int)from, (long int)to);
    if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);

    for (i=0; i<igraph_vector_size(&etypes); i++) {
      if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
        igraph_strvector_get(&enames, i, &name);
        IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name,
                     igraph_ess_1(edge), &numv));
        if (!isnan(VECTOR(numv)[0])) {
          IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
          ret=fprintf(outstream, "      <data key=\"e_%s\">%g</data>\n",
                      name_escaped, VECTOR(numv)[0]);
          igraph_Free(name_escaped);
          if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        }
      } else if (VECTOR(etypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
        char *s, *s_escaped;
        igraph_strvector_get(&enames, i, &name);
        IGRAPH_CHECK(igraph_i_xml_escape(name, &name_escaped));
        ret=fprintf(outstream, "      <data key=\"e_%s\">", name_escaped);
        igraph_Free(name_escaped);
        IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name,
                     igraph_ess_1(edge), &strv));
        igraph_strvector_get(&strv, 0, &s);
        IGRAPH_CHECK(igraph_i_xml_escape(s, &s_escaped));
        ret=fprintf(outstream, "%s", s_escaped);
        igraph_Free(s_escaped);
        if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
        ret=fprintf(outstream, "</data>\n");
        if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
    }

    ret=fprintf(outstream, "    </edge>\n");
    if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    IGRAPH_EIT_NEXT(it);
  }
  igraph_eit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  
  ret=fprintf(outstream, "  </graph>\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  fprintf(outstream, "</graphml>\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  
  igraph_strvector_destroy(&gnames);
  igraph_strvector_destroy(&vnames);
  igraph_strvector_destroy(&enames);
  igraph_vector_destroy(&gtypes);
  igraph_vector_destroy(&vtypes);
  igraph_vector_destroy(&etypes);
  igraph_vector_destroy(&numv);
  igraph_strvector_destroy(&strv);
  IGRAPH_FINALLY_CLEAN(8);

  return 0;
}
