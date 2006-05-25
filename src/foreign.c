/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2005 Gabor Csardi <csardi@rmki.kfki.hu>
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
#include <config.h>

#include <ctype.h>		/* isspace */
#include "memory.h"

#ifdef HAVE_LIBXML
#include <libxml/encoding.h>
#include <libxml/parser.h>
#endif

/**
 * \section about_loadsave 
 * 
 * <para>These functions can write a graph to a file, or read a graph
 * from a file.</para>
 * 
 * <para>Note that as \a igraph uses the traditional C streams, it is
 * possible to read/write files from/to memory, at least on GNU
 * operating systems supporting \quote non-standard\endquote streams.</para>
 */

/**
 * \ingroup loadsave
 * \function igraph_read_graph_edgelist
 * \brief Reads an edge list from a file and creates a graph.
 * 
 * This format is simply a series of even number integers separated by
 * whitespace. The one edge (ie. two integers) per line format is thus
 * not required (but recommended for readability). Edges of directed
 * graphs are assumed to be in from, to order.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream Pointer to a stream, it should be readable.
 * \param n The number of vertices in the graph. If smaller than the
 *        largest integer in the file it will be ignored. It is thus
 *        safe to supply zero here.
 * \param directed Logical, if true the graph is directed, if false it
 *        will be undirected.
 * \return Error code:
 *         \c IGRAPH_PARSEERROR: if there is a
 *         problem reading the file, or the file is syntactically
 *         incorrect. 
 * 
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges. It is assumed that
 * reading an integer requires O(1)
 * time. 
 */

int igraph_read_graph_edgelist(igraph_t *graph, FILE *instream, 
			       igraph_integer_t n, igraph_bool_t directed) {

  igraph_vector_t edges=IGRAPH_VECTOR_NULL;
  long int from, to;
  int c;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_CHECK(igraph_vector_reserve(&edges, 100));

  /* skip all whitespace */
  do {
    c = getc (instream);
  } while (isspace (c));
  ungetc (c, instream);
  
  while (!feof(instream)) {
    int read;
    read=fscanf(instream, "%li", &from);
    if (read != 1) { 
      IGRAPH_ERROR("parsing edgelist file failed", IGRAPH_PARSEERROR); 
    }
    read=fscanf(instream, "%li", &to);
    if (read != 1) { 
      IGRAPH_ERROR("parsing edgelist file failed", IGRAPH_PARSEERROR); 
    }
    IGRAPH_CHECK(igraph_vector_push_back(&edges, from));
    IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
    
    /* skip all whitespace */
    do {
      c = getc (instream);
    } while (isspace (c));
    ungetc (c, instream);    
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, n, directed));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

extern int igraph_ncol_yyparse();
extern FILE *igraph_ncol_yyin;
igraph_vector_t *igraph_ncol_vector=0;
igraph_vector_t *igraph_ncol_weights=0;
igraph_trie_t *igraph_ncol_trie=0;

/**
 * \ingroup loadsave
 * \function igraph_read_graph_ncol
 * \brief Reads a <code>.ncol</code> file used by LGL, also
 * useful for creating graphs from \quote named\endquote (and
 * optionally weighted) edge lists. 
 * 
 * This format is used by the Large Graph Layout program
 * (http://bioinformatics.icmb.utexas.edu/lgl/), and it is simply a
 * symbolic weighted edge list. It is a simple text file with one edge
 * per line. An edge is defined by two symbolic vertex names separated
 * by whitespace. (The symbolic vertex names themselves cannot contain 
 * whitespace. They might follow by an optional number, this will be
 * the weight of the edge; the number can be negative and can be in
 * scientific notation. If there is no weight specified to an edge it
 * is assumed to be zero.
 *
 * The resulting graph is always undirected.
 * LGL cannot deal with files which contain multiple or loop edges, 
 * this is however not checked here, as \a igraph is happy with
 * these.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream Pointer to a stream, it should be readable.
 * \param predefnames Pointer to the symbolic names of the vertices in
 *        the file. If \c NULL is given here then vertex ids will be
 *        assigned to vertex names in the order of their appearence in
 *        the \c .ncol file. If it is not \c NULL and some unknown
 *        vertex names are found in the \c .ncol file then new vertex
 *        ids will be assigned to them. 
 * \param names Logical value, if TRUE the symbolic names of the
 *        vertices will be added to the graph as a vertex attribute
 *        called \quote name\endquote.
 * \param weights Logical value, if TRUE the weights of the
 *        edges is added to the graph as an edge attribute called
 *        \quote weight\endquote.
 * \param directed Whether to create a directed graph. As this format
 *        was originally used only for undirected graphs there is no
 *        information in the file about the directedness of the graph.
 *        Set this parameter to \c IGRAPH_DIRECTED or \c
 *        IGRAPH_UNDIRECTED to create a directed or undirected graph.
 * \return Error code:
 *         \c IGRAPH_PARSEERROR: if there is a
 *          problem reading 
 *         the file, or the file is syntactically incorrect.
 *
 * Time complexity:
 * O(|V|+|E|log(|V|)) if we neglect
 * the time required by the parsing. As usual
 * |V| is the number of vertices,
 * while |E| is the number of edges. 
 * 
 * \sa \ref igraph_read_graph_lgl(), \ref igraph_write_graph_ncol()
 */

int igraph_read_graph_ncol(igraph_t *graph, FILE *instream, 
			   igraph_strvector_t *predefnames,
			   igraph_bool_t names, igraph_bool_t weights, igraph_bool_t directed) {
  
  igraph_vector_t edges, ws;
  igraph_trie_t trie=IGRAPH_TRIE_NULL;
  long int no_predefined=0;
  igraph_vector_ptr_t name, weight;
  igraph_vector_ptr_t *pname=0, *pweight=0;
  igraph_i_attribute_record_t namerec, weightrec;
  const char *namestr="name", *weightstr="weight";

  IGRAPH_CHECK(igraph_empty(graph, 0, directed));
  IGRAPH_FINALLY(igraph_destroy, graph);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  IGRAPH_TRIE_INIT_FINALLY(&trie, names);
  IGRAPH_VECTOR_INIT_FINALLY(&ws, 0);

  /* Add the predefined names, if any */
  if (predefnames != 0) {
    long int i, id, n;
    char *key;
    n=no_predefined=igraph_strvector_size(predefnames);
    for (i=0; i<n; i++) {
      igraph_strvector_get(predefnames, i, &key);
      igraph_trie_get(&trie, key, &id);
      if (id != i) {
	IGRAPH_WARNING("reading NCOL file, duplicate entry in predefnames");
	no_predefined--;
      }
    }
  }
  
  igraph_ncol_vector=&edges;
  igraph_ncol_weights=&ws;
  igraph_ncol_trie=&trie;
  igraph_ncol_yyin=instream;

  igraph_ncol_yyparse();

  if (predefnames != 0 && 
      igraph_trie_size(&trie) != no_predefined) {
    IGRAPH_WARNING("unknown vertex/vertices found, predefnames extended");    
  }

  if (names) {
    igraph_strvector_t *namevec;
    IGRAPH_CHECK(igraph_vector_ptr_init(&name, 1)); 
    pname=&name;
    igraph_trie_getkeys(&trie, &namevec); /* dirty */
    namerec.name=namestr;
    namerec.type=IGRAPH_ATTRIBUTE_STRING;
    namerec.value=namevec;
    VECTOR(name)[0]=&namerec;
  }

  if (weights) {
    IGRAPH_CHECK(igraph_vector_ptr_init(&weight, 1)); 
    pweight=&weight;
    weightrec.name=weightstr;
    weightrec.type=IGRAPH_ATTRIBUTE_NUMERIC;
    weightrec.value=&ws;
    VECTOR(weight)[0]=&weightrec;
  }

  IGRAPH_CHECK(igraph_add_vertices(graph, igraph_vector_max(&edges)+1, pname));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, pweight));

  if (pname) {
    igraph_vector_ptr_destroy(pname);
  }
  if (pweight) {
    igraph_vector_ptr_destroy(pweight);
  }
  igraph_vector_destroy(&ws); 
  igraph_trie_destroy(&trie);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}

extern int igraph_lgl_yyparse();
extern FILE *igraph_lgl_yyin;
igraph_vector_t *igraph_lgl_vector=0;
igraph_vector_t *igraph_lgl_weights=0;
igraph_trie_t *igraph_lgl_trie=0;

/**
 * \ingroup loadsave
 * \function igraph_read_graph_lgl
 * \brief Reads a graph from an <code>.lgl</code> file
 * 
 * The <code>.lgl</code> format is used by the Large Graph
 * Layout visualization software
 * (http://bioinformatics.icmb.utexas.edu/lgl/), it can 
 * describe undirected optionally weighted graphs. From the LGL
 * manual: 
 * 
 * \blockquote <para>The second format is the LGL file format
 * (<code>.lgl</code> file 
 * suffix). This is yet another graph file format that tries to be as
 * stingy as possible with space, yet keeping the edge file in a human
 * readable (not binary) format. The format itself is like the
 * following:
 * \verbatim # vertex1name
vertex2name [optionalWeight]
vertex3name [optionalWeight] \endverbatim
 * Here, the first vertex of an edge is preceded with a pound sign
 * '#'.  Then each vertex that shares an edge with that vertex is
 * listed one per line on subsequent lines.</para> \endblockquote
 * 
 * LGL cannot handle loop and multiple edges or directed graphs, but
 * in \a igraph it is not an error to have multiple and loop edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream A stream, it should be readable.
 * \param names Logical value, if TRUE the symbolic names of the
 *        vertices will be added to the graph as a vertex attribute
 *        called \quote name\endquote.
 * \param weights Logical value, if TRUE the weights of the
 *        edges is added to the graph as an edge attribute called
 *        \quote weight\endquote.
 * \return Error code:
 *         \c IGRAPH_PARSEERROR: if there is a
 *         problem reading the file, or the file is syntactically
 *         incorrect. 
 *
 * Time complexity:
 * O(|V|+|E|log(|V|)) if we neglect
 * the time required by the parsing. As usual
 * |V| is the number of vertices,
 * while |E| is the number of edges. 
 * 
 * \sa \ref igraph_read_graph_ncol(), \ref igraph_write_graph_lgl()
 */

int igraph_read_graph_lgl(igraph_t *graph, FILE *instream,
			  igraph_bool_t names, igraph_bool_t weights) {

  igraph_vector_t edges=IGRAPH_VECTOR_NULL, ws=IGRAPH_VECTOR_NULL;
  igraph_trie_t trie=IGRAPH_TRIE_NULL;
  igraph_vector_ptr_t name, weight;
  igraph_vector_ptr_t *pname=0, *pweight=0;
  igraph_i_attribute_record_t namerec, weightrec;
  const char *namestr="name", *weightstr="weight";
  
  IGRAPH_VECTOR_INIT_FINALLY(&ws, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  IGRAPH_TRIE_INIT_FINALLY(&trie, names);
  
  igraph_lgl_vector=&edges;
  igraph_lgl_weights=&ws;
  igraph_lgl_trie=&trie;
  igraph_lgl_yyin=instream;

  igraph_lgl_yyparse();
  
  IGRAPH_CHECK(igraph_empty(graph, 0, IGRAPH_UNDIRECTED));
  IGRAPH_FINALLY(igraph_destroy, graph);

  if (names) {
    igraph_strvector_t *namevec;
    IGRAPH_CHECK(igraph_vector_ptr_init(&name, 1)); 
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &name);
    pname=&name;
    igraph_trie_getkeys(&trie, &namevec); /* dirty */
    namerec.name=namestr;
    namerec.type=IGRAPH_ATTRIBUTE_STRING;
    namerec.value=namevec;
    VECTOR(name)[0]=&namerec;
  }

  if (weights) {
    IGRAPH_CHECK(igraph_vector_ptr_init(&weight, 1)); 
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &weight);
    pweight=&weight;
    weightrec.name=weightstr;
    weightrec.type=IGRAPH_ATTRIBUTE_NUMERIC;
    weightrec.value=&ws;
    VECTOR(weight)[0]=&weightrec;
  }

  IGRAPH_CHECK(igraph_add_vertices(graph, igraph_vector_max(&edges)+1, pname));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, pweight));
  
  if (pweight) {
    igraph_vector_ptr_destroy(pweight);
    IGRAPH_FINALLY_CLEAN(1);
  }
  if (pname) {
    igraph_vector_ptr_destroy(pname);
    IGRAPH_FINALLY_CLEAN(1);
  }
  igraph_trie_destroy(&trie);
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&ws);
  IGRAPH_FINALLY_CLEAN(4);
  
  return 0;
}

extern int igraph_pajek_yyparse();
extern FILE *igraph_pajek_yyin;
igraph_vector_t *igraph_pajek_vector=0;
igraph_stack_t *igraph_pajek_coords=0;
igraph_bool_t igraph_pajek_directed;
long int igraph_pajek_vcount=0;
long int igraph_pajek_actfrom, igraph_pajek_actto;
int igraph_pajek_mode=0;	/* 0 - general, 1 - vertex, 2 - edge */

int igraph_read_graph_pajek(igraph_t *graph, FILE *instream) {
  igraph_vector_t edges;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  igraph_pajek_vector=&edges;
  igraph_pajek_yyin=instream;
  
  igraph_pajek_yyparse();
  
  IGRAPH_CHECK(igraph_create(graph, &edges, igraph_pajek_vcount, 
			     igraph_pajek_directed));
  igraph_vector_destroy(&edges);

  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

#ifdef HAVE_LIBXML
int igraph_i_libxml2_read_callback(void *instream, char* buffer, int len) {
  int res;  
  res=fread(buffer, 1, len, (FILE*)instream);
  if (res) return res;
  if (feof((FILE*)instream)) return 0;
  return -1;
}

int igraph_i_libxml2_close_callback(void *instream) { return 0; }

struct igraph_i_graphml_key_data {
  enum { I_GRAPHML_GRAPH, I_GRAPHML_NODE, I_GRAPHML_EDGE,
      I_GRAPHML_UNKNOWN_TARGET } what_for;
  enum { I_GRAPHML_BOOLEAN, I_GRAPHML_INTEGER, I_GRAPHML_LONG,
      I_GRAPHML_FLOAT, I_GRAPHML_DOUBLE, I_GRAPHML_STRING,
      I_GRAPHML_UNKNOWN_TYPE } type;
  char* name;
  char* default_value;
  long id;
};

struct igraph_i_graphml_parser_state {
  enum { START, INSIDE_GRAPHML, INSIDE_GRAPH, INSIDE_NODE, INSIDE_EDGE,
      INSIDE_KEY, INSIDE_DEFAULT, INSIDE_DATA, FINISH, UNKNOWN, ERROR } st;
  igraph_t *g;
  igraph_trie_t node_trie;
  igraph_trie_t key_trie;
  igraph_vector_t edgelist;
  igraph_vector_ptr_t key_list;
  struct igraph_i_graphml_key_data *current_key;
  unsigned int prev_state;
  unsigned int unknown_depth;
  int index;
  igraph_bool_t successful, edges_directed, directed_default;
};

void igraph_i_graphml_sax_handler_error(void *state0, const char* msg, ...) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  state->successful=0;
  state->st=ERROR;
  /* TODO: use the message */
}

xmlEntityPtr igraph_i_graphml_sax_handler_get_entity(void *state0,
						     const xmlChar* name) {
  return xmlGetPredefinedEntity(name);
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
  
  state->st=START;
  state->successful=1;
  state->edges_directed=0;
  state->current_key=NULL;
  igraph_vector_init(&state->edgelist, 0);
  igraph_vector_ptr_init(&state->key_list, 0);
  igraph_trie_init(&state->node_trie, 1);
  /*igraph_trie_init(&state->key_trie, 1);*/
}

void igraph_i_graphml_sax_handler_end_document(void *state0) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  struct igraph_i_graphml_key_data *kd;
  long i, l;

  if (state->index<0) {
    igraph_empty(state->g, igraph_trie_size(&state->node_trie),
		 state->directed_default);
    igraph_add_edges(state->g, &state->edgelist, 0);
  }

  l=igraph_vector_ptr_size(&state->key_list);
  for (i=0; i<l; i++) {
    kd=(struct igraph_i_graphml_key_data*)VECTOR(state->key_list)[i];
    if (kd->name) Free(kd->name);
    if (kd->default_value) Free(kd->default_value);
  }
  
  igraph_trie_destroy(&state->node_trie);
  /*igraph_trie_destroy(&state->key_trie);*/
  igraph_vector_ptr_destroy_all(&state->key_list);
  igraph_vector_destroy(&state->edgelist);
}

void igraph_i_graphml_sax_handler_start_element(void *state0,
						const xmlChar* name,
						const xmlChar** attrs) {
  struct igraph_i_graphml_parser_state *state=
    (struct igraph_i_graphml_parser_state*)state0;
  xmlChar** it;
  long int id1, id2;
  unsigned short int directed;

  switch (state->st) {
  case START:
    /* If we are in the START state and received a graphml tag,
     * change to INSIDE_GRAPHML state. Otherwise, change to UNKNOWN. */
    if (!strcmp(name, "graphml"))
      state->st=INSIDE_GRAPHML;
    else
      igraph_i_graphml_handle_unknown_start_tag(state);
    break;
    
  case INSIDE_GRAPHML:
    /* If we are in the INSIDE_GRAPHML state and received a graph tag,
     * change to INSIDE_GRAPH state if the state->index counter reached
     * zero (this is to handle multiple graphs in the same file).
     * Otherwise, change to UNKNOWN. */
    if (!strcmp(name, "graph")) {
      if (state->index==0) {
	state->st=INSIDE_GRAPH;
	for (it=(xmlChar**)attrs; *it; it+=2) {
	  if (!strcmp(*it, "edgedefault")) {
	    if (!strcmp(*(it+1), "directed")) state->edges_directed=1;
	    else if (!strcmp(*(it+1), "undirected")) state->edges_directed=0;
	  }
	}
      }
      state->index--;
    } else if (!strcmp(name, "key")) {
      if (!state->current_key)
	state->current_key=Calloc(1, struct igraph_i_graphml_key_data);
      state->current_key->name=NULL;
      state->current_key->id=-1;
      state->current_key->default_value=NULL;
      state->current_key->what_for=I_GRAPHML_UNKNOWN_TARGET;
      state->current_key->type=I_GRAPHML_UNKNOWN_TYPE;
      for (it=(xmlChar**)attrs; *it; it+=2) {
	if (!strcmp(*it, "id")) {
	  /*igraph_trie_get(&state->key_trie, (char*)*(it+1), &state->current_key->id);*/
	} else if (!strcmp(*it, "attr.name")) {
	  state->current_key->name=strdup(*(it+1));
	} else if (!strcmp(*it, "attr.type")) {
	  if (!strcmp(*(it+1), "boolean")) state->current_key->type=I_GRAPHML_BOOLEAN;
	  else if (!strcmp(*(it+1), "string")) state->current_key->type=I_GRAPHML_STRING;
	  else if (!strcmp(*(it+1), "float")) state->current_key->type=I_GRAPHML_FLOAT;
	  else if (!strcmp(*(it+1), "double")) state->current_key->type=I_GRAPHML_DOUBLE;
	  else if (!strcmp(*(it+1), "int")) state->current_key->type=I_GRAPHML_INTEGER;
	  else if (!strcmp(*(it+1), "long")) state->current_key->type=I_GRAPHML_LONG;
	  else {
	    // some warning here?
	  }
	} else if (!strcmp(*it, "for")) {
	  if (!strcmp(*(it+1), "graph")) state->current_key->what_for=I_GRAPHML_GRAPH;
	  else if (!strcmp(*(it+1), "node")) state->current_key->what_for=I_GRAPHML_NODE;
	  else if (!strcmp(*(it+1), "edge")) state->current_key->what_for=I_GRAPHML_EDGE;
	  else {
	    // some warning here?
	  }
	}
      }
      state->st=INSIDE_KEY;
    } else
      igraph_i_graphml_handle_unknown_start_tag(state);
    break;

  case INSIDE_KEY:
    /* If we are in the INSIDE_KEY state, check for default tag */
    if (!strcmp(name, "default")) state->st=INSIDE_DEFAULT;
    else igraph_i_graphml_handle_unknown_start_tag(state);
    break;

  case INSIDE_DEFAULT:
    /* If we are in the INSIDE_DEFAULT state, every further tag will be unknown */
    igraph_i_graphml_handle_unknown_start_tag(state);
    break;
    
  case INSIDE_GRAPH:
    /* If we are in the INSIDE_GRAPH state, check for node and edge tags */
    if (!strcmp(name, "edge")) {
      id1=-1; id2=-1; directed=state->edges_directed;
      for (it=(xmlChar**)attrs; *it; it+=2) {
	if (!strcmp(*it, "source")) {
	  /* TODO: convert from xmlChar to char instead of just typecasting */
	  igraph_trie_get(&state->node_trie, (char*)*(it+1), &id1);
	}
	if (!strcmp(*it, "target")) {
	  igraph_trie_get(&state->node_trie, (char*)*(it+1), &id2);
	}
	if (!strcmp(*it, "directed")) {
	  if (!strcmp((char*)*(it+1), "true")) directed=1; else directed=0;
	}
      }
      if (id1>=0 && id2>=0) {
	igraph_vector_push_back(&state->edgelist, id1);
	igraph_vector_push_back(&state->edgelist, id2);
	if (!directed && state->directed_default) {
	  igraph_vector_push_back(&state->edgelist, id2);
	  igraph_vector_push_back(&state->edgelist, id1);
	}
      } else {
	igraph_i_graphml_sax_handler_error(state, "Edge with missing source or target encountered");
	return;
      }
      state->st=INSIDE_EDGE;
    } else if (!strcmp(name, "node")) {
      for (it=(xmlChar**)attrs; *it; it+=2) {
	if (!strcmp(*it, "id")) {
	  it++;
	  igraph_trie_get(&state->node_trie, (char*)*it, &id1);
	  break;
	}
      }
      state->st=INSIDE_NODE;
    } else if (!strcmp(name, "data")) {
      state->prev_state=state->st;
      state->st=INSIDE_DATA;
    } else
      igraph_i_graphml_handle_unknown_start_tag(state);
    break;
    
  case INSIDE_NODE:
    /* We don't expect any tags inside node */
    igraph_i_graphml_handle_unknown_start_tag(state);
    break;
    
  case INSIDE_EDGE:
    /* We don't expect any tags inside edge */
    igraph_i_graphml_handle_unknown_start_tag(state);
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
    if (!state->current_key->name ||
	(state->current_key->what_for == I_GRAPHML_UNKNOWN_TARGET) ||
	(state->current_key->type == I_GRAPHML_UNKNOWN_TYPE)) {
      /* some warning here? */
    } else {
      /* printf("Encountered new key: %s\n", state->current_key->name);
      printf("Target: %d\n", state->current_key->what_for);
      printf("Type: %d\n", state->current_key->type);
      printf("Default value: %s\n", state->current_key->default_value); */
      igraph_vector_ptr_push_back(&state->key_list, state->current_key);
      state->current_key=NULL;
    }
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
    if (state->current_key->default_value) {
      // continuation of character data, so append to the previous
      state->current_key->default_value=realloc(state->current_key->default_value,
						strlen(state->current_key->default_value)+len+1);
      if (state->current_key->default_value)
	strncat(state->current_key->default_value, (const char*)ch, (size_t)len);
      // TODO: this is not really efficient, since it reallocs every time
      // Maybe it's time to implement a string type in igraph?
    } else 
      state->current_key->default_value=strdup((const char*)ch);
    break;
    
  case INSIDE_DATA:
    break;
    
  default:
    // just ignore it
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

/**
 * \ingroup loadsave
 * \function igraph_read_graph_graphml
 * \brief Reads a graph from a GraphML file.
 * 
 * GraphML is an XML-based file format for representing various types of
 * graphs. Currently only the most basic import functionality is implemented
 * in igraph: it can read GraphML files without nested graphs and hyperedges.
 * Attributes of the graph are not loaded yet.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream A stream, it should be readable.
 * \param directed Whether the imported graph should be directed. Please
 *              note that if you ask for a directed graph, but the
 *              GraphML file to be read contains an undirected graph,
 *              the resulting \c igraph_t graph will be undirected as
 *              well, but it will appear as directed. If you ask for an
 *              undirected graph, the result will be undirected even if
 *              the original GraphML file contained a directed graph.
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
 */
int igraph_read_graph_graphml(igraph_t *graph, FILE *instream,
			      igraph_bool_t directed, int index) {
#ifdef HAVE_LIBXML
  xmlParserCtxtPtr ctxt;
  struct igraph_i_graphml_parser_state state;
  int res;
  char buffer[4096];

  if (index<0)
    IGRAPH_ERROR("Graph index must be non-negative", IGRAPH_EINVAL);
  
  // Create a progressive parser context
  state.g=graph;
  state.index=index<0?0:index;
  state.directed_default=directed;
  ctxt=xmlCreateIOParserCtxt(&igraph_i_graphml_sax_handler, &state,
			     igraph_i_libxml2_read_callback,
			     igraph_i_libxml2_close_callback,
			     instream, XML_CHAR_ENCODING_NONE);
  if (ctxt==NULL)
    IGRAPH_ERROR("Can't create progressive parser context", IGRAPH_PARSEERROR);

  // Parse the file
  while ((res=fread(buffer, 1, 4096, instream))>0) {
    xmlParseChunk(ctxt, buffer, res, 0);
    if (!state.successful) break;
  }
  xmlParseChunk(ctxt, buffer, res, 1);
  
  // Free the context
  xmlFreeParserCtxt(ctxt);
  if (!state.successful)
    IGRAPH_ERROR("Malformed GraphML file", IGRAPH_PARSEERROR);
  if (state.index>=0)
    IGRAPH_ERROR("Graph index was too large", IGRAPH_EINVAL);
  
  return 0;
#else
  IGRAPH_ERROR("GraphML support is disabled", IGRAPH_UNIMPLEMENTED);
#endif
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_edgelist
 * \brief Writes the edge list of a graph to a file.
 * 
 * One edge is written per line, separated by a single space.
 * For directed graphs edges are written in from, to order.
 * \param graph The graph object to write.
 * \param outstream Pointer to a stream, it should be writable.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error writing the
 *         file. 
 * 
 * Time complexity: O(|E|), the
 * number of edges in the  graph. It is assumed that writing an
 * integer to the file requires O(1)
 * time. 
 */

int igraph_write_graph_edgelist(const igraph_t *graph, FILE *outstream) {

  igraph_eit_t it;
  
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_FROM), 
				 &it));
  IGRAPH_FINALLY(igraph_eit_destroy, &it);

  while (!IGRAPH_EIT_END(it)) {
    igraph_integer_t from, to;
    int ret;
    igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
    ret=fprintf(outstream, "%li %li\n", 
		(long int) from,
		(long int) to);
    if (ret < 0) {
      IGRAPH_ERROR("Write error", IGRAPH_EFILE);
    }
    IGRAPH_EIT_NEXT(it);
  }
  
  igraph_eit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/** 
 * \ingroup loadsave
 * \function igraph_write_graph_ncol
 * \brief Writes the graph to a file in <code>.ncol</code> format
 * 
 * <code>.ncol</code> is a format used by LGL, see \ref
 * igraph_read_graph_ncol() for details. 
 * 
 * Note that having multiple or loop edges in an
 * <code>.ncol</code> file breaks the  LGL software but 
 * \a igraph does not check for this condition. 
 * \param graph The graph to write.
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param names The name of the vertex attribute, if symbolic names
 *        are written to the file. If not, supply 0 here.
 * \param weights The name of the edge attribute, if they are also
 *        written to the file. If you don't want weights, supply 0
 *        here.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error writing the
 *         file. 
 * 
 * Time complexity: O(|E|), the
 * number of edges. All file operations are expected to have time
 * complexity O(1). 
 *
 * \sa \ref igraph_read_graph_ncol(), \ref igraph_write_graph_lgl()
 */

int igraph_write_graph_ncol(const igraph_t *graph, FILE *outstream, 
			    const char *names, const char *weights) {
  igraph_eit_t it;
  igraph_attribute_type_t nametype, weighttype;
  
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_FROM), 
				 &it));
  IGRAPH_FINALLY(igraph_eit_destroy, &it);

  /* Check if we have the names attribute */
  if (names && !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX,
					    names)) {
    names=0;
    IGRAPH_WARNING("names attribute does not exists");
  } 
  if (names) {
    IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &nametype,
					    IGRAPH_ATTRIBUTE_VERTEX, names));
  }
  if (names && nametype != IGRAPH_ATTRIBUTE_NUMERIC && 
      nametype != IGRAPH_ATTRIBUTE_STRING) {
    IGRAPH_WARNING("ignoring names attribute, unknown attribute type");
    names=0;
  }

  /* Check the weights as well */
  if (weights && !igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
					       weights)) {
    weights=0;
    IGRAPH_WARNING("weights attribute does not exists");
  }
  if (weights) {
    IGRAPH_CHECK(igraph_i_attribute_gettype(graph, &weighttype, 
					    IGRAPH_ATTRIBUTE_EDGE, weights));
  }
  if (weights && weighttype != IGRAPH_ATTRIBUTE_NUMERIC) {
    IGRAPH_WARNING("ignoring weights attribute, unknown attribute type");
    weights=0;
  }

  if (names==0 && weights ==0) {
    /* No names, no weights */
    while (!IGRAPH_EIT_END(it)) {
      igraph_integer_t from, to;
      int ret;
      igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
      ret=fprintf(outstream, "%li %li\n",
		  (long int) from,
		  (long int) to);
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
  } else if (weights==0) {
    /* No weights, but use names */
    igraph_strvector_t nvec;
    IGRAPH_CHECK(igraph_strvector_init(&nvec, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
    IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names, 
							   igraph_vss_all(),
							   &nvec));
    while (!IGRAPH_EIT_END(it)) {
      igraph_integer_t edge=IGRAPH_EIT_GET(it);
      igraph_integer_t from, to;
      int ret=0;
      char *str1, *str2;
      igraph_edge(graph, edge, &from, &to);
      igraph_strvector_get(&nvec, from, &str1);
      igraph_strvector_get(&nvec, to, &str2);
      ret=fprintf(outstream, "%s %s\n", str1, str2);
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
    IGRAPH_FINALLY_CLEAN(1);
  } else if (names==0) {
    /* No names but weights */
    igraph_strvector_t wvec;
    IGRAPH_CHECK(igraph_strvector_init(&wvec, igraph_ecount(graph)));
    IGRAPH_FINALLY(igraph_strvector_destroy, &wvec);
    IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, weights, 
							 igraph_ess_all(IGRAPH_EDGEORDER_FROM), 
							 &wvec));
    while (!IGRAPH_EIT_END(it)) {
      igraph_integer_t edge=IGRAPH_EIT_GET(it);
      igraph_integer_t from, to;
      int ret=0;
      char *str1;
      igraph_edge(graph, edge, &from, &to);
      igraph_strvector_get(&wvec, edge, &str1);
      ret=fprintf(outstream, "%li %li %s\n", 
		  (long int)from, (long int)to, str1);
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);      
    }
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    /* Both names and weights */
    igraph_strvector_t nvec, wvec;
    IGRAPH_CHECK(igraph_strvector_init(&wvec, igraph_ecount(graph)));
    IGRAPH_FINALLY(igraph_strvector_destroy, &wvec);
    IGRAPH_CHECK(igraph_strvector_init(&nvec, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
    IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, weights, 
							 igraph_ess_all(IGRAPH_EDGEORDER_FROM), 
							 &wvec));
    IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names, 
							   igraph_vss_all(),
							   &nvec));
    while (!IGRAPH_EIT_END(it)) {
      igraph_integer_t edge=IGRAPH_EIT_GET(it);
      igraph_integer_t from, to;
      int ret=0;
      char *str1, *str2, *str3;
      igraph_edge(graph, edge, &from, &to);
      igraph_strvector_get(&nvec, from, &str1);
      igraph_strvector_get(&nvec, to, &str2);
      igraph_strvector_get(&wvec, edge, &str3);
      ret=fprintf(outstream, "%s %s ", str1, str2);
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      ret=fprintf(outstream, "%s\n", str3);
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
    IGRAPH_FINALLY_CLEAN(2);
  }
  
  igraph_eit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_lgl
 * \brief Writes the graph to a file in <code>.lgl</code> format
 *
 * <code>.lgl</code> is a format used by LGL, see \ref
 * igraph_read_graph_lgl() for details.
 *
 * Note that having multiple or loop edges in an
 * <code>.lgl</code> file breaks the  LGL software but \a igraph
 * does not check for this condition. 
 * \param graph The graph to write. 
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param names The name of the vertex attribute, if symbolic names
 *        are written to the file. If not supply 0 here.
 * \param weights The name of the edge attribute, if they are also
 *        written to the file. If you don't want weights supply 0
 *        here.
 * \param isolates Logical, if TRUE isolated vertices are also written
 *        to the file. If FALSE they will be omitted.
 * \return Error code:
 *         \c IGRAPH_EFILE if there is an error
 *         writing the file. 
 *
 * Time complexity: O(|E|), the
 * number of edges if \p isolates is
 * FALSE, O(|V|+|E|) otherwise. All
 * file operations are expected to have time complexity 
 * O(1). 
 *
 * \sa \ref igraph_read_graph_ncol(), \ref igraph_write_graph_lgl()
 */

int igraph_write_graph_lgl(const igraph_t *graph, FILE *outstream,
			   const char *names, const char *weights,
			   igraph_bool_t isolates) {
  igraph_eit_t it;
  long int actvertex=-1;
  
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_FROM),
				 &it));
  IGRAPH_FINALLY(igraph_eit_destroy, &it);
/*   if (names==0 && weights==0) { */
    /* No names, no weights */
    while (!IGRAPH_EIT_END(it)) {
      igraph_integer_t from, to;
      int ret;
      igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
      if (from==actvertex) {
	ret=fprintf(outstream, "%li\n", (long int)to);
      } else {
	actvertex=from;
	ret=fprintf(outstream, "# %li\n%li\n", (long int)from, (long int)to);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
/*   } else if (weights==0) { */
/*     /\* No weights but use names *\/ */
/*     while (!igraph_es_end(graph, &it)) { */
/*       igraph_attribute_type_t type; */
/*       long int from=igraph_es_from(graph, &it); */
/*       long int to=igraph_es_to(graph, &it);       */
/*       void *ptr1, *ptr2; */
/*       int ret; */
/*       IGRAPH_CHECK(igraph_get_vertex_attribute(graph, names, to, &ptr2, &type)); */
/*       if (from==actvertex) { */
/* 	if (type==IGRAPH_ATTRIBUTE_NUM) { */
/* 	  ret=fprintf(outstream, "%f\n", *(igraph_real_t*)ptr2); */
/* 	} else { */
/* 	  ret=fprintf(outstream, "%s\n", (char*)ptr2); */
/* 	} */
/*       } else { */
/* 	actvertex=from; */
/* 	IGRAPH_CHECK(igraph_get_vertex_attribute(graph, names, from, &ptr1, &type)); */
/* 	if (type==IGRAPH_ATTRIBUTE_NUM) { */
/* 	  ret=fprintf(outstream, "# %f\n%f\n", *(igraph_real_t*)ptr1, *(igraph_real_t*)ptr2); */
/* 	} else { */
/* 	  ret=fprintf(outstream, "# %s\n%s\n", (char*)ptr1, (char*)ptr2); */
/* 	} */
/*       } */
/*       if (ret<0) { */
/* 	IGRAPH_ERROR("Write failed", IGRAPH_EFILE); */
/*       } */
/*       igraph_es_next(graph, &it); */
/*     } */
/*   } else if (names==0) { */
/*     /\* No names but weights *\/ */
/*     while (!igraph_es_end(graph, &it)) { */
/*       igraph_attribute_type_t type; */
/*       long int from=igraph_es_from(graph, &it); */
/*       long int to=igraph_es_to(graph, &it); */
/*       long int edge=igraph_es_get(graph, &it); */
/*       void *ptr; */
/*       int ret; */
/*       IGRAPH_CHECK(igraph_get_edge_attribute(graph, weights, edge, &ptr, &type)); */
/*       if (from==actvertex) { */
/* 	if (type==IGRAPH_ATTRIBUTE_NUM) { */
/* 	  ret=fprintf(outstream, "%li %f\n", to, *(igraph_real_t*)ptr); */
/* 	} else { */
/* 	  ret=fprintf(outstream, "%li %s\n", to, (char*)ptr); */
/* 	} */
/*       } else { */
/* 	actvertex=from; */
/* 	if (type==IGRAPH_ATTRIBUTE_NUM) { */
/* 	  ret=fprintf(outstream, "# %li\n%li %f\n", from, to, *(igraph_real_t*)ptr); */
/* 	} else { */
/* 	  ret=fprintf(outstream, "# %li\n%li %s\n", from, to, (char*)ptr); */
/* 	} */
/*       } */
/*       if (ret<0) { */
/* 	IGRAPH_ERROR("Write failed", IGRAPH_EFILE); */
/*       } */
/*       igraph_es_next(graph, &it); */
/*     } */
/*   } else { */
/*     /\* Both names and weights *\/ */
/*     while (!igraph_es_end(graph, &it)) { */
/*       igraph_attribute_type_t vtype, etype; */
/*       long int from=igraph_es_from(graph, &it); */
/*       long int to=igraph_es_to(graph, &it); */
/*       long int edge=igraph_es_get(graph, &it); */
/*       void *ptr1, *ptr2, *ptr; */
/*       int ret; */
/*       IGRAPH_CHECK(igraph_get_vertex_attribute(graph, names, to, &ptr2, &vtype)); */
/*       IGRAPH_CHECK(igraph_get_edge_attribute(graph, weights, edge, &ptr, &etype)); */
/*       if (from==actvertex) { */
/* 	if (vtype==IGRAPH_ATTRIBUTE_NUM) { */
/* 	  ret=fprintf(outstream, "%f ", *(igraph_real_t*)ptr2); */
/* 	} else { */
/* 	  ret=fprintf(outstream, "%s ", (char*)ptr2); */
/* 	}  */
/*       } else { */
/* 	actvertex=from; */
/* 	IGRAPH_CHECK(igraph_get_vertex_attribute(graph, names, from, &ptr1, &vtype)); */
/* 	if (vtype==IGRAPH_ATTRIBUTE_NUM) { */
/* 	  ret=fprintf(outstream, "# %f\n%f ", *(igraph_real_t*)ptr1, *(igraph_real_t*)ptr2); */
/* 	} else { */
/* 	  ret=fprintf(outstream, "# %s\n%s ", (char*)ptr1, (char*)ptr2); */
/* 	} */
/*       } */
/*       if (ret<0) { */
/* 	IGRAPH_ERROR("Write failed", IGRAPH_EFILE); */
/*       } */
/*       if (etype==IGRAPH_ATTRIBUTE_NUM) { */
/* 	ret=fprintf(outstream, "%f\n", *(igraph_real_t*)ptr); */
/*       } else { */
/* 	ret=fprintf(outstream, "%s\n", (char*)ptr); */
/*       } */
/*       if (ret<0) { */
/* 	IGRAPH_ERROR("Write failed", IGRAPH_EFILE); */
/*       } */
/*       igraph_es_next(graph, &it); */
/*     } */
/*   } */

  if (isolates) {
    long int nov=igraph_vcount(graph);
    long int i;
/*     igraph_attribute_type_t type; */
    void *ptr;
    int ret=0;
    igraph_vector_t deg;

    IGRAPH_VECTOR_INIT_FINALLY(&deg, 1);
    for (i=0; i<nov; i++) {
      igraph_degree(graph, &deg, igraph_vss_1(i), IGRAPH_ALL, IGRAPH_LOOPS);
      if (VECTOR(deg)[0]==0) {
/* 	if (names==0) { */
	  ret=fprintf(outstream, "# %li\n", i);
/* 	} else { */
/* 	  IGRAPH_CHECK(igraph_get_vertex_attribute(graph, names, i, &ptr, &type)); */
/* 	  if (type==IGRAPH_ATTRIBUTE_NUM) { */
/* 	    ret=fprintf(outstream, "# %f\n", *(igraph_real_t*)ptr); */
/* 	  } else { */
/* 	    ret=fprintf(outstream, "# %s\n", (char*)ptr); */
/* 	  } */
/* 	} */
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
    }
    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(1);
  }  
  
  igraph_eit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_graphml
 * \brief Writes the graph to a file in GraphML format
 *
 * GraphML is an XML-based file format for representing various types of
 * graphs. See the GraphML Primer (<ulink>http://graphml.graphdrawing.org/primer/graphml-primer.html</ulink>)
 * for detailed format description.
 * 
 * No attributes are written at present.
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
 */
int igraph_write_graph_graphml(const igraph_t *graph, FILE *outstream) {
  int ret;
  igraph_integer_t l, vc;
  igraph_eit_t it;
  
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
  ret=fprintf(outstream, "  <graph id=\"G\" edgedefault=\"%s\">\n", (igraph_is_directed(graph)?"directed":"undirected"));
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  
  /* Let's dump the nodes first */
  vc=igraph_vcount(graph);
  for (l=0; l<vc; l++) {
    ret=fprintf(outstream, "    <node id=\"n%ld\" />\n", (long)l);
    if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  }
  
  /* Now the edges */
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(0), &it));
  IGRAPH_FINALLY(igraph_eit_destroy, &it);
  while (!IGRAPH_EIT_END(it)) {
    igraph_integer_t from, to;
    igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
    ret=fprintf(outstream, "    <edge source=\"n%ld\" target=\"n%ld\" />\n", 
		(long int)from, (long int)to);
    if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
    IGRAPH_EIT_NEXT(it);
  }
  igraph_eit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  
  ret=fprintf(outstream, "  </graph>\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  fprintf(outstream, "</graphml>\n");
  if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
  
  return 0;
}
