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
#include "config.h"

#include <ctype.h>		/* isspace */
#include <string.h>
#include "memory.h"

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
 * </para><para>
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

    IGRAPH_ALLOW_INTERRUPTION();
    
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

extern int igraph_ncol_yyparse(void);
extern FILE *igraph_ncol_yyin;
extern int igraph_i_ncol_eof;
long int igraph_ncol_mylineno;
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
 * </para><para>
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
 * </para><para>
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
  igraph_ncol_mylineno=1;
  igraph_i_ncol_eof=0;

  igraph_ncol_yyparse();

  if (predefnames != 0 && 
      igraph_trie_size(&trie) != no_predefined) {
    IGRAPH_WARNING("unknown vertex/vertices found, predefnames extended");    
  }

  if (names) {
    const igraph_strvector_t *namevec;
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

extern int igraph_lgl_yyparse(void);
extern FILE *igraph_lgl_yyin;
extern int igraph_i_lgl_eof;
long int igraph_lgl_mylineno;
igraph_vector_t *igraph_lgl_vector=0;
igraph_vector_t *igraph_lgl_weights=0;
igraph_trie_t *igraph_lgl_trie=0;

/**
 * \ingroup loadsave
 * \function igraph_read_graph_lgl
 * \brief Reads a graph from an <code>.lgl</code> file
 * 
 * </para><para>
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
 * </para><para>
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
  igraph_lgl_mylineno=1;
  igraph_i_lgl_eof=0;

  igraph_lgl_yyparse();
  
  IGRAPH_CHECK(igraph_empty(graph, 0, IGRAPH_UNDIRECTED));
  IGRAPH_FINALLY(igraph_destroy, graph);

  if (names) {
    const igraph_strvector_t *namevec;
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

  IGRAPH_CHECK(igraph_add_vertices(graph, igraph_trie_size(&trie), pname));
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

extern int igraph_pajek_yyparse(void);
extern FILE *igraph_pajek_yyin;
extern int igraph_i_pajek_eof;
long int igraph_pajek_mylineno;
igraph_vector_t *igraph_pajek_vector=0;
igraph_bool_t igraph_pajek_directed;
long int igraph_pajek_vcount=0;
long int igraph_pajek_actfrom, igraph_pajek_actto;
int igraph_pajek_mode=0;	/* 0 - general, 1 - vertex, 2 - edge */
igraph_trie_t *igraph_i_pajek_vertex_attribute_names;
igraph_vector_ptr_t *igraph_i_pajek_vertex_attributes;
igraph_trie_t *igraph_i_pajek_edge_attribute_names;
igraph_vector_ptr_t *igraph_i_pajek_edge_attributes;
long int igraph_i_pajek_vertexid=0;
long int igraph_i_pajek_actvertex=0;
long int igraph_i_pajek_actedge=0;

/* int vector_print(igraph_vector_t *v) { */
/*   long int i, size=igraph_vector_size(v); */
/*   for (i=0; i<size; i++) { */
/*     printf("%f|", VECTOR(*v)[i]); */
/*   } */
/*   printf("\n"); */
/*   return 0; */
/* } */

/* int strvector_print(igraph_strvector_t *sv) { */
/*   long int i, size=igraph_strvector_size(sv); */
/*   char *str; */
/*   for (i=0; i<size; i++) { */
/*     igraph_strvector_get(sv, i, &str); */
/*     printf("%s|", str); */
/*   } */
/*   printf("\n"); */
/*   return 0; */
/* } */

/**
 * \function igraph_read_graph_pajek
 * \brief Reads a file in Pajek format
 * 
 * \param graph Pointer to an uninitialized graph object.
 * \param file An already opened file handler.
 * \return Error code.
 * 
 * </para><para>
 * Only a subset of the Pajek format is implemented. This is partially
 * because this format is not very well documented, but also because
 * <command>igraph</command> does not support some Pajek features, like 
 * multigraphs.
 * 
 * </para><para> 
 * The list of the current limitations:
 * \olist
 * \oli Only <filename>.net</filename> files are supported, Pajek
 * project files (<filename>.paj</filename>) are not. These might be
 * supported in the future if there is need for it.
 * \oli Time events networks are not supported.
 * \oli Hypergraphs (ie. graphs with non-binary edges) are not
 * supported.
 * \oli Graphs with both directed and non-directed edges are not
 * supported, are they cannot be represented in
 * <command>igraph</command>.
 * \oli Bipartite or affiliation networks are not supported. They can
 * be imported but the vertex type information is omitted.
 * \oli Only Pajek networks are supported, permutations, hierarchies,
 * clusters and vectors are not.
 * \oli Graphs with multiple edge sets are not supported.
 * \endolist
 * 
 * </para><para>
 * If there are attribute handlers installed,
 * <command>igraph</command> also reads the vertex and edge attributes
 * from the file. Most attributes are renamed to be more informative: 
 * `\c color' instead of `\c c', `\c xfact' instead of `\c x_fact',
 * `\c yfact' instead of `y_fact', `\c labeldist' instead of `\c lr',
 * `\c labeldegree2' instead of `\c lphi', `\c framewidth' instead of `\c bw',
 * `\c fontsize'
 * instead of `\c fos', `\c rotation' instead of `\c phi', `\c radius' instead
 * of `\c r',
 * `\c diamondratio' instead of `\c q', `\c labeldegree' instead of `\c la',
 * `\c vertexsize'
 * instead of `\c size', `\c color' instead of `\c ic', `\c framecolor' instead of
 * `\c bc', `\c labelcolor' instead of `\c lc', these belong to vertices. 
 * 
 * </para><para>
 * Edge attributes are also renamed, `\c s' to `\c arrowsize', `\c w'
 * to `\c edgewidth', `\c h1' to `\c hook1', `\c h2' to `\c hook2',
 * `\c a1' to `\c angle1', `\c a2' to `\c angle2', `\c k1' to 
 * `\c velocity1', `\c k2' to `\c velocity2', `\c ap' to `\c
 * arrowpos', `\c lp' to `\c labelpos', `\c lr' to 
 * `\c labelangle', `\c lphi' to `\c labelangle2', `\c la' to `\c
 * labeldegree', `\c fos' to 
 * `\c fontsize', `\c a' to `\c arrowtype', `\c p' to `\c
 * linepattern', `\c l' to `\c label', `\c lc' to 
 * `\c labelcolor', `\c c' to `\c color'.
 * 
 * </para><para>
 * In addition the following vertex attributes might be added: `\c id'
 * if there are vertex ids in the file, `\c x' and `\c y' or `\c x'
 * and `\c y' and `\c z' if there are vertex coordinates in the file,
 * `\c color-red', `\c color-green' and `\c color-blue' if the vertex
 * color is given in RGB notation, `\c framecolor-red', `\c
 * framecolor-green' and `\c framecolor-blue` if the frame color is
 * given in RGB notation and finally `\c labelcolor-red', `\c
 * labelcolor-green' and `\c labelcolor-blue' if the label color is
 * given in RGB notation.
 * 
 * </para><para>The following additional edge attributes might be
 * added: `\c weight' if there are edge weights present, `\c
 * color-red', `\c color-green' and `\c color-blue' if the edge color
 * is given in RGB notation. 
 * 
 * </para><para>
 * See the pajek homepage:
 * http://vlado.fmf.uni-lj.si/pub/networks/pajek/ for more info on
 * Pajek and the Pajek manual:
 * http://vlado.fmf.uni-lj.si/pub/networks/pajek/doc/pajekman.pdf for
 * information on the Pajek file format.
 *
 * </para><para>
 * Time complexity: O(|V|+|E|+|A|), |V| is the number of vertices, |E|
 * the number of edges, |A| the number of attributes (vertex + edge)
 * in the graph if there are attribute handlers installed.
 * 
 * \sa \ref igraph_write_graph_pajek() for writing Pajek files, \ref
 * igraph_read_graph_graphml() for reading GraphML files.
 */

int igraph_read_graph_pajek(igraph_t *graph, FILE *instream) {

  igraph_vector_t edges;
  igraph_trie_t vattrnames;
  igraph_vector_ptr_t vattrs;
  igraph_trie_t eattrnames;
  igraph_vector_ptr_t eattrs;
  /* igraph_hashtable_t vattrhash; */
  /* igraph_hashtable_t eattrhash; */
  long int i, j;
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  IGRAPH_TRIE_INIT_FINALLY(&vattrnames, 1);
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&vattrs, 0);
  IGRAPH_TRIE_INIT_FINALLY(&eattrnames, 1);
  IGRAPH_VECTOR_PTR_INIT_FINALLY(&eattrs, 0);

  igraph_pajek_vector=&edges;
  igraph_pajek_yyin=instream;

  igraph_pajek_mode=0;
  igraph_pajek_vcount=0;
  igraph_i_pajek_vertexid=0;
  igraph_i_pajek_vertex_attribute_names=&vattrnames;
  igraph_i_pajek_vertex_attributes=&vattrs;
  igraph_i_pajek_edge_attribute_names=&eattrnames;
  igraph_i_pajek_edge_attributes=&eattrs;
  igraph_i_pajek_actedge=0;
  igraph_pajek_mylineno=1;
  igraph_i_pajek_eof=0;

  igraph_pajek_yyparse();

  for (i=0; i<igraph_vector_ptr_size(&eattrs); i++) {
    igraph_i_attribute_record_t *rec=VECTOR(eattrs)[i];
    if (rec->type==IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_vector_t *vec=(igraph_vector_t*)rec->value;
      long int origsize=igraph_vector_size(vec);
      igraph_vector_resize(vec, igraph_i_pajek_actedge);
      for (j=origsize; j<igraph_i_pajek_actedge; j++) {
	VECTOR(*vec)[j] = IGRAPH_NAN;
      }
    } else if (rec->type==IGRAPH_ATTRIBUTE_STRING) {
      igraph_strvector_t *strvec=(igraph_strvector_t*)rec->value;
      long int origsize=igraph_strvector_size(strvec);
      igraph_strvector_resize(strvec, igraph_i_pajek_actedge);
      for (j=origsize; j<igraph_i_pajek_actedge; j++) {
	igraph_strvector_set(strvec, j, "");
      }
    }
  }

  IGRAPH_CHECK(igraph_empty(graph, 0, igraph_pajek_directed));
  IGRAPH_FINALLY(igraph_destroy, graph);
  IGRAPH_CHECK(igraph_add_vertices(graph, igraph_pajek_vcount, &vattrs));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, &eattrs));

  for (i=0; i<igraph_vector_ptr_size(&vattrs); i++) {
    igraph_i_attribute_record_t *rec=VECTOR(vattrs)[i];
    if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_vector_t *vec=(igraph_vector_t*) rec->value;
      igraph_vector_destroy(vec);
      Free(vec);
    } else if (rec->type==IGRAPH_ATTRIBUTE_STRING) {
      igraph_strvector_t *strvec=(igraph_strvector_t *)rec->value;
      igraph_strvector_destroy(strvec);
      Free(strvec);
    }
    Free(rec);
  }

  for (i=0; i<igraph_vector_ptr_size(&eattrs); i++) {
    igraph_i_attribute_record_t *rec=VECTOR(eattrs)[i];
    if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_vector_t *vec=(igraph_vector_t*) rec->value;
      igraph_vector_destroy(vec);
      Free(vec);
    } else if (rec->type==IGRAPH_ATTRIBUTE_STRING) {
      igraph_strvector_t *strvec=(igraph_strvector_t *)rec->value;
      igraph_strvector_destroy(strvec);
      Free(strvec);
    }
    Free(rec);
  }

  igraph_vector_destroy(&edges);  
  igraph_vector_ptr_destroy(&eattrs);
  igraph_trie_destroy(&eattrnames);
  igraph_vector_ptr_destroy(&vattrs);
  igraph_trie_destroy(&vattrnames);

  IGRAPH_FINALLY_CLEAN(6);
  return 0;
}

/**
 * \function igraph_read_graph_dimacs
 */

int igraph_read_graph_dimacs(igraph_t *graph, FILE *instream,
			     igraph_integer_t *source, 
			     igraph_integer_t *target, 
			     igraph_vector_t *capacity, 
			     igraph_bool_t directed) {
  
  igraph_vector_t edges;
  long int no_of_nodes=-1;
  long int no_of_edges=-1;
  long int tsource=-1;
  long int ttarget=-1;
  char problem[6];
  char c;      
  
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  if (capacity) {
    igraph_vector_clear(capacity);
  }
  
  while (!feof(instream)) {
    int read;
    char str[3];
    
    IGRAPH_ALLOW_INTERRUPTION();
    
    read=fscanf(instream, "%2c", str);
    if (feof(instream)) {
      break;
    }
    if (read != 1) {
      IGRAPH_ERROR("parsing dimacs file failed", IGRAPH_PARSEERROR);
    }
    switch (str[0]) {
      long int tmp;
      long int from, to;
      igraph_real_t cap;

    case 'c':
      /* comment */
      break;

    case 'p':
      if (no_of_nodes != -1) {
	IGRAPH_ERROR("reading dimacs file failed, double 'p' line", 
		     IGRAPH_PARSEERROR);
      }
      read=fscanf(instream, "%5s %li %li", problem, 
		  &no_of_nodes, &no_of_edges);
      if (read != 3) {
	IGRAPH_ERROR("reading dimacs file failed", IGRAPH_PARSEERROR);
      }
      IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges*2));
      if (capacity) {
	IGRAPH_CHECK(igraph_vector_reserve(capacity, no_of_edges));
      }
      break;

    case 'n':
      read=fscanf(instream, "%li %1c", &tmp, str);
      if (str[0]=='s') {
	if (tsource != -1) {
	  IGRAPH_ERROR("reading dimacsfile: multiple source vertex line", 
		       IGRAPH_PARSEERROR);
	} else {
	  tsource=tmp;
	}
      } else if (str[0]=='t') {
	if (ttarget != -1) {
	  IGRAPH_ERROR("reading dimacsfile: multiple source vertex line", 
		       IGRAPH_PARSEERROR);
	} else {
	  ttarget=tmp;
	}
      } else {
	IGRAPH_ERROR("invalid node descriptor line in dimacs file",
		     IGRAPH_PARSEERROR);
      }
      
      break;
      
    case 'a':
      read=fscanf(instream, "%li %li %lf", &from, &to, &cap);
      if (read != 3) {
	IGRAPH_ERROR("reading dimacs file", IGRAPH_PARSEERROR);
      }
      IGRAPH_CHECK(igraph_vector_push_back(&edges, from-1));
      IGRAPH_CHECK(igraph_vector_push_back(&edges, to-1));
      if (capacity) {
	IGRAPH_CHECK(igraph_vector_push_back(capacity, cap));
      }
      break;           

    default:
      IGRAPH_ERROR("unknown line type in dimacs file", IGRAPH_PARSEERROR);
    }

    /* Go to next line */
    while (!feof(instream) && (c=getc(instream)) != '\n') ;      
  }

  if (source) {
    *source=tsource-1;
  }
  if (target) {
    *target=ttarget-1;
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, directed));
  igraph_vector_destroy(&edges);

  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

int igraph_i_read_graph_graphdb_getword(FILE *instream) {
  int b1, b2;
  unsigned char c1, c2;
  b1 = fgetc(instream);
  b2 = fgetc(instream);
  if (b1 != EOF) {
    c1=b1; c2=b2;
    return c1 | (c2<<8);
  } else {
    return -1;
  }
}

int igraph_read_graph_graphdb(igraph_t *graph, FILE *instream, 
			      igraph_bool_t directed) {
  
  igraph_vector_t edges;
  long int nodes;
  long int i, j;
  igraph_bool_t end=0;

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  nodes=igraph_i_read_graph_graphdb_getword(instream);
  if (nodes<0) {
    IGRAPH_ERROR("Can't read from file", IGRAPH_EFILE);
  }
  for (i=0; !end && i<nodes; i++) {
    long int len=igraph_i_read_graph_graphdb_getword(instream);
    if (len<0) {
      end=1;
      break;
    }
    for (j=0; ! end && j<len; j++) {
      long int to=igraph_i_read_graph_graphdb_getword(instream);
      if (to<0) {
	end=1; 
	break;
      }
      IGRAPH_CHECK(igraph_vector_push_back(&edges, i));
      IGRAPH_CHECK(igraph_vector_push_back(&edges, to));
    }
  }

  if (end) {
    IGRAPH_ERROR("Truncated graphdb file", IGRAPH_EFILE);
  }
  
  IGRAPH_CHECK(igraph_create(graph, &edges, nodes, directed));
  igraph_vector_destroy(&edges);
  
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}
		             

/**
 * \ingroup loadsave
 * \function igraph_write_graph_edgelist
 * \brief Writes the edge list of a graph to a file.
 * 
 * </para><para>
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
 * </para><para>
 * <code>.ncol</code> is a format used by LGL, see \ref
 * igraph_read_graph_ncol() for details. 
 * 
 * </para><para>
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
    igraph_strvector_destroy(&nvec);
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
    igraph_strvector_destroy(&wvec);
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
    igraph_strvector_destroy(&nvec);
    igraph_strvector_destroy(&wvec);
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
 * </para><para>
 * <code>.lgl</code> is a format used by LGL, see \ref
 * igraph_read_graph_lgl() for details.
 *
 * </para><para>
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
  
  if (names==0 && weights==0) {
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
  } else if (weights==0) {
    /* No weights but use names */
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
      igraph_strvector_get(&nvec, to, &str2);

      if (from==actvertex) {
	ret=fprintf(outstream, "%s\n", str2);
      } else {
	actvertex=from;
	igraph_strvector_get(&nvec, from, &str1);
	ret=fprintf(outstream, "# %s\n%s\n", str1, str2);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
    IGRAPH_FINALLY_CLEAN(1);
  } else if (names==0) {
    igraph_strvector_t wvec;
    IGRAPH_CHECK(igraph_strvector_init(&wvec, igraph_ecount(graph)));
    IGRAPH_FINALLY(igraph_strvector_destroy, &wvec);
    IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, weights,
							 igraph_ess_all(IGRAPH_EDGEORDER_FROM),
							 &wvec));
    /* No names but weights */
    while (!IGRAPH_EIT_END(it)) {
      igraph_integer_t edge=IGRAPH_EIT_GET(it);
      igraph_integer_t from, to;
      int ret=0;
      char *str1;
      igraph_edge(graph, edge, &from, &to);
      igraph_strvector_get(&wvec, edge, &str1);
      if (from==actvertex) {
	ret=fprintf(outstream, "%li %s\n", (long)to, str1);
      } else {
	actvertex=from;
	ret=fprintf(outstream, "# %li\n%li %s\n", (long)from, (long)to, str1);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
    igraph_strvector_destroy(&wvec);
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
      igraph_strvector_get(&nvec, to, &str2);
      igraph_strvector_get(&wvec, edge, &str3);
      if (from==actvertex) {
	ret=fprintf(outstream, "%s ", str2);
      } else {
	actvertex=from;
	igraph_strvector_get(&nvec, from, &str1);
	ret=fprintf(outstream, "# %s\n%s ", str1, str2);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      ret=fprintf(outstream, "%s\n", str3);
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
    igraph_strvector_destroy(&nvec);
    igraph_strvector_destroy(&wvec);
    IGRAPH_FINALLY_CLEAN(2);
  }

  if (isolates) {
    long int nov=igraph_vcount(graph);
    long int i;
    int ret=0;
    igraph_vector_t deg;
    igraph_strvector_t nvec;
    char *str;

    IGRAPH_VECTOR_INIT_FINALLY(&deg, 1);
    IGRAPH_CHECK(igraph_strvector_init(&nvec, 1));
    IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
    for (i=0; i<nov; i++) {
      igraph_degree(graph, &deg, igraph_vss_1(i), IGRAPH_ALL, IGRAPH_LOOPS);
      if (VECTOR(deg)[0]==0) {
 	if (names==0) { 
	  ret=fprintf(outstream, "# %li\n", i);
	} else {
	  IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, names,
								 igraph_vss_1(i),
								 &nvec));
	  igraph_strvector_get(&nvec, 0, &str);
	  ret=fprintf(outstream, "# %s\n", str);
	}
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
    }
    igraph_strvector_destroy(&nvec);
    igraph_vector_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(2);
  }  
  
  igraph_eit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/* Order matters here! */
#define V_ID                0
#define V_X                 1
#define V_Y                 2
#define V_Z                 3
#define V_SHAPE             4
#define V_XFACT             5
#define V_YFACT             6
#define V_COLOR_RED         7
#define V_COLOR_GREEN       8
#define V_COLOR_BLUE        9
#define V_FRAMECOLOR_RED   10
#define V_FRAMECOLOR_GREEN 11
#define V_FRAMECOLOR_BLUE  12
#define V_LABELCOLOR_RED   13
#define V_LABELCOLOR_GREEN 14
#define V_LABELCOLOR_BLUE  15
#define V_LABELDIST        16
#define V_LABELDEGREE2     17
#define V_FRAMEWIDTH       18
#define V_FONTSIZE         19
#define V_ROTATION         20
#define V_RADIUS           21
#define V_DIAMONDRATIO     22
#define V_LABELDEGREE      23
#define V_VERTEXSIZE       24
#define V_FONT             25
#define V_URL              26
#define V_COLOR            27
#define V_FRAMECOLOR       28
#define V_LABELCOLOR       29
#define V_LAST             30

#define E_WEIGHT            0
#define E_COLOR_RED         1
#define E_COLOR_GREEN       2
#define E_COLOR_BLUE        3
#define E_ARROWSIZE         4
#define E_EDGEWIDTH         5
#define E_HOOK1             6
#define E_HOOK2             7
#define E_ANGLE1            8
#define E_ANGLE2            9
#define E_VELOCITY1        10
#define E_VELOCITY2        11
#define E_ARROWPOS         12
#define E_LABELPOS         13
#define E_LABELANGLE       14
#define E_LABELANGLE2      15
#define E_LABELDEGREE      16
#define E_FONTSIZE         17
#define E_ARROWTYPE        18
#define E_LINEPATTERN      19
#define E_LABEL            20
#define E_LABELCOLOR       21
#define E_COLOR            22
#define E_LAST             23

/**
 * \function igraph_write_graph_pajek
 * \brief Writes a graph to a file in Pajek format.
 * 
 * </para><para>
 * The Pajek vertex and edge parameters (like color) are determined by
 * the attributes of the vertices and edges, of course this requires
 * an attribute handler to be installed. The names of the
 * corresponding vertex and edge attributes are listed at \ref
 * igraph_read_graph_pajek(), eg. the `\c color' vertex attributes
 * determines the color (`\c c' in Pajek) parameter.
 * \param graph The graph object to write.
 * \param outstream The file to write to. It should be opened and
 * writable.
 * \return Error code.
 *
 * Time complexity: O(|V|+|E|+|A|), |V| is the number of vertices, |E|
 * is the number of edges, |A| the number of attributes (vertex +
 * edge) in the graph if there are attribute handlers installed. 
 *
 * \sa \ref igraph_read_graph_pajek() for reading Pajek graphs, \ref
 * igraph_write_graph_graphml() for writing a graph in GraphML format,
 * this suites <command>igraph</command> graphs better.
 */

int igraph_write_graph_pajek(const igraph_t *graph, FILE *outstream) {
  long int no_of_nodes=igraph_vcount(graph);
  long int i, j;

  igraph_attribute_type_t vtypes[V_LAST], etypes[E_LAST];
  igraph_bool_t write_vertex_attrs=0;  

  /* Same order as the #define's */
  const char *vnames[] = { "id", "x", "y", "z", "shape", "xfact", "yfact",
			   "color-red", "color-green", "color-blue",
			   "framecolor-red", "framecolor-green", 
			   "framecolor-blue", "labelcolor-red", 
			   "labelcolor-green", "labelcolor-blue",
			   "labeldist", "labeldegree2", "framewidth",
			   "fontsize", "rotation", "radius", 
			   "diamondratio", "labeldegree", "vertexsize", 
			   "font", "url", "color", "framecolor",
			   "labelcolor" };

  const char *vnumnames[] = { "xfact", "yfact", "labeldist", 
			      "labeldegree2", "framewidth", "fontsize",
			      "rotation", "radius", "diamondratio",
			      "labeldegree", "vertexsize" };
  const char *vnumnames2[]= { "x_fact", "y_fact", "lr", "lphi", "bw",
			      "fos", "phi", "r", "q", "la", "size" };
  const char *vstrnames[] = { "font", "url", "color", "framecolor", 
			      "labelcolor" };
  const char *vstrnames2[]= { "font", "url", "ic", "bc", "lc" };  
  
  const char *enames[] = { "weight", "color-red", "color-green", "color-blue", 
			   "arrowsize", "edgewidth", "hook1", "hook2", 
			   "angle1", "angle2", "velocity1", "velocity2",
			   "arrowpos", "labelpos", "labelangle",
			   "labelangle2", "labeldegree", "fontsize",
			   "arrowtype", "linepattern", "label", "labelcolor",
			   "color" };
  const char *enumnames[] = { "arrowsize", "edgewidth", "hook1", "hook2",
			      "angle1", "angle2", "velocity1", "velocity2", 
			      "arrowpos", "labelpos", "labelangle",
			      "labelangle2", "labeldegree", "fontsize" };
  const char *enumnames2[]= { "s", "w", "h1", "h2", "a1", "a2", "k1", "k2", 
			      "ap", "lp", "lr", "lphi", "la", "fos" };
  const char *estrnames[] = { "arrowtype", "linepattern", "label",
			      "labelcolor", "color" };
  const char *estrnames2[]= { "a", "p", "l", "lc", "c" };

  const char *newline="\n\r";
  
  igraph_es_t es;
  igraph_eit_t eit;

  igraph_vector_t numv;
  igraph_strvector_t strv;

  igraph_vector_t ex_numa;
  igraph_vector_t ex_stra;
  igraph_vector_t vx_numa;
  igraph_vector_t vx_stra;

  IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
  IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);

  IGRAPH_VECTOR_INIT_FINALLY(&ex_numa, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&ex_stra, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vx_numa, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vx_stra, 0);

  /* Write header */
  if (fprintf(outstream, "*Vertices %li%s", no_of_nodes, newline) < 0) {
    IGRAPH_ERROR("Cannot write pajek file", IGRAPH_EFILE);
  }

  /* Check the vertex attributes */
  memset(vtypes, 0, sizeof(vtypes[0])*V_LAST);
  for (i=0; i<V_LAST; i++) {
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, 
				    vnames[i])) { 
      igraph_i_attribute_gettype(graph, &vtypes[i], IGRAPH_ATTRIBUTE_VERTEX, 
				 vnames[i]);
      write_vertex_attrs=1;
    } else {
      vtypes[i]=-1;
    }
  }
  for (i=0; i<sizeof(vnumnames)/sizeof(const char*); i++) {
    igraph_attribute_type_t type;
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, 
				    vnumnames[i])) {
      igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_VERTEX, 
				 vnumnames[i]);
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	IGRAPH_CHECK(igraph_vector_push_back(&vx_numa, i));
      }
    }
  }
  for (i=0; i<sizeof(vstrnames)/sizeof(const char*); i++) {
    igraph_attribute_type_t type;
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_VERTEX, 
				    vstrnames[i])) {
      igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_VERTEX, 
				 vstrnames[i]);
      if (type==IGRAPH_ATTRIBUTE_STRING) {
	IGRAPH_CHECK(igraph_vector_push_back(&vx_stra, i));
      }
    }
  }

  /* Write vertices */
  if (write_vertex_attrs) {
    for (i=0; i<no_of_nodes; i++) {

      /* vertex id */
      fprintf(outstream, "%li", i+1);
      if (vtypes[V_ID] == IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_ID], 
						   igraph_vss_1(i), &numv);
	fprintf(outstream, " \"%g\"", VECTOR(numv)[0]);
      } else if (vtypes[V_ID] == IGRAPH_ATTRIBUTE_STRING) {
	char *s;
	igraph_i_attribute_get_string_vertex_attr(graph, vnames[V_ID],
						  igraph_vss_1(i), &strv);
	igraph_strvector_get(&strv, 0, &s);
	fprintf(outstream, " \"%s\"", s);
      }
      
      /* coordinates */
      if (vtypes[V_X] == IGRAPH_ATTRIBUTE_NUMERIC &&
	  vtypes[V_Y] == IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_X], 
						   igraph_vss_1(i), &numv);
	fprintf(outstream, " %g", VECTOR(numv)[0]);
	igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_Y], 
						   igraph_vss_1(i), &numv);
	fprintf(outstream, " %g", VECTOR(numv)[0]);
	if (vtypes[V_Z] == IGRAPH_ATTRIBUTE_NUMERIC) {
	  igraph_i_attribute_get_numeric_vertex_attr(graph, vnames[V_Z], 
						     igraph_vss_1(i), &numv);
	  fprintf(outstream, " %g", VECTOR(numv)[0]);
	}
      }
      
      /* shape */
      if (vtypes[V_SHAPE] == IGRAPH_ATTRIBUTE_STRING) {
	char *s;
	igraph_i_attribute_get_string_vertex_attr(graph, vnames[V_SHAPE],
						  igraph_vss_1(i), &strv);
	igraph_strvector_get(&strv, 0, &s);
	fprintf(outstream, " %s", s);
      }
      
      /* numeric parameters */
      for (j=0; j<igraph_vector_size(&vx_numa); j++) {
	int idx=VECTOR(vx_numa)[j];
	igraph_i_attribute_get_numeric_vertex_attr(graph, vnumnames[idx],
						   igraph_vss_1(i), &numv);
	fprintf(outstream, " %s %g", vnumnames2[idx], VECTOR(numv)[0]);
      }

      /* string parameters */
      for (j=0; j<igraph_vector_size(&vx_stra); j++) {
	int idx=VECTOR(vx_numa)[j];
	char *s;
	igraph_i_attribute_get_string_vertex_attr(graph, vstrnames[idx],
						  igraph_vss_1(i), &strv);
	igraph_strvector_get(&strv, 0, &s);
	fprintf(outstream, " %s \"%s\"", vstrnames2[idx], s);
      }      
      
      /* trailing newline */
      fprintf(outstream, "%s", newline);
    }
  }

  /* edges header */
  if (igraph_is_directed(graph)) {
    fprintf(outstream, "*Arcs%s", newline);
  } else {
    fprintf(outstream, "*Edges%s", newline);
  }
  
  IGRAPH_CHECK(igraph_es_all(&es, IGRAPH_EDGEORDER_ID));
  IGRAPH_FINALLY(igraph_es_destroy, &es);
  IGRAPH_CHECK(igraph_eit_create(graph, es, &eit));
  IGRAPH_FINALLY(igraph_eit_destroy, &eit);

  /* Check edge attributes */
  for (i=0; i<E_LAST; i++) {
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE,
				    enames[i])) {
      igraph_i_attribute_gettype(graph, &etypes[i], IGRAPH_ATTRIBUTE_EDGE,
				 enames[i]);
    } else {
      etypes[i]=-1;
    }
  }
  for (i=0; i<sizeof(enumnames)/sizeof(const char*); i++) {
    igraph_attribute_type_t type;
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, 
				    enumnames[i])) {
      igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_EDGE, 
				 enumnames[i]);
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	IGRAPH_CHECK(igraph_vector_push_back(&ex_numa, i));
      }
    }
  }
  for (i=0; i<sizeof(estrnames)/sizeof(const char*); i++) {
    igraph_attribute_type_t type;
    if (igraph_i_attribute_has_attr(graph, IGRAPH_ATTRIBUTE_EDGE, 
				    estrnames[i])) {
      igraph_i_attribute_gettype(graph, &type, IGRAPH_ATTRIBUTE_EDGE, 
				 estrnames[i]);
      if (type==IGRAPH_ATTRIBUTE_STRING) {
	IGRAPH_CHECK(igraph_vector_push_back(&ex_stra, i));
      }
    }
  }
  
  for (i=0; !IGRAPH_EIT_END(eit); IGRAPH_EIT_NEXT(eit), i++) {
    long int edge=IGRAPH_EIT_GET(eit);
    igraph_integer_t from, to;
    igraph_edge(graph, edge, &from,  &to);
    fprintf(outstream, "%li %li", (long int) from+1, (long int) to+1);
    
    /* Weights */
    if (etypes[E_WEIGHT] == IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_i_attribute_get_numeric_edge_attr(graph, enames[E_WEIGHT],
					       igraph_ess_1(edge), &numv);
      fprintf(outstream, " %g", VECTOR(numv)[0]);
    }
    
    /* numeric parameters */
    for (j=0; j<igraph_vector_size(&ex_numa); j++) {
      int idx=VECTOR(ex_numa)[j];
      igraph_i_attribute_get_numeric_edge_attr(graph, enumnames[idx],
					       igraph_ess_1(edge), &numv);
      fprintf(outstream, " %s %g", enumnames2[idx], VECTOR(numv)[0]);
    }
    
    /* string parameters */
    for (j=0; j<igraph_vector_size(&ex_stra); j++) {
      int idx=VECTOR(ex_stra)[j];
      char *s;
      igraph_i_attribute_get_string_edge_attr(graph, estrnames[idx],
					      igraph_ess_1(edge), &strv);
      igraph_strvector_get(&strv, 0, &s);
      fprintf(outstream, " %s \"%s\"", estrnames2[idx], s);
    }

    /* trailing newline */
    fprintf(outstream, "%s", newline);
  }

  igraph_eit_destroy(&eit);
  igraph_es_destroy(&es);
  igraph_vector_destroy(&ex_numa);
  igraph_vector_destroy(&ex_stra);
  igraph_vector_destroy(&vx_numa);
  igraph_vector_destroy(&vx_stra);
  igraph_strvector_destroy(&strv);
  igraph_vector_destroy(&numv);
  IGRAPH_FINALLY_CLEAN(8);
  return 0;
}

int igraph_write_graph_dimacs(const igraph_t *graph, FILE *outstream,
			      long int source, long int target,
			      const igraph_vector_t *capacity) {
  
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  igraph_eit_t it;
  long int i=0;
  int ret;

  if (igraph_vector_size(capacity) != no_of_edges) {
    IGRAPH_ERROR("invalid capacity vector length", IGRAPH_EINVAL);
  }
  
  IGRAPH_CHECK(igraph_eit_create(graph, igraph_ess_all(IGRAPH_EDGEORDER_ID),
				 &it));
  IGRAPH_FINALLY(igraph_eit_destroy, &it);
  
  ret=fprintf(outstream, 
	      "c created by igraph\np max %li %li\nn %li s\nn %li t\n",
	      no_of_nodes, no_of_edges, source+1, target+1);
  if (ret < 0) {
    IGRAPH_ERROR("Write error", IGRAPH_EFILE);
  }
  

  while (!IGRAPH_EIT_END(it)) {
    igraph_integer_t from, to, cap;
    igraph_edge(graph, IGRAPH_EIT_GET(it), &from, &to);
    cap=VECTOR(*capacity)[i++];
    ret=fprintf(outstream, "a %li %li %g\n",
		(long int) from+1, (long int) to+1, cap);
    if (ret < 0) {
      IGRAPH_ERROR("Write error", IGRAPH_EFILE);
    }
    IGRAPH_EIT_NEXT(it);
  }
  
  igraph_eit_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;  
}
