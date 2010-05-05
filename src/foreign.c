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
#include "igraph_math.h"
#include "gml_tree.h"
#include "memory.h"

#include <ctype.h>		/* isspace */
#include <string.h>
#include <time.h>

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
extern char *igraph_i_ncol_errmsg;

/**
 * \ingroup loadsave
 * \function igraph_read_graph_ncol
 * \brief Reads a <code>.ncol</code> file used by LGL.
 *
 * Also useful for creating graphs from \quote named\endquote (and
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
  igraph_i_ncol_errmsg=0;

  if (igraph_ncol_yyparse()) {
    if (igraph_i_ncol_errmsg) {
      IGRAPH_ERROR(igraph_i_ncol_errmsg, IGRAPH_PARSEERROR);
    } else {
      IGRAPH_ERROR("Cannot read NCOL file", IGRAPH_PARSEERROR);
    }
  }

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
extern char *igraph_i_lgl_errmsg;

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
  igraph_i_lgl_errmsg=0;

  if (igraph_lgl_yyparse()) {
    if (igraph_i_lgl_errmsg) {
      IGRAPH_ERROR(igraph_i_lgl_errmsg, IGRAPH_PARSEERROR);
    } else {
      IGRAPH_ERROR("Cannot read LGL file", IGRAPH_PARSEERROR);
    }
  }

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
long int igraph_pajek_vcount;
long int igraph_pajek_actfrom, igraph_pajek_actto;
int igraph_pajek_mode=0;	/* 0 - general, 1 - vertex, 2 - edge */
igraph_trie_t *igraph_i_pajek_vertex_attribute_names;
igraph_vector_ptr_t *igraph_i_pajek_vertex_attributes;
igraph_trie_t *igraph_i_pajek_edge_attribute_names;
igraph_vector_ptr_t *igraph_i_pajek_edge_attributes;
long int igraph_i_pajek_vertexid=0;
long int igraph_i_pajek_actvertex=0;
long int igraph_i_pajek_actedge=0;
extern char *igraph_i_pajek_errmsg;

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
  igraph_pajek_vcount=-1;
  igraph_i_pajek_vertexid=0;
  igraph_i_pajek_vertex_attribute_names=&vattrnames;
  igraph_i_pajek_vertex_attributes=&vattrs;
  igraph_i_pajek_edge_attribute_names=&eattrnames;
  igraph_i_pajek_edge_attributes=&eattrs;
  igraph_i_pajek_actedge=0;
  igraph_pajek_mylineno=1;
  igraph_i_pajek_eof=0;
  igraph_i_pajek_errmsg=0;

  if (igraph_pajek_yyparse()) {
    if (igraph_i_pajek_errmsg) {
      IGRAPH_ERROR(igraph_i_pajek_errmsg, IGRAPH_PARSEERROR);
    } else {
      IGRAPH_ERROR("Cannot read Pajek file", IGRAPH_PARSEERROR);
    }
  }

  if (igraph_pajek_vcount < 0)
    IGRAPH_ERROR("invalid vertex count in Pajek file", IGRAPH_EINVAL);

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
      igraph_Free(vec);
    } else if (rec->type==IGRAPH_ATTRIBUTE_STRING) {
      igraph_strvector_t *strvec=(igraph_strvector_t *)rec->value;
      igraph_strvector_destroy(strvec);
      igraph_Free(strvec);
    }
    igraph_free( (char*)(rec->name));
    igraph_Free(rec);
  }

  for (i=0; i<igraph_vector_ptr_size(&eattrs); i++) {
    igraph_i_attribute_record_t *rec=VECTOR(eattrs)[i];
    if (rec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_vector_t *vec=(igraph_vector_t*) rec->value;
      igraph_vector_destroy(vec);
      igraph_Free(vec);
    } else if (rec->type==IGRAPH_ATTRIBUTE_STRING) {
      igraph_strvector_t *strvec=(igraph_strvector_t *)rec->value;
      igraph_strvector_destroy(strvec);
      igraph_Free(strvec);
    }
    igraph_free( (char*)(rec->name));
    igraph_Free(rec);
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
 * \brief Read a graph in DIMACS format.
 * 
 * This function reads the DIMACS file format, more specifically the 
 * version for network flow problems, see the files at
 * ftp://dimacs.rutgers.edu/pub/netflow/general-info/
 * 
 * </para><para>
 * This is a line-oriented text file (ASCII) format. The first
 * character of each line defines the type of the line. If the first
 * character is <code>c</code> the line is a comment line and it is
 * ignored. There is one problem line (<code>p</code> in the file, it
 * must appear before any node and arc descriptor lines. The problem
 * line has three fields separated by spaces: the problem type
 * (<code>min</code>, <code>max</code> or <code>asn</code>), the
 * number of vertices and number of edges in the graph.
 * Exactly two node identification lines are expected
 * (<code>n</code>), one for the source, one for the target vertex.
 * These have two fields: the id of the vertex and the type of the
 * vertex, either <code>s</code> (=source) or <code>t</code>
 * (=target). Arc lines start with <code>a</code> and have three
 * fields: the source vertex, the target vertex and the edge capacity.
 * 
 * </para><para>
 * Vertex ids are numbered from 1.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The file to read from.
 * \param source Pointer to an integer, the id of the source node will
 *    be stored here. (The igraph vertex id, which is one less than
 *    the actual number in the file.) It is ignored if
 *    <code>NULL</code>.
 * \param target Pointer to an integer, the (igraph) id of the target
 *    node will be stored here. It is ignored if <code>NULL</code>.
 * \param capacity Pointer to an initialized vector, the capacity of
 *    the edges will be stored here if not <code>NULL</code>.
 * \param directed Boolean, whether to create a directed graph.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|+c), the number of vertices plus the
 * number of edges, plus the size of the file in characters.
 * 
 * \sa \ref igraph_write_graph_dimacs()
 */

int igraph_read_graph_dimacs(igraph_t *graph, FILE *instream,
			     igraph_strvector_t *problem,
			     igraph_vector_t *label,
			     igraph_integer_t *source, 
			     igraph_integer_t *target,
			     igraph_vector_t *capacity,			     
			     igraph_bool_t directed) {
  
  igraph_vector_t edges;
  long int no_of_nodes=-1;
  long int no_of_edges=-1;
  long int tsource=-1;
  long int ttarget=-1;
  char prob[21];
  char c;      
  int problem_type=0;

#define PROBLEM_EDGE  1
#define PROBLEM_MAX   2  

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
      long int tmp, tmp2;
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
      read=fscanf(instream, "%20s %li %li", prob, 
		  &no_of_nodes, &no_of_edges);
      if (read != 3) {
	IGRAPH_ERROR("reading dimacs file failed", IGRAPH_PARSEERROR);
      }
      if (!strcmp(prob, "edge")) {
	/* edge list */
	problem_type=PROBLEM_EDGE;
	if (label) {
	  long int i;
	  IGRAPH_CHECK(igraph_vector_resize(label, no_of_nodes));
	  for (i=0; i<no_of_nodes; i++) {
	    VECTOR(*label)[i]=i+1;
	  }
	}
      } else if (!strcmp(prob, "max")) {
	/* maximum flow problem */
	problem_type=PROBLEM_MAX;
	if (capacity) {
	  IGRAPH_CHECK(igraph_vector_reserve(capacity, no_of_edges));
	}
      } else {
	IGRAPH_ERROR("Unknown problem type, should be 'edge' or 'max'",
		     IGRAPH_PARSEERROR);
      }
      if (problem) {
	igraph_strvector_clear(problem);
	IGRAPH_CHECK(igraph_strvector_add(problem, prob));
      }
      IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges*2));
      break;

    case 'n':
      /* for MAX this is either the source or target vertex,
	 for EDGE this is a vertex label */
      if (problem_type == PROBLEM_MAX) {
	str[0]='x';
	read=fscanf(instream, "%li %1s", &tmp, str);
	if (str[0]=='s') {
	  if (tsource != -1) {
	    IGRAPH_ERROR("reading dimacsfile: multiple source vertex line", 
			 IGRAPH_PARSEERROR);
	  } else {
	    tsource=tmp;
	  }
	} else if (str[0]=='t') {
	  if (ttarget != -1) {
	    IGRAPH_ERROR("reading dimacsfile: multiple target vertex line", 
			 IGRAPH_PARSEERROR);
	  } else {
	    ttarget=tmp;
	  }
	} else {
	  IGRAPH_ERROR("invalid node descriptor line in dimacs file",
		       IGRAPH_PARSEERROR);
	}
      } else {
	read=fscanf(instream, "%li %li", &tmp, &tmp2);
	if (label) {
	  VECTOR(*label)[tmp]=tmp2;
	}
      }
      
      break;
      
    case 'a':
      /* This is valid only for MAX, a weighted edge */
      if (problem_type != PROBLEM_MAX) { 
	IGRAPH_ERROR("'a' lines are allowed only in MAX problem files", 
		     IGRAPH_PARSEERROR);
      }
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
      
    case 'e':
      /* Edge line, only in EDGE */
      if (problem_type != PROBLEM_EDGE) {
	IGRAPH_ERROR("'e' lines are allowed only in EDGE problem files",
		     IGRAPH_PARSEERROR);
      }
      read=fscanf(instream, "%li %li", &from, &to);
      if (read != 2) {
	IGRAPH_ERROR("reading dimacs file", IGRAPH_PARSEERROR);
      }
      IGRAPH_CHECK(igraph_vector_push_back(&edges, from-1));
      IGRAPH_CHECK(igraph_vector_push_back(&edges, to-1));
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

/**
 * \function igraph_read_graph_graphdb
 * \brief Read a graph in the binary graph database format.
 * 
 * This is a binary format, used in the graph database
 * for isomorphism testing (http://amalfi.dis.unina.it/graph/)
 * From the graph database homepage
 * (http://amalfi.dis.unina.it/graph/db/doc/graphdbat-2.html):
 * </para>
 * 
 * \blockquote <para>
 * The graphs are stored in a compact binary format, one graph per
 * file. The file is composed of 16 bit words, which are represented
 * using the so-called little-endian convention, i.e. the least
 * significant byte of the word is stored first.</para>
 * 
 * <para>
 * Then, for each node, the file contains the list of edges coming
 * out of the node itself. The list is represented by a word encoding
 * its length, followed by a word for each edge, representing the
 * destination node of the edge. Node numeration is 0-based, so the
 * first node of the graph has index 0.</para> \endblockquote
 * 
 * <para>
 * Only unlabelled graphs are implemented.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read from.
 * \param directed Logical scalar, whether to create a directed graph.
 * \return Error code.
 * 
 * Time complexity: O(|V|+|E|), the number of vertices plus the 
 * number of edges.
 */

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

extern int igraph_gml_yyparse(void);
extern FILE *igraph_gml_yyin;
extern int igraph_gml_eof;
extern igraph_gml_tree_t *igraph_i_gml_parsed_tree;
long int igraph_gml_mylineno;
extern char *igraph_i_gml_errmsg;

void igraph_i_gml_destroy_attrs(igraph_vector_ptr_t **ptr) {  
  long int i;
  igraph_vector_ptr_t *vec;
  for (i=0; i<3; i++) {
    long int j;
    vec=ptr[i];
    for (j=0; j<igraph_vector_ptr_size(vec); j++) {
      igraph_i_attribute_record_t *atrec=VECTOR(*vec)[j];
      if (atrec->type == IGRAPH_ATTRIBUTE_NUMERIC) {
	igraph_vector_t *value=(igraph_vector_t*)atrec->value;
	igraph_vector_destroy(value);
	igraph_Free(value);
      } else {
	igraph_strvector_t *value=(igraph_strvector_t*)atrec->value;
	igraph_strvector_destroy(value);
	igraph_Free(value);
      }
      igraph_Free(atrec->name);
      igraph_Free(atrec);
    }
    igraph_vector_ptr_destroy(vec);
  }
}

igraph_real_t igraph_i_gml_toreal(igraph_gml_tree_t *node, long int pos) {

  igraph_real_t value=0.0;
  int type=igraph_gml_tree_type(node, pos);
  
  switch (type) {
  case IGRAPH_I_GML_TREE_INTEGER:
    value=igraph_gml_tree_get_integer(node, pos);
    break;
  case IGRAPH_I_GML_TREE_REAL:
    value=igraph_gml_tree_get_real(node, pos);
    break;
  default:
    IGRAPH_ERROR("Internal error while parsing GML file", IGRAPH_FAILURE);
    break;
  }
  
  return value;
}

const char *igraph_i_gml_tostring(igraph_gml_tree_t *node, long int pos) {
  
  int type=igraph_gml_tree_type(node, pos);
  static char tmp[256];
  const char *p=tmp;
  long int i;
  igraph_real_t d;

  switch (type) {
  case IGRAPH_I_GML_TREE_INTEGER:
    i=igraph_gml_tree_get_integer(node, pos);
    snprintf(tmp, sizeof(tmp)/sizeof(char), "%li", i);
    break;
  case IGRAPH_I_GML_TREE_REAL:
    d=igraph_gml_tree_get_real(node, pos);
    snprintf(tmp, sizeof(tmp)/sizeof(char), "%g", d);
    break;
  case IGRAPH_I_GML_TREE_STRING:
    p=igraph_gml_tree_get_string(node, pos);
    break;
  default:
    break;
  }

  return p;
}

/**
 * \function igraph_read_graph_gml
 * \brief Read a graph in GML format.
 * 
 * GML is a simple textual format, see
 * http://www.infosun.fim.uni-passau.de/Graphlet/GML/ for details.
 * 
 * </para><para>
 * Although all syntactically correct GML can be parsed, 
 * we implement only a subset of this format, some attributes might be
 * ignored. Here is a list of all the differences:
 * \olist
 * \oli Only <code>node</code> and <code>edge</code> attributes are 
 *      used, and only if they have a simple type: integer, real or
 *      string. So if an attribute is an array or a record, then it is 
 *      ignored. This is also true if only some values of the
 *      attribute are complex. 
 * \oli Top level attributes except for <code>Version</code> and the
 *      first <code>graph</code> attribute are completely ignored.
 * \oli Graph attributes except for <code>node</code> and
 *      <code>edge</code> are completely ignored.
 * \oli There is no maximum line length. 
 * \oli There is no maximum keyword length.
 * \oli Character entities in strings are not interpreted.
 * \oli We allow <code>inf</code> (infinity) and <code>nan</code>
 *      (not a number) as a real number. This is case insensitive, so
 *      <code>nan</code>, <code>NaN</code> and <code>NAN</code> are equal.
 * \endolist
 * 
 * </para><para> Please contact us if you cannot live with these
 * limitations of the GML parser.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream The stream to read the GML file from. 
 * \return Error code.
 * 
 * Time complexity: should be proportional to the length of the file.
 * 
 * \sa \ref igraph_read_graph_graphml() for a more modern format, 
 * \ref igraph_write_graph_gml() for writing GML files.
 */

int igraph_read_graph_gml(igraph_t *graph, FILE *instream) {
  
  long int i, p;
  long int no_of_nodes=0, no_of_edges=0;
  igraph_trie_t trie;
  igraph_vector_t edges;
  igraph_bool_t directed=IGRAPH_UNDIRECTED;
  igraph_gml_tree_t *gtree;
  long int gidx;
  igraph_trie_t vattrnames;
  igraph_trie_t eattrnames;
  igraph_trie_t gattrnames;
  igraph_vector_ptr_t gattrs=IGRAPH_VECTOR_PTR_NULL, 
    vattrs=IGRAPH_VECTOR_PTR_NULL, eattrs=IGRAPH_VECTOR_PTR_NULL;
  igraph_vector_ptr_t *attrs[3];
  long int edgeptr=0;
  
  attrs[0]=&gattrs; attrs[1]=&vattrs; attrs[2]=&eattrs;
  
  igraph_gml_yyin=instream;
  igraph_gml_mylineno=1;
  igraph_gml_eof=0;
  igraph_i_gml_errmsg=0;

  i=igraph_gml_yyparse();
  if (i != 0) {
    if (igraph_i_gml_errmsg) {
      IGRAPH_ERROR(igraph_i_gml_errmsg, IGRAPH_PARSEERROR);
    } else {
      IGRAPH_ERROR("Cannot read GML file", IGRAPH_PARSEERROR);
    }
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  /* Check version, if present, integer and not '1' then ignored */
  i=igraph_gml_tree_find(igraph_i_gml_parsed_tree, "Version", 0);
  if (i>=0 &&
      igraph_gml_tree_type(igraph_i_gml_parsed_tree, i)==IGRAPH_I_GML_TREE_INTEGER &&
      igraph_gml_tree_get_integer(igraph_i_gml_parsed_tree, i) != 1) {
    igraph_gml_tree_destroy(igraph_i_gml_parsed_tree);
    IGRAPH_ERROR("Unknown GML version", IGRAPH_UNIMPLEMENTED);
    /* RETURN HERE!!!! */
  }
  
  /* get the graph */
  gidx=igraph_gml_tree_find(igraph_i_gml_parsed_tree, "graph", 0);
  if (gidx==-1) {
    IGRAPH_ERROR("No 'graph' object in GML file", IGRAPH_PARSEERROR);
  }
  if (igraph_gml_tree_type(igraph_i_gml_parsed_tree, gidx) !=
      IGRAPH_I_GML_TREE_TREE) {
    IGRAPH_ERROR("Invalid type for 'graph' object in GML file", IGRAPH_PARSEERROR);
  }
  gtree=igraph_gml_tree_get_tree(igraph_i_gml_parsed_tree, gidx);

  IGRAPH_FINALLY(igraph_i_gml_destroy_attrs, &attrs);
  igraph_vector_ptr_init(&gattrs, 0);
  igraph_vector_ptr_init(&vattrs, 0);
  igraph_vector_ptr_init(&eattrs, 0);

  IGRAPH_TRIE_INIT_FINALLY(&trie, 0);
  IGRAPH_TRIE_INIT_FINALLY(&vattrnames, 0);
  IGRAPH_TRIE_INIT_FINALLY(&eattrnames, 0);
  IGRAPH_TRIE_INIT_FINALLY(&gattrnames, 0);

  /* Is is directed? */
  i=igraph_gml_tree_find(gtree, "directed", 0);
  if (i>=0 && igraph_gml_tree_type(gtree, i)==IGRAPH_I_GML_TREE_INTEGER) {
    if (igraph_gml_tree_get_integer(gtree, i) == 1) {
      directed=IGRAPH_DIRECTED;
    }
  }

  /* Now we go over all objects in the graph and collect the attribute names and
     types. Plus we collect node ids. We also do some checks. */
  for (i=0; i<igraph_gml_tree_length(gtree); i++) {
    long int j;
    char cname[100];
    const char *name=igraph_gml_tree_name(gtree, i);
    if (!strcmp(name, "node")) {
      igraph_gml_tree_t *node;
      igraph_bool_t hasid;
      no_of_nodes++;
      if (igraph_gml_tree_type(gtree, i) != IGRAPH_I_GML_TREE_TREE) {
	IGRAPH_ERROR("'node' is not a list", IGRAPH_PARSEERROR);
      }
      node=igraph_gml_tree_get_tree(gtree, i);
      hasid=0;
      for (j=0; j<igraph_gml_tree_length(node); j++) {
	const char *name=igraph_gml_tree_name(node, j);
	long int trieid, triesize=igraph_trie_size(&vattrnames);
	IGRAPH_CHECK(igraph_trie_get(&vattrnames, name, &trieid));
	if (trieid==triesize) {
	  /* new attribute */
	  igraph_i_attribute_record_t *atrec=igraph_Calloc(1, igraph_i_attribute_record_t);
	  int type=igraph_gml_tree_type(node, j);
	  if (!atrec) {
	    IGRAPH_ERROR("Cannot read GML file", IGRAPH_ENOMEM);
	  }
	  IGRAPH_CHECK(igraph_vector_ptr_push_back(&vattrs, atrec));
	  atrec->name=strdup(name);
	  if (type==IGRAPH_I_GML_TREE_INTEGER || type==IGRAPH_I_GML_TREE_REAL) {
	    atrec->type=IGRAPH_ATTRIBUTE_NUMERIC;
	  } else {
	    atrec->type=IGRAPH_ATTRIBUTE_STRING;
	  }
	} else {
	  /* already seen, should we update type? */
	  igraph_i_attribute_record_t *atrec=VECTOR(vattrs)[trieid];
	  int type1=atrec->type;
	  int type2=igraph_gml_tree_type(node, j);
	  if (type1==IGRAPH_ATTRIBUTE_NUMERIC && type2==IGRAPH_I_GML_TREE_STRING) {
	    atrec->type=IGRAPH_ATTRIBUTE_STRING;
	  }
	}
	/* check id */
	if (!hasid && !strcmp(name, "id")) {
	  long int id;
	  if (igraph_gml_tree_type(node, j) != IGRAPH_I_GML_TREE_INTEGER) {
	    IGRAPH_ERROR("Non-integer node id in GML file", IGRAPH_PARSEERROR);
	  }
	  id=igraph_gml_tree_get_integer(node, j);
	  snprintf(cname, sizeof(cname)/sizeof(char)-1, "%li", id);
	  IGRAPH_CHECK(igraph_trie_get(&trie, cname, &id));
	  hasid=1;
	}
      }
      if (!hasid) {
	IGRAPH_ERROR("Node without 'id' while parsing GML file", IGRAPH_PARSEERROR);
      }
    } else if (!strcmp(name, "edge")) {
      igraph_gml_tree_t *edge;
      igraph_bool_t has_source=0, has_target=0;
      no_of_edges++;
      if (igraph_gml_tree_type(gtree, i) != IGRAPH_I_GML_TREE_TREE) {
	IGRAPH_ERROR("'edge' is not a list", IGRAPH_PARSEERROR);
      }
      edge=igraph_gml_tree_get_tree(gtree, i);
      has_source=has_target=0;
      for (j=0; j<igraph_gml_tree_length(edge); j++) {
	const char *name=igraph_gml_tree_name(edge, j);
	if (!strcmp(name, "source")) {
	  has_source=1;
	  if (igraph_gml_tree_type(edge, j) != IGRAPH_I_GML_TREE_INTEGER) {
	    IGRAPH_ERROR("Non-integer 'source' for an edge in GML file",
			 IGRAPH_PARSEERROR);
	  }
	} else if (!strcmp(name, "target")) {
	  has_target=1;
	  if (igraph_gml_tree_type(edge, j) != IGRAPH_I_GML_TREE_INTEGER) {
	    IGRAPH_ERROR("Non-integer 'source' for an edge in GML file",
			 IGRAPH_PARSEERROR);
	  }
	} else {
	  long int trieid, triesize=igraph_trie_size(&eattrnames);
	  IGRAPH_CHECK(igraph_trie_get(&eattrnames, name, &trieid));
	  if (trieid==triesize) {
	    /* new attribute */
	    igraph_i_attribute_record_t *atrec=igraph_Calloc(1, igraph_i_attribute_record_t);
	    int type=igraph_gml_tree_type(edge, j);
	    if (!atrec) {
	      IGRAPH_ERROR("Cannot read GML file", IGRAPH_ENOMEM);
	    }
	    IGRAPH_CHECK(igraph_vector_ptr_push_back(&eattrs, atrec));
	    atrec->name=strdup(name);
	    if (type==IGRAPH_I_GML_TREE_INTEGER || type==IGRAPH_I_GML_TREE_REAL) {
	      atrec->type=IGRAPH_ATTRIBUTE_NUMERIC;
	    } else {
	      atrec->type=IGRAPH_ATTRIBUTE_STRING;
	    }
	  } else {
	    /* already seen, should we update type? */
	    igraph_i_attribute_record_t *atrec=VECTOR(eattrs)[trieid];
	    int type1=atrec->type;
	    int type2=igraph_gml_tree_type(edge, j);
	    if (type1==IGRAPH_ATTRIBUTE_NUMERIC && type2==IGRAPH_I_GML_TREE_STRING) {
	      atrec->type=IGRAPH_ATTRIBUTE_STRING;
	    }
	  }
	}
      }	/* for */
      if (!has_source) {
	IGRAPH_ERROR("No 'source' for edge in GML file", IGRAPH_PARSEERROR);
      }
      if (!has_target) {
	IGRAPH_ERROR("No 'target' for edge in GML file", IGRAPH_PARSEERROR);
      }
    } else {
      /* anything to do? Maybe add as graph attribute.... */
    }
  }

  /* check vertex id uniqueness */
  if (igraph_trie_size(&trie) != no_of_nodes) {
    IGRAPH_ERROR("Node 'id' not unique", IGRAPH_PARSEERROR);
  }
  
  /* now we allocate the vectors and strvectors for the attributes */
  for (i=0; i<igraph_vector_ptr_size(&vattrs); i++) {
    igraph_i_attribute_record_t *atrec=VECTOR(vattrs)[i];
    int type=atrec->type;
    if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_vector_t *p=igraph_Calloc(1, igraph_vector_t);
      atrec->value=p;
      IGRAPH_CHECK(igraph_vector_init(p, no_of_nodes));
    } else if (type == IGRAPH_ATTRIBUTE_STRING) {
      igraph_strvector_t *p=igraph_Calloc(1, igraph_strvector_t);
      atrec->value=p;
      IGRAPH_CHECK(igraph_strvector_init(p, no_of_nodes));
    } else {
      IGRAPH_WARNING("A composite attribute ignored");
    }
  }

  for (i=0; i<igraph_vector_ptr_size(&eattrs); i++) {
    igraph_i_attribute_record_t *atrec=VECTOR(eattrs)[i];
    int type=atrec->type;
    if (type == IGRAPH_ATTRIBUTE_NUMERIC) {
      igraph_vector_t *p=igraph_Calloc(1, igraph_vector_t);
      atrec->value=p;
      IGRAPH_CHECK(igraph_vector_init(p, no_of_edges));
    } else if (type == IGRAPH_ATTRIBUTE_STRING) {
      igraph_strvector_t *p=igraph_Calloc(1, igraph_strvector_t);
      atrec->value=p;
      IGRAPH_CHECK(igraph_strvector_init(p, no_of_edges));
    } else {
      IGRAPH_WARNING("A composite attribute ignored");
    }
  }

  /* Ok, now the edges, attributes too */
  IGRAPH_CHECK(igraph_vector_resize(&edges, no_of_edges*2));
  p=-1;
  while ( (p=igraph_gml_tree_find(gtree, "edge", p+1)) != -1) {
    igraph_gml_tree_t *edge;
    long int from, to, fromidx=0, toidx=0;
    char name[100];
    long int j;
    edge=igraph_gml_tree_get_tree(gtree, p);
    for (j=0; j<igraph_gml_tree_length(edge); j++) {
      const char *n=igraph_gml_tree_name(edge, j);
      if (!strcmp(n, "source")) {
	fromidx=igraph_gml_tree_find(edge, "source", 0);
      } else if (!strcmp(n, "target")) {
	toidx=igraph_gml_tree_find(edge, "target", 0);
      } else {
	long int edgeid=edgeptr/2;
	long int trieidx;
	igraph_i_attribute_record_t *atrec;
	int type;
	igraph_trie_get(&eattrnames, n, &trieidx);
	atrec=VECTOR(eattrs)[trieidx];
	type=atrec->type;
	if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	  igraph_vector_t *v=(igraph_vector_t *)atrec->value;
	  VECTOR(*v)[edgeid]=igraph_i_gml_toreal(edge, j);
	} else if (type==IGRAPH_ATTRIBUTE_STRING) {
	  igraph_strvector_t *v=(igraph_strvector_t *)atrec->value;
	  const char *value=igraph_i_gml_tostring(edge, j);
	  IGRAPH_CHECK(igraph_strvector_set(v, edgeid, value));
	}
      }
    }
    from=igraph_gml_tree_get_integer(edge, fromidx);
    to=igraph_gml_tree_get_integer(edge, toidx);
    snprintf(name, sizeof(name)/sizeof(char)-1, "%li", from);
    IGRAPH_CHECK(igraph_trie_get(&trie, name, &from));
    snprintf(name, sizeof(name)/sizeof(char)-1, "%li", to);
    IGRAPH_CHECK(igraph_trie_get(&trie, name, &to));
    if (igraph_trie_size(&trie) != no_of_nodes) {
      IGRAPH_ERROR("Unkown node id found at an edge", IGRAPH_PARSEERROR);
    }
    VECTOR(edges)[edgeptr++]=from;
    VECTOR(edges)[edgeptr++]=to;
  }

  /* and add vertex attributes */
  for (i=0; i<igraph_gml_tree_length(gtree); i++) {
    const char *n;
    char name[100];
    long int j, k;
    n=igraph_gml_tree_name(gtree, i);
    if (!strcmp(n, "node")) {
      igraph_gml_tree_t *node=igraph_gml_tree_get_tree(gtree, i);
      long int iidx=igraph_gml_tree_find(node, "id", 0);
      long int id=igraph_gml_tree_get_integer(node, iidx);
      snprintf(name, sizeof(name)/sizeof(char)-1, "%li", id);
      igraph_trie_get(&trie, name, &id);
      for (j=0; j<igraph_gml_tree_length(node); j++) {
	const char *aname=igraph_gml_tree_name(node, j);
	igraph_i_attribute_record_t *atrec;
	int type;
	igraph_trie_get(&vattrnames, aname, &k);
	atrec=VECTOR(vattrs)[k];
	type=atrec->type;
	if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	  igraph_vector_t *v=(igraph_vector_t *)atrec->value;
	  VECTOR(*v)[id]=igraph_i_gml_toreal(node, j);
	} else if (type==IGRAPH_ATTRIBUTE_STRING) {
	  igraph_strvector_t *v=(igraph_strvector_t *)atrec->value;
	  const char *value=igraph_i_gml_tostring(node, j);
	  IGRAPH_CHECK(igraph_strvector_set(v, id, value));
	}
      }
    }
  }
  
  igraph_gml_tree_destroy(igraph_i_gml_parsed_tree);
  
  igraph_trie_destroy(&trie);
  igraph_trie_destroy(&gattrnames);
  igraph_trie_destroy(&vattrnames);
  igraph_trie_destroy(&eattrnames);
  IGRAPH_FINALLY_CLEAN(4);

  IGRAPH_CHECK(igraph_empty_attrs(graph, 0, directed, 0)); /* TODO */
  IGRAPH_CHECK(igraph_add_vertices(graph, no_of_nodes, &vattrs));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, &eattrs));

  igraph_i_gml_destroy_attrs(attrs);
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(2);
  
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
    igraph_vector_t wvec;
    IGRAPH_VECTOR_INIT_FINALLY(&wvec, igraph_ecount(graph));
    IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, weights, 
							 igraph_ess_all(IGRAPH_EDGEORDER_ID), 
							 &wvec));
    while (!IGRAPH_EIT_END(it)) {
      igraph_integer_t edge=IGRAPH_EIT_GET(it);
      igraph_integer_t from, to;
      int ret=0;
      igraph_edge(graph, edge, &from, &to);
      ret=fprintf(outstream, "%li %li %f\n", 
		  (long int)from, (long int)to, VECTOR(wvec)[(long int)edge]);
      if (ret<0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);      
    }
    igraph_vector_destroy(&wvec);
    IGRAPH_FINALLY_CLEAN(1);
  } else {
    /* Both names and weights */
    igraph_strvector_t nvec;
	igraph_vector_t wvec;
    IGRAPH_VECTOR_INIT_FINALLY(&wvec, igraph_ecount(graph));
    IGRAPH_CHECK(igraph_strvector_init(&nvec, igraph_vcount(graph)));
    IGRAPH_FINALLY(igraph_strvector_destroy, &nvec);
    IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, weights, 
							 igraph_ess_all(IGRAPH_EDGEORDER_ID), 
							 &wvec));
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
      ret=fprintf(outstream, "%s %s ", str1, str2);
      if (ret<0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      ret=fprintf(outstream, "%f\n", VECTOR(wvec)[(long int)edge]);
      if (ret<0) {
        IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      IGRAPH_EIT_NEXT(it);
    }
    igraph_strvector_destroy(&nvec);
    igraph_vector_destroy(&wvec);
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
							 igraph_ess_all(IGRAPH_EDGEORDER_ID),
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
							 igraph_ess_all(IGRAPH_EDGEORDER_ID),
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
 * writable. Make sure that you open the file in binary format if you use MS Windows,
 * otherwise end of line characters will be messed up. (igraph will be able
 * to read back these messed up files, but Pajek won't.)
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

  const char *newline="\x0d\x0a";
  
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
      } else {
	fprintf(outstream, " \"%li\"", i+1);
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

/**
 * \function igraph_write_graph_dimacs
 * \brief Write a graph in DIMACS format.
 * 
 * This function writes a graph to an output stream in DIMACS format,
 * describing a maximum flow problem.
 * See ftp://dimacs.rutgers.edu/pub/netflow/general-info/
 * 
 * </para><para>
 * This file format is discussed in the documentation of \ref
 * igraph_read_graph_dimacs(), see that for more information.
 * 
 * \param graph The graph to write to the stream.
 * \param outstream The stream.
 * \param source Integer, the id of the source vertex for the maximum
 *     flow. 
 * \param target Integer, the id of the target vertex.
 * \param capacity Pointer to an initialized vector containing the
 *     edge capacity values.
 * \return Error code.
 * 
 * Time complexity: O(|E|), the number of edges in the graph.
 * 
 * \sa igraph_read_graph_dimacs()
 */

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

int igraph_i_gml_convert_to_key(const char *orig, char **key) {
  static int no=1;
  char strno[50];
  long int i, len=strlen(orig), newlen=0, plen=0;
  igraph_bool_t pref=0;  
  /* do we need a prefix? */
  if (len==0 || !isalpha(orig[0])) { 
    pref=1; no++;
    snprintf(strno, sizeof(strno)-1, "igraph");
    plen=newlen=strlen(strno);
  }
  for (i=0; i<len; i++) {
    if (isalnum(orig[i])) { newlen++; }
  }
  *key=igraph_Calloc(newlen+1, char);
  if (! *key) {
    IGRAPH_ERROR("Writing GML file failed", IGRAPH_ENOMEM);
  }
  memcpy(*key, strno, plen*sizeof(char));
  for (i=0; i<len; i++) {
    if (isalnum(orig[i])) { 
      (*key)[plen++] = orig[i];
    }
  }
  (*key)[newlen]='\0';
  
  return 0;
}

#define CHECK(cmd) do { ret=cmd; if (ret<0) IGRAPH_ERROR("Write failed", IGRAPH_EFILE); } while (0)

/** 
 * \function igraph_write_graph_gml
 * \brief Write the graph to a stream in GML format 
 * 
 * GML is a quite general textual format, see 
 * http://www.infosun.fim.uni-passau.de/Graphlet/GML/ for details.
 * 
 * </para><para> The graph, vertex and edges attributes are written to the
 * file as well, if they are numeric of string.
 * 
 * </para><para> As igraph is more forgiving about attribute names, it might 
 * be neccessary to simplify the them before writing to the GML file.
 * This way we'll have a syntactically correct GML file. The following 
 * simple procedure is performed on each attribute name: first the alphanumeric 
 * characters are extracted, the others are ignored. Then if the first character
 * is not a letter then the attribute name is prefixed with <quote>igraph</quote>.
 * Note that this might result identical names for two attributes, igraph 
 * does not check this. 
 * 
 * </para><para> The <quote>id</quote> vertex attribute is treated specially. 
 * If the <parameter>id</parameter> argument is not 0 then it should be a numeric 
 * vector with the vertex ids and the <quote>id</quote> vertex attribute is 
 * ignored (if there is one). If <parameter>id</parameter> is 0 and there is a 
 * numeric <quote>id</quote> vertex attribute that is used instead. If ids
 * are not specified in either way then the regular igraph vertex ids are used.
 * 
 * </para><para> Note that whichever way vertex ids are specified, their 
 * uniqueness is not checked.
 * 
 * </para><para> If the graph has edge attributes named <quote>source</quote> 
 * or <quote>target</quote> they're silently ignored. GML uses these attributes
 * to specify the edges, so we cannot write them to the file. Rename them 
 * before calling this function if you want to preserve them.
 * \param graph The graph to write to the stream.
 * \param outstream The stream to write the file to.
 * \param id Either <code>NULL</code> or a numeric vector with the vertex ids.
 *        See details above.
 * \param creator An optional string to write to the stream in the creator line.
 *        If this is 0 then the current date and time is added.
 * \return Error code.
 * 
 * Time complexity: should be proportional to the number of characters written
 * to the file.
 * 
 * \sa \ref igraph_read_graph_gml() for reading GML files, 
 * \ref igraph_read_graph_graphml() for a more modern format.
 */

int igraph_write_graph_gml(const igraph_t *graph, FILE *outstream, 
			   const igraph_vector_t *id, const char *creator) {
  int ret;
  igraph_strvector_t gnames, vnames, enames;
  igraph_vector_t gtypes, vtypes, etypes;
  igraph_vector_t numv;
  igraph_strvector_t strv;
  long int i;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);

  igraph_vector_t v_myid;
  const igraph_vector_t *myid=id;

  time_t curtime=time(0);
  char *timestr=ctime(&curtime);
  timestr[strlen(timestr)-1]='\0'; /* nicely remove \n */
  
  CHECK(fprintf(outstream, 
		"Creator \"igraph version %s %s\"\nVersion 1\ngraph\n[\n", 
		PACKAGE_VERSION, creator ? creator : timestr));
  
  IGRAPH_STRVECTOR_INIT_FINALLY(&gnames, 0);
  IGRAPH_STRVECTOR_INIT_FINALLY(&vnames, 0);
  IGRAPH_STRVECTOR_INIT_FINALLY(&enames, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&gtypes, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vtypes, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&etypes, 0);
  IGRAPH_CHECK(igraph_i_attribute_get_info(graph, 
					   &gnames, &gtypes,
					   &vnames, &vtypes,
					   &enames, &etypes));  

  IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
  IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);

  /* Check whether there is an 'id' node attribute if the supplied is 0 */
  if (!id) {
    igraph_bool_t found=0; 
    for (i=0; i<igraph_vector_size(&vtypes); i++) {
      char *n;
      igraph_strvector_get(&vnames, i, &n);
      if (!strcmp(n, "id") && VECTOR(vtypes)[i]==IGRAPH_ATTRIBUTE_NUMERIC) { 
	found=1; break; 
      }
    }
    if (found) {
      IGRAPH_VECTOR_INIT_FINALLY(&v_myid, no_of_nodes);
      IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, "id", 
							      igraph_vss_all(),
							      &v_myid));
      myid=&v_myid;
    }
  }      

  /* directedness */
  CHECK(fprintf(outstream, "  directed %i\n", igraph_is_directed(graph) ? 1 : 0));

  /* Graph attributes first */
  for (i=0; i<igraph_vector_size(&gtypes); i++) {
    char *name, *newname;
    igraph_strvector_get(&gnames, i, &name);
    IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
    if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
      IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
      CHECK(fprintf(outstream, "  %s %g\n", newname, VECTOR(numv)[0]));
      igraph_Free(newname);
    } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
      char *s;
      IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
      igraph_strvector_get(&strv, 0, &s);
      CHECK(fprintf(outstream, "  %s \"%s\"\n", newname, s));
      igraph_Free(newname);
    } else {
      igraph_Free(newname);
      IGRAPH_WARNING("A non-numeric, non-string graph attribute ignored");
    }
  } 
  
  /* Now come the vertices */
  for (i=0; i<no_of_nodes; i++) {
    long int j;
    CHECK(fprintf(outstream, "  node\n  [\n"));
    /* id */
    CHECK(fprintf(outstream, "    id %li\n", myid ? (long int)VECTOR(*myid)[i] : i));
    /* other attributes */
    for (j=0; j<igraph_vector_size(&vtypes); j++) {
      int type=VECTOR(vtypes)[j];
      char *name, *newname;
      igraph_strvector_get(&vnames, j, &name);
      if (!strcmp(name, "id")) { continue; }	
      IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name, 
								igraph_vss_1(i),
								&numv));
	CHECK(fprintf(outstream, "    %s %g\n", newname, VECTOR(numv)[0]));
      } else if (type==IGRAPH_ATTRIBUTE_STRING) { 
	char *s;
	IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name,
							       igraph_vss_1(i),
							       &strv));
	igraph_strvector_get(&strv, 0, &s);
	CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
      }
      igraph_Free(newname);
    }
    CHECK(fprintf(outstream, "  ]\n"));
  }

  /* The edges too */
  for (i=0; i<no_of_edges; i++) {
    long int from=IGRAPH_FROM(graph, i);
    long int to=IGRAPH_TO(graph, i);
    long int j;
    CHECK(fprintf(outstream, "  edge\n  [\n"));
    /* source and target */
    CHECK(fprintf(outstream, "    source %li\n", 
		  myid ? (long int)VECTOR(*myid)[from] : from));
    CHECK(fprintf(outstream, "    target %li\n", 
		  myid ? (long int)VECTOR(*myid)[to] : to));

    /* other attributes */
    for (j=0; j<igraph_vector_size(&etypes); j++) {
      int type=VECTOR(etypes)[j];
      char *name, *newname;
      igraph_strvector_get(&enames, j, &name);
      if (!strcmp(name, "source") || !strcmp(name, "target")) { continue; }	
      IGRAPH_CHECK(igraph_i_gml_convert_to_key(name, &newname));
      if (type==IGRAPH_ATTRIBUTE_NUMERIC) {
	IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name, 
							      igraph_ess_1(i),
							      &numv));
	CHECK(fprintf(outstream, "    %s %g\n", newname, VECTOR(numv)[0]));
      } else if (type==IGRAPH_ATTRIBUTE_STRING) { 
	char *s;
	IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name,
							     igraph_ess_1(i),
							     &strv));
	igraph_strvector_get(&strv, 0, &s);
	CHECK(fprintf(outstream, "    %s \"%s\"\n", newname, s));
      }
      igraph_Free(newname);
    }
    CHECK(fprintf(outstream, "  ]\n"));
  }

  CHECK(fprintf(outstream, "]\n"));

  if (&v_myid == myid) { 
    igraph_vector_destroy(&v_myid);
    IGRAPH_FINALLY_CLEAN(1);
  }

  igraph_strvector_destroy(&strv);
  igraph_vector_destroy(&numv);
  igraph_vector_destroy(&etypes);
  igraph_vector_destroy(&vtypes);
  igraph_vector_destroy(&gtypes);
  igraph_strvector_destroy(&enames);
  igraph_strvector_destroy(&vnames);
  igraph_strvector_destroy(&gnames);
  IGRAPH_FINALLY_CLEAN(8);
  
  return 0;
}

int igraph_i_dot_escape(const char *orig, char **result) {
  /* do we have to escape the string at all? */
  long int i, j, len=strlen(orig), newlen=0;
  igraph_bool_t need_quote=0, is_number=1;
  for (i=0; i<len; i++) {
	if (isdigit(orig[i])) { newlen++; }
	else if (orig[i] == '-' && i==0) { newlen++; }
	else if (orig[i] == '.') {
	  if (is_number) { newlen++; }
	  else { need_quote=1; newlen++; }
	} else if (orig[i] == '_') {
	  is_number=0; newlen++;
	} else if (orig[i] == '\\') {
	  need_quote=1; is_number=0; newlen+=2; /* will be escaped */
	} else if (orig[i] == '"') {
	  need_quote=1; is_number=0; newlen+=2; /* will be escaped */
	} else if (isalpha(orig[i])) {
	  is_number=0; newlen++;
	} else {
	  is_number=0; need_quote=1; newlen++;
	}
  }
  if (is_number && orig[len-1] == '.') is_number=0;
  if (!is_number && isdigit(orig[0])) need_quote=1;

  if (is_number || !need_quote) {
	*result=strdup(orig);
	if (!*result) IGRAPH_ERROR("Writing DOT file failed", IGRAPH_ENOMEM);
  } else {
	*result=igraph_Calloc(newlen+3, char);
	(*result)[0]='"';
	(*result)[newlen+1]='"';
	(*result)[newlen+2]='\0';
	for (i=0, j=1; i<len; i++) {
	  if (orig[i] == '\\' || orig[i] == '"') (*result)[j++] = '\\';
	  (*result)[j++] = orig[i];
	}
  }
  
  return 0;
}

/**
 * \function igraph_write_graph_dot
 * \brief Write the graph to a stream in DOT format
 *
 * DOT is the format used by the widely known GraphViz software, see
 * http://www.graphviz.org for details. The grammar of the DOT format
 * can be found here: http://www.graphviz.org/doc/info/lang.html
 *
 * </para><para>This is only a preliminary implementation, only the vertices
 * and the edges are written but not the attributes or any visualization
 * information.
 *
 * \param graph The graph to write to the stream.
 * \param outstream The stream to write the file to.
 *
 * Time complexity: should be proportional to the number of characters written
 * to the file.
 * 
 * \sa \ref igraph_write_graph_graphml() for a more modern format.
 */
int igraph_write_graph_dot(const igraph_t *graph, FILE* outstream) {
  int ret;
  long int i, j;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  char edgeop[3];
  igraph_strvector_t gnames, vnames, enames;
  igraph_vector_t gtypes, vtypes, etypes;
  igraph_vector_t numv;
  igraph_strvector_t strv;

  IGRAPH_STRVECTOR_INIT_FINALLY(&gnames, 0);
  IGRAPH_STRVECTOR_INIT_FINALLY(&vnames, 0);
  IGRAPH_STRVECTOR_INIT_FINALLY(&enames, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&gtypes, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&vtypes, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&etypes, 0);
  IGRAPH_CHECK(igraph_i_attribute_get_info(graph, 
					   &gnames, &gtypes,
					   &vnames, &vtypes,
					   &enames, &etypes));  

  IGRAPH_VECTOR_INIT_FINALLY(&numv, 1);
  IGRAPH_STRVECTOR_INIT_FINALLY(&strv, 1);

  CHECK(fprintf(outstream, "/* Created by igraph %s */\n",
	PACKAGE_VERSION));

  if (igraph_is_directed(graph)) {
	CHECK(fprintf(outstream, "digraph {\n"));
	strcpy(edgeop, "->");
  } else {
	CHECK(fprintf(outstream, "graph {\n"));
	strcpy(edgeop, "--");
  }

  /* Write the graph attributes */
  if (igraph_vector_size(&gtypes)>0) {
	CHECK(fprintf(outstream, "  graph [\n"));
	for (i=0; i<igraph_vector_size(&gtypes); i++) {
	  char *name, *newname;
	  igraph_strvector_get(&gnames, i, &name);
	  IGRAPH_CHECK(igraph_i_dot_escape(name, &newname));
	  if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_NUMERIC) {
		IGRAPH_CHECK(igraph_i_attribute_get_numeric_graph_attr(graph, name, &numv));
		if (VECTOR(numv)[0] == (long)VECTOR(numv)[0])
		  CHECK(fprintf(outstream, "    %s=%ld\n", newname, (long)VECTOR(numv)[0]));
		else
		  CHECK(fprintf(outstream, "    %s=%f\n", newname, VECTOR(numv)[0]));
		igraph_Free(newname);
	  } else if (VECTOR(gtypes)[i] == IGRAPH_ATTRIBUTE_STRING) {
		char *s, *news;
		IGRAPH_CHECK(igraph_i_attribute_get_string_graph_attr(graph, name, &strv));
		igraph_strvector_get(&strv, 0, &s);
		IGRAPH_CHECK(igraph_i_dot_escape(s, &news));
		CHECK(fprintf(outstream, "    %s=%s\n", newname, news));
		igraph_Free(newname);
		igraph_Free(news);
	  } else {
		IGRAPH_WARNING("A non-numeric, non-string graph attribute ignored");
	  }
	}
	CHECK(fprintf(outstream, "  ];\n"));
  }

  /* Write the vertices */
  if (igraph_vector_size(&vtypes) > 0) {
	for (i=0; i<no_of_nodes; i++) {
	  CHECK(fprintf(outstream, "  %ld [\n", i));
	  for (j=0; j<igraph_vector_size(&vtypes); j++) {
		char *name, *newname;
		igraph_strvector_get(&vnames, j, &name);
		IGRAPH_CHECK(igraph_i_dot_escape(name, &newname));
		if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
		  IGRAPH_CHECK(igraph_i_attribute_get_numeric_vertex_attr(graph, name, igraph_vss_1(i), &numv));
		  if (VECTOR(numv)[0] == (long)VECTOR(numv)[0])
			CHECK(fprintf(outstream, "    %s=%ld\n", newname, (long)VECTOR(numv)[0]));
		  else
			CHECK(fprintf(outstream, "    %s=%f\n", newname, VECTOR(numv)[0]));
		  igraph_Free(newname);
		} else if (VECTOR(vtypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
		  char *s, *news;
		  IGRAPH_CHECK(igraph_i_attribute_get_string_vertex_attr(graph, name, igraph_vss_1(i), &strv));
		  igraph_strvector_get(&strv, 0, &s);
		  IGRAPH_CHECK(igraph_i_dot_escape(s, &news));
		  CHECK(fprintf(outstream, "    %s=%s\n", newname, news));
		  igraph_Free(newname);
		  igraph_Free(news);
		} else {
		  IGRAPH_WARNING("A non-numeric, non-string graph attribute ignored");
		}
	  }
	  CHECK(fprintf(outstream, "  ];\n"));
	}
  } else {
	for (i=0; i<no_of_nodes; i++)
	  CHECK(fprintf(outstream, "  %ld;\n", i));
  }
  CHECK(fprintf(outstream, "\n"));

  /* Write the edges */
  if (igraph_vector_size(&etypes) > 0) {
	for (i=0; i<no_of_edges; i++) {
	  long int from=IGRAPH_FROM(graph, i);
	  long int to=IGRAPH_TO(graph, i);
	  CHECK(fprintf(outstream, "  %ld %s %ld [\n", from, edgeop, to));
	  for (j=0; j<igraph_vector_size(&etypes); j++) {
		char *name, *newname;
		igraph_strvector_get(&enames, j, &name);
		IGRAPH_CHECK(igraph_i_dot_escape(name, &newname));
		if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_NUMERIC) {
		  IGRAPH_CHECK(igraph_i_attribute_get_numeric_edge_attr(graph, name, igraph_ess_1(i), &numv));
		  if (VECTOR(numv)[0] == (long)VECTOR(numv)[0])
			CHECK(fprintf(outstream, "    %s=%ld\n", newname, (long)VECTOR(numv)[0]));
		  else
			CHECK(fprintf(outstream, "    %s=%f\n", newname, VECTOR(numv)[0]));
		  igraph_Free(newname);
		} else if (VECTOR(etypes)[j] == IGRAPH_ATTRIBUTE_STRING) {
		  char *s, *news;
		  IGRAPH_CHECK(igraph_i_attribute_get_string_edge_attr(graph, name, igraph_ess_1(i), &strv));
		  igraph_strvector_get(&strv, 0, &s);
		  IGRAPH_CHECK(igraph_i_dot_escape(s, &news));
		  CHECK(fprintf(outstream, "    %s=%s\n", newname, news));
		  igraph_Free(newname);
		  igraph_Free(news);
		} else {
		  IGRAPH_WARNING("A non-numeric, non-string graph attribute ignored");
		}
	  }
	  CHECK(fprintf(outstream, "  ];\n"));
	}
  } else {
	for (i=0; i<no_of_edges; i++) {
	  long int from=IGRAPH_FROM(graph, i);
	  long int to=IGRAPH_TO(graph, i);
	  CHECK(fprintf(outstream, "  %ld %s %ld;\n", from, edgeop, to));
	}
  }
  CHECK(fprintf(outstream, "}\n"));
  
  igraph_strvector_destroy(&strv);
  igraph_vector_destroy(&numv);
  igraph_vector_destroy(&etypes);
  igraph_vector_destroy(&vtypes);
  igraph_vector_destroy(&gtypes);
  igraph_strvector_destroy(&enames);
  igraph_strvector_destroy(&vnames);
  igraph_strvector_destroy(&gnames);
  IGRAPH_FINALLY_CLEAN(8);
  
  return 0;
}

#undef CHECK
