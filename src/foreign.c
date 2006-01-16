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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"

#include <ctype.h>		/* isspace */

/**
 * \section about_loadsave 
 * 
 * <para>These functions can write a graph to a file, or read a graph
 * from a file.</para>
 * 
 * <para>Note that as &igraph; uses the traditional C streams it is
 * possible to read/write files from/to memory, at least on GNU
 * operating systems supporting <quote>non-standard</quote> streams.</para>
 */

/**
 * \ingroup loadsave
 * \function igraph_read_graph_edgelist
 * \brief Reads an edge list from a file and creates a graph
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
 *         <constant>IGRAPH_PARSEERROR</constant>: if there is a
 *         problem reading the file, or the file is syntactically
 *         incorrect. 
 * 
 * Time complexity: O(|V|+|E|), the
 * number of vertices plus the number of edges. It is assumed that
 * reading an integer requires O(1)
 * time. 
 */

int igraph_read_graph_edgelist(igraph_t *graph, FILE *instream, 
			       integer_t n, bool_t directed) {

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
 * \brief Reads a <literal>.ncol</literal> file used by LGL, also
 * useful for creating graphs from <quote>named</quote> (and
 * optionally weighted) edge lists. 
 * 
 * This format is used by the Large Graph Layout program
 * (<ulink url="http://bioinformatics.icmb.utexas.edu/lgl/">http://bioinformatics.icmb.utexas.edu/lgl/</ulink>), and it is simply a
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
 * this is however not checked here, as &igraph; is happy with
 * these.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream Pointer to a stream, it should be readable.
 * \param names Logical value, if TRUE the symbolic names of the
 *        vertices will be added to the graph as a vertex attribute
 *        called <quote>name</quote>.
 * \param weights Logical value, if TRUE the weights of the
 *        edges is added to the graph as an edge attribute called
 *        <quote>weight</quote>.
 * \return Error code:
 *         <constant>IGRAPH_PARSEERROR</constant>: if there is a
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
			  bool_t names, bool_t weights) {
  
  igraph_vector_t edges, ws;
  igraph_trie_t trie=IGRAPH_TRIE_NULL;

  IGRAPH_TRIE_INIT_FINALLY(&trie, names);
  IGRAPH_VECTOR_INIT_FINALLY(&ws, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);

  igraph_ncol_vector=&edges;
  igraph_ncol_weights=&ws;
  igraph_ncol_trie=&trie;
  igraph_ncol_yyin=instream;

  igraph_ncol_yyparse();

  IGRAPH_CHECK(igraph_create(graph, &edges, 0, 0));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  if (weights) {
    /* TODO: we cannot really assume that the edge ids are the same
       as the order they were added. TODO: add edges together with
       attributes */
    long int i;
    IGRAPH_CHECK(igraph_add_edge_attribute(graph, "weight", 
					   IGRAPH_ATTRIBUTE_NUM));
    for (i=0; i<igraph_ecount(graph); i++) {
      IGRAPH_CHECK(igraph_set_edge_attribute(graph, "weight", i, 
					     &VECTOR(ws)[i])); 
    }
  }
  igraph_vector_destroy(&ws);
  IGRAPH_FINALLY_CLEAN(1);
  
  if (names) {
    long int i;
    IGRAPH_CHECK(igraph_add_vertex_attribute(graph, "name", 
					     IGRAPH_ATTRIBUTE_STR));
    for (i=0; i<igraph_vcount(graph); i++) {
      char *str;
      igraph_trie_idx(&trie, i, &str);
      IGRAPH_CHECK(igraph_set_vertex_attribute(graph, "name", i, str));
    }
  }
  igraph_trie_destroy(&trie);
  IGRAPH_FINALLY_CLEAN(1);

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
 * \brief Reads a graph from an <literal>.lgl</literal> file
 * 
 * The <literal>.lgl</literal> format is used by the Large Graph
 * Layout visualization software
 * (<ulink url="http://bioinformatics.icmb.utexas.edu/lgl/">http://bioinformatics.icmb.utexas.edu/lgl/</ulink>), it can 
 * describe undirected optionally weighted graphs. From the LGL
 * manual: 
 * 
 * <blockquote><para>The second format is the LGL file format
 * (<literal>.lgl</literal> file 
 * suffix). This is yet another graph file format that tries to be as
 * stingy as possible with space, yet keeping the edge file in a human
 * readable (not binary) format. The format itself is like the
 * following:
 * <programlisting>
 * # vertex1name
 * vertex2name [optionalWeight]
 * vertex3name [optionalWeight]
 * </programlisting>
 * Here, the first vertex of an edge is preceded with a pound sign
 * '#'.  Then each vertex that shares an edge with that vertex is
 * listed one per line on subsequent lines.</para></blockquote>
 * 
 * LGL cannot handle loop and multiple edges or directed graphs, but
 * in &igraph; it is not an error to have multiple and loop edges.
 * \param graph Pointer to an uninitialized graph object.
 * \param instream A stream, it should be readable.
 * \param names Logical value, if TRUE the symbolic names of the
 *        vertices will be added to the graph as a vertex attribute
 *        called <quote>name</quote>.
 * \param weights Logical value, if TRUE the weights of the
 *        edges is added to the graph as an edge attribute called
 *        <quote>weight</quote>.
 * \return Error code:
 *         <constant>IGRAPH_PARSEERROR</constant>: if there is a
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
			  bool_t names, bool_t weights) {

  igraph_vector_t edges=IGRAPH_VECTOR_NULL, ws=IGRAPH_VECTOR_NULL;
  igraph_trie_t trie=IGRAPH_TRIE_NULL;
  
  IGRAPH_TRIE_INIT_FINALLY(&trie, names);
  IGRAPH_VECTOR_INIT_FINALLY(&ws, 0);
  IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
  
  igraph_lgl_vector=&edges;
  igraph_lgl_weights=&ws;
  igraph_lgl_trie=&trie;
  igraph_lgl_yyin=instream;

  igraph_lgl_yyparse();
  
  IGRAPH_CHECK(igraph_create(graph, &edges, igraph_trie_size(&trie), 0));
  igraph_vector_destroy(&edges);
  IGRAPH_FINALLY_CLEAN(1);
  
  if (weights) {
    /* TODO: we cannot really assume that the edge ids are the same
       as the order they were added. TODO: add edges together with
       attributes */
    long int i;
    IGRAPH_CHECK(igraph_add_edge_attribute(graph, "weight", 
					   IGRAPH_ATTRIBUTE_NUM));
    for (i=0; i<igraph_ecount(graph); i++) {
      IGRAPH_CHECK(igraph_set_edge_attribute(graph, "weight", i, 
					     &VECTOR(ws)[i]));
    }
  }
  igraph_vector_destroy(&ws);
  IGRAPH_FINALLY_CLEAN(1);
  
  if (names) {
    long int i;
    IGRAPH_CHECK(igraph_add_vertex_attribute(graph, "name", 
					     IGRAPH_ATTRIBUTE_STR));
    for (i=0; i<igraph_vcount(graph); i++) {
      char *str;
      igraph_trie_idx(&trie, i, &str);
      IGRAPH_CHECK(igraph_set_vertex_attribute(graph, "name", i, str));
    }
  }
  igraph_trie_destroy(&trie);
  IGRAPH_FINALLY_CLEAN(1);

  return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_edgelist
 * \brief Writes the edge list of a graph to a file
 * 
 * One edge is written per line, separated by a single space.
 * For directed graphs edges are written in from, to order.
 * \param graph The graph object to write.
 * \param outstream Pointer to a stream, it should be writable.
 * \return Error code:
 *         <constant>IGRAPH_EFILE</constant> if there is an error writing the
 *         file. 
 * 
 * Time complexity: O(|E|), the
 * number of edges in the  graph. It is assumed that writing an
 * integer to the file requires O(1)
 * time. 
 */

int igraph_write_graph_edgelist(const igraph_t *graph, FILE *outstream) {

  igraph_es_t it;
  
  IGRAPH_CHECK(igraph_es_fromorder(graph, &it));
  IGRAPH_FINALLY(igraph_es_destroy, &it);

  while (!igraph_es_end(graph, &it)) {
    int ret=fprintf(outstream, "%li %li\n", 
		    (long int) igraph_es_from(graph, &it),
		    (long int) igraph_es_to(graph, &it));
    if (ret < 0) {
      IGRAPH_ERROR("Write error", IGRAPH_EFILE);
    }
    igraph_es_next(graph, &it);
  }

  igraph_es_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/** 
 * \ingroup loadsave
 * \function igraph_write_graph_ncol
 * \brief Writes the graph to a file in <literal>.ncol</literal> format
 * 
 * <literal>.ncol</literal> is a format used by LGL, see \ref
 * igraph_read_graph_ncol() for details. 
 * 
 * Note that having multiple or loop edges in an
 * <literal>.ncol</literal> file breaks the  LGL software but 
 * &igraph; does not check for this condition. 
 * \param graph The graph to write.
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param names The name of the vertex attribute, if symbolic names
 *        are written to the file. If not supply 0 here.
 * \param weights The name of the edge attribute, if they are also
 *        written to the file. If you don't want weights supply 0
 *        here.
 * \return Error code:
 *         <constant>IGRAPH_EFILE</constant> if there is an error writing the
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
  igraph_es_t it;
  
  IGRAPH_CHECK(igraph_es_fromorder(graph, &it));
  IGRAPH_FINALLY(igraph_es_destroy, &it);
  if (names==0 && weights ==0) {
    /* No names, no weights */
    while (!igraph_es_end(graph, &it)) {
      int ret=fprintf(outstream, "%li %li\n",
		      (long int) igraph_es_from(graph, &it),
		      (long int) igraph_es_to(graph, &it));
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  } else if (weights==0) {
    /* No weights, but use names */
    while (!igraph_es_end(graph, &it)) {
      igraph_attribute_type_t type;
      long int from=igraph_es_from(graph, &it);
      long int to  =igraph_es_to  (graph, &it);
      void *ptr1, *ptr2;
      int ret=0;
      igraph_get_vertex_attribute(graph, names, from, &ptr1, &type);
      igraph_get_vertex_attribute(graph, names, to,   &ptr2, &type);
      if (type==IGRAPH_ATTRIBUTE_NUM) {
	ret=fprintf(outstream, "%f %f\n", *(real_t*)ptr1, *(real_t*)ptr2);
      } else if (type==IGRAPH_ATTRIBUTE_STR) {
	ret=fprintf(outstream, "%s %s\n", (char*)ptr1, (char*)ptr2);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  } else if (names==0) {
    /* No names but weights */
    while (!igraph_es_end(graph, &it)) {
      igraph_attribute_type_t type;
      long int from=igraph_es_from(graph, &it);
      long int to  =igraph_es_to  (graph, &it);
      long int edge=igraph_es_get(graph, &it);
      void *ptr;
      int ret=0;
      igraph_get_edge_attribute(graph, weights, edge, &ptr, &type);
      if (type==IGRAPH_ATTRIBUTE_NUM) {
	ret=fprintf(outstream, "%li %li %f\n", from, to, *(real_t*)ptr);
      } else if (type==IGRAPH_ATTRIBUTE_STR) {
	ret=fprintf(outstream, "%li %li %s\n", from, to, (char*)ptr);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  } else {
    /* Both names and weights */
    while (!igraph_es_end(graph, &it)) {
      igraph_attribute_type_t vtype, etype;
      long int from=igraph_es_from(graph, &it);
      long int to  =igraph_es_to  (graph, &it);
      long int edge=igraph_es_get(graph, &it);
      void *ptr, *ptr1, *ptr2;
      int ret=0;
      igraph_get_vertex_attribute(graph, names, from, &ptr1, &vtype);
      igraph_get_vertex_attribute(graph, names, to,   &ptr2, &vtype);
      igraph_get_edge_attribute(graph, weights, edge, &ptr, &etype);
      if (vtype==IGRAPH_ATTRIBUTE_NUM) {
	ret=fprintf(outstream, "%f %f ", *(real_t*)ptr1, *(real_t*)ptr2);
      } else if (vtype==IGRAPH_ATTRIBUTE_STR) {
	ret=fprintf(outstream, "%s %s ", (char*)ptr1, (char*)ptr2);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      if (etype==IGRAPH_ATTRIBUTE_NUM) {
	ret=fprintf(outstream, "%f\n", *(real_t*)ptr);
      } else if (etype==IGRAPH_ATTRIBUTE_STR) {
	ret=fprintf(outstream, "%s\n", (char*)ptr);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  }
  
  igraph_es_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

/**
 * \ingroup loadsave
 * \function igraph_write_graph_lgl
 * \brief Writes the graph to a file in <literal>.lgl</literal> format
 *
 * <literal>.lgl</literal> is a format used by LGL, see \ref
 * igraph_read_graph_lgl() for details.
 *
 * Note that having multiple or loop edges in an
 * <literal>.lgl</literal> file breaks the  LGL software but &igraph;
 * does not check for this condition. 
 * \param graph The graph to write. 
 * \param outstream The stream object to write to, it should be
 *        writable.
 * \param names The name of the vertex attribute, if symbolic names
 *        are written to the file. If not supply 0 here.
 * \param weights The name of the edge attribute, if they are also
 *        written to the file. If you don't want weights supply 0
 *        here.
 * \param isolates Logical, if TRUE isolate vertices are also written
 *        to the file. If FALSE they will be omitted.
 * \return Error code:
 *         <constant>IGRAPH_EFILE</constant> if there is an error
 *         writing the file. 
 *
 * Time complexity: O(|E|), the
 * number of edges if isolates is
 * FALSE, O(|V|+|E|) otherwise. All
 * file operations are expected to have time complexity 
 * O(1). 
 *
 * \sa \ref igraph_read_graph_ncol(), \ref igraph_write_graph_lgl()
 */

int igraph_write_graph_lgl(const igraph_t *graph, FILE *outstream,
			   const char *names, const char *weights,
			   bool_t isolates) {
  igraph_es_t it;
  long int actvertex=-1;
  
  IGRAPH_CHECK(igraph_es_fromorder(graph, &it));
  IGRAPH_FINALLY(igraph_es_destroy, &it);
  if (names==0 && weights==0) {
    /* No names, no weights */
    while (!igraph_es_end(graph, &it)) {
      long int from=igraph_es_from(graph, &it);
      long int to=igraph_es_to(graph, &it);
      int ret;
      if (from==actvertex) {
	ret=fprintf(outstream, "%li\n", to);
      } else {
	actvertex=from;
	ret=fprintf(outstream, "# %li\n%li\n", from, to);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  } else if (weights==0) {
    /* No weights but use names */
    while (!igraph_es_end(graph, &it)) {
      igraph_attribute_type_t type;
      long int from=igraph_es_from(graph, &it);
      long int to=igraph_es_to(graph, &it);      
      void *ptr1, *ptr2;
      int ret;
      igraph_get_vertex_attribute(graph, names, to, &ptr2, &type);
      if (from==actvertex) {
	if (type==IGRAPH_ATTRIBUTE_NUM) {
	  ret=fprintf(outstream, "%f\n", *(real_t*)ptr2);
	} else {
	  ret=fprintf(outstream, "%s\n", (char*)ptr2);
	}
      } else {
	actvertex=from;
	igraph_get_vertex_attribute(graph, names, from, &ptr1, &type);
	if (type==IGRAPH_ATTRIBUTE_NUM) {
	  ret=fprintf(outstream, "# %f\n%f\n", *(real_t*)ptr1, *(real_t*)ptr2);
	} else {
	  ret=fprintf(outstream, "# %s\n%s\n", (char*)ptr1, (char*)ptr2);
	}
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  } else if (names==0) {
    /* No names but weights */
    while (!igraph_es_end(graph, &it)) {
      igraph_attribute_type_t type;
      long int from=igraph_es_from(graph, &it);
      long int to=igraph_es_to(graph, &it);
      long int edge=igraph_es_get(graph, &it);
      void *ptr;
      int ret;
      igraph_get_edge_attribute(graph, weights, edge, &ptr, &type);
      if (from==actvertex) {
	if (type==IGRAPH_ATTRIBUTE_NUM) {
	  ret=fprintf(outstream, "%li %f\n", to, *(real_t*)ptr);
	} else {
	  ret=fprintf(outstream, "%li %s\n", to, (char*)ptr);
	}
      } else {
	actvertex=from;
	if (type==IGRAPH_ATTRIBUTE_NUM) {
	  ret=fprintf(outstream, "# %li\n%li %f\n", from, to, *(real_t*)ptr);
	} else {
	  ret=fprintf(outstream, "# %li\n%li %s\n", from, to, (char*)ptr);
	}
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  } else {
    /* Both names and weights */
    while (!igraph_es_end(graph, &it)) {
      igraph_attribute_type_t vtype, etype;
      long int from=igraph_es_from(graph, &it);
      long int to=igraph_es_to(graph, &it);
      long int edge=igraph_es_get(graph, &it);
      void *ptr1, *ptr2, *ptr;
      int ret;
      igraph_get_vertex_attribute(graph, names, to, &ptr2, &vtype);
      igraph_get_edge_attribute(graph, weights, edge, &ptr, &etype);
      if (from==actvertex) {
	if (vtype==IGRAPH_ATTRIBUTE_NUM) {
	  ret=fprintf(outstream, "%f ", *(real_t*)ptr2);
	} else {
	  ret=fprintf(outstream, "%s ", (char*)ptr2);
	} 
      } else {
	actvertex=from;
	igraph_get_vertex_attribute(graph, names, from, &ptr1, &vtype);
	if (vtype==IGRAPH_ATTRIBUTE_NUM) {
	  ret=fprintf(outstream, "# %f\n%f ", *(real_t*)ptr1, *(real_t*)ptr2);
	} else {
	  ret=fprintf(outstream, "# %s\n%s ", (char*)ptr1, (char*)ptr2);
	}
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      if (etype==IGRAPH_ATTRIBUTE_NUM) {
	ret=fprintf(outstream, "%f\n", *(real_t*)ptr);
      } else {
	ret=fprintf(outstream, "%s\n", (char*)ptr);
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
      igraph_es_next(graph, &it);
    }
  }

  if (isolates) {
    long int nov=igraph_vcount(graph);
    long int i;
    igraph_attribute_type_t type;
    void *ptr;
    int ret=0;
    igraph_vs_t it;

    /* TODO: eliminate iterators from here, dirty... */
    IGRAPH_CHECK(igraph_vs_adj(graph, &it, 0, IGRAPH_ALL));
    for (i=0; i<nov; i++) {
      igraph_vs_adj_set(graph, &it, i, IGRAPH_ALL);
      if (igraph_vs_end(graph, &it)) {
	if (names==0) {
	  ret=fprintf(outstream, "# %li\n", i);
	} else {
	  igraph_get_vertex_attribute(graph, names, i, &ptr, &type);
	  if (type==IGRAPH_ATTRIBUTE_NUM) {
	    ret=fprintf(outstream, "# %f\n", *(real_t*)ptr);
	  } else {
	    ret=fprintf(outstream, "# %s\n", (char*)ptr);
	  }
	}
      }
      if (ret<0) {
	IGRAPH_ERROR("Write failed", IGRAPH_EFILE);
      }
    }
  }  
  
  igraph_es_destroy(&it);
  IGRAPH_FINALLY_CLEAN(1);
  return 0;
}

