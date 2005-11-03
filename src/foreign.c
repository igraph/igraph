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
 * \ingroup loadsave
 * \brief Reads an edge list from a file and creates a graph
 * 
 * This format is simply a series of even number integers separated by
 * whitespace. The one edge (ie. two integers) per line format is thus
 * not required (but recommended for readability). Edges of directed
 * graphs are assumed to be in from, to order.
 * @param graph Pointer to an uninitialized graph object.
 * @param instream Pointer to a stream, it should be readable.
 * @param n The number of vertices in the graph. If smaller than the
 *        largest integer in the file it will be ignored. It is thus
 *        safe to supply zero here.
 * @param directed Logical, if true the graph is directed, if false it
 *        will be undirected.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|+|E|)</code>, the number of vertices
 * plus the number of edges. It is assumed that reading an integer
 * requires <code>O(1)</code> time.
 */

int igraph_read_graph_edgelist(igraph_t *graph, FILE *instream, 
			       integer_t n, bool_t directed) {

  vector_t edges;
  vector_init(&edges, 0);
  vector_reserve(&edges, 100);
  long int from, to;
  int c;

  /* skip all whitespace */
  do {
    c = getc (instream);
  } while (isspace (c));
  ungetc (c, instream);
  
  while (!feof(instream)) {
    int read;
    read=fscanf(instream, "%li", &from);
    if (read != 1) { igraph_error("Parse error in file\n"); }
    read=fscanf(instream, "%li", &to);
    if (read != 1) { igraph_error("Parse error in file\n"); }
    vector_push_back(&edges, from);
    vector_push_back(&edges, to);
    
    /* skip all whitespace */
    do {
      c = getc (instream);
    } while (isspace (c));
    ungetc (c, instream);    
  }
  
  igraph_create(graph, &edges, n, directed);
  vector_destroy(&edges);
  return 0;
}

extern int igraph_ncol_yyparse();
extern FILE *igraph_ncol_yyin;
vector_t *igraph_ncol_vector=0;
vector_t *igraph_ncol_weights=0;
igraph_trie_t *igraph_ncol_trie=0;

/**
 * \ingroup loadsave
 * \brief Reads a \a .ncol file used by LGL, also useful for creating
 * graphs from "named" (and optionally weighted) edge lists.
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
 * @param graph Pointer to an uninitialized graph object.
 * @param instream Pointer to a stream, it should be readable.
 * @param names Logical value, if TRUE the symbolic names of the
 *        vertices will be added to the graph as a vertex attribute
 *        called "name".
 * @param weights Logical value, if TRUE the weights of the
 *        edges is added to the graph as an edge attribute called
 *        "weight".
 * @return Error code.
 *
 * Time complexity: <code>O(|V|+|E|)</code> if we neglect the time
 * required by the parsing. As usual <code>|V|</code> is the number of
 * vertices, while <code>|E|</code> is the number of edges.
 */

int igraph_read_graph_ncol(igraph_t *graph, FILE *instream, 
			  bool_t names, bool_t weights) {
  
  vector_t edges, ws;
  igraph_trie_t trie;

  igraph_ncol_vector=&edges;
  vector_init(&edges, 0);
  igraph_ncol_weights=&ws;
  vector_init(&ws, 0);
  igraph_ncol_trie=&trie;
  igraph_trie_init(&trie, names);
  igraph_ncol_yyin=instream;
  igraph_ncol_yyparse();

  igraph_create(graph, &edges, 0, 0);
  vector_destroy(&edges);
  
  if (weights) {
    /* TODO: we cannot really assume that the edge ids are the same
       as the order they were added. TODO: add edges together with
       attributes */
    long int i;
    igraph_add_edge_attribute(graph, "weight", IGRAPH_ATTRIBUTE_NUM);
    for (i=0; i<igraph_ecount(graph); i++) {
      igraph_set_edge_attribute(graph, "weight", i, &VECTOR(ws)[i]);
    }
  }
  vector_destroy(&ws);
  
  if (names) {
    long int i;
    igraph_add_vertex_attribute(graph, "name", IGRAPH_ATTRIBUTE_STR);
    for (i=0; i<igraph_vcount(graph); i++) {
      char *str;
      igraph_trie_idx(&trie, i, &str);
      igraph_set_vertex_attribute(graph, "name", i, str);
    }
  }
  igraph_trie_destroy(&trie);

  return 0;
}

int igraph_read_graph_pajek(igraph_t *graph, FILE *instream) {
  /* TODO */
  return 0;
}

/**
 * \ingroup loadsave
 * \brief Writes the edge list of a graph to a file
 * 
 * One edge is written per line, separated by a single space.
 * For directed graphs edges are written in from, to order.
 * @param graph The graph object to write.
 * @param outstream Pointer to a stream, it should be writable.
 * @return Error code.
 * 
 * Time complexity: <code>O(|E|)</code>, the number of edges in the
 * graph. It is assumed that writing an integer to the file requires
 * <code>O(1)</code> time.
 */

int igraph_write_graph_edgelist(igraph_t *graph, FILE *outstream) {

  igraph_iterator_t it;
  
  igraph_iterator_efromorder(graph, &it);
  while (!igraph_end(graph, &it)) {
    fprintf(outstream, "%li %li\n", 
	    (long int) igraph_get_vertex_from(graph, &it),
	    (long int) igraph_get_vertex_to(graph, &it));
    igraph_next(graph, &it);
  }
  
  return 0;
}

/** 
 * \ingroup loadsave
 * \brief Writes the graph to a file in \a .ncol format
 * 
 * \a .ncol is a format used by LGL, see igraph_read_graph_ncol() for 
 * details. 
 * 
 * Note that having multiple or loop edges in an \a .ncol file breaks
 * the  LGL software but \a igraph does not check for this condition.
 * @param graph The graph to write.
 * @param outstream The stream object to write to, it should be
 *        writable.
 * @param names The name of the vertex attribute, if symbolic names
 *        are written to the file. If not supply 0 here.
 * @param weights The name of the edge attribute, if they are also
 *        written to the file. If you don't want weights supply 0
 *        here.
 * @return Error code.
 * 
 * Time complexity: <code>O(|E|)</code>, the number of edges. All file
 * operations are expected to have time complexity <code>O(1)</code>.
 */

int igraph_write_graph_ncol(igraph_t *graph, FILE *outstream, 
			    const char *names, const char *weights) {
  igraph_iterator_t it;
  
  igraph_iterator_efromorder(graph, &it);
  if (names==0 && weights ==0) {
    /* No names, no weights */
    while (!igraph_end(graph, &it)) {
      fprintf(outstream, "%li %li\n",
	      (long int) igraph_get_vertex_from(graph, &it),
	      (long int) igraph_get_vertex_to(graph, &it));
      igraph_next(graph, &it);
    }
  } else if (weights==0) {
    /* No weights, but use names */
    while (!igraph_end(graph, &it)) {
      igraph_attribute_type_t type;
      long int from=igraph_get_vertex_from(graph, &it);
      long int to  =igraph_get_vertex_to  (graph, &it);
      void *ptr1, *ptr2;
      igraph_get_vertex_attribute(graph, names, from, &ptr1, &type);
      igraph_get_vertex_attribute(graph, names, to,   &ptr2, &type);
      if (type==IGRAPH_ATTRIBUTE_NUM) {
	fprintf(outstream, "%f %f\n", *(real_t*)ptr1, *(real_t*)ptr2);
      } else if (type==IGRAPH_ATTRIBUTE_STR) {
	fprintf(outstream, "%s %s\n", (char*)ptr1, (char*)ptr2);
      }
      igraph_next(graph, &it);
    }
  } else if (names==0) {
    /* No names but weights */
    while (!igraph_end(graph, &it)) {
      igraph_attribute_type_t type;
      long int from=igraph_get_vertex_from(graph, &it);
      long int to  =igraph_get_vertex_to  (graph, &it);
      long int edge=igraph_get_edge(graph, &it);
      void *ptr;
      igraph_get_edge_attribute(graph, weights, edge, &ptr, &type);
      if (type==IGRAPH_ATTRIBUTE_NUM) {
	fprintf(outstream, "%li %li %f\n", from, to, *(real_t*)ptr);
      } else if (type==IGRAPH_ATTRIBUTE_STR) {
	fprintf(outstream, "%li %li %s\n", from, to, (char*)ptr);
      }
      igraph_next(graph, &it);
    }
  } else {
    /* Both names and weights */
    while (!igraph_end(graph, &it)) {
      igraph_attribute_type_t vtype, etype;
      long int from=igraph_get_vertex_from(graph, &it);
      long int to  =igraph_get_vertex_to  (graph, &it);
      long int edge=igraph_get_edge(graph, &it);
      void *ptr, *ptr1, *ptr2;
      igraph_get_vertex_attribute(graph, names, from, &ptr1, &vtype);
      igraph_get_vertex_attribute(graph, names, to,   &ptr2, &vtype);
      igraph_get_edge_attribute(graph, weights, edge, &ptr, &etype);
      if (vtype==IGRAPH_ATTRIBUTE_NUM) {
	fprintf(outstream, "%f %f ", *(real_t*)ptr1, *(real_t*)ptr2);
      } else if (vtype==IGRAPH_ATTRIBUTE_STR) {
	fprintf(outstream, "%s %s ", (char*)ptr1, (char*)ptr2);
      }
      if (etype==IGRAPH_ATTRIBUTE_NUM) {
	fprintf(outstream, "%f\n", *(real_t*)ptr);
      } else if (etype==IGRAPH_ATTRIBUTE_STR) {
	fprintf(outstream, "%s\n", (char*)ptr);
      }
      igraph_next(graph, &it);
    }
  }

  return 0;
}

int igraph_write_graph_pajek(igraph_t *graph, FILE *outstream) {
  /* TODO */
  return 0;
}


/* /\* Skip space in the string *\/ */
/* char *REST_i_skip_space(char *p) { */
/*   while (isspace(*p)) p++; */
/*   return p; */
/* } */

/* /\* The next non-empty line *\/ */
/* char *REST_i_next_line(SEXP lines, long int *currentline) { */
/*   char *result; */
/*   Rboolean l=FALSE; */
/*   (*currentline)=(*currentline)-1; */
/*   while (!l && (*currentline)+1<GET_LENGTH(lines)){ */
/*     result=CHAR(STRING_ELT(lines, (*currentline)+1))-1; */
/*     while (!l && (*(result+1) != '\0')) { */
/*       l=(!isspace(*(result+1))); */
/*       result++; */
/*     } */
/*     (*currentline)++; */
/*   } */
/*   if (l) */
/*     { return result; } */
/*   else */
/*     { return 0; } */
/* } */

/* /\* The next line _starting_ with a star *\/ */
/* char *REST_i_next_star_line(SEXP lines, long int *currentline) { */
/*   char *result; */
/*   Rboolean l=FALSE; */
/*   (*currentline)=(*currentline)-1; */
/*   while (!l && (*currentline)+1<GET_LENGTH(lines)) { */
/*     result=CHAR(STRING_ELT(lines, (*currentline)+1)); */
/*     while (isspace(*result)) { result++; } */
/*     l= (*result=='*'); */
/*     (*currentline)++; */
/*   } */
/*   if (l) */
/*     { return result; } */
/*   else */
/*     { return 0; } */
/* } */

/* /\* Read a (possibly) quoted string *\/ */
/* /\* TODO: nested " *\/ */
/* SEXP REST_i_get_quoted_string(char **posp) { */
/*   char *tmp; char *start; */
/*   int len=0; */
/*   *posp=REST_i_skip_space(*posp); */
/*   if ((**posp)=='"') { */
/*     (*posp)++; */
/*     start=*posp; */
/*     while ((**posp) != '"' && (**posp) != '\0') { (*posp)++; len++; } */
/*     if ((**posp) == '"') { (*posp)++; } */
/*   } else { */
/*     start=*posp; */
/*     while (!isspace(**posp) && (**posp) != '\0') { (*posp)++; len++; } */
/*   } */
/*   tmp=Calloc(len+1, unsigned char); */
/*   memcpy(tmp, start, len*sizeof(unsigned char)); */
/*   tmp[len]='\0'; */
/*   return CREATE_STRING_VECTOR(tmp); */
/* } */

/* /\* the color convertation matrix *\/ */
/* static long int REST_i_color_table_length=97; */
/* static char *REST_i_color_table[] = { "lightyellow", "yellow", */
/* 				      "GreenYellow","GreenYellow", */
/* 				      "Yellow", "Yellow", */
/* 				      "Goldenrod", "Goldenrod", */
/* 				      "Orange", "Orange", */
/* 				      "Maroon", "Maroon", */
/* 				      "Red", "Red", */
/* 				      "OrangeRed", "OrangeRed", */
/* 				      "Salmon", "Salmon", */
/* 				      "Magenta", "Magenta", */
/* 				      "VioletRed", "VioletRed", */
/* 				      "Gray20", "Gray20", */
/* 				      "Gray35", "Gray35", */
/* 				      "Gray55", "Gray55", */
/* 				      "Gray70", "Gray70", */
/* 				      "Gray85", "Gray85", */
/* 				      "Lavender", "Lavender", */
/* 				      "Thistle", "Thistle", */
/* 				      "Orchid", "Orchid", */
/* 				      "DarkOrchid", "DarkOrchid", */
/* 				      "Purple", "Purple", */
/* 				      "Plum", "Plum", */
/* 				      "Violet", "Violet", */
/* 				      "BlueViolet", "BlueViolet", */
/* 				      "CadetBlue", "CadetBlue", */
/* 				      "CornflowerBlue", "CornflowerBlue", */
/* 				      "MidnightBlue", "MidnightBlue", */
/* 				      "NavyBlue", "NavyBlue", */
/* 				      "RoyalBlue", "RoyalBlue", */
/* 				      "Blue", "Blue", */
/* 				      "Cyan", "Cyan", */
/* 				      "SkyBlue", "SkyBlue", */
/* 				      "Turquoise", "Turquoise", */
/* 				      "Aquamarine", "Aquamarine", */
/* 				      "Gray10", "Gray10", */
/* 				      "Gray25", "Gray25", */
/* 				      "Gray40", "Gray40", */
/* 				      "Gray60", "Gray60", */
/* 				      "Gray75", "Gray75", */
/* 				      "Gray90", "Gray90", */
/* 				      "SeaGreen", "SeaGreen", */
/* 				      "Green", "Green", */
/* 				      "ForestGreen", "ForestGreen", */
/* 				      "LimeGreen", "LimeGreen", */
/* 				      "YellowGreen", "YellowGreen", */
/* 				      "SpringGreen", "SpringGreen", */
/* 				      "Brown", "Brown", */
/* 				      "Tan", "Tan", */
/* 				      "Gray", "Gray", */
/* 				      "Black", "Black", */
/* 				      "White", "White", */
/* 				      "LightYellow", "LightYellow", */
/* 				      "LightCyan", "LightCyan", */
/* 				      "LightGreen", "LightGreen", */
/* 				      "Pink", "Pink", */
/* 				      "Gray15", "Gray15", */
/* 				      "Gray30", "Gray30", */
/* 				      "Gray45", "Gray45", */
/* 				      "Gray65", "Gray65", */
/* 				      "Gray80", "Gray80", */
/* 				      "Gray95", "Gray95", */
/* 	      /\* These are made up according to the Pajek manual *\/ */
/* 				      "Dandelion", "#ffb528", */
/* 				      "Apricot", "#ffae7a", */
/* 				      "Peach", "#ff7f4c", */
/* 				      "Melon", "#ff8a7f", */
/* 				      "YellowOrange", "#ff9400", */
/* 				      "BurntOrange", "#ff7d00", */
/* 				      "Bittersweet", "#c20200", */
/* 				      "RedOrange", "#ff3a21", */
/* 				      "Mahogany", "#a60000", */
/* 				      "BrickRed", "#b80000", */
/* 				      "RubineRed", "#ff00de", */
/* 				      "WildStrawberry", "#ff0a9c", */
/* 				      "CarnationPink", "#ff5eff", */
/* 				      "Rhodamine", "#ff2eff", */
/* 				      "Mulberry", "#a314fa", */
/* 				      "RedViolet", "#9700a8", */
/* 				      "Gray05", "Gray5", */
/* 				      "Fuchsia", "#7302eb", */
/* 				      "RoyalPurple", "#3f19ff", */
/* 				      "Periwinkle", "#6e73ff", */
/* 				      "Cerulean", "#0fe3ff", */
/* 				      "ProcessBlue", "#0affff", */
/* 				      "TealBlue", "#1efaa3", */
/* 				      "BlueGreen", "#26ffab", */
/* 				      "Emerald", "#00ff7f", */
/* 				      "JungleGreen", "#02ff7a", */
/* 				      "PineGreen", "#00bf28", */
/* 				      "OliveGreen", "#009900", */
/* 				      "RawSienna", "#8c0000", */
/* 				      "Sepia", "#4c0000", */
/* 				      "LightMagenta", "#ffccff", */
/* 				      "LightPurple", "#ccccff", */
/* 				      "LightOrange", "#ffccb3", */
/* 				      "Canary", "#ffff7f", */
/* 				      "LFadedGreen", "#e6ffcc", */
/* 				      "LSkyBlue", "#d9f3ff" }; */

/* /\* Convert Pajek color to string *\/ */
/* SEXP REST_i_pajek_color(char **posp) { */
/*   SEXP pcolor=REST_i_get_quoted_string(posp); */
/*   char *pcolorstr=CHAR(pcolor); */
/*   int l=0; */
/*   int idx=-1; */
/*   while (!l && idx+1 < REST_i_color_table_length) { */
/*     l=!strcasecmp(pcolorstr, REST_i_color_table[2*(idx+1)]); */
/*     idx += 1; */
/*   } */
/*   if (!l) { */
/*     return NA_STRING; */
/*   } else { */
/*     return CREATE_STRING_VECTOR(REST_i_color_table[2*idx+1]); */
/*   } */
/* } */
    
/* SEXP REST_import_pajek(SEXP interface, SEXP lines, SEXP other, */
/* 		       SEXP pattributes) { */

/*   SEXP result; */
/*   long int currentline; */
/*   char *currentpos; */
/*   long int no_of_lines; */
/*   long int no_of_nodes; */
  
/*   dqueue_t edges; */
/*   SEXP args; */
/*   int directed=0; /\* 0-not known, 1-directed, 2-undirected *\/ */
/*   SEXP dirattr; */
/*   long int unprot=0; */
/*   SEXP edgevec; */
/*   long int i, j; */

/*   SEXP graph; */
/*   REST_i_ptrtable_t ptrtable; */

/*   Rboolean attributes; */
/*   /\* do we have these attributes already? *\/ */
/*   SEXP attr_coords=0, attr_name=0, attr_ic=0; */

/*   SEXP names=0; */

/*   no_of_lines=GET_LENGTH(lines); */
/*   attributes=LOGICAL(pattributes)[0]; */

/*   if (no_of_lines==0) { */
/*     error("invalid file format: empty file"); */
/*   } */

/*   currentline=0; */
/*   currentpos=REST_i_next_line(lines, &currentline); */
  
/*   /\* Optional "*network" *\/ */
/*   if (currentpos != 0 && !strncasecmp(currentpos, "*network", 8)) { */
/*     currentline++; */
/*     currentpos=REST_i_next_line(lines, &currentline); */
/*   } */
  
/*   /\* Number of vertices *\/ */
/*   if (currentpos == 0 || */
/*       strncasecmp(currentpos, "*vertices", 9)) { */
/*     error("invalid file format, can't find '*Vertices'"); */
/*   } */
  
/*   currentpos+=9; */
/*   no_of_nodes=strtol(currentpos, 0, 0); */

/*   /\* Vertices *\/ */
/*   if (!attributes) { */
/*     /\* Read the vertices, no attributes, fast *\/ */
/*     currentline++; */
/*     currentpos=REST_i_next_star_line(lines, &currentline); */
/*   } else { */
/*     /\* Read the vertices *\/ */
/*     PROTECT(attr_name=NEW_LIST(no_of_nodes)); unprot++; */
/*     currentline++; */
/*     currentpos=REST_i_next_line(lines, &currentline); */
/*     while (currentpos != 0 && *currentpos != '*') { */
/*       char *tail; */
/*       long int vid=strtol(currentpos, &tail, 0); */
/*       currentpos=REST_i_skip_space(tail); */
/*       if (*currentpos != '\0') { */
/* 	/\* label *\/ */
/* 	SET_VECTOR_ELT(attr_name, vid-1, */
/* 		       ScalarString(REST_i_get_quoted_string(&currentpos))); */
/* 	currentpos=REST_i_skip_space(currentpos); */
/*       } */
/*       while (*currentpos != '\0') { */
/* 	/\* we found an(other) attribute *\/ */
/* 	if (isdigit(*currentpos)) { */
/* 	  /\* coordinates (2 or 3) *\/ */
/* 	  char *tail; */
/* 	  SEXP tmp=NEW_NUMERIC(3);    /\* no protect, risky??? *\/ */
/* 	  REAL(tmp)[0]=strtod(currentpos, &currentpos); */
/* 	  REAL(tmp)[1]=strtod(currentpos, &currentpos); */
/* 	  REAL(tmp)[2]=strtod(currentpos, &currentpos); */
/* 	  if (attr_coords==0) { */
/* 	    PROTECT(attr_coords=NEW_LIST(no_of_nodes)); unprot++; */
/* 	  } */
/* 	  SET_VECTOR_ELT(attr_coords, vid-1, tmp); */
/* 	} else if (!strncasecmp(currentpos, "triangle", 8)) { */
/* 	  /\* shape *\/ */
/* 	  currentpos+=8; */
/* 	} else if (!strncasecmp(currentpos, "ellipse", 7)) { */
/* 	  /\* shape *\/ */
/* 	  currentpos+=7; */
/* 	} else if (!strncasecmp(currentpos, "diamond", 7)) { */
/* 	  /\* shape *\/ */
/* 	  currentpos+=7; */
/* 	} else if (!strncasecmp(currentpos, "s_size", 6)) { */
/* 	  /\* vertex size *\/ */
/* 	  currentpos+=6; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "x_fact", 6)) { */
/* 	  /\* magnification in x direction ??? *\/ */
/* 	  currentpos+=6; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "y_fact", 6)) { */
/* 	  /\* magnification in y direction ??? *\/ */
/* 	  currentpos+=6; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "cross", 5)) { */
/* 	  /\* shape *\/ */
/* 	  currentpos+=5; */
/* 	} else if (!strncasecmp(currentpos, "empty", 5)) { */
/* 	  /\* shape *\/ */
/* 	  currentpos+=5; */
/* 	} else if (!strncasecmp(currentpos, "SHAPE", 5)) { */
/* 	  /\* external name of vertex ??? *\/ */
/* 	  currentpos+=5; */
/* 	  REST_i_get_quoted_string(&currentpos); */
/* 	} else if (!strncasecmp(currentpos, "polar", 5)) { */
/* 	  /\* where edges can join the shape, polar coords *\/ */
/* 	  currentpos+=5; */
/* 	  strtod(currentpos, &currentpos); */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "lphi", 4)) { */
/* 	  /\* position of label in degrees *\/ */
/* 	  currentpos+=4; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "font", 4)) { */
/* 	  /\* postscript font for labels *\/ */
/* 	  currentpos+=4; */
/* 	  REST_i_get_quoted_string(&currentpos); */
/* 	} else if (!strncasecmp(currentpos, "cart", 4)) { */
/* 	  /\* where edges can join the shape, cartesian coords *\/ */
/* 	  currentpos+=4; */
/* 	  strtod(currentpos, &currentpos); */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "circ", 4)) { */
/* 	  /\* where edges can join the shape, ??? *\/ */
/* 	  currentpos+=4; */
/* 	  strtod(currentpos, &currentpos); */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "box", 3)) { */
/* 	  /\* shape *\/ */
/* 	  currentpos+=3; */
/* 	} else if (!strncasecmp(currentpos, "phi", 3)) { */
/* 	  /\* rotation of object in degrees *\/ */
/* 	  currentpos+=3; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "fos", 3)) { */
/* 	  /\* font size for label *\/ */
/* 	  currentpos+=3; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "sh", 2)) { */
/* 	  /\* shape (internal?) *\/ */
/* 	  currentpos+=2; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "ic", 2)) { */
/* 	  /\* interior color of vertex *\/ */
/* 	  SEXP color; */
/* 	  currentpos+=2; */
/* 	  color=REST_i_pajek_color(&currentpos); */
/* 	  if (attr_ic==0) { */
/* 	    PROTECT(attr_ic=NEW_LIST(no_of_nodes)); unprot++; */
/* 	  } */
/* 	  SET_VECTOR_ELT(attr_ic, vid-1, ScalarString(color)); */
/* 	} else if (!strncasecmp(currentpos, "bc", 2)) { */
/* 	  /\* boundary color of vertex *\/ */
/* 	  currentpos+=2; */
/* 	  REST_i_pajek_color(&currentpos); */
/* 	} else if (!strncasecmp(currentpos, "bw", 2)) { */
/* 	  /\* boundary width of vertex *\/ */
/* 	  currentpos+=2; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "lc", 2)) { */
/* 	  /\* label color *\/ */
/* 	  currentpos+=2; */
/* 	  REST_i_pajek_color(&currentpos); */
/* 	} else if (!strncasecmp(currentpos, "la", 2)) { */
/* 	  /\* label angle in degrees ??? *\/ */
/* 	  currentpos+=2; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "lr", 2)) { */
/* 	  /\* label distance from vertex *\/ */
/* 	  currentpos+=2; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "r", 1)) { */
/* 	  /\* for rectangle and diamond: radius of corners *\/ */
/* 	  currentpos+=1; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else if (!strncasecmp(currentpos, "q", 1)) { */
/* 	  /\* diamonds: ratio between top and middle side *\/ */
/* 	  currentpos+=1; */
/* 	  strtod(currentpos, &currentpos); */
/* 	} else { */
/* 	  UNPROTECT(unprot); */
/* 	  error("Pajek file invalid in line %d", currentline+1); */
/* 	} */
/* 	currentpos=REST_i_skip_space(currentpos); */
/*       } */

/*       currentline++; */
/*       currentpos=REST_i_next_line(lines, &currentline); */
/*     } */
/*   } /\* if (! attributes) *\/ */

/*   /\* Edges *\/ */

/*   dqueue_init(&edges, 1000); */
  
/*   /\* inv: on '*' or end of file *\/ */
/*   while (currentpos != 0) { */
/*     if (!strncasecmp(currentpos, "*arcslist", 9) */
/* 	|| !strncasecmp(currentpos, "*edgeslist", 10)) { */
      
/*       Rboolean arcs; */
/*       long int startsize=dqueue_size(&edges); */
/*       if (!strncasecmp(currentpos, "*arcslist", 9)) { */
/* 	arcs=TRUE; */
/*       } else { */
/* 	arcs=FALSE; */
/*       } */
      
/*       currentline++; */
/*       currentpos=REST_i_next_line(lines, &currentline); */
/*       while (currentpos != 0 && *currentpos != '*') { */
/* 	long int from=strtol(currentpos, &currentpos, 0); */
/* 	currentpos=REST_i_skip_space(currentpos); */
/* 	while (*currentpos != '\0') { */
/* 	  char *tail; */
/* 	  dqueue_push(&edges, from); */
/* 	  dqueue_push(&edges, strtol(currentpos, &tail, 0)); */
/* 	  if (tail==currentpos) { */
/* 	    UNPROTECT(unprot); */
/* 	    error("Error, numeric expected in line %d", currentline); */
/* 	  } */
/* 	  currentpos=REST_i_skip_space(tail); */
/* 	} */
/* 	currentline++; */
/* 	currentpos=REST_i_next_line(lines, &currentline); */
/*       } */

/*       if (dqueue_size(&edges) != startsize && directed==0) { */
/* 	if (arcs) { directed=1; } else { directed=2; } */
/*       } */

/*     } else if (!strncasecmp(currentpos, "*arcs", 5) || */
/* 	       !strncasecmp(currentpos, "*edges", 6)) { */
/*       Rboolean arcs; */
/*       long int startsize=dqueue_size(&edges); */
/*       if (!strncasecmp(currentpos, "*arcs", 5)) */
/* 	{ arcs=TRUE; } else { arcs=FALSE; } */
      
/*       currentline++; */
/*       currentpos=REST_i_next_line(lines, &currentline); */
/*       while (currentpos != 0 && *currentpos != '*') { */
/* 	char *tail; */
/* 	dqueue_push(&edges, strtol(currentpos, &tail, 0)); */
/* 	dqueue_push(&edges, strtol(tail, 0, 0)); */
/* 	currentline++; */
/* 	currentpos=REST_i_next_line(lines, &currentline); */
/*       } */
      
/*       if (dqueue_size(&edges) != startsize && directed==0) { */
/* 	if (arcs) { directed=1; } else { directed=2; } */
/*       } */
      
/*     } else if (!strncasecmp(currentpos, "*matrix", 7)) { */
/*       int i, j; */

/*       currentline++; */
/*       currentpos=REST_i_next_line(lines, &currentline); */
/*       for (i=0; i<no_of_nodes; i++) { */
/* 	for (j=0; j<no_of_nodes; j++) { */
/* 	  double elem=strtod(currentpos, &currentpos); */
/* 	  if (elem>0) { */
/* 	    dqueue_push(&edges, i+1); */
/* 	    dqueue_push(&edges, j+1); */
/* 	  } */
/* 	} */
/* 	currentline++; */
/* 	currentpos=REST_i_next_line(lines, &currentline); */
/*       } */
      
/*     } else { */
/*       error("Cannot interpret pajek file, error in line %d.", currentline+1); */
/*     } */
/*   } */
  
/*   i=dqueue_size(&edges); */
/*   PROTECT(edgevec=NEW_NUMERIC(i)); unprot++; */
/*   for (j=0; j<i; j++) { */
/*     REAL(edgevec)[j]=dqueue_pop(&edges); */
/*   } */

/*   if (isNull(REST_i_get_list_element(other, "directed"))) { */
/*     SEXP tmp, names; */
/*     PROTECT(tmp=NEW_LIST(1)); */
/*     SET_VECTOR_ELT(tmp, 0, ScalarLogical(directed==1 ? TRUE : FALSE)); */
/*     PROTECT(names=ScalarString(CREATE_STRING_VECTOR("directed"))); */
/*     SET_NAMES(tmp, names); */
/*     UNPROTECT(2); */
/*     /\* TODO: dirty Hack: calling 'interface' explicitly *\/ */
/*     PROTECT(other=EVAL(lang3(interface, */
/* 			     ScalarString(CREATE_STRING_VECTOR("c")), */
/* 			     AS_LIST(list2(tmp, other))))); */
/*     unprot++; */
/*   } */

/*   ptrtable = REST_i_getptrtable(0); /\* default *\/ */
/*   graph=ptrtable.graph_empty(interface, other); */
/*   ptrtable = REST_i_getptrtable(graph); */
/*   PROTECT(result= */
/* 	  ptrtable.add_edges(interface,  */
/* 			     ptrtable.add_vertices(interface, graph,  */
/* 						   no_of_nodes), */
/* 			     edgevec)); unprot++; */
  
/*   /\* add the attributes *\/ */
/*   if (attr_coords != 0) { */
/*     UNPROTECT(1); */
/*     PROTECT(result=ptrtable.set_vertex_attribute */
/* 	    (interface, ptrtable.add_vertex_attribute(interface, result,  */
/* 						      "coords",  */
/* 						      ScalarReal(NA_REAL)), */
/* 	     "coords", NULL_USER_OBJECT, attr_coords)); */
/*   } */
/*   if (attr_name != 0) { */
/*     UNPROTECT(1); */
/*     PROTECT(result=ptrtable.set_vertex_attribute */
/* 	    (interface, ptrtable.add_vertex_attribute(interface,  */
/* 						      result, "name",  */
/* 						      ScalarReal(NA_REAL)), */
/* 	    "name", NULL_USER_OBJECT, attr_name)); */
/*   } */
/*   if (attr_ic != 0) { */
/*     UNPROTECT(1); */
/*     PROTECT(result=ptrtable.set_vertex_attribute */
/* 	    (interface, ptrtable.add_vertex_attribute(interface, result,  */
/* 						      "color",  */
/* 						      ScalarReal(NA_REAL)), */
/* 	     "color", NULL_USER_OBJECT, attr_ic)); */
/*   } */
  
/*   dqueue_destroy(&edges); */
/*   UNPROTECT(unprot); */
/*   return result; */
/* } */
