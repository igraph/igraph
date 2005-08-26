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

int igraph_read_graph(const char *filename, integer_t format) {
  /* TODO */
  return 0;
}

int igraph_write_graph(igraph_t *graph, const char *filename, 
		       integer_t format) {
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
