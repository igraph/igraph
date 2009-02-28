/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

%{

/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include <stdio.h>
#include <string.h>
#include "types.h"
#include "memory.h"
#include "error.h"
#include "attributes.h"
#include "config.h"
#include "igraph_math.h"
#include <math.h>
extern int igraph_pajek_yylex(void);
extern int igraph_pajek_mylineno;
extern char *igraph_pajek_yytext;
extern int igraph_pajek_yyleng;
char *igraph_i_pajek_errmsg;
int igraph_pajek_yyerror(char *s);
int igraph_i_pajek_add_string_vertex_attribute(const char *name, 
					       const char *value,
					       int len);
int igraph_i_pajek_add_string_edge_attribute(const char *name, 
					     const char *value,
					     int len);
int igraph_i_pajek_add_numeric_vertex_attribute(const char *name, 
						igraph_real_t value);
int igraph_i_pajek_add_numeric_edge_attribute(const char *name, 
					      igraph_real_t value);
int igraph_i_pajek_add_numeric_attribute(igraph_trie_t *names,
					 igraph_vector_ptr_t *attrs,
					 long int count,
					 const char *attrname,
					 igraph_integer_t vid,
					 igraph_real_t number);
int igraph_i_pajek_add_string_attribute(igraph_trie_t *names,
					igraph_vector_ptr_t *attrs,
					long int count,
					const char *attrname,
					igraph_integer_t vid,
					const char *str);
void igraph_i_pajek_reset_scanner(void);
extern igraph_vector_t *igraph_pajek_vector;
extern igraph_bool_t igraph_pajek_directed;
extern long int igraph_pajek_vcount;
extern int igraph_pajek_mode;
extern long int igraph_pajek_actfrom, igraph_pajek_actto;
extern igraph_trie_t *igraph_i_pajek_vertex_attribute_names;
extern igraph_vector_t *igraph_i_pajek_vertex_attribute_types;
extern igraph_vector_ptr_t *igraph_i_pajek_vertex_attributes;
extern igraph_trie_t *igraph_i_pajek_edge_attribute_names;
extern igraph_vector_t *igraph_i_pajek_edge_attribute_types;
extern igraph_vector_ptr_t *igraph_i_pajek_edge_attributes;
extern igraph_real_t igraph_pajek_get_number(const char *str, long int len);
extern long int igraph_i_pajek_actvertex;
extern long int igraph_i_pajek_actedge;
%}

%output="y.tab.c"
%name-prefix="igraph_pajek_yy"
%defines

%union {
  long int intnum;
  double   realnum;  
  struct {
    char *str;
    int len;
  } string;  
}

%type <intnum>   longint;
%type <intnum>   arcfrom;
%type <intnum>   arcto;
%type <intnum>   edgefrom;
%type <intnum>   edgeto;
%type <realnum>  number;
%type <string>   word;
%type <string>   vpwordpar;
%type <string>   epwordpar;
%type <intnum>   vertex;

%token NEWLINE
%token NUM
%token ALNUM
%token QSTR
%token PSTR
%token NETWORKLINE
%token VERTICESLINE
%token ARCSLINE
%token EDGESLINE
%token ARCSLISTLINE
%token EDGESLISTLINE
%token MATRIXLINE

%token VP_X_FACT
%token VP_Y_FACT
%token VP_IC
%token VP_BC
%token VP_LC
%token VP_LR
%token VP_LPHI
%token VP_BW
%token VP_FOS
%token VP_PHI
%token VP_R
%token VP_Q
%token VP_LA
%token VP_LR
%token VP_FONT
%token VP_URL
%token VP_SIZE

%token EP_C
%token EP_S
%token EP_A
%token EP_W
%token EP_H1
%token EP_H2
%token EP_A1
%token EP_A2
%token EP_K1
%token EP_K2
%token EP_AP
%token EP_P
%token EP_L
%token EP_LP
%token EP_LR
%token EP_LPHI
%token EP_LC
%token EP_LA
%token EP_SIZE
%token EP_FOS

%%

input: nethead vertices edgeblock;

nethead: /* empty */ | NETWORKLINE words NEWLINE;

vertices: verticeshead NEWLINE vertdefs;

verticeshead: VERTICESLINE longint {
  igraph_pajek_vcount=$2;
};

vertdefs: /* empty */  | vertdefs vertexline;

vertexline: NEWLINE |
            vertex NEWLINE |
            vertex { igraph_i_pajek_actvertex=$1; } vertexid vertexcoords shape params NEWLINE { }
;

vertex: longint { $$=$1; igraph_pajek_mode=1; };

vertexid: word {
  igraph_i_pajek_add_string_vertex_attribute("id", $1.str, $1.len);
};

vertexcoords: /* empty */ 
            | number number { 
  igraph_i_pajek_add_numeric_vertex_attribute("x", $1);
  igraph_i_pajek_add_numeric_vertex_attribute("y", $2);
	    }
            | number number number { 
  igraph_i_pajek_add_numeric_vertex_attribute("x", $1);
  igraph_i_pajek_add_numeric_vertex_attribute("y", $2);
  igraph_i_pajek_add_numeric_vertex_attribute("z", $3);
	    };

shape: /* empty */ | word { 
  igraph_i_pajek_add_string_vertex_attribute("shape", $1.str, $1.len);	      
};

params: /* empty */ | params param;

param:
       vpword
     | VP_X_FACT number {
	 igraph_i_pajek_add_numeric_vertex_attribute("xfact", $2);
       }
     | VP_Y_FACT number {
	 igraph_i_pajek_add_numeric_vertex_attribute("yfact", $2);
       }
     | VP_IC number number number { /* RGB color */
         igraph_i_pajek_add_numeric_vertex_attribute("color-red", $2);
	 igraph_i_pajek_add_numeric_vertex_attribute("color-green", $3);
	 igraph_i_pajek_add_numeric_vertex_attribute("color-blue", $4);
       }
     | VP_BC number number number {
         igraph_i_pajek_add_numeric_vertex_attribute("framecolor-red", $2);
	 igraph_i_pajek_add_numeric_vertex_attribute("framecolor-green", $3);
	 igraph_i_pajek_add_numeric_vertex_attribute("framecolor-blue", $4);
       }
     | VP_LC number number number {
         igraph_i_pajek_add_numeric_vertex_attribute("labelcolor-red", $2);
	 igraph_i_pajek_add_numeric_vertex_attribute("labelcolor-green", $3);
	 igraph_i_pajek_add_numeric_vertex_attribute("labelcolor-blue", $4);
       }
     | VP_LR number {
         igraph_i_pajek_add_numeric_vertex_attribute("labeldist", $2);
     }  
     | VP_LPHI number {
         igraph_i_pajek_add_numeric_vertex_attribute("labeldegree2", $2);
     }  
     | VP_BW number {
         igraph_i_pajek_add_numeric_vertex_attribute("framewidth", $2);
     }         
     | VP_FOS number {
         igraph_i_pajek_add_numeric_vertex_attribute("fontsize", $2);
     }          
     | VP_PHI number {       
         igraph_i_pajek_add_numeric_vertex_attribute("rotation", $2);
     }         
     | VP_R number {
         igraph_i_pajek_add_numeric_vertex_attribute("radius", $2);
     }         
     | VP_Q number {
         igraph_i_pajek_add_numeric_vertex_attribute("diamondratio", $2);
     }         
     | VP_LA number {
         igraph_i_pajek_add_numeric_vertex_attribute("labeldegree", $2);
     }         
     | VP_SIZE number {
         igraph_i_pajek_add_numeric_vertex_attribute("vertexsize", $2);
     } 
;

vpword: VP_FONT { igraph_pajek_mode=3; } vpwordpar { 
         igraph_pajek_mode=1;
	 igraph_i_pajek_add_string_vertex_attribute("font", $3.str, $3.len);
     } 
     | VP_URL { igraph_pajek_mode=3; } vpwordpar {
         igraph_pajek_mode=1;
	 igraph_i_pajek_add_string_vertex_attribute("url", $3.str, $3.len);
     } 
     | VP_IC { igraph_pajek_mode=3; } vpwordpar {
         igraph_pajek_mode=1;
	 igraph_i_pajek_add_string_vertex_attribute("color", $3.str, $3.len);
     } 
     | VP_BC { igraph_pajek_mode=3; } vpwordpar {
         igraph_pajek_mode=1;
	 igraph_i_pajek_add_string_vertex_attribute("framecolor", 
						    $3.str, $3.len);
     } 
     | VP_LC { igraph_pajek_mode=3; } vpwordpar {
         igraph_pajek_mode=1;
	 igraph_i_pajek_add_string_vertex_attribute("labelcolor", 
						    $3.str, $3.len);
     } 
;

vpwordpar: word { $$=$1; };

edgeblock: /* empty */ | edgeblock arcs | edgeblock edges | edgeblock arcslist | edgeblock edgeslist | edgeblock adjmatrix;

arcs: ARCSLINE NEWLINE arcsdefs { igraph_pajek_directed=1; };

arcsdefs: /* empty */ | arcsdefs arcsline;

arcsline: NEWLINE | 
          arcfrom arcto { igraph_i_pajek_actedge++;
	                  igraph_pajek_mode=2; } weight edgeparams NEWLINE  { 
  igraph_vector_push_back(igraph_pajek_vector, $1-1);
  igraph_vector_push_back(igraph_pajek_vector, $2-1); }
;

arcfrom: longint;

arcto: longint;

edges: EDGESLINE NEWLINE edgesdefs { igraph_pajek_directed=0; };

edgesdefs: /* empty */ | edgesdefs edgesline;

edgesline: NEWLINE | 
          edgefrom edgeto { igraph_i_pajek_actedge++; 
	                    igraph_pajek_mode=2; } weight edgeparams NEWLINE { 
  igraph_vector_push_back(igraph_pajek_vector, $1-1);
  igraph_vector_push_back(igraph_pajek_vector, $2-1); }
;

edgefrom: longint;

edgeto: longint;

weight: /* empty */ | number {
  igraph_i_pajek_add_numeric_edge_attribute("weight", $1);
};

edgeparams: /* empty */ | edgeparams edgeparam;

edgeparam:
     epword
   | EP_C number number number {
       igraph_i_pajek_add_numeric_edge_attribute("color-red", $2);
       igraph_i_pajek_add_numeric_edge_attribute("color-green", $3);
       igraph_i_pajek_add_numeric_edge_attribute("color-blue", $4);
   }
   | EP_S number { 
       igraph_i_pajek_add_numeric_edge_attribute("arrowsize", $2);
   }
   | EP_W number {
       igraph_i_pajek_add_numeric_edge_attribute("edgewidth", $2);
   }
   | EP_H1 number {
       igraph_i_pajek_add_numeric_edge_attribute("hook1", $2);
   }
   | EP_H2 number {
       igraph_i_pajek_add_numeric_edge_attribute("hook2", $2);
   }
   | EP_A1 number {
       igraph_i_pajek_add_numeric_edge_attribute("angle1", $2);
   }
   | EP_A2 number {
       igraph_i_pajek_add_numeric_edge_attribute("angle2", $2);
   }
   | EP_K1 number {
       igraph_i_pajek_add_numeric_edge_attribute("velocity1", $2);
   }
   | EP_K2 number {
       igraph_i_pajek_add_numeric_edge_attribute("velocity2", $2);
   }
   | EP_AP number {
       igraph_i_pajek_add_numeric_edge_attribute("arrowpos", $2);
   }
   | EP_LP number {
       igraph_i_pajek_add_numeric_edge_attribute("labelpos", $2);
   }
   | EP_LR number {
       igraph_i_pajek_add_numeric_edge_attribute("labelangle", $2);
   }
   | EP_LPHI number {
       igraph_i_pajek_add_numeric_edge_attribute("labelangle2", $2);
   }
   | EP_LA number {
       igraph_i_pajek_add_numeric_edge_attribute("labeldegree", $2);
   }
   | EP_SIZE number {		/* what is this??? */
       igraph_i_pajek_add_numeric_edge_attribute("arrowsize", $2);
   }
   | EP_FOS number {
       igraph_i_pajek_add_numeric_edge_attribute("fontsize", $2);
   }
;

epword: EP_A { igraph_pajek_mode=4; } epwordpar {
      igraph_pajek_mode=2;
      igraph_i_pajek_add_string_edge_attribute("arrowtype", $3.str, $3.len);
    }
    | EP_P { igraph_pajek_mode=4; } epwordpar {
      igraph_pajek_mode=2;
      igraph_i_pajek_add_string_edge_attribute("linepattern", $3.str, $3.len);
    }
    | EP_L { igraph_pajek_mode=4; } epwordpar {
      igraph_pajek_mode=2;
      igraph_i_pajek_add_string_edge_attribute("label", $3.str, $3.len);
    }
    | EP_LC { igraph_pajek_mode=4; } epwordpar {
      igraph_pajek_mode=2;
      igraph_i_pajek_add_string_edge_attribute("labelcolor", $3.str, $3.len);
    }
    | EP_C { igraph_pajek_mode=4; } epwordpar {
      igraph_pajek_mode=2;
      igraph_i_pajek_add_string_edge_attribute("color", $3.str, $3.len);
    }
;

epwordpar: word { igraph_pajek_mode=2; $$=$1; };

arcslist: ARCSLISTLINE NEWLINE arcslistlines { igraph_pajek_directed=1; };

arcslistlines: /* empty */ | arcslistlines arclistline;

arclistline: NEWLINE | arclistfrom arctolist NEWLINE;

arctolist: /* empty */ | arctolist arclistto;

arclistfrom: longint { igraph_pajek_mode=0; igraph_pajek_actfrom=fabs($1)-1; };

arclistto: longint { 
  igraph_vector_push_back(igraph_pajek_vector, igraph_pajek_actfrom); 
  igraph_vector_push_back(igraph_pajek_vector, fabs($1)-1); 
};

edgeslist: EDGESLISTLINE NEWLINE edgelistlines { igraph_pajek_directed=0; };

edgelistlines: /* empty */ | edgelistlines edgelistline;

edgelistline: NEWLINE | edgelistfrom edgetolist NEWLINE;

edgetolist: /* empty */ | edgetolist edgelistto;

edgelistfrom: longint { igraph_pajek_mode=0; igraph_pajek_actfrom=fabs($1)-1; };

edgelistto: longint { 
  igraph_vector_push_back(igraph_pajek_vector, igraph_pajek_actfrom); 
  igraph_vector_push_back(igraph_pajek_vector, fabs($1)-1); 
};

/* -----------------------------------------------------*/

adjmatrix: matrixline NEWLINE adjmatrixlines;

matrixline: MATRIXLINE { igraph_pajek_actfrom=0; igraph_pajek_actto=0; };

adjmatrixlines: /* empty */ | adjmatrixlines adjmatrixline;

adjmatrixline: adjmatrixnumbers NEWLINE { igraph_pajek_actfrom++; igraph_pajek_actto=0; };

adjmatrixnumbers: /* empty */ | adjmatrixentry adjmatrixnumbers;

adjmatrixentry: longint {
  if ($1>0) { 
    igraph_vector_push_back(igraph_pajek_vector, igraph_pajek_actfrom);
    igraph_vector_push_back(igraph_pajek_vector, igraph_pajek_actto);
  }
  igraph_pajek_actto++;
};

/* -----------------------------------------------------*/

longint: NUM { $$=igraph_pajek_get_number(igraph_pajek_yytext,
					  igraph_pajek_yyleng); };

number: NUM  { $$=igraph_pajek_get_number(igraph_pajek_yytext,
					  igraph_pajek_yyleng); };

words: /* empty */ | words word;

word: ALNUM { $$.str=igraph_pajek_yytext; $$.len=igraph_pajek_yyleng; }
      | NUM { $$.str=igraph_pajek_yytext; $$.len=igraph_pajek_yyleng; }
      | QSTR { $$.str=igraph_pajek_yytext+1; $$.len=igraph_pajek_yyleng-2; };

%%

int igraph_pajek_yyerror(char *s)
{
  static char str[300];  
  igraph_i_pajek_reset_scanner();
  snprintf(str, sizeof(str), "Parse error in Pajek file, line %li (%s)", 
	   (long)igraph_pajek_mylineno, s);
  igraph_i_pajek_errmsg = str;
  return 0;
}

igraph_real_t igraph_pajek_get_number(const char *str, long int length) {
  igraph_real_t num;
  char *tmp=igraph_Calloc(length+1, char);
  
  strncpy(tmp, str, length);
  tmp[length]='\0';
  sscanf(tmp, "%lf", &num);
  igraph_Free(tmp);
  return num;
} 

/* TODO: NA's */

int igraph_i_pajek_add_numeric_attribute(igraph_trie_t *names,
					 igraph_vector_ptr_t *attrs,
					 long int count,
					 const char *attrname,
					 igraph_integer_t vid,
					 igraph_real_t number) {
  long int attrsize=igraph_trie_size(names);
  long int id;
  igraph_vector_t *na;
  igraph_i_attribute_record_t *rec;

  igraph_trie_get(names, attrname, &id);
  if (id == attrsize) {
    /* add a new attribute */
    rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    na=igraph_Calloc(1, igraph_vector_t);
    igraph_vector_init(na, count);
    rec->name=strdup(attrname);
    rec->type=IGRAPH_ATTRIBUTE_NUMERIC;
    rec->value=na;
    igraph_vector_ptr_push_back(attrs, rec);
  }
  rec=VECTOR(*attrs)[id];
  na=(igraph_vector_t*)rec->value;
  if (igraph_vector_size(na) == vid) {
    IGRAPH_CHECK(igraph_vector_push_back(na, number));
  } else if (igraph_vector_size(na) < vid) {
    long int origsize=igraph_vector_size(na);
    IGRAPH_CHECK(igraph_vector_resize(na, (long int)vid+1));
    for (;origsize<count; origsize++) {
      VECTOR(*na)[origsize] = IGRAPH_NAN;
    }
    VECTOR(*na)[(long int) vid] = number;
  } else { 
    VECTOR(*na)[(long int) vid] = number;
  }    

  return 0;
}

/* TODO: NA's */

int igraph_i_pajek_add_string_attribute(igraph_trie_t *names,
					igraph_vector_ptr_t *attrs,
					long int count,
					const char *attrname,
					igraph_integer_t vid,
					const char *str) {
  long int attrsize=igraph_trie_size(names);
  long int id;
  igraph_strvector_t *na;
  igraph_i_attribute_record_t *rec;
  long int i;

  igraph_trie_get(names, attrname, &id);
  if (id == attrsize) {
    /* add a new attribute */
    rec=igraph_Calloc(1, igraph_i_attribute_record_t);
    na=igraph_Calloc(1, igraph_strvector_t);
    igraph_strvector_init(na, count);
    for (i=0; i<count; i++) {
      igraph_strvector_set(na, i, "");
    }
    rec->name=strdup(attrname);
    rec->type=IGRAPH_ATTRIBUTE_STRING;
    rec->value=na;
    igraph_vector_ptr_push_back(attrs, rec);
  }
  rec=VECTOR(*attrs)[id];
  na=(igraph_strvector_t*)rec->value;
  if (igraph_strvector_size(na) <= vid) { 
    long int origsize=igraph_strvector_size(na);
    IGRAPH_CHECK(igraph_strvector_resize(na, vid+1));
    for (;origsize<count; origsize++) {
      igraph_strvector_set(na, origsize, "");
    }
  }
  igraph_strvector_set(na, vid, str);

  return 0;
}

int igraph_i_pajek_add_string_vertex_attribute(const char *name, 
					       const char *value,
					       int len) {
  char *tmp;
  int ret;

  tmp=igraph_Calloc(len+1, char);
  if (tmp==0) {
    IGRAPH_ERROR("cannot add element to hash table", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp);
  strncpy(tmp, value, len);
  tmp[len]='\0';

  ret=igraph_i_pajek_add_string_attribute(igraph_i_pajek_vertex_attribute_names,
					  igraph_i_pajek_vertex_attributes,
					  igraph_pajek_vcount,
					  name, igraph_i_pajek_actvertex-1,
					  tmp);
  
  igraph_Free(tmp);
  IGRAPH_FINALLY_CLEAN(1);
  
  return ret;
}

int igraph_i_pajek_add_string_edge_attribute(const char *name, 
					     const char *value,
					     int len) {
  char *tmp;
  int ret;

  tmp=igraph_Calloc(len+1, char);
  if (tmp==0) {
    IGRAPH_ERROR("cannot add element to hash table", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(free, tmp);
  strncpy(tmp, value, len);
  tmp[len]='\0';
  
  ret=igraph_i_pajek_add_string_attribute(igraph_i_pajek_edge_attribute_names,
					  igraph_i_pajek_edge_attributes,
					  igraph_i_pajek_actedge,
					  name, igraph_i_pajek_actedge-1,
					  tmp);

  igraph_Free(tmp);
  IGRAPH_FINALLY_CLEAN(1);
  
  return ret;
}

int igraph_i_pajek_add_numeric_vertex_attribute(const char *name, 
						igraph_real_t value) {
  
  return
    igraph_i_pajek_add_numeric_attribute(igraph_i_pajek_vertex_attribute_names,
					 igraph_i_pajek_vertex_attributes,
					 igraph_pajek_vcount,
					 name, igraph_i_pajek_actvertex-1,
					 value);
}

int igraph_i_pajek_add_numeric_edge_attribute(const char *name, 
					      igraph_real_t value) {

  return
    igraph_i_pajek_add_numeric_attribute(igraph_i_pajek_edge_attribute_names,
					 igraph_i_pajek_edge_attributes,
					 igraph_i_pajek_actedge,
					 name, igraph_i_pajek_actedge-1,
					 value);
}
