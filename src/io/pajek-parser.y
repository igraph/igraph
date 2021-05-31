/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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
#include <math.h>

#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"
#include "igraph_attributes.h"
#include "config.h"

#include "core/math.h"
#include "io/pajek-header.h"
#include "io/parsers/pajek-parser.h" /* it must come first because of YYSTYPE */
#include "io/parsers/pajek-lexer.h"
#include "internal/hacks.h"

int igraph_pajek_yyerror(YYLTYPE* locp,
                         igraph_i_pajek_parsedata_t *context,
                         const char *s);

int igraph_i_pajek_add_string_vertex_attribute(const char *name,
                                               const char *value,
                                               int len,
                                               igraph_i_pajek_parsedata_t *context);
int igraph_i_pajek_add_string_edge_attribute(const char *name,
                                             const char *value,
                                             int len,
                                             igraph_i_pajek_parsedata_t *context);
int igraph_i_pajek_add_numeric_vertex_attribute(const char *name,
                                                igraph_real_t value,
                                                igraph_i_pajek_parsedata_t *context);
int igraph_i_pajek_add_numeric_edge_attribute(const char *name,
                                              igraph_real_t value,
                                              igraph_i_pajek_parsedata_t *context);
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

int igraph_i_pajek_add_bipartite_type(igraph_i_pajek_parsedata_t *context);
int igraph_i_pajek_check_bipartite(igraph_i_pajek_parsedata_t *context);

extern igraph_real_t igraph_pajek_get_number(const char *str, long int len);
extern long int igraph_i_pajek_actvertex;
extern long int igraph_i_pajek_actedge;

#define scanner context->scanner

%}

%pure-parser
/* bison: do not remove the equals sign; macOS XCode ships with bison 2.3, which
 * needs the equals sign */
%name-prefix="igraph_pajek_yy"
%defines
%locations
%error-verbose
%parse-param { igraph_i_pajek_parsedata_t* context }
%lex-param { void *scanner }

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
%token ERROR

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

input: nethead vertices edgeblock {
  if (context->vcount2 > 0) { igraph_i_pajek_check_bipartite(context); }
 };

nethead: /* empty */ | NETWORKLINE words NEWLINE;

vertices: verticeshead NEWLINE vertdefs;

verticeshead: VERTICESLINE longint {
  context->vcount=$2;
  context->vcount2=0;
            }
            | VERTICESLINE longint longint {
  context->vcount=$2;
  context->vcount2=$3;
  igraph_i_pajek_add_bipartite_type(context);
};

vertdefs: /* empty */  | vertdefs vertexline;

vertexline: NEWLINE |
            vertex NEWLINE |
            vertex { context->actvertex=$1; } vertexid vertexcoords shape params NEWLINE { }
;

vertex: longint { $$=$1; context->mode=1; };

vertexid: word {
  igraph_i_pajek_add_string_vertex_attribute("id", $1.str, $1.len, context);
  igraph_i_pajek_add_string_vertex_attribute("name", $1.str, $1.len, context);
};

vertexcoords: /* empty */
            | number number {
  igraph_i_pajek_add_numeric_vertex_attribute("x", $1, context);
  igraph_i_pajek_add_numeric_vertex_attribute("y", $2, context);
            }
            | number number number {
  igraph_i_pajek_add_numeric_vertex_attribute("x", $1, context);
  igraph_i_pajek_add_numeric_vertex_attribute("y", $2, context);
  igraph_i_pajek_add_numeric_vertex_attribute("z", $3, context);
            };

shape: /* empty */ | word {
  igraph_i_pajek_add_string_vertex_attribute("shape", $1.str, $1.len, context);
};

params: /* empty */ | params param;

param:
       vpword
     | VP_X_FACT number {
         igraph_i_pajek_add_numeric_vertex_attribute("xfact", $2, context);
       }
     | VP_Y_FACT number {
         igraph_i_pajek_add_numeric_vertex_attribute("yfact", $2, context);
       }
     | VP_IC number number number { /* RGB color */
         igraph_i_pajek_add_numeric_vertex_attribute("color-red", $2, context);
         igraph_i_pajek_add_numeric_vertex_attribute("color-green", $3, context);
         igraph_i_pajek_add_numeric_vertex_attribute("color-blue", $4, context);
       }
     | VP_BC number number number {
         igraph_i_pajek_add_numeric_vertex_attribute("framecolor-red", $2, context);
         igraph_i_pajek_add_numeric_vertex_attribute("framecolor-green", $3, context);
         igraph_i_pajek_add_numeric_vertex_attribute("framecolor-blue", $4, context);
       }
     | VP_LC number number number {
         igraph_i_pajek_add_numeric_vertex_attribute("labelcolor-red", $2, context);
         igraph_i_pajek_add_numeric_vertex_attribute("labelcolor-green", $3, context);
         igraph_i_pajek_add_numeric_vertex_attribute("labelcolor-blue", $4, context);
       }
     | VP_LR number {
         igraph_i_pajek_add_numeric_vertex_attribute("labeldist", $2, context);
     }
     | VP_LPHI number {
         igraph_i_pajek_add_numeric_vertex_attribute("labeldegree2", $2, context);
     }
     | VP_BW number {
         igraph_i_pajek_add_numeric_vertex_attribute("framewidth", $2, context);
     }
     | VP_FOS number {
         igraph_i_pajek_add_numeric_vertex_attribute("fontsize", $2, context);
     }
     | VP_PHI number {
         igraph_i_pajek_add_numeric_vertex_attribute("rotation", $2, context);
     }
     | VP_R number {
         igraph_i_pajek_add_numeric_vertex_attribute("radius", $2, context);
     }
     | VP_Q number {
         igraph_i_pajek_add_numeric_vertex_attribute("diamondratio", $2, context);
     }
     | VP_LA number {
         igraph_i_pajek_add_numeric_vertex_attribute("labeldegree", $2, context);
     }
     | VP_SIZE number {
         igraph_i_pajek_add_numeric_vertex_attribute("vertexsize", $2, context);
     }
;

vpword: VP_FONT { context->mode=3; } vpwordpar {
         context->mode=1;
         igraph_i_pajek_add_string_vertex_attribute("font", $3.str, $3.len, context);
     }
     | VP_URL { context->mode=3; } vpwordpar {
         context->mode=1;
         igraph_i_pajek_add_string_vertex_attribute("url", $3.str, $3.len, context);
     }
     | VP_IC { context->mode=3; } vpwordpar {
         context->mode=1;
         igraph_i_pajek_add_string_vertex_attribute("color", $3.str, $3.len, context);
     }
     | VP_BC { context->mode=3; } vpwordpar {
         context->mode=1;
         igraph_i_pajek_add_string_vertex_attribute("framecolor",
                                                    $3.str, $3.len, context);
     }
     | VP_LC { context->mode=3; } vpwordpar {
         context->mode=1;
         igraph_i_pajek_add_string_vertex_attribute("labelcolor",
                                                    $3.str, $3.len, context);
     }
;

vpwordpar: word { $$=$1; };

edgeblock: /* empty */ | edgeblock arcs | edgeblock edges | edgeblock arcslist | edgeblock edgeslist | edgeblock adjmatrix;

arcs:   ARCSLINE NEWLINE arcsdefs        { context->directed=1; }
      | ARCSLINE number NEWLINE arcsdefs { context->directed=1; };

arcsdefs: /* empty */ | arcsdefs arcsline;

arcsline: NEWLINE |
          arcfrom arcto { context->actedge++;
                          context->mode=2; } weight edgeparams NEWLINE  {
  igraph_vector_push_back(context->vector, $1-1);
  igraph_vector_push_back(context->vector, $2-1); }
;

arcfrom: longint;

arcto: longint;

edges:   EDGESLINE NEWLINE edgesdefs { context->directed=0; }
       | EDGESLINE number NEWLINE edgesdefs { context->directed=0; }

edgesdefs: /* empty */ | edgesdefs edgesline;

edgesline: NEWLINE |
          edgefrom edgeto { context->actedge++;
                            context->mode=2; } weight edgeparams NEWLINE {
  igraph_vector_push_back(context->vector, $1-1);
  igraph_vector_push_back(context->vector, $2-1); }
;

edgefrom: longint;

edgeto: longint;

weight: /* empty */ | number {
  igraph_i_pajek_add_numeric_edge_attribute("weight", $1, context);
};

edgeparams: /* empty */ | edgeparams edgeparam;

edgeparam:
     epword
   | EP_C number number number {
       igraph_i_pajek_add_numeric_edge_attribute("color-red", $2, context);
       igraph_i_pajek_add_numeric_edge_attribute("color-green", $3, context);
       igraph_i_pajek_add_numeric_edge_attribute("color-blue", $4, context);
   }
   | EP_S number {
       igraph_i_pajek_add_numeric_edge_attribute("arrowsize", $2, context);
   }
   | EP_W number {
       igraph_i_pajek_add_numeric_edge_attribute("edgewidth", $2, context);
   }
   | EP_H1 number {
       igraph_i_pajek_add_numeric_edge_attribute("hook1", $2, context);
   }
   | EP_H2 number {
       igraph_i_pajek_add_numeric_edge_attribute("hook2", $2, context);
   }
   | EP_A1 number {
       igraph_i_pajek_add_numeric_edge_attribute("angle1", $2, context);
   }
   | EP_A2 number {
       igraph_i_pajek_add_numeric_edge_attribute("angle2", $2, context);
   }
   | EP_K1 number {
       igraph_i_pajek_add_numeric_edge_attribute("velocity1", $2, context);
   }
   | EP_K2 number {
       igraph_i_pajek_add_numeric_edge_attribute("velocity2", $2, context);
   }
   | EP_AP number {
       igraph_i_pajek_add_numeric_edge_attribute("arrowpos", $2, context);
   }
   | EP_LP number {
       igraph_i_pajek_add_numeric_edge_attribute("labelpos", $2, context);
   }
   | EP_LR number {
       igraph_i_pajek_add_numeric_edge_attribute("labelangle", $2, context);
   }
   | EP_LPHI number {
       igraph_i_pajek_add_numeric_edge_attribute("labelangle2", $2, context);
   }
   | EP_LA number {
       igraph_i_pajek_add_numeric_edge_attribute("labeldegree", $2, context);
   }
   | EP_SIZE number { /* what is this??? */
       igraph_i_pajek_add_numeric_edge_attribute("arrowsize", $2, context);
   }
   | EP_FOS number {
       igraph_i_pajek_add_numeric_edge_attribute("fontsize", $2, context);
   }
;

epword: EP_A { context->mode=4; } epwordpar {
      context->mode=2;
      igraph_i_pajek_add_string_edge_attribute("arrowtype", $3.str, $3.len, context);
    }
    | EP_P { context->mode=4; } epwordpar {
      context->mode=2;
      igraph_i_pajek_add_string_edge_attribute("linepattern", $3.str, $3.len, context);
    }
    | EP_L { context->mode=4; } epwordpar {
      context->mode=2;
      igraph_i_pajek_add_string_edge_attribute("label", $3.str, $3.len, context);
    }
    | EP_LC { context->mode=4; } epwordpar {
      context->mode=2;
      igraph_i_pajek_add_string_edge_attribute("labelcolor", $3.str, $3.len, context);
    }
    | EP_C { context->mode=4; } epwordpar {
      context->mode=2;
      igraph_i_pajek_add_string_edge_attribute("color", $3.str, $3.len, context);
    }
;

epwordpar: word { context->mode=2; $$=$1; };

arcslist: ARCSLISTLINE NEWLINE arcslistlines { context->directed=1; };

arcslistlines: /* empty */ | arcslistlines arclistline;

arclistline: NEWLINE | arclistfrom arctolist NEWLINE;

arctolist: /* empty */ | arctolist arclistto;

arclistfrom: longint { context->mode=0; context->actfrom=labs($1)-1; };

arclistto: longint {
  igraph_vector_push_back(context->vector, context->actfrom);
  igraph_vector_push_back(context->vector, labs($1)-1);
};

edgeslist: EDGESLISTLINE NEWLINE edgelistlines { context->directed=0; };

edgelistlines: /* empty */ | edgelistlines edgelistline;

edgelistline: NEWLINE | edgelistfrom edgetolist NEWLINE;

edgetolist: /* empty */ | edgetolist edgelistto;

edgelistfrom: longint { context->mode=0; context->actfrom=labs($1)-1; };

edgelistto: longint {
  igraph_vector_push_back(context->vector, context->actfrom);
  igraph_vector_push_back(context->vector, labs($1)-1);
};

/* -----------------------------------------------------*/

adjmatrix: matrixline NEWLINE adjmatrixlines;

matrixline: MATRIXLINE { context->actfrom=0;
                         context->actto=0;
                         context->directed=(context->vcount2==0);
                       };

adjmatrixlines: /* empty */ | adjmatrixlines adjmatrixline;

adjmatrixline: adjmatrixnumbers NEWLINE { context->actfrom++; context->actto=0; };

adjmatrixnumbers: /* empty */ | adjmatrixentry adjmatrixnumbers;

adjmatrixentry: number {
  if ($1 != 0) {
    if (context->vcount2==0) {
      context->actedge++;
      igraph_i_pajek_add_numeric_edge_attribute("weight", $1, context);
      igraph_vector_push_back(context->vector, context->actfrom);
      igraph_vector_push_back(context->vector, context->actto);
    } else if (context->vcount2 + context->actto < context->vcount) {
      context->actedge++;
      igraph_i_pajek_add_numeric_edge_attribute("weight", $1, context);
      igraph_vector_push_back(context->vector, context->actfrom);
      igraph_vector_push_back(context->vector,
                              context->vcount2+context->actto);
    }
  }
  context->actto++;
};

/* -----------------------------------------------------*/

longint: NUM { $$=igraph_pajek_get_number(igraph_pajek_yyget_text(scanner),
                                          igraph_pajek_yyget_leng(scanner)); };

number: NUM  { $$=igraph_pajek_get_number(igraph_pajek_yyget_text(scanner),
                                          igraph_pajek_yyget_leng(scanner)); };

words: /* empty */ | words word;

word: ALNUM { $$.str=igraph_pajek_yyget_text(scanner);
              $$.len=igraph_pajek_yyget_leng(scanner); }
      | NUM { $$.str=igraph_pajek_yyget_text(scanner);
              $$.len=igraph_pajek_yyget_leng(scanner); }
      | QSTR { $$.str=igraph_pajek_yyget_text(scanner)+1;
               $$.len=igraph_pajek_yyget_leng(scanner)-2; };

%%

int igraph_pajek_yyerror(YYLTYPE* locp,
                         igraph_i_pajek_parsedata_t *context,
                         const char *s) {
  snprintf(context->errmsg, sizeof(context->errmsg)/sizeof(char)-1,
           "Parse error in Pajek file, line %i (%s)",
           locp->first_line, s);
  return 0;
}

igraph_real_t igraph_pajek_get_number(const char *str, long int length) {
  igraph_real_t num;
  char *tmp=IGRAPH_CALLOC(length+1, char);

  strncpy(tmp, str, length);
  tmp[length]='\0';
  sscanf(tmp, "%lf", &num);
  IGRAPH_FREE(tmp);
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
  igraph_attribute_record_t *rec;

  igraph_trie_get(names, attrname, &id);
  if (id == attrsize) {
    /* add a new attribute */
    rec=IGRAPH_CALLOC(1, igraph_attribute_record_t);
    na=IGRAPH_CALLOC(1, igraph_vector_t);
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
  igraph_attribute_record_t *rec;
  long int i;

  igraph_trie_get(names, attrname, &id);
  if (id == attrsize) {
    /* add a new attribute */
    rec=IGRAPH_CALLOC(1, igraph_attribute_record_t);
    na=IGRAPH_CALLOC(1, igraph_strvector_t);
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
                                               int len,
                                               igraph_i_pajek_parsedata_t *context) {
  char *tmp;
  int ret;

  tmp=IGRAPH_CALLOC(len+1, char);
  if (tmp==0) {
    IGRAPH_ERROR("cannot add element to hash table", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmp);
  strncpy(tmp, value, len);
  tmp[len]='\0';

  ret=igraph_i_pajek_add_string_attribute(context->vertex_attribute_names,
                                          context->vertex_attributes,
                                          context->vcount,
                                          name, context->actvertex-1,
                                          tmp);

  IGRAPH_FREE(tmp);
  IGRAPH_FINALLY_CLEAN(1);

  return ret;
}

int igraph_i_pajek_add_string_edge_attribute(const char *name,
                                             const char *value,
                                             int len,
                                             igraph_i_pajek_parsedata_t *context) {
  char *tmp;
  int ret;

  tmp=IGRAPH_CALLOC(len+1, char);
  if (tmp==0) {
    IGRAPH_ERROR("cannot add element to hash table", IGRAPH_ENOMEM);
  }
  IGRAPH_FINALLY(igraph_free, tmp);
  strncpy(tmp, value, len);
  tmp[len]='\0';

  ret=igraph_i_pajek_add_string_attribute(context->edge_attribute_names,
                                          context->edge_attributes,
                                          context->actedge,
                                          name, context->actedge-1,
                                          tmp);

  IGRAPH_FREE(tmp);
  IGRAPH_FINALLY_CLEAN(1);

  return ret;
}

int igraph_i_pajek_add_numeric_vertex_attribute(const char *name,
                                                igraph_real_t value,
                                                igraph_i_pajek_parsedata_t *context) {

  return
    igraph_i_pajek_add_numeric_attribute(context->vertex_attribute_names,
                                         context->vertex_attributes,
                                         context->vcount,
                                         name, context->actvertex-1,
                                         value);
}

int igraph_i_pajek_add_numeric_edge_attribute(const char *name,
                                              igraph_real_t value,
                                              igraph_i_pajek_parsedata_t *context) {

  return
    igraph_i_pajek_add_numeric_attribute(context->edge_attribute_names,
                                         context->edge_attributes,
                                         context->actedge,
                                         name, context->actedge-1,
                                         value);
}

int igraph_i_pajek_add_bipartite_type(igraph_i_pajek_parsedata_t *context) {

  const char *attrname="type";
  igraph_trie_t *names=context->vertex_attribute_names;
  igraph_vector_ptr_t *attrs=context->vertex_attributes;
  int i, n=context->vcount, n1=context->vcount2;
  long int attrid, attrsize=igraph_trie_size(names);
  igraph_attribute_record_t *rec;
  igraph_vector_t *na;

  if (n1 > n) {
    IGRAPH_ERROR("Invalid number of vertices in bipartite Pajek file",
                 IGRAPH_PARSEERROR);
  }

  igraph_trie_get(names, attrname, &attrid);
  if (attrid != attrsize) {
    IGRAPH_ERROR("Duplicate 'type' attribute in Pajek file, "
                 "this should not happen", IGRAPH_EINTERNAL);
  }

  /* add a new attribute */
  rec=IGRAPH_CALLOC(1, igraph_attribute_record_t);
  na=IGRAPH_CALLOC(1, igraph_vector_t);
  igraph_vector_init(na, n);
  rec->name=strdup(attrname);
  rec->type=IGRAPH_ATTRIBUTE_NUMERIC;
  rec->value=na;
  igraph_vector_ptr_push_back(attrs, rec);

  for (i=0; i<n1; i++) {
    VECTOR(*na)[i] = 0;
  }
  for (i=n1; i<n; i++) {
    VECTOR(*na)[i] = 1;
  }

  return 0;
}

int igraph_i_pajek_check_bipartite(igraph_i_pajek_parsedata_t *context) {
  const igraph_vector_t *edges=context->vector;
  int i, n1=context->vcount2;
  int ne=igraph_vector_size(edges);

  for (i=0; i<ne; i+=2) {
    int v1=VECTOR(*edges)[i];
    int v2=VECTOR(*edges)[i+1];
    if ( (v1 < n1 && v2 < n1) || (v1 > n1 && v2 > n1) ) {
      IGRAPH_WARNING("Invalid edge in bipartite graph");
    }
  }

  return 0;
}
