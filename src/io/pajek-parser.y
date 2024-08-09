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

#include "igraph_attributes.h"
#include "igraph_error.h"
#include "igraph_memory.h"
#include "igraph_types.h"

#include "io/pajek-header.h"
#include "io/parsers/pajek-parser.h" /* it must come first because of YYSTYPE */
#include "io/parsers/pajek-lexer.h"
#include "io/parse_utils.h"
#include "internal/hacks.h" /* strdup */

#include <stdio.h>
#include <string.h>
#include <math.h>

int igraph_pajek_yyerror(YYLTYPE* locp,
                         igraph_i_pajek_parsedata_t *context,
                         const char *s);

static igraph_error_t add_string_vertex_attribute(const char *name,
                                               const char *value,
                                               size_t len,
                                               igraph_i_pajek_parsedata_t *context);
static igraph_error_t add_string_edge_attribute(const char *name,
                                             const char *value,
                                             size_t len,
                                             igraph_i_pajek_parsedata_t *context);
static igraph_error_t add_numeric_vertex_attribute(const char *name,
                                                igraph_real_t value,
                                                igraph_i_pajek_parsedata_t *context);
static igraph_error_t add_numeric_edge_attribute(const char *name,
                                              igraph_real_t value,
                                              igraph_i_pajek_parsedata_t *context);
static igraph_error_t add_numeric_attribute(igraph_trie_t *names,
                                         igraph_vector_ptr_t *attrs,
                                         igraph_integer_t count,
                                         const char *attrname,
                                         igraph_integer_t vid,
                                         igraph_real_t number);
static igraph_error_t add_string_attribute(igraph_trie_t *names,
                                        igraph_vector_ptr_t *attrs,
                                        igraph_integer_t count,
                                        const char *attrname,
                                        igraph_integer_t vid,
                                        const char *str,
                                        igraph_integer_t str_len);

static igraph_error_t add_bipartite_type(igraph_i_pajek_parsedata_t *context);
static igraph_error_t check_bipartite(igraph_i_pajek_parsedata_t *context);

static igraph_error_t make_dynstr(const char *src, size_t len, char **res);
static igraph_bool_t is_standard_vattr(const char *attrname);
static igraph_bool_t is_standard_eattr(const char *attrname);
static igraph_error_t deconflict_attrname(char **attrname);

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
  igraph_integer_t intnum;
  igraph_real_t    realnum;
  struct {
    char *str;
    size_t len;
  } string;
  char *dynstr;
}

%type <intnum>   integer;
%type <intnum>   vertex;
%type <realnum>  number;
%type <string>   word;
%type <string>   parstrval;
%type <dynstr>   parname;

%destructor { free($$); } parname;

%token NEWLINE       "end of line"
%token NUM           "number"
%token ALNUM         "word"
%token QSTR          "quoted string"
%token NETWORKLINE   "*Network line"
%token VERTICESLINE  "*Vertices line"
%token ARCSLINE      "*Arcs line"
%token EDGESLINE     "*Edges line"
%token ARCSLISTLINE  "*Arcslist line"
%token EDGESLISTLINE "*Edgeslist line"
%token MATRIXLINE    "*Matrix line"
%token END 0         "end of file" /* friendly name for $end */
%token ERROR

%token VP_X_FACT
%token VP_Y_FACT
%token VP_PHI
%token VP_R
%token VP_Q
%token VP_IC
%token VP_BC
%token VP_BW
%token VP_LC
%token VP_LA
%token VP_LR
%token VP_LPHI
%token VP_FOS
%token VP_FONT
%token VP_URL

%token EP_H1
%token EP_H2
%token EP_W
%token EP_C
%token EP_P
%token EP_A
%token EP_S
%token EP_A1
%token EP_K1
%token EP_A2
%token EP_K2
%token EP_AP
%token EP_L
%token EP_LP
%token EP_LR
%token EP_LPHI
%token EP_LC
%token EP_LA
%token EP_FOS
%token EP_FONT

%%

input: nethead vertices edgeblock final_newlines {
  if (context->vcount2 > 0) { check_bipartite(context); }
  if (! context->eof) {
    /* In Pajek files, an empty line after *Vertices signifies the end of the network data.
     * If there is more data after one or more empty lines, we warn the user, as this
     * may indicate file corruption, for example a stray empty lines before *Edges. */
    IGRAPH_WARNINGF("Empty line encountered, ignoring rest of file after line %d.", @4.first_line);
  }
  YYACCEPT; /* stop parsing even if there is more data in the file. */
 };

final_newlines: /* empty */ | NEWLINE final_newlines ;

nethead: /* empty */ | NETWORKLINE ;

vertices: verticeshead NEWLINE vertdefs ;

verticeshead: VERTICESLINE integer {
  context->vcount=$2;
  context->vcount2=0;
  if (context->vcount < 0) {
    IGRAPH_YY_ERRORF("Invalid vertex count in Pajek file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->vcount);
  }
  if (context->vcount > IGRAPH_PAJEK_MAX_VERTEX_COUNT) {
    IGRAPH_YY_ERRORF("Vertex count too large in Pajek file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->vcount);
  }
            }
            | VERTICESLINE integer integer {
  context->vcount=$2;
  context->vcount2=$3;
  if (context->vcount < 0) {
    IGRAPH_YY_ERRORF("Invalid vertex count in Pajek file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->vcount);
  }
  if (context->vcount > IGRAPH_PAJEK_MAX_VERTEX_COUNT) {
    IGRAPH_YY_ERRORF("Vertex count too large in Pajek file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->vcount);
  }
  if (context->vcount2 < 0) {
    IGRAPH_YY_ERRORF("Invalid two-mode vertex count in Pajek file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->vcount2);
  }
  if (context->vcount2 > IGRAPH_PAJEK_MAX_VERTEX_COUNT) {
    IGRAPH_YY_ERRORF("2-mode vertex count too large in Pajek file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->vcount2);
  }
  IGRAPH_YY_CHECK(add_bipartite_type(context));
};

vertdefs: /* empty */  | vertdefs vertexline;

vertexline: vertex NEWLINE |
            vertex { context->actvertex=$1; } vertexid vertexcoords shape vertparams NEWLINE { }
;

vertex: integer {
  igraph_integer_t v = $1;
  /* Per feedback from Pajek's authors, negative signs should be ignored for vertex IDs.
   * See https://nascol.discourse.group/t/pajek-arcslist-edgelist-format/44/2
   * This applies to all of *Edges, *Arcs, *Edgeslist, *Arcslist and *Vertices section.
   * IGRAPH_INTEGER_MIN cannot be negated on typical platforms so we keep it as-is.
   */
  if (v < 0 && v > IGRAPH_INTEGER_MIN) {
    v = -v;
  }
  if (v < 1 || v > context->vcount) {
      IGRAPH_YY_ERRORF(
                  "Invalid vertex ID (%" IGRAPH_PRId ") in Pajek file. "
                  "The number of vertices is %" IGRAPH_PRId ".",
                  IGRAPH_EINVAL, v, context->vcount);
  }
  $$ = v;
};

vertexid: word {
  IGRAPH_YY_CHECK(add_string_vertex_attribute("id", $1.str, $1.len, context));
  IGRAPH_YY_CHECK(add_string_vertex_attribute("name", $1.str, $1.len, context));
};

vertexcoords: /* empty */
            | number number {
  IGRAPH_YY_CHECK(add_numeric_vertex_attribute("x", $1, context));
  IGRAPH_YY_CHECK(add_numeric_vertex_attribute("y", $2, context));
            }
            | number number number {
  IGRAPH_YY_CHECK(add_numeric_vertex_attribute("x", $1, context));
  IGRAPH_YY_CHECK(add_numeric_vertex_attribute("y", $2, context));
  IGRAPH_YY_CHECK(add_numeric_vertex_attribute("z", $3, context));
            };

shape: /* empty */ | word {
  IGRAPH_YY_CHECK(add_string_vertex_attribute("shape", $1.str, $1.len, context));
};

vertparams: /* empty */ | vertparams vertparam;

vertparam:
       vpword
     | VP_X_FACT number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("xfact", $2, context));
       }
     | VP_Y_FACT number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("yfact", $2, context));
       }
     | VP_LR number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("labeldist", $2, context));
     }
     | VP_LPHI number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("labeldegree2", $2, context));
     }
     | VP_BW number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("framewidth", $2, context));
     }
     | VP_FOS number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("fontsize", $2, context));
     }
     | VP_PHI number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("rotation", $2, context));
     }
     | VP_R number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("radius", $2, context));
     }
     | VP_Q number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("diamondratio", $2, context));
     }
     | VP_LA number {
         IGRAPH_YY_CHECK(add_numeric_vertex_attribute("labeldegree", $2, context));
     }
;

vpword: VP_FONT parstrval {
         IGRAPH_YY_CHECK(add_string_vertex_attribute("font", $2.str, $2.len, context));
     }
     | VP_URL parstrval {
         IGRAPH_YY_CHECK(add_string_vertex_attribute("url", $2.str, $2.len, context));
     }
     | VP_IC parstrval {
         IGRAPH_YY_CHECK(add_string_vertex_attribute("color", $2.str, $2.len, context));
     }
     | VP_BC parstrval {
         IGRAPH_YY_CHECK(add_string_vertex_attribute("framecolor", $2.str, $2.len, context));
     }
     | VP_LC parstrval {
         IGRAPH_YY_CHECK(add_string_vertex_attribute("labelcolor", $2.str, $2.len, context));
     }
     | parname parstrval {
         IGRAPH_FINALLY(igraph_free, $1);
         if (is_standard_vattr($1)) {
          IGRAPH_YY_CHECK(deconflict_attrname(&$1));
          /* update address on finally stack */
          IGRAPH_FINALLY_CLEAN(1);
          IGRAPH_FINALLY(igraph_free, $1);
         }
         IGRAPH_YY_CHECK(add_string_vertex_attribute(
           $1, $2.str, $2.len, context));
         IGRAPH_FREE($1);
         IGRAPH_FINALLY_CLEAN(1);
     }
;

edgeblock: /* empty */ | edgeblock arcs | edgeblock edges | edgeblock arcslist | edgeblock edgeslist | edgeblock adjmatrix;

arcs:   ARCSLINE NEWLINE arcsdefs        { context->directed=true; }
      | ARCSLINE number NEWLINE arcsdefs { context->directed=true; };

arcsdefs: /* empty */ | arcsdefs arcsline;

arcsline: vertex vertex { context->actedge++; } weight edgeparams NEWLINE  {
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $1-1));
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $2-1)); }
;

edges:   EDGESLINE NEWLINE edgesdefs { context->directed=0; }
       | EDGESLINE number NEWLINE edgesdefs { context->directed=0; }

edgesdefs: /* empty */ | edgesdefs edgesline;

edgesline: vertex vertex { context->actedge++; } weight edgeparams NEWLINE {
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $1-1));
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $2-1)); }
;

weight: /* empty */ | number {
  IGRAPH_YY_CHECK(add_numeric_edge_attribute("weight", $1, context));
};

edgeparams: /* empty */ | edgeparams edgeparam;

edgeparam:
     epword
   | EP_S number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("arrowsize", $2, context));
   }
   | EP_W number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("edgewidth", $2, context));
   }
   | EP_H1 number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("hook1", $2, context));
   }
   | EP_H2 number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("hook2", $2, context));
   }
   | EP_A1 number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("angle1", $2, context));
   }
   | EP_A2 number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("angle2", $2, context));
   }
   | EP_K1 number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("velocity1", $2, context));
   }
   | EP_K2 number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("velocity2", $2, context));
   }
   | EP_AP number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("arrowpos", $2, context));
   }
   | EP_LP number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("labelpos", $2, context));
   }
   | EP_LR number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("labelangle", $2, context));
   }
   | EP_LPHI number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("labelangle2", $2, context));
   }
   | EP_LA number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("labeldegree", $2, context));
   }
   | EP_FOS number {
       IGRAPH_YY_CHECK(add_numeric_edge_attribute("fontsize", $2, context));
   }
;

epword: EP_A parstrval {
      IGRAPH_YY_CHECK(add_string_edge_attribute("arrowtype", $2.str, $2.len, context));
    }
    | EP_P parstrval {
      IGRAPH_YY_CHECK(add_string_edge_attribute("linepattern", $2.str, $2.len, context));
    }
    | EP_L parstrval {
      IGRAPH_YY_CHECK(add_string_edge_attribute("label", $2.str, $2.len, context));
    }
    | EP_LC parstrval {
      IGRAPH_YY_CHECK(add_string_edge_attribute("labelcolor", $2.str, $2.len, context));
    }
    | EP_C parstrval {
      IGRAPH_YY_CHECK(add_string_edge_attribute("color", $2.str, $2.len, context));
    }
    | EP_FONT parstrval {
      IGRAPH_YY_CHECK(add_string_edge_attribute("font", $2.str, $2.len, context));
    }
    | parname parstrval {
        IGRAPH_FINALLY(igraph_free, $1);
        if (is_standard_eattr($1)) {
          IGRAPH_YY_CHECK(deconflict_attrname(&$1));
          /* update address on finally stack */
          IGRAPH_FINALLY_CLEAN(1);
          IGRAPH_FINALLY(igraph_free, $1);
        }
        IGRAPH_YY_CHECK(add_string_edge_attribute(
           $1, $2.str, $2.len, context));
        IGRAPH_FREE($1);
        IGRAPH_FINALLY_CLEAN(1);
     }
;

arcslist: ARCSLISTLINE NEWLINE arcslistlines { context->directed=true; };

arcslistlines: /* empty */ | arcslistlines arclistline;

arclistline: arclistfrom arctolist NEWLINE;

arctolist: /* empty */ | arctolist arclistto;

arclistfrom: vertex { context->actfrom=$1-1; };

arclistto: vertex {
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, context->actfrom));
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $1-1));
};

edgeslist: EDGESLISTLINE NEWLINE edgelistlines { context->directed=0; };

edgelistlines: /* empty */ | edgelistlines edgelistline;

edgelistline: edgelistfrom edgetolist NEWLINE;

edgetolist: /* empty */ | edgetolist edgelistto;

edgelistfrom: vertex { context->actfrom=$1-1; };

edgelistto: vertex {
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, context->actfrom));
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, $1-1));
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
      IGRAPH_YY_CHECK(add_numeric_edge_attribute("weight", $1, context));
      IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, context->actfrom));
      IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, context->actto));
    } else if (context->vcount2 + context->actto < context->vcount) {
      context->actedge++;
      IGRAPH_YY_CHECK(add_numeric_edge_attribute("weight", $1, context));
      IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector, context->actfrom));
      IGRAPH_YY_CHECK(igraph_vector_int_push_back(context->vector,
                              context->vcount2+context->actto));
    }
  }
  context->actto++;
};

/* -----------------------------------------------------*/

integer: NUM {
  igraph_integer_t val;
  IGRAPH_YY_CHECK(igraph_i_parse_integer(igraph_pajek_yyget_text(scanner),
                                         igraph_pajek_yyget_leng(scanner),
                                         &val));
  $$=val;
};

number: NUM  {
  igraph_real_t val;
  IGRAPH_YY_CHECK(igraph_i_parse_real(igraph_pajek_yyget_text(scanner),
                                      igraph_pajek_yyget_leng(scanner),
                                      &val));
  $$=val;
};

parname: word {
  IGRAPH_YY_CHECK(make_dynstr($1.str, $1.len, &$$));
};

parstrval: word { $$=$1; };

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

/* TODO: NA's */

static igraph_error_t add_numeric_attribute(igraph_trie_t *names,
                                            igraph_vector_ptr_t *attrs,
                                            igraph_integer_t count,
                                            const char *attrname,
                                            igraph_integer_t elem_id,
                                            igraph_real_t number) {
  igraph_integer_t attrsize = igraph_trie_size(names);
  igraph_integer_t id;
  igraph_vector_t *na;
  igraph_attribute_record_t *rec;

  IGRAPH_CHECK(igraph_trie_get(names, attrname, &id));
  if (id == attrsize) {
    /* add a new attribute */
    rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
    CHECK_OOM_RP(rec);
    IGRAPH_FINALLY(igraph_free, rec);

    na = IGRAPH_CALLOC(1, igraph_vector_t);
    CHECK_OOM_RP(na);
    IGRAPH_FINALLY(igraph_free, na);
    IGRAPH_VECTOR_INIT_FINALLY(na, count);

    rec->name = strdup(attrname);
    CHECK_OOM_RP(rec->name);
    IGRAPH_FINALLY(igraph_free, (void *) rec->name);

    rec->type = IGRAPH_ATTRIBUTE_NUMERIC;
    rec->value = na;

    IGRAPH_CHECK(igraph_vector_ptr_push_back(attrs, rec));
    IGRAPH_FINALLY_CLEAN(4); /* ownership of rec transferred to attrs */
  }

  rec = VECTOR(*attrs)[id];
  na = (igraph_vector_t *) rec->value;
  if (igraph_vector_size(na) == elem_id) {
    IGRAPH_CHECK(igraph_vector_push_back(na, number));
  } else if (igraph_vector_size(na) < elem_id) {
    igraph_integer_t origsize=igraph_vector_size(na);
    IGRAPH_CHECK(igraph_vector_resize(na, elem_id+1));
    for (;origsize<count; origsize++) {
      VECTOR(*na)[origsize] = IGRAPH_NAN;
    }
    VECTOR(*na)[elem_id] = number;
  } else {
    VECTOR(*na)[elem_id] = number;
  }

  return IGRAPH_SUCCESS;
}

/* TODO: NA's */

static igraph_error_t make_dynstr(const char *src, size_t len, char **res) {
  *res = strndup(src, len);
  CHECK_OOM_RP(*res);
  return IGRAPH_SUCCESS;
}

static igraph_error_t add_string_attribute(igraph_trie_t *names,
                                           igraph_vector_ptr_t *attrs,
                                           igraph_integer_t count,
                                           const char *attrname,
                                           igraph_integer_t elem_id,
                                           const char *str,
                                           igraph_integer_t str_len) {
  igraph_integer_t attrsize=igraph_trie_size(names);
  igraph_integer_t id;
  igraph_strvector_t *na;
  igraph_attribute_record_t *rec;

  if (attrname[0] == '\0') {
    /* This is relevant only for custom attributes, which are always of string type.
       No need to add the same check for numerical attributes. */
    IGRAPH_ERROR("\"\" is not allowed as a parameter name in Pajek files.", IGRAPH_PARSEERROR);
  }

  IGRAPH_CHECK(igraph_trie_get(names, attrname, &id));
  if (id == attrsize) {

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    /* There are 21 standard vertex attributes and 21 standard edge attributes.
     * We refuse to allow more to reduce memory usage when fuzzing. */
    if (attrsize > 21) {
      IGRAPH_ERROR("Too many attributes in Pajek file.", IGRAPH_PARSEERROR);
    }
#endif

    /* add a new attribute */
    rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
    CHECK_OOM_RP(rec);
    IGRAPH_FINALLY(igraph_free, rec);

    na = IGRAPH_CALLOC(1, igraph_strvector_t);
    CHECK_OOM_RP(na);
    IGRAPH_FINALLY(igraph_free, na);
    IGRAPH_STRVECTOR_INIT_FINALLY(na, count);

    rec->name = strdup(attrname);
    CHECK_OOM_RP(rec->name);
    IGRAPH_FINALLY(igraph_free, (char *) rec->name);

    rec->type = IGRAPH_ATTRIBUTE_STRING;
    rec->value = na;

    IGRAPH_CHECK(igraph_vector_ptr_push_back(attrs, rec));
    IGRAPH_FINALLY_CLEAN(4); /* ownership of rec transferred to attrs */
  }

  rec = VECTOR(*attrs)[id];
  na = (igraph_strvector_t *) rec->value;
  if (igraph_strvector_size(na) <= elem_id) {
    IGRAPH_CHECK(igraph_strvector_resize(na, elem_id+1));
  }
  IGRAPH_CHECK(igraph_strvector_set_len(na, elem_id, str, str_len));

  return IGRAPH_SUCCESS;
}

static igraph_error_t add_string_vertex_attribute(const char *name,
                                                  const char *value,
                                                  size_t len,
                                                  igraph_i_pajek_parsedata_t *context) {

  return add_string_attribute(context->vertex_attribute_names,
                              context->vertex_attributes,
                              context->vcount,
                              name, context->actvertex-1,
                              value, len);
}

static igraph_error_t add_string_edge_attribute(const char *name,
                                                const char *value,
                                                size_t len,
                                                igraph_i_pajek_parsedata_t *context) {

  return add_string_attribute(context->edge_attribute_names,
                              context->edge_attributes,
                              context->actedge,
                              name, context->actedge-1,
                              value, len);
}

static igraph_error_t add_numeric_vertex_attribute(const char *name,
                                                   igraph_real_t value,
                                                   igraph_i_pajek_parsedata_t *context) {

  return add_numeric_attribute(context->vertex_attribute_names,
                               context->vertex_attributes,
                               context->vcount,
                               name, context->actvertex-1,
                               value);
}

static igraph_error_t add_numeric_edge_attribute(const char *name,
                                                 igraph_real_t value,
                                                 igraph_i_pajek_parsedata_t *context) {

  return add_numeric_attribute(context->edge_attribute_names,
                               context->edge_attributes,
                               context->actedge,
                               name, context->actedge-1,
                               value);
}

static igraph_error_t add_bipartite_type(igraph_i_pajek_parsedata_t *context) {

  const char *attrname="type";
  igraph_trie_t *names=context->vertex_attribute_names;
  igraph_vector_ptr_t *attrs=context->vertex_attributes;
  igraph_integer_t n=context->vcount, n1=context->vcount2;
  igraph_integer_t attrid, attrsize = igraph_trie_size(names);
  igraph_attribute_record_t *rec;
  igraph_vector_bool_t *na;

  if (n1 > n) {
    IGRAPH_ERROR("Invalid number of vertices in bipartite Pajek file.",
                 IGRAPH_PARSEERROR);
  }

  IGRAPH_CHECK(igraph_trie_get(names, attrname, &attrid));

  /* It should not be possible for the "type" attribute to be already
   * present at this point. */
  IGRAPH_ASSERT(attrid == attrsize);

  /* add a new attribute */
  rec = IGRAPH_CALLOC(1, igraph_attribute_record_t);
  CHECK_OOM_RP(rec);
  IGRAPH_FINALLY(igraph_free, rec);

  na = IGRAPH_CALLOC(1, igraph_vector_bool_t);
  CHECK_OOM_RP(na);
  IGRAPH_FINALLY(igraph_free, na);
  IGRAPH_VECTOR_BOOL_INIT_FINALLY(na, n);

  rec->name = strdup(attrname);
  CHECK_OOM_RP(rec->name);
  IGRAPH_FINALLY(igraph_free, (char *) rec->name);

  rec->type = IGRAPH_ATTRIBUTE_BOOLEAN;
  rec->value = na;

  IGRAPH_CHECK(igraph_vector_ptr_push_back(attrs, rec));
  IGRAPH_FINALLY_CLEAN(4); /* ownership of 'rec' transferred to 'attrs' */

  for (igraph_integer_t i=0; i<n1; i++) {
    VECTOR(*na)[i] = false;
  }
  for (igraph_integer_t i=n1; i<n; i++) {
    VECTOR(*na)[i] = true;
  }

  return IGRAPH_SUCCESS;
}

static igraph_error_t check_bipartite(igraph_i_pajek_parsedata_t *context) {
  const igraph_vector_int_t *edges=context->vector;
  igraph_integer_t n1=context->vcount2;
  igraph_integer_t ne=igraph_vector_int_size(edges);

  for (igraph_integer_t i=0; i<ne; i+=2) {
    igraph_integer_t v1 = VECTOR(*edges)[i];
    igraph_integer_t v2 = VECTOR(*edges)[i+1];
    if ( (v1 < n1 && v2 < n1) || (v1 > n1 && v2 > n1) ) {
      IGRAPH_WARNING("Invalid edge in bipartite graph.");
    }
  }

  return IGRAPH_SUCCESS;
}

/* Check if attrname is a standard vertex attribute name used by igraph
   for Pajek data. All of these must be listed here to prevent overwriting
   standard attributes, or crashes due to incompatible attribute types. */
static igraph_bool_t is_standard_vattr(const char *attrname) {
  const char *names[] = {
    /* vertex names: */
    "id", /* TODO: remove for 0.11 */ "name",
    /* other vertex attributes: */
    "type", "x", "y", "z",
    /* vertex parameters: */
    "xfact", "yfact",
    "labeldist", "labeldegree2", "framewidth",
    "fontsize", "rotation", "radius",
    "diamondratio", "labeldegree",
    "font", "url", "color", "framecolor",
    "labelcolor"
  };
  for (size_t i=0; i < sizeof(names) / sizeof(names[0]); i++) {
    if (strcmp(attrname, names[i]) == 0) {
      return true;
    }
  }
  return false;
}

/* Check if attrname is a standard edge attribute name used by igraph
   for Pajek data. All of these must be listed here to prevent overwriting
   standard attributes, or crashes due to incompatible attribute types. */
static igraph_bool_t is_standard_eattr(const char *attrname) {
  const char *names[] = {
    /* other edge attributes: */
    "weight",
    /* edge parameters: */
    "arrowsize", "edgewidth", "hook1", "hook2",
    "angle1", "angle2", "velocity1", "velocity2",
    "arrowpos", "labelpos", "labelangle",
    "labelangle2", "labeldegree", "fontsize", "font",
    "arrowtype", "linepattern", "label", "labelcolor",
    "color"
  };
  for (size_t i=0; i < sizeof(names) / sizeof(names[0]); i++) {
    if (strcmp(attrname, names[i]) == 0) {
      return true;
    }
  }
  return false;
}

/* Add a _ character at the end of an attribute name to avoid conflict
 * with standard Pajek attributes. */
static igraph_error_t deconflict_attrname(char **attrname) {
  size_t len = strlen(*attrname);
  char *tmp = IGRAPH_REALLOC(*attrname, len+2, char);
  CHECK_OOM_RP(tmp);
  tmp[len] = '_';
  tmp[len+1] = '\0';
  *attrname = tmp;
  return IGRAPH_SUCCESS;
}
