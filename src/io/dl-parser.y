/*
   IGraph library.
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "internal/hacks.h"
#include "io/dl-header.h"
#include "io/parsers/dl-parser.h"
#include "io/parsers/dl-lexer.h"
#include "io/parse_utils.h"

int igraph_dl_yyerror(YYLTYPE* locp, igraph_i_dl_parsedata_t* context,
                      const char *s);
static igraph_error_t igraph_i_dl_add_str(char *newstr, yy_size_t length,
                        igraph_i_dl_parsedata_t *context);
static igraph_error_t igraph_i_dl_add_edge(igraph_integer_t from, igraph_integer_t to,
                         igraph_i_dl_parsedata_t *context);
static igraph_error_t igraph_i_dl_add_edge_w(igraph_integer_t from, igraph_integer_t to,
                           igraph_real_t weight,
                           igraph_i_dl_parsedata_t *context);
static igraph_error_t igraph_i_dl_check_vid(igraph_integer_t dl_vid);

#define scanner context->scanner

%}

%pure-parser
/* bison: do not remove the equals sign; macOS XCode ships with bison 2.3, which
 * needs the equals sign */
%name-prefix="igraph_dl_yy"
%defines
%locations
%error-verbose
%parse-param { igraph_i_dl_parsedata_t* context }
%lex-param { void* scanner }

%union {
  igraph_integer_t integer;
  igraph_real_t real;
};

%type <integer> integer elabel;
%type <real> weight;

%token NUM              "number"
%token NEWLINE          "end of line"
%token DL               "DL"
%token NEQ              "n=vertexcount"
%token DATA             "data:"
%token LABELS           "labels:"
%token LABELSEMBEDDED   "labels embedded:"
%token FORMATFULLMATRIX
%token FORMATEDGELIST1
%token FORMATNODELIST1
%token DIGIT            "binary digit"
%token LABEL            "label"
%token EOFF
%token END 0            "end of file" /* friendly name for $end */
%token ERROR

%%

input: DL NEQ integer NEWLINE rest trail eof {
  context->n=$3;
  if (context->n < 0) {
    IGRAPH_YY_ERRORF("Invalid vertex count in DL file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->n);
  }
  if (context->n > IGRAPH_DL_MAX_VERTEX_COUNT) {
    IGRAPH_YY_ERRORF("Vertex count too large in DL file (%" IGRAPH_PRId ").", IGRAPH_EINVAL, context->n);
  }
};

trail: | trail newline;

eof: | EOFF;

rest:    formfullmatrix { context->type=IGRAPH_DL_MATRIX; }
      |  edgelist1      { context->type=IGRAPH_DL_EDGELIST1; }
      |  nodelist1      { context->type=IGRAPH_DL_NODELIST1; }
;

formfullmatrix:  FORMATFULLMATRIX newline fullmatrix {} | fullmatrix {} ;

newline: | NEWLINE ;

fullmatrix:   DATA newline fullmatrixdata { }
            | LABELS newline labels newline DATA newline fullmatrixdata { }
            | LABELSEMBEDDED newline DATA newline labeledfullmatrixdata { }
;

labels:       {} /* nothing, empty matrix */
            | labels newline LABEL {
              IGRAPH_YY_CHECK(igraph_i_dl_add_str(igraph_dl_yyget_text(scanner),
                                                  igraph_dl_yyget_leng(scanner),
                                                  context)); }
;

fullmatrixdata: {} | fullmatrixdata zerooneseq NEWLINE {
  context->from += 1;
  context->to = 0;
 } ;

zerooneseq: | zerooneseq zeroone { } ;

zeroone: DIGIT {
  /* TODO: What if the digit is neither 0 or 1? Are multigraphs allowed? */
  char c = igraph_dl_yyget_text(scanner)[0];
  if (c == '1') {
    IGRAPH_YY_CHECK(igraph_vector_int_push_back(&context->edges,
                                         context->from));
    IGRAPH_YY_CHECK(igraph_vector_int_push_back(&context->edges,
                                         context->to));
  } else if (c != '0') {
      IGRAPH_YY_ERRORF("Unexpected digit '%c' in adjacency matrix in DL file.",
                    IGRAPH_EINVAL, c);
  }
  context->to += 1;
} ;

labeledfullmatrixdata: reallabeledfullmatrixdata {} ;

reallabeledfullmatrixdata: labelseq NEWLINE labeledmatrixlines {} ;

labelseq: | labelseq newline label ;

label: LABEL { IGRAPH_YY_CHECK(igraph_i_dl_add_str(igraph_dl_yyget_text(scanner),
                                                   igraph_dl_yyget_leng(scanner),
                                                   context)); };

labeledmatrixlines: labeledmatrixline {
                 context->from += 1;
                 context->to = 0;
               }
             | labeledmatrixlines labeledmatrixline {
                 context->from += 1;
                 context->to = 0;
               };

labeledmatrixline: LABEL zerooneseq NEWLINE { } ;

/*-----------------------------------------------------------*/

edgelist1: FORMATEDGELIST1 newline edgelist1rest {} ;

edgelist1rest:   DATA newline edgelist1data {}
             | LABELS newline labels newline DATA newline edgelist1data {}
             | LABELSEMBEDDED newline DATA newline labelededgelist1data {}
             | LABELS newline labels newline LABELSEMBEDDED newline DATA newline labelededgelist1data {}
             | LABELSEMBEDDED newline LABELS newline labels newline DATA newline labelededgelist1data {}
;

edgelist1data: {} /* nothing, empty graph */
             | edgelist1data edgelist1dataline {}
;

edgelist1dataline: integer integer weight NEWLINE {
                    igraph_integer_t from = $1, to = $2;
                    IGRAPH_YY_CHECK(igraph_i_dl_check_vid(from));
                    IGRAPH_YY_CHECK(igraph_i_dl_check_vid(to));
                    IGRAPH_YY_CHECK(igraph_i_dl_add_edge_w(from-1, to-1, $3, context)); }
                 | integer integer NEWLINE {
                    igraph_integer_t from = $1, to = $2;
                    IGRAPH_YY_CHECK(igraph_i_dl_check_vid(from));
                    IGRAPH_YY_CHECK(igraph_i_dl_check_vid(to));
                    IGRAPH_YY_CHECK(igraph_i_dl_add_edge(from-1, to-1, context));
} ;

integer: NUM {
    igraph_integer_t val;
    IGRAPH_YY_CHECK(igraph_i_parse_integer(igraph_dl_yyget_text(scanner),
                                           igraph_dl_yyget_leng(scanner),
                                           &val));
    $$=val;
};

labelededgelist1data: {} /* nothing, empty graph */
                    | labelededgelist1data labelededgelist1dataline {}
;

labelededgelist1dataline: elabel elabel weight NEWLINE {
                          IGRAPH_YY_CHECK(igraph_i_dl_add_edge_w($1, $2, $3, context)); }
                        | elabel elabel NEWLINE {
                          IGRAPH_YY_CHECK(igraph_i_dl_add_edge($1, $2, context));
 };

weight: NUM {
    igraph_real_t val;
    IGRAPH_YY_CHECK(igraph_i_parse_real(igraph_dl_yyget_text(scanner),
                                        igraph_dl_yyget_leng(scanner),
                                        &val));
    $$=val;
};

elabel: LABEL {
  igraph_integer_t trie_id;

  /* Copy label list to trie, if needed */
  if (igraph_strvector_size(&context->labels) != 0) {
    igraph_integer_t i, id, n=igraph_strvector_size(&context->labels);
    for (i=0; i<n; i++) {
      IGRAPH_YY_CHECK(igraph_trie_get(&context->trie, igraph_strvector_get(&context->labels, i), &id));
    }
    igraph_strvector_clear(&context->labels);
  }
  IGRAPH_YY_CHECK(igraph_trie_get_len(&context->trie, igraph_dl_yyget_text(scanner),
                                   igraph_dl_yyget_leng(scanner), &trie_id));
  IGRAPH_ASSERT(0 <= trie_id && trie_id < IGRAPH_DL_MAX_VERTEX_COUNT);
  $$ = trie_id;
 };

/*-----------------------------------------------------------*/

nodelist1: FORMATNODELIST1 newline nodelist1rest {} ;

nodelist1rest:   DATA nodelist1data {}
             | LABELS newline labels newline DATA newline nodelist1data {}
             | LABELSEMBEDDED newline DATA newline labelednodelist1data {}
             | LABELS newline labels newline LABELSEMBEDDED newline DATA newline labelednodelist1data {}
             | LABELSEMBEDDED newline LABELS newline labels newline DATA newline labelednodelist1data {}
;

nodelist1data: {} /* nothing, empty graph */
             | nodelist1data nodelist1dataline {}
;

nodelist1dataline: from tolist NEWLINE {} ;

from: NUM {
  IGRAPH_YY_CHECK(igraph_i_parse_integer(igraph_dl_yyget_text(scanner),
                  igraph_dl_yyget_leng(scanner),
                  &context->from));
  IGRAPH_YY_CHECK(igraph_i_dl_check_vid(context->from));
} ;

tolist: {} | tolist integer {
  igraph_integer_t to = $2;
  IGRAPH_YY_CHECK(igraph_i_dl_check_vid(to));
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(&context->edges,
                                              context->from-1));
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(&context->edges, to-1));
 } ;

labelednodelist1data: {} /* nothing, empty graph */
                    | labelednodelist1data labelednodelist1dataline {}
;

labelednodelist1dataline: fromelabel labeltolist NEWLINE { } ;

fromelabel: elabel {
  context->from=$1;
 };

labeltolist: | labeltolist elabel {
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(&context->edges,
                                              context->from));
  IGRAPH_YY_CHECK(igraph_vector_int_push_back(&context->edges, $2));
 } ;

%%

int igraph_dl_yyerror(YYLTYPE* locp,
                      igraph_i_dl_parsedata_t* context,
                      const char *s) {
  snprintf(context->errmsg, sizeof(context->errmsg)/sizeof(char)-1,
           "Parse error in DL file, line %i (%s)",
           locp->first_line, s);
  return 0;
}

static igraph_error_t igraph_i_dl_add_str(char *newstr, yy_size_t length,
                        igraph_i_dl_parsedata_t *context) {
  IGRAPH_CHECK(igraph_strvector_push_back_len(&context->labels, newstr, length));
  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_dl_add_edge(igraph_integer_t from, igraph_integer_t to,
                         igraph_i_dl_parsedata_t *context) {
  //IGRAPH_CHECK(igraph_i_dl_check_vid(from+1));
  //IGRAPH_CHECK(igraph_i_dl_check_vid(to+1));
  IGRAPH_CHECK(igraph_vector_int_push_back(&context->edges, from));
  IGRAPH_CHECK(igraph_vector_int_push_back(&context->edges, to));
  return IGRAPH_SUCCESS;
}

static igraph_error_t igraph_i_dl_add_edge_w(igraph_integer_t from, igraph_integer_t to,
                           igraph_real_t weight,
                           igraph_i_dl_parsedata_t *context) {
  igraph_integer_t n=igraph_vector_size(&context->weights);
  igraph_integer_t n2=igraph_vector_int_size(&context->edges)/2;
  if (n != n2) {
    IGRAPH_CHECK(igraph_vector_resize(&context->weights, n2));
    for (; n<n2; n++) {
      VECTOR(context->weights)[n]=IGRAPH_NAN;
    }
  }
  IGRAPH_CHECK(igraph_i_dl_add_edge(from, to, context));
  IGRAPH_CHECK(igraph_vector_push_back(&context->weights, weight));
  return IGRAPH_SUCCESS;
}

/* Raise an error if the vertex index is invalid in the DL file.
 * DL files use 1-based vertex indices. */
static igraph_error_t igraph_i_dl_check_vid(igraph_integer_t dl_vid) {
    if (dl_vid < 1) {
        IGRAPH_ERRORF("Invalid vertex index in DL file: %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, dl_vid);
    }
    if (dl_vid > IGRAPH_DL_MAX_VERTEX_COUNT) {
        IGRAPH_ERRORF("Vertex index too large in DL file: %" IGRAPH_PRId ".",
                      IGRAPH_EINVAL, dl_vid);
    }
    return IGRAPH_SUCCESS;
}
