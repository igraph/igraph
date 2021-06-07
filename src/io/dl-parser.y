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

#include "config.h"

#include "core/math.h"
#include "internal/hacks.h"
#include "io/dl-header.h"
#include "io/parsers/dl-parser.h"
#include "io/parsers/dl-lexer.h"

#include <stdio.h>

int igraph_dl_yyerror(YYLTYPE* locp, igraph_i_dl_parsedata_t* context,
                      const char *s);
int igraph_i_dl_add_str(char *newstr, int length,
                        igraph_i_dl_parsedata_t *context);
int igraph_i_dl_add_edge(long int from, long int to,
                         igraph_i_dl_parsedata_t *context);
int igraph_i_dl_add_edge_w(long int from, long int to,
                           igraph_real_t weight,
                           igraph_i_dl_parsedata_t *context);

extern igraph_real_t igraph_pajek_get_number(const char *str, long int len);

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
  long int integer;
  igraph_real_t real;
};

%type <integer> integer elabel;
%type <real> weight;

%token NUM
%token NEWLINE
%token DL
%token NEQ
%token DATA
%token LABELS
%token LABELSEMBEDDED
%token FORMATFULLMATRIX
%token FORMATEDGELIST1
%token FORMATNODELIST1
%token DIGIT
%token LABEL
%token EOFF
%token ERROR

%%

input: DL NEQ integer NEWLINE rest trail eof { context->n=$3; };

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
              igraph_i_dl_add_str(igraph_dl_yyget_text(scanner),
                                  igraph_dl_yyget_leng(scanner),
                                  context); }
;

fullmatrixdata: {} | fullmatrixdata zerooneseq NEWLINE {
  context->from += 1;
  context->to = 0;
 } ;

zerooneseq: | zerooneseq zeroone { } ;

zeroone: DIGIT {
  if (igraph_dl_yyget_text(scanner)[0]=='1') {
    IGRAPH_CHECK(igraph_vector_push_back(&context->edges,
                                         context->from));
    IGRAPH_CHECK(igraph_vector_push_back(&context->edges,
                                         context->to));
  }
  context->to += 1;
} ;

labeledfullmatrixdata: reallabeledfullmatrixdata {} ;

reallabeledfullmatrixdata: labelseq NEWLINE labeledmatrixlines {} ;

labelseq: | labelseq newline label ;

label: LABEL { igraph_i_dl_add_str(igraph_dl_yyget_text(scanner),
                                   igraph_dl_yyget_leng(scanner),
                                   context); };

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
                   igraph_i_dl_add_edge_w($1-1, $2-1, $3, context); }
                 | integer integer NEWLINE {
                   igraph_i_dl_add_edge($1-1, $2-1, context);
} ;

integer: NUM { $$=igraph_pajek_get_number(igraph_dl_yyget_text(scanner),
                                          igraph_dl_yyget_leng(scanner)); };

labelededgelist1data: {} /* nothing, empty graph */
                    | labelededgelist1data labelededgelist1dataline {}
;

labelededgelist1dataline: elabel elabel weight NEWLINE {
                          igraph_i_dl_add_edge_w($1, $2, $3, context); }
                        | elabel elabel NEWLINE {
                          igraph_i_dl_add_edge($1, $2, context);
 };

weight: NUM { $$=igraph_pajek_get_number(igraph_dl_yyget_text(scanner),
                                         igraph_dl_yyget_leng(scanner)); };

elabel: LABEL {
  /* Copy label list to trie, if needed */
  if (igraph_strvector_size(&context->labels) != 0) {
    long int i, id, n=igraph_strvector_size(&context->labels);
    for (i=0; i<n; i++) {
      igraph_trie_get(&context->trie,
                      STR(context->labels, i), &id);
    }
    igraph_strvector_clear(&context->labels);
  }
  igraph_trie_get2(&context->trie, igraph_dl_yyget_text(scanner),
                   igraph_dl_yyget_leng(scanner), &$$);
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

from: NUM { context->from=igraph_pajek_get_number(igraph_dl_yyget_text(scanner),
                                                          igraph_dl_yyget_leng(scanner)); } ;

tolist: {} | tolist integer {
  IGRAPH_CHECK(igraph_vector_push_back(&context->edges,
                                       context->from-1));
  IGRAPH_CHECK(igraph_vector_push_back(&context->edges, $2-1));
 } ;

labelednodelist1data: {} /* nothing, empty graph */
                    | labelednodelist1data labelednodelist1dataline {}
;

labelednodelist1dataline: fromelabel labeltolist NEWLINE { } ;

fromelabel: elabel {
  context->from=$1;
 };

labeltolist: | labeltolist elabel {
  IGRAPH_CHECK(igraph_vector_push_back(&context->edges,
                                       context->from));
  IGRAPH_CHECK(igraph_vector_push_back(&context->edges, $2));
 } ;

%%

int igraph_dl_yyerror(YYLTYPE* locp, igraph_i_dl_parsedata_t* context,
                      const char *s) {
  snprintf(context->errmsg,
           sizeof(context->errmsg)/sizeof(char)-1,
           "%s in line %i", s, locp->first_line);
  return 0;
}

int igraph_i_dl_add_str(char *newstr, int length,
                        igraph_i_dl_parsedata_t *context) {
  int tmp=newstr[length];
  newstr[length]='\0';
  IGRAPH_CHECK(igraph_strvector_add(&context->labels, newstr));
  newstr[length]=tmp;
  return 0;
}

int igraph_i_dl_add_edge(long int from, long int to,
                         igraph_i_dl_parsedata_t *context) {
  IGRAPH_CHECK(igraph_vector_push_back(&context->edges, from));
  IGRAPH_CHECK(igraph_vector_push_back(&context->edges, to));
  return 0;
}

int igraph_i_dl_add_edge_w(long int from, long int to,
                           igraph_real_t weight,
                           igraph_i_dl_parsedata_t *context) {
  long int n=igraph_vector_size(&context->weights);
  long int n2=igraph_vector_size(&context->edges)/2;
  if (n != n2) {
    igraph_vector_resize(&context->weights, n2);
    for (; n<n2; n++) {
      VECTOR(context->weights)[n]=IGRAPH_NAN;
    }
  }
  IGRAPH_CHECK(igraph_i_dl_add_edge(from, to, context));
  IGRAPH_CHECK(igraph_vector_push_back(&context->weights, weight));
  return 0;
}
