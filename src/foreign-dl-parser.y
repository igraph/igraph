/* 
   IGraph library.
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Copyright (C) 2009  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include "config.h"
#include "igraph.h"
#include "memory.h"
#include "foreign-dl-header.h"
#include "foreign-dl-parser.h"
#include <stdio.h>

extern int igraph_dl_eof;
extern long int igraph_dl_mylineno;
extern char *igraph_i_dl_errmsg;
int igraph_dl_yyerror(char *s);
extern igraph_i_dl_parsedata_t igraph_i_dl_data;

int igraph_i_dl_add_str(char *newstr, int length);
int igraph_i_dl_add_edge(long int from, long int to);
int igraph_i_dl_add_edge_w(long int from, long int to, igraph_real_t weight);

extern char *igraph_dl_yytext;
extern int igraph_dl_yyleng;

extern igraph_real_t igraph_pajek_get_number(const char *str, long int len);
 
%}

%output="y.tab.c"
%name-prefix="igraph_dl_yy"
%defines

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

%%

input: DL NEQ integer NEWLINE rest trail eof { igraph_i_dl_data.n=$3; };

trail: | trail newline;

eof: | EOFF;

rest:    formfullmatrix { igraph_i_dl_data.type=IGRAPH_DL_MATRIX; }
      |  edgelist1      { igraph_i_dl_data.type=IGRAPH_DL_EDGELIST1; }
      |  nodelist1      { igraph_i_dl_data.type=IGRAPH_DL_NODELIST1; }
; 

formfullmatrix:  FORMATFULLMATRIX newline fullmatrix {} | fullmatrix {} ;

newline: | NEWLINE ;

fullmatrix:   DATA newline fullmatrixdata { }
            | LABELS newline labels newline DATA newline fullmatrixdata { }
            | LABELSEMBEDDED newline DATA newline labeledfullmatrixdata { }
;

labels: 	    {}		/* nothing, empty matrix */
            | labels newline LABEL { 
	      igraph_i_dl_add_str(igraph_dl_yytext, igraph_dl_yyleng); }
;

fullmatrixdata: {} | fullmatrixdata zerooneseq NEWLINE {
  igraph_i_dl_data.from += 1;
  igraph_i_dl_data.to = 0;
 } ;

zerooneseq: | zerooneseq zeroone { } ;

zeroone: DIGIT {
  if (igraph_dl_yytext[0]=='1') {
    IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, 
					 igraph_i_dl_data.from));
    IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, 
					 igraph_i_dl_data.to));
  }
  igraph_i_dl_data.to += 1;
} ;

labeledfullmatrixdata: reallabeledfullmatrixdata {} ;

reallabeledfullmatrixdata: labelseq NEWLINE labeledmatrixlines {} ;

labelseq: | labelseq newline label ;

label: LABEL { igraph_i_dl_add_str(igraph_dl_yytext, igraph_dl_yyleng); };

labeledmatrixlines: labeledmatrixline {
	         igraph_i_dl_data.from += 1; 
		 igraph_i_dl_data.to = 0;
               } 
             | labeledmatrixlines labeledmatrixline { 
	         igraph_i_dl_data.from += 1; 
		 igraph_i_dl_data.to = 0;
               };

labeledmatrixline: LABEL zerooneseq NEWLINE { } ;

/*-----------------------------------------------------------*/

edgelist1: FORMATEDGELIST1 newline edgelist1rest {} ;

edgelist1rest:   DATA edgelist1data {}
             | LABELS newline labels newline DATA newline edgelist1data {}
             | LABELSEMBEDDED newline DATA newline labelededgelist1data {}
             | LABELS newline labels newline LABELSEMBEDDED newline DATA newline labelededgelist1data {}
             | LABELSEMBEDDED newline LABELS newline labels newline DATA newline labelededgelist1data {}
;

edgelist1data: 		{}	/* nothing, empty graph */
             | edgelist1data edgelist1dataline {}
;

edgelist1dataline: integer integer weight NEWLINE {
                     igraph_i_dl_add_edge_w($1-1, $2-1, $3); }
                 | integer integer NEWLINE {
		     igraph_i_dl_add_edge($1-1, $2-1);
} ;

integer: NUM { $$=igraph_pajek_get_number(igraph_dl_yytext, 
					  igraph_dl_yyleng); };

labelededgelist1data: 	{}	/* nothing, empty graph */
             | labelededgelist1data labelededgelist1dataline {}
;

labelededgelist1dataline: elabel elabel weight NEWLINE {
                            igraph_i_dl_add_edge_w($1, $2, $3); }
                        | elabel elabel NEWLINE {
			  igraph_i_dl_add_edge($1, $2);
 };

weight: NUM { $$=igraph_pajek_get_number(igraph_dl_yytext, 
					  igraph_dl_yyleng); };

elabel: LABEL {
  /* Copy label list to trie, if needed */
  if (igraph_strvector_size(&igraph_i_dl_data.labels) != 0) {
    long int i, id, n=igraph_strvector_size(&igraph_i_dl_data.labels);
    for (i=0; i<n; i++) {
      igraph_trie_get(&igraph_i_dl_data.trie,
		      STR(igraph_i_dl_data.labels, i), &id);
    }
    igraph_strvector_clear(&igraph_i_dl_data.labels);
  }
  igraph_trie_get2(&igraph_i_dl_data.trie, igraph_dl_yytext, 
		   igraph_dl_yyleng, &$$);
 };

/*-----------------------------------------------------------*/

nodelist1: FORMATNODELIST1 newline nodelist1rest {} ;

nodelist1rest:   DATA nodelist1data {}
             | LABELS newline labels newline DATA newline nodelist1data {}
             | LABELSEMBEDDED newline DATA newline labelednodelist1data {}
             | LABELS newline labels newline LABELSEMBEDDED newline DATA newline labelednodelist1data {}
             | LABELSEMBEDDED newline LABELS newline labels newline DATA newline labelednodelist1data {}
;

nodelist1data: 		{}	/* nothing, empty graph */
             | nodelist1data nodelist1dataline {}
;

nodelist1dataline: from tolist NEWLINE {} ;

from: NUM { igraph_i_dl_data.from=igraph_pajek_get_number(igraph_dl_yytext,
							  igraph_dl_yyleng); } ;

tolist: {} | tolist integer { 
  IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, 
				       igraph_i_dl_data.from-1)); 
  IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, $2-1));
 } ;

labelednodelist1data: 		{}	/* nothing, empty graph */
           | labelednodelist1data labelednodelist1dataline {}
;

labelednodelist1dataline: fromelabel labeltolist NEWLINE { } ;

fromelabel: elabel {
  igraph_i_dl_data.from=$1;
 };

labeltolist: | labeltolist elabel {
  IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, 
				       igraph_i_dl_data.from));
  IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, $2));
 } ;

%%

int igraph_dl_yyerror(char *s) {
  static char str[300];

  snprintf(str, sizeof(str)-1, "Parse error in DL file, line %li (%s)", 
	   (long)igraph_dl_mylineno, s);
  igraph_i_dl_errmsg=str;
  return 0;
}

int igraph_i_dl_add_str(char *newstr, int length) {
  int tmp=newstr[length];
  newstr[length]='\0';
  IGRAPH_CHECK(igraph_strvector_add(&igraph_i_dl_data.labels, newstr));
  newstr[length]=tmp;
  return 0;
}

int igraph_i_dl_add_edge(long int from, long int to) {
  IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, from));
  IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.edges, to));
  return 0;
}

int igraph_i_dl_add_edge_w(long int from, long int to, igraph_real_t weight) {
  long int n=igraph_vector_size(&igraph_i_dl_data.weights);
  long int n2=igraph_vector_size(&igraph_i_dl_data.edges)/2;
  if (n != n2) {
    igraph_vector_resize(&igraph_i_dl_data.weights, n2);
    for (; n<n2; n++) {
      VECTOR(igraph_i_dl_data.weights)[n]=IGRAPH_NAN;
    }
  }
  IGRAPH_CHECK(igraph_i_dl_add_edge(from, to));
  IGRAPH_CHECK(igraph_vector_push_back(&igraph_i_dl_data.weights, weight));
  return 0;
}
