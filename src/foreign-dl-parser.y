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
int igraph_i_dl_add_matcol();

extern char *igraph_dl_yytext;
extern int igraph_dl_yyleng;

extern igraph_real_t igraph_pajek_get_number(const char *str, long int len);
 
%}

%output="y.tab.c"
%name-prefix="igraph_dl_yy"
%defines

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

input: DL NEQ NUM NEWLINE rest eof { 
  igraph_i_dl_data.n=igraph_pajek_get_number(igraph_dl_yytext, 
					      igraph_dl_yyleng);  } ;

eof: | EOFF;

rest:    formfullmatrix { igraph_i_dl_data.type=IGRAPH_DL_MATRIX; }
      |  edgelist1      { igraph_i_dl_data.type=IGRAPH_DL_EDGELIST1; }
      |  nodelist1      { igraph_i_dl_data.type=IGRAPH_DL_NODELIST1; }
; 

formfullmatrix:  FORMATFULLMATRIX newline fullmatrix {} | fullmatrix {} ;

newline: | NEWLINE ;

fullmatrix:   DATA newline fullmatrixdata { }
            | LABELS newline labels newline DATA newline fullmatrixdata {
              }
            | LABELSEMBEDDED newline DATA newline labeledfullmatrixdata {
              }
;

labels: 	    {}		/* nothing, empty matrix */
            | labels newline LABEL { 
	      igraph_i_dl_add_str(igraph_dl_yytext, igraph_dl_yyleng); }
;

fullmatrixdata: {} | fullmatrixdata zerooneseq NEWLINE {
  igraph_i_dl_add_matcol();
 } ;

zerooneseq: | zerooneseq zeroone {
  igraph_vector_bool_push_back(&igraph_i_dl_data.zerooneseq,
			       igraph_i_dl_data.zeroone);
 } ;

zeroone: DIGIT { igraph_i_dl_data.zeroone=igraph_dl_yytext[0]-'0'; } ;

labeledfullmatrixdata: {} | reallabeledfullmatrixdata {} ;

reallabeledfullmatrixdata: labelseq newline labeledmatrixlines {} ;

labelseq: {} | labelseq newline LABEL {} ;

labeledmatrixlines: labeledmatrixline {} 
             | labeledmatrixlines labeledmatrixline {}
;

labeledmatrixline: LABEL zerooneseq NEWLINE {} ;

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

edgelist1dataline: NUM NUM NEWLINE {} ;

labelededgelist1data: 	{}	/* nothing, empty graph */
             | labelededgelist1data labelededgelist1dataline {}
;

labelededgelist1dataline: LABEL LABEL NEWLINE {} ;

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

from: NUM {} ;

tolist: {} | tolist NUM {} ;

labelednodelist1data: 		{}	/* nothing, empty graph */
           | labelednodelist1data labelednodelist1dataline {}
;

labelednodelist1dataline: LABEL labeltolist NEWLINE {} ;

labeltolist: | labeltolist LABEL {} ;

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

int igraph_i_dl_add_matcol() {
  igraph_vector_bool_t *vec=&igraph_i_dl_data.zerooneseq;
  igraph_matrix_bool_t *mat=&igraph_i_dl_data.matrix;
  long int nrow=igraph_matrix_bool_nrow(mat);
  long int ncol=igraph_matrix_bool_ncol(mat);
  long int len=igraph_vector_bool_size(vec);
  if (nrow==0) {
    /* Trick to allocate all memory at once */
    IGRAPH_CHECK(igraph_matrix_bool_resize(&igraph_i_dl_data.matrix, len, len));
    IGRAPH_CHECK(igraph_matrix_bool_resize(&igraph_i_dl_data.matrix, len, 0));
  }
  IGRAPH_CHECK(igraph_matrix_bool_add_cols(mat, 1));
  IGRAPH_CHECK(igraph_matrix_bool_set_col(mat, vec, ncol));
  igraph_vector_bool_clear(vec);
  return 0;
}
