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
#include "foreign-dl-parser.h"
#include <stdio.h>

extern FILE *igraph_dl_yyin;
extern int igraph_dl_eof;
extern long int igraph_dl_mylineno;
extern char *igraph_i_dl_errmsg;
int igraph_dl_yyerror(char *s);

%}

%output="y.tab.c"
%name-prefix="igraph_dl_yy"
%defines
%token-table

%union {
}

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

input: DL NEQ NUM NEWLINE rest eof {} ;

eof: | EOFF;

rest:    formfullmatrix {}
      |  edgelist1  {} 
      |  nodelist1  {}
; 

formfullmatrix:  FORMATFULLMATRIX newline fullmatrix {} | fullmatrix {} ;

newline: | NEWLINE ;

fullmatrix:   DATA newline fullmatrixdata {}
            | LABELS newline labels newline DATA newline fullmatrixdata {}
            | LABELSEMBEDDED newline DATA newline labeledfullmatrixdata {}
;

labels: 	    {}		/* nothing, empty matrix */
            | labels newline LABEL {}
;

fullmatrixdata: {} | fullmatrixdata zerooneseq NEWLINE {} ;

zerooneseq: | zerooneseq zeroone {} ;

zeroone: DIGIT { } ;

labeledfullmatrixdata: {} | reallabeledfullmatrixdata ;

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

