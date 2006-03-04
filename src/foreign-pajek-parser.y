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
extern int igraph_pajek_yylex();
extern int igraph_pajek_mylineno;
extern char *igraph_pajek_yytext;
extern int igraph_pajek_yyleng;
int igraph_pajek_yyerror(char *s);
#include "types.h"
#include "memory.h"
#include "error.h"
#include <math.h>
extern igraph_vector_t *igraph_pajek_vector;
extern igraph_stack_t *igraph_pajek_coords;
extern bool_t igraph_pajek_directed;
extern long int igraph_pajek_vcount;
extern int igraph_pajek_mode;
extern long int igraph_pajek_actfrom, igraph_pajek_actto;
real_t igraph_pajek_get_number(const char *str, long int len);
%}

%output="y.tab.c"
%name-prefix="igraph_pajek_yy"
%defines

%union {
  long int intnum;
  double   realnum;
  char *string;
}

%type <intnum>   longint;
%type <intnum>   arcfrom;
%type <intnum>   arcto;
%type <intnum>   edgefrom;
%type <intnum>   edgeto;
%type <realnum>  number;
%type <string>   word;

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

vertices: VERTICESLINE longint NEWLINE vertdefs {
  igraph_pajek_vcount=$2; 
};

vertdefs: /* empty */  | vertdefs vertexline;

vertexline: NEWLINE |
            vertex NEWLINE |
            vertex vertexid vertexcoords shape params NEWLINE;

vertex: longint { igraph_pajek_mode=1; } 

vertexid: word ;

vertexcoords: /* empty */ 
            | number number {
  if (igraph_pajek_coords) {
    igraph_stack_push(igraph_pajek_coords, $1);
    igraph_stack_push(igraph_pajek_coords, $2);
  }
}
            | number number number {
  if (igraph_pajek_coords) {
    igraph_stack_push(igraph_pajek_coords, $1);
    igraph_stack_push(igraph_pajek_coords, $2);
    igraph_stack_push(igraph_pajek_coords, $3);
  }
}

shape: /* empty */ | word;

params: /* empty */ | params param;

param:
       vpword vpwordpar  
     | VP_X_FACT number
     | VP_Y_FACT number
     | VP_IC number number number
     | VP_BC number number number
     | VP_LC number number number
     | VP_LR number
     | VP_LPHI number
     | VP_BW number
     | VP_FOS number
     | VP_PHI number
     | VP_R number
     | VP_Q number
     | VP_LA number
     | VP_SIZE number
;

vpword: vpword2 { igraph_pajek_mode=3; };

vpword2: VP_FONT | VP_URL | VP_IC | VP_BC | VP_LC;

vpwordpar: word { igraph_pajek_mode=1; }

edgeblock: /* empty */ | edgeblock arcs | edgeblock edges | edgeblock arcslist | edgeblock edgeslist | edgeblock adjmatrix;

arcs: ARCSLINE NEWLINE arcsdefs { igraph_pajek_directed=0; };

arcsdefs: /* empty */ | arcsdefs arcsline;

arcsline: NEWLINE | 
          arcfrom arcto weight edgeparams NEWLINE  { 
  igraph_vector_push_back(igraph_pajek_vector, $1-1);
  igraph_vector_push_back(igraph_pajek_vector, $2-1); }
;

arcfrom: longint { igraph_pajek_mode=2; };

arcto: longint;

edges: EDGESLINE NEWLINE edgesdefs { igraph_pajek_directed=0; }

edgesdefs: /* empty */ | edgesdefs edgesline;

edgesline: NEWLINE | 
           edgefrom edgeto weight edgeparams NEWLINE { 
  igraph_vector_push_back(igraph_pajek_vector, $1-1);
  igraph_vector_push_back(igraph_pajek_vector, $2-1); }
;

edgefrom: longint { igraph_pajek_mode=2; }

edgeto: longint;

weight: /* empty */ | number;

edgeparams: /* empty */ | edgeparams edgeparam;

edgeparam:
     epword epwordpar
   | EP_C number number number
   | EP_S number
   | EP_W number
   | EP_H1 number
   | EP_H2 number
   | EP_A1 number
   | EP_A2 number
   | EP_K1 number
   | EP_K2 number
   | EP_AP number
   | EP_LP number
   | EP_LR number
   | EP_LPHI number
   | EP_LA number
   | EP_SIZE number
   | EP_FOS number
;

epword: epword2 { igraph_pajek_mode=3; };

epword2: EP_A | EP_P | EP_L | EP_LC | EP_C;

epwordpar: word { igraph_pajek_mode=2; };

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

matrixline: MATRIXLINE { igraph_pajek_actfrom=0; igraph_pajek_actto=0; }

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

word: ALNUM { igraph_pajek_yytext+=igraph_pajek_yyleng; }
      | NUM { igraph_pajek_get_number(igraph_pajek_yytext, 
				      igraph_pajek_yyleng); }
      | QSTR { 
  igraph_pajek_yytext++; 
  while (*igraph_pajek_yytext != '"') {
    igraph_pajek_yytext++;
  }
  igraph_pajek_yytext++; }
      | PSTR { 
  igraph_pajek_yytext++; 
  while (*igraph_pajek_yytext != ')') {
    igraph_pajek_yytext++;
  }
  igraph_pajek_yytext++;
};

%%

int igraph_pajek_yyerror(char *s)
{
  fprintf(stderr, "hiba: %i\n", igraph_pajek_mylineno);
  IGRAPH_ERROR("Cannot read pajek file", IGRAPH_PARSEERROR);
}

real_t igraph_pajek_get_number(const char *str, long int length) {
  real_t num;
  char *tmp=Calloc(length+1, char);
  
  strncpy(tmp, str, length);
  tmp[length]='\0';
  sscanf(tmp, "%lf", &num);
  Free(tmp);
  return num;
} 
