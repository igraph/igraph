/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

static igraph_error_handler_t *igraph_i_error_handler=0;

static char *igraph_i_error_strings[]={ "No error",            /* 0 */
					"Failed",              /* 1 */
					"Out of memory",       /* 2 */
					"Parse error",         /* 3 */
					"Invalid value",       /* 4 */
					"Already exists",      /* 5 */
					"Invalid edge vector", /* 6 */
					"Invalid vertex id",   /* 7 */
					"Non-square matrix",   /* 8 */
					"Invalid mode",        /* 9 */
					"File operation error" /* 10 */
};

const char* igraph_strerror(const int igraph_errno) {
  return igraph_i_error_strings[igraph_errno];
}

int igraph_error(const char *reason, const char *file, int line,
		 int igraph_errno) {

  if (igraph_i_error_handler) {
    igraph_i_error_handler(reason, file, line, igraph_errno);
    return igraph_errno;
  }  else {
    fprintf(stderr, "Error at %s:%i :%s, %s\n", file, line, reason,
	    igraph_strerror(igraph_errno));
    abort();
  }
}

void igraph_error_handler_abort (const char *reason, const char *file,
				 int line, int igraph_errno) {
  fprintf(stderr, "Error at %s:%i :%s, %s\n", file, line, reason,
	  igraph_strerror(igraph_errno));
  abort();
}

void igraph_error_handler_ignore (const char *reason, const char *file,
				  int line, int igraph_errno) {
  fprintf(stderr, "Error at %s:%i :%s, %s\n", file, line, reason,
	  igraph_strerror(igraph_errno));
}

igraph_error_handler_t *
igraph_set_error_handler (igraph_error_handler_t * new_handler)
{
  igraph_error_handler_t * previous_handler = igraph_i_error_handler;
  igraph_i_error_handler = new_handler;
  return previous_handler;
}

struct igraph_i_protectedPtr igraph_i_finally_stack[100];

/*
 * Adds another element to the free list
 */

void IGRAPH_FINALLY_REAL(void (*func)(void*), void* ptr) {
  int no=igraph_i_finally_stack[0].all;
  assert (no<100);
  igraph_i_finally_stack[no].ptr=ptr;
  igraph_i_finally_stack[no].func=func;
  igraph_i_finally_stack[0].all ++;
}

void IGRAPH_FINALLY_CLEAN(int minus) { 
  igraph_i_finally_stack[0].all -= minus;
  if (igraph_i_finally_stack[0].all < 0) {
    fprintf(stderr, "corrupt finally stack\n");
    igraph_i_finally_stack[0].all = 0;
  }
}

void IGRAPH_FINALLY_FREE() {
  int p;
  for (p=igraph_i_finally_stack[0].all-1; p>=0; p++) {
    igraph_i_finally_stack[p].func(igraph_i_finally_stack[p].ptr);
  }
  igraph_i_finally_stack[0].all=0;
}
