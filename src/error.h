/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef IGRAPH_ERROR_H
#define IGRAPH_ERROR_H

/* This file contains the igraph error handling.
 * Most bits are taken literally from the GSL library (with the GSL_
 * prefix renamed to IGRAPH_), as i couldn't find a better way to do
 * them. */

enum {
  IGRAPH_SUCCESS    = 0,
  IGRAPH_FAILURE    = 1,
  IGRAPH_ENOMEM     = 2,
  IGRAPH_PARSEERROR = 3,
  IGRAPH_EINVAL     = 4,
  IGRAPH_EXISTS     = 5,
  IGRAPH_EINVEVECTOR= 6,
  IGRAPH_EINVVID    = 7,
  IGRAPH_NONSQUARE  = 8,
  IGRAPH_EINVMODE   = 9,
  IGRAPH_EFILE      = 10
};

/**
 */

#define IGRAPH_ERROR(reason, igraph_errno) \
       do { \
       igraph_error (reason, __FILE__, __LINE__, igraph_errno) ; \
       return igraph_errno ; \
       } while (0)

/**
 */

int igraph_error(const char *reason, const char *file, int line,
		 int igraph_errno);

/**
 */

const char* igraph_strerror(const int igraph_errno);

/**
 */

typedef void igraph_error_handler_t (const char * reason, const char * file,
				     int line, int igraph_errno);

igraph_error_handler_t igraph_error_handler_abort;
igraph_error_handler_t igraph_error_handler_ignore;
igraph_error_handler_t igraph_error_handler_printignore;

/**
 */

igraph_error_handler_t *
igraph_set_error_handler(igraph_error_handler_t new_handler);

#define IGRAPH_ERROR_SELECT_2(a,b)       ((a) != IGRAPH_SUCCESS ? (a) : ((b) != IGRAPH_SUCCESS ? (b) : IGRAPH_SUCCESS))
#define IGRAPH_ERROR_SELECT_3(a,b,c)     ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_2(b,c))
#define IGRAPH_ERROR_SELECT_4(a,b,c,d)   ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_3(b,c,d))
#define IGRAPH_ERROR_SELECT_5(a,b,c,d,e) ((a) != IGRAPH_SUCCESS ? (a) : IGRAPH_ERROR_SELECT_4(b,c,d,e))

/* Now comes the more conveninent error handling macro arsenal.
 * Ideas taken from exception.{h,c} by Laurent Deniau see
 * http://cern.ch/Laurent.Deniau/html/oopc/exception.html for more 
 * information. We don't use the exception handling code though.  */

struct igraph_i_protectedPtr {
  int all;
  void *ptr;
  void (*func)(void*);
};

typedef void igraph_finally_func_t (void*);

void IGRAPH_FINALLY_REAL(void (*func)(void*), void* ptr);
void IGRAPH_FINALLY_CLEAN(int); 
void IGRAPH_FINALLY_FREE();

#define IGRAPH_FINALLY(func, ptr) \
  IGRAPH_FINALLY_REAL((igraph_finally_func_t*)(func), (ptr))

#define IGRAPH_FERROR(reason, errno) do { \
  IGRAPH_FINALLY_FREE(); \
  IGRAPH_ERROR(reason, errno); \
  } while (0)

#define IGRAPH_CHECK(a) do { \
                 int igraph_i_ret=(a); \
                 if (igraph_i_ret != 0) { \
                    IGRAPH_FERROR("", igraph_i_ret); \
                 } } while(0)

#endif
