/* -*- mode: C -*-  */
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

#ifndef IGRAPH_PROGRESS_H
#define IGRAPH_PROGRESS_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

#include "igraph_types.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Progress handlers                                  */
/* -------------------------------------------------- */

typedef int igraph_progress_handler_t(const char *message, igraph_real_t percent,
				      void *data);

extern igraph_progress_handler_t igraph_progress_handler_stderr;

igraph_progress_handler_t *
igraph_set_progress_handler(igraph_progress_handler_t new_handler);

int igraph_progress(const char *message, igraph_real_t percent, void *data);

int igraph_progressf(const char *message, igraph_real_t percent, void *data, 
		     ...);

#define IGRAPH_PROGRESS(message, percent, data) \
  do { \
    if (igraph_progress((message), (percent), (data)) != IGRAPH_SUCCESS) { \
      IGRAPH_FINALLY_FREE(); \
      return IGRAPH_INTERRUPTED; \
    } \
  } while (0)

__END_DECLS

#endif
