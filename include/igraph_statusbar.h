/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#ifndef IGRAPH_STATUSBAR
#define IGRAPH_STATUSBAR

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

typedef int igraph_status_handler_t(const char *message, void *data);
extern igraph_status_handler_t igraph_status_handler_stderr;

igraph_status_handler_t *
igraph_set_status_handler(igraph_status_handler_t new_handler);

int igraph_status(const char *message, void *data);

#define IGRAPH_STATUS(message, data) \
  do { \
    if (igraph_status((message), (data)) != IGRAPH_SUCCESS) { \
      IGRAPH_FINALLY_FREE(); \
      return IGRAPH_INTERRUPTED; \
    } \
  } while (0)

int igraph_statusf(const char *message, void *data, ...);

#define IGRAPH_STATUSF(args) \
  do { \
    if (igraph_statusf args != IGRAPH_SUCCESS) { \
      IGRAPH_FINALLY_FREE(); \
      return IGRAPH_INTERRUPTED; \
    } \
  } while (0)

__END_DECLS

#endif
