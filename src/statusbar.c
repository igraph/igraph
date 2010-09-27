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

#include "igraph_statusbar.h"
#include "igraph_error.h"
#include <stdio.h>
#include <stdarg.h>

static igraph_status_handler_t *igraph_i_status_handler=0;

int igraph_status(const char *message, void *data) {
  if (igraph_i_status_handler) {
    if (igraph_i_status_handler(message, data) != IGRAPH_SUCCESS) { 
      return IGRAPH_INTERRUPTED;
    }
  }
  return IGRAPH_SUCCESS;
}

int igraph_statusf(const char *message, void *data, ...) {
  char buffer[300];
  va_list ap;
  va_start(ap, data);
  vsnprintf(buffer, sizeof(buffer)-1, message, ap);
  return igraph_status(buffer, data);
}

int igraph_status_handler_stderr(const char *message, void *data) {
  fputs(message, stderr);
  return 0;
}

igraph_status_handler_t *
igraph_set_status_handler(igraph_status_handler_t new_handler) {
  igraph_status_handler_t *previous_handler=igraph_i_status_handler;
  igraph_i_status_handler = new_handler;
  return previous_handler;
}

