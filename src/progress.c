/* -*- mode: C -*-  */
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "config.h"

static igraph_progress_handler_t *igraph_i_progress_handler=0;

int igraph_progress(const char *message, igraph_real_t percent, void *data) {
  if (igraph_i_progress_handler) {
    if (igraph_i_progress_handler(message, percent, data) != IGRAPH_SUCCESS)
      return IGRAPH_INTERRUPTED;
  }
  return IGRAPH_SUCCESS;
}

#ifndef USING_R
int igraph_progress_handler_stderr(const char *message, igraph_real_t percent,
				  void* data) {
  fputs(message, stderr);
  fprintf(stderr, "%.1f percent ready\n", (double)percent);
  return 0;
}
#endif

igraph_progress_handler_t *
igraph_set_progress_handler(igraph_progress_handler_t new_handler) {
  igraph_progress_handler_t *previous_handler=igraph_i_progress_handler;
  igraph_i_progress_handler = new_handler;
  return previous_handler;
}
