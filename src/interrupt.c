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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

#include "igraph.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

igraph_interruption_handler_t *igraph_i_interruption_handler=0;

int igraph_allow_interruption(void* data) {
  if (igraph_i_interruption_handler) {
    return igraph_i_interruption_handler(data);
  }
  return IGRAPH_SUCCESS;
}

igraph_interruption_handler_t *
igraph_set_interruption_handler (igraph_interruption_handler_t * new_handler)
{
  igraph_interruption_handler_t * previous_handler = igraph_i_interruption_handler;
  igraph_i_interruption_handler = new_handler;
  return previous_handler;
}
