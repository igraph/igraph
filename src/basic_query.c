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

bool_t igraph_are_connected(igraph_t *graph, integer_t v1, integer_t v2) {

  igraph_iterator_t it;
  bool_t res=0;

  igraph_iterator_vneis(graph, &it, v1, 1);
  
  while (!res && !igraph_end(graph, &it)) {
    res= (igraph_get_vertex(graph, &it) == v2);
    igraph_next(graph, &it);
  }
  
  igraph_iterator_destroy(graph, &it);  
  
  return res;
}
