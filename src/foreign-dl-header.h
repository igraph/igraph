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

#include "igraph_types.h"

typedef enum { IGRAPH_DL_MATRIX, 
	       IGRAPH_DL_EDGELIST1, IGRAPH_DL_NODELIST1 } igraph_i_dl_type_t;

typedef struct {
  long int n;
  long int from, to;
  igraph_vector_t edges;
  igraph_vector_t weights;
  igraph_strvector_t labels;
  igraph_trie_t trie;
  igraph_i_dl_type_t type;
} igraph_i_dl_parsedata_t;
