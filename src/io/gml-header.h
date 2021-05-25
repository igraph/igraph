/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi@rmki.kfki.hu>
   334 Harvard street, Cambridge MA, 02139 USA

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

#include "igraph_error.h"

#include "io/gml-tree.h"

typedef struct {
    void *scanner;
    int eof;
    int depth;
    char errmsg[300];
    igraph_gml_tree_t *tree;
} igraph_i_gml_parsedata_t;

/**
 * Initializes a GML parser context.
 */
int igraph_i_gml_parsedata_init(igraph_i_gml_parsedata_t* context);

/**
 * Destroys a GML parser context, freeing all memory currently used by the
 * context.
 */
void igraph_i_gml_parsedata_destroy(igraph_i_gml_parsedata_t* context);
