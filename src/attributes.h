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

#ifndef REST_ATTRIBUTES_H
#define REST_ATTRIBUTES_H

#include "types.h"

typedef enum { IGRAPH_ATTRIBUTE_NUM=0, 
	       IGRAPH_ATTRIBUTE_STR } igraph_attribute_type_t;

typedef struct s_igraph_attribute_list {
  long int len;
  igraph_strvector_t names;
  igraph_vector_t types;
  vector_ptr_t data;
} igraph_attribute_list_t;

int igraph_attribute_list_init(igraph_attribute_list_t *al, long int len);
void igraph_attribute_list_destroy(igraph_attribute_list_t *al);
int igraph_attribute_list_add(igraph_attribute_list_t *al,
			      const char *name, igraph_attribute_type_t type);
int igraph_attribute_list_remove(igraph_attribute_list_t *al, 
				 const char *name);
int igraph_attribute_list_get(const igraph_attribute_list_t *al, 
			      const char *name, 
			      long int idx, void **value, 
			      igraph_attribute_type_t *type);
int igraph_attribute_list_set(igraph_attribute_list_t *al, const char *name,
			      long int idx, const void *value);
int igraph_attribute_list_get_many(const igraph_attribute_list_t *al, 
				   const char *name,
				   const igraph_vector_t *idx, void **value);
int igraph_attribute_list_set_many(igraph_attribute_list_t *al, 
				   const char *name,
				   const igraph_vector_t *idx, const void *value);
int igraph_attribute_list_get_all(const igraph_attribute_list_t *al, 
				  const char *name,
				  void **value, igraph_attribute_type_t *type);
long int igraph_attribute_list_size(const igraph_attribute_list_t *al);
int igraph_attribute_list_add_elem(igraph_attribute_list_t *al, long int ne);
int igraph_attribute_list_names(const igraph_attribute_list_t *al,
				igraph_strvector_t *names, igraph_vector_t *types);
int igraph_attribute_list_copy(igraph_attribute_list_t *to,
			       const igraph_attribute_list_t *from);
int igraph_attribute_list_get_type(const igraph_attribute_list_t *al, 
				   const char *name,
				   igraph_attribute_type_t *type);
bool_t igraph_attribute_list_has(const igraph_attribute_list_t *graph, 
				 const char *name);

void igraph_attribute_list_remove_elem_idx(igraph_attribute_list_t *al, 
					   long int *index, long int nremove);
void igraph_attribute_list_remove_elem_neg(igraph_attribute_list_t *al,
					   const igraph_vector_t *neg, 
					   long int nremove);

#endif
