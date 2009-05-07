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

#ifndef IGRAPH_INTERFACE_H
#define IGRAPH_INTERFACE_H

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
#include "igraph_datatype.h"
#include "igraph_iterators.h"

__BEGIN_DECLS

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

int igraph_empty(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed);
int igraph_empty_attrs(igraph_t *graph, igraph_integer_t n, igraph_bool_t directed, void *attr);
int igraph_destroy(igraph_t *graph);
int igraph_copy(igraph_t *to, const igraph_t *from);
int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges, 
		     void *attr);
int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, 
			void *attr);
int igraph_delete_edges(igraph_t *graph, igraph_es_t edges);
int igraph_delete_vertices(igraph_t *graph, const igraph_vs_t vertices);
igraph_integer_t igraph_vcount(const igraph_t *graph);
igraph_integer_t igraph_ecount(const igraph_t *graph);
int igraph_neighbors(const igraph_t *graph, igraph_vector_t *neis, igraph_integer_t vid, 
		     igraph_neimode_t mode); 
igraph_bool_t igraph_is_directed(const igraph_t *graph);
int igraph_degree(const igraph_t *graph, igraph_vector_t *res, 
		  const igraph_vs_t vids, igraph_neimode_t mode, 
		  igraph_bool_t loops);
int igraph_edge(const igraph_t *graph, igraph_integer_t eid, 
		igraph_integer_t *from, igraph_integer_t *to);		
int igraph_edges(const igraph_t *graph, igraph_es_t eids,
		 igraph_vector_t *edges);
int igraph_get_eid(const igraph_t *graph, igraph_integer_t *eid,
		   igraph_integer_t from, igraph_integer_t to,
		   igraph_bool_t directed);
int igraph_get_eid2(const igraph_t *graph, igraph_integer_t *eid,
		   igraph_integer_t pfrom, igraph_integer_t pto,
		    igraph_bool_t directed);
int igraph_get_eids(const igraph_t *graph, igraph_vector_t *eids,
		    const igraph_vector_t *pairs, igraph_bool_t directed);
int igraph_adjacent(const igraph_t *graph, igraph_vector_t *eids, igraph_integer_t vid,
		    igraph_neimode_t mode);

#define IGRAPH_FROM(g,e) (VECTOR((g)->from)[(long int)(e)])
#define IGRAPH_TO(g,e)   (VECTOR((g)->to)  [(long int)(e)])
#define IGRAPH_OTHER(g,e,v) (IGRAPH_TO(g,(e))==(v) ? IGRAPH_FROM((g),(e)) : IGRAPH_TO((g),(e)))

__END_DECLS

#endif
