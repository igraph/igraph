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

#ifndef IGRAPH_PATHS_H
#define IGRAPH_PATHS_H

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

int igraph_diameter(const igraph_t *graph, igraph_integer_t *res, 
		    igraph_integer_t *from, igraph_integer_t *to,
		    igraph_vector_t *path,
		    igraph_bool_t directed, igraph_bool_t unconn);
int igraph_diameter_dijkstra(const igraph_t *graph,
			     const igraph_vector_t *weights,
			     igraph_real_t *pres,
			     igraph_integer_t *pfrom,
			     igraph_integer_t *pto,
			     igraph_vector_t *path,
			     igraph_bool_t directed,
			     igraph_bool_t unconn);

int igraph_shortest_paths(const igraph_t *graph, igraph_matrix_t *res, 
			  const igraph_vs_t from, const igraph_vs_t to, 
			  igraph_neimode_t mode);
int igraph_get_shortest_paths(const igraph_t *graph, igraph_vector_ptr_t *res,
			      igraph_integer_t from, const igraph_vs_t to, 
			      igraph_neimode_t mode);
int igraph_get_all_shortest_paths(const igraph_t *graph,
				  igraph_vector_ptr_t *res, 
				  igraph_vector_t *nrgeo,
				  igraph_integer_t from, const igraph_vs_t to, 
				  igraph_neimode_t mode);
int igraph_shortest_paths_dijkstra(const igraph_t *graph,
				   igraph_matrix_t *res,
				   const igraph_vs_t from,
				   const igraph_vs_t to,
				   const igraph_vector_t *weights, 
				   igraph_neimode_t mode);
int igraph_shortest_paths_bellman_ford(const igraph_t *graph,
				   igraph_matrix_t *res,
				   const igraph_vs_t from,
				   const igraph_vs_t to,
				   const igraph_vector_t *weights, 
				   igraph_neimode_t mode);
int igraph_get_shortest_paths_dijkstra(const igraph_t *graph,
                                       igraph_vector_ptr_t *res,
				       igraph_integer_t from,
				       igraph_vs_t to,
				       const igraph_vector_t *weights,
				       igraph_neimode_t mode); 
int igraph_shortest_paths_johnson(const igraph_t *graph,
				  igraph_matrix_t *res,
				  const igraph_vs_t from,
				  const igraph_vs_t to,
				  const igraph_vector_t *weights);

int igraph_average_path_length(const igraph_t *graph, igraph_real_t *res,
			       igraph_bool_t directed, igraph_bool_t unconn);
int igraph_path_length_hist(const igraph_t *graph, igraph_vector_t *res,
			    igraph_real_t *unconnected, igraph_bool_t directed);


__END_DECLS

#endif
