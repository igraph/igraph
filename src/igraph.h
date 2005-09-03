/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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

#ifndef RESTGAME_H
#define RESTGAME_H

#include "types.h"

typedef struct igraph_s {
  integer_t n;			/* Number of vertices           */
  bool_t directed;		/* Is it directed?              */
  vector_t from;		/* The edge list, first column  */
  vector_t to;                  /* The edge list, second column */
  vector_t oi;			/* Out index                    */
  vector_t ii;			/* In index                     */
  vector_t os;			/* Out index start              */
  vector_t is;			/* In index start               */
} igraph_t;

/* -------------------------------------------------- */
/* Interface                                          */
/* -------------------------------------------------- */

int igraph_empty(igraph_t *graph, integer_t n, bool_t directed);
int igraph_destroy(igraph_t *graph);
int igraph_add_edges(igraph_t *graph, vector_t *edges);
int igraph_add_vertices(igraph_t *graph, integer_t nv);
int igraph_delete_edges(igraph_t *graph, vector_t *edges);
int igraph_delete_vertices(igraph_t *graph, vector_t *vertices);
integer_t igraph_vcount(igraph_t *graph);
integer_t igraph_ecount(igraph_t *graph);
int igraph_neighbors(igraph_t *graph, vector_t *neis, integer_t vid, 
		   integer_t mode); 
bool_t igraph_is_directed(igraph_t *graph);
int igraph_degree(igraph_t *graph, vector_t *res, vector_t *vids, 
		  integer_t mode, bool_t loops);

/* -------------------------------------------------- */
/* Iterators                                          */
/* -------------------------------------------------- */

struct igraph_iterator_t;

typedef struct igraph_iterator_t {
  integer_t type;
  void *data;
  int (*next)(igraph_t *graph, struct igraph_iterator_t *it);
  int (*prev)(igraph_t *graph, struct igraph_iterator_t *it);
  bool_t (*end)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertex)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertexfrom)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertexto)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getvertexnei)(igraph_t *graph, struct igraph_iterator_t *it);
  integer_t (*getedge)(igraph_t *graph, struct igraph_iterator_t *it);
} igraph_iterator_t;

/* The constructors & destructor */
int igraph_iterator_vid(igraph_t *graph, igraph_iterator_t *it);
int igraph_iterator_eid(igraph_t *graph, igraph_iterator_t *it);
int igraph_iterator_efromorder(igraph_t *graph, igraph_iterator_t *it);
int igraph_iterator_eneis(igraph_t *graph, igraph_iterator_t *it, 
			  integer_t vid, integer_t mode);

int igraph_iterator_destroy(igraph_t *graph, igraph_iterator_t *it);

/* Generics, common functions */
int igraph_next(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_nei(igraph_t *graph, igraph_iterator_t *it);

/* Generics, vertices */
integer_t igraph_get_vertex(igraph_t *graph, igraph_iterator_t *it);

/* Generics, edges */
integer_t igraph_get_vertex_from(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_to(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_edge(igraph_t *graph, igraph_iterator_t *it);

/* Specifics, simple vertex iterator */
int igraph_next_vid(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_vid(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_vid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_vid(igraph_t *graph, igraph_iterator_t *it);

/* Specifics, simple edge iterator */
int igraph_next_eid(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_eid(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_eid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_from_eid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_to_eid(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_edge_eid(igraph_t *graph, igraph_iterator_t *it); 

/* Specifics, edge iterator according to the 'from' order */
int igraph_next_efromorder(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_efromorder(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_efromorder(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_from_efromorder(igraph_t *graph, 
					    igraph_iterator_t *it);
integer_t igraph_get_vertex_to_efromorder(igraph_t *graph, 
					  igraph_iterator_t *it);
integer_t igraph_get_edge_efromorder(igraph_t *graph, igraph_iterator_t *it);


/* Iterates over the edges to and/or from a vertex */
int igraph_next_eneis(igraph_t *graph, igraph_iterator_t *it);
int igraph_prev_eneis(igraph_t *graph, igraph_iterator_t *it);
bool_t igraph_end_eneis(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_from_eneis(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_vertex_to_eneis(igraph_t *graph, igraph_iterator_t *it);
integer_t igraph_get_edge_eneis(igraph_t *graph, igraph_iterator_t *it); 
int igraph_iterator_eneis_set(igraph_t *graph, igraph_iterator_t *it, 
			      integer_t vid, integer_t mode);
integer_t igraph_get_vertex_nei_eneis(igraph_t *graph, igraph_iterator_t *it);

/* TODO: attributes */

/* -------------------------------------------------- */
/* Error handling                                     */
/* -------------------------------------------------- */

int igraph_error(const char *msg);

/* -------------------------------------------------- */
/* Constructors, deterministic                        */
/* -------------------------------------------------- */

int igraph_create(igraph_t *graph, vector_t *edges, integer_t n, 
		  bool_t directed);
int igraph_adjacency(igraph_t *graph, vector_t *adjmatrix,
		     integer_t dirmode);
int igraph_star(igraph_t *graph, integer_t n, integer_t mode, 
		integer_t center, bool_t directed);
int igraph_lattice(igraph_t *graph, vector_t *dimvector, integer_t nei, 
		   bool_t directed, bool_t mutual, bool_t circular);
int igraph_ring(igraph_t *graph, integer_t n, bool_t directed, bool_t mutual);
int igraph_tree(igraph_t *graph, integer_t n, integer_t children);

/* -------------------------------------------------- */
/* Constructors, games (=stochastic)                  */
/* -------------------------------------------------- */

int igraph_barabasi_game(igraph_t *graph, integer_t n, integer_t m, 
			 vector_t *outseq, bool_t outpref, bool_t directed);
int igraph_erdos_renyi_game(igraph_t *graph, integer_t n, real_t p,
			    bool_t directed, bool_t loops);
int igraph_degree_sequence_game(igraph_t *graph, vector_t *out_deg,
				vector_t *in_deg, integer_t method);
int igraph_growing_random_game(igraph_t *graph, integer_t n, 
			       integer_t m, bool_t directed, bool_t citation);
int igraph_aging_prefatt_game(igraph_t *graph, integer_t n, integer_t m,
			      integer_t aging_type, real_t aging_exp);

/* -------------------------------------------------- */
/* Basic query functions                              */
/* -------------------------------------------------- */

bool_t igraph_are_connected(igraph_t *graph, integer_t v1, integer_t v2);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

int igraph_diameter(igraph_t *graph, integer_t *res, 
		    bool_t directed, bool_t unconn);
int igraph_minimum_spanning_tree_unweighted(igraph_t *graph, igraph_t *mst);
int igraph_minimum_spanning_tree_prim(igraph_t *graph, igraph_t *mst,
				      vector_t *weights);
int igraph_closeness(igraph_t *graph, vector_t *res, vector_t *vids, 
		     integer_t mode);
int igraph_shortest_paths(igraph_t *graph, matrix_t *res, 
			  vector_t *from, integer_t mode);
int igraph_get_shortest_paths(igraph_t *graph, vector_t *res,
			      integer_t from, integer_t mode);
int igraph_subcomponent(igraph_t *graph, vector_t *res, real_t vid, 
			integer_t mode);
int igraph_subgraph(igraph_t *graph, igraph_t *res, integer_t  mode);
int igraph_simplify(igraph_t *graph, bool_t remove_loops, 
		    bool_t remove_multiple); 
int igraph_betweenness (igraph_t *graph, vector_t *res, vector_t *vids, 
			bool_t directed);
int igraph_edge_betweenness (igraph_t *graph, vector_t *result, 
			     bool_t directed);
/* TODO: degree.distribution (?) */

/* -------------------------------------------------- */
/* Components                                         */
/* -------------------------------------------------- */

int igraph_clusters(igraph_t *graph, vector_t *membership, vector_t *csize, 
		    integer_t mode);
int igraph_is_connected(igraph_t *graph, bool_t *res, integer_t mode);
/* TODO: cluster.distribution (?) */

/* -------------------------------------------------- */
/* Layouts                                            */
/* -------------------------------------------------- */

int igraph_layout_random(igraph_t *graph, vector_t *res);
int igraph_layout_circle(igraph_t *graph, vector_t *res);
int igraph_layout_fruchterman_reingold(igraph_t *graph, vector_t *res, 
				       integer_t niter, real_t coolexp,
				       integer_t frame, vector_t *initial,
				       real_t initemp);
int igraph_layout_kamada_kawai(igraph_t *graph, matrix_t *res,
			       integer_t niter, real_t sigma, 
			       real_t initemp, real_t coolexp,
			       real_t kkconst);
int igraph_layout_springs(igraph_t *graph, vector_t *res,
			  real_t mass, real_t equil, real_t k,
			  real_t repeqdis, real_t kfr, bool_t repulse);

/* -------------------------------------------------- */
/* Centrality                                         */
/* -------------------------------------------------- */

/* TODO: evcent */

/* -------------------------------------------------- */
/* Cocitation                                         */
/* -------------------------------------------------- */

int igraph_cocitation(igraph_t *graph, matrix_t *res, vector_t *vids);
int igraph_bibcoupling(igraph_t *graph, matrix_t *res, vector_t *vids);

/* -------------------------------------------------- */
/* Community Structure                                */
/* -------------------------------------------------- */

/* TODO: eb.community */
/* TODO: cut.community */
/* TODO: edge.type.matrix */
/* TODO: modularity */
/* TODO:  */

/* -------------------------------------------------- */
/* Conversion                                         */
/* -------------------------------------------------- */

int igraph_get_adjacency(igraph_t *graph, vector_t *res);
int igraph_get_edgelist(igraph_t *graph, vector_t *res);

/* -------------------------------------------------- */
/* Read and write foreign formats                     */
/* -------------------------------------------------- */

int igraph_read_graph(const char *filename, integer_t format);
int igraph_write_graph(igraph_t *graph, const char *filename, 
		       integer_t format);

/* -------------------------------------------------- */
/* Dynamics measurement                               */
/* -------------------------------------------------- */

int igraph_measure_dynamics_idage(igraph_t *graph, matrix_t *akl, 
				  matrix_t *sd,
				  vector_t *st, integer_t agebins,
				  integer_t maxind, bool_t lsd);
int igraph_measure_dynamics_idage_st(igraph_t *graph, vector_t *res,
				     matrix_t *akl);

/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

int igraph_running_mean(vector_t *data, vector_t *res, integer_t binwidth);

#endif
