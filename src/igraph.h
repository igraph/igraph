/* -*- mode: C -*-  */
/* 
   IGraph R package.
   Copyright (C) 2003, 2004  Gabor Csardi <csardi@rmki.kfki.hu>
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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

/* Atmenetileg amig megoldodik a memoriafelszabaditas */
#define R_CheckUserInterrupt(); ;

#define I(x) (INTEGER(x)[0])
#define R(x) (REAL(x)[0])

#define RMATRIX(x, i, j) (REAL(x) [i-1 + \
        INTEGER(GET_DIM(x))[0] * (j-1)])

#define IMATRIX(x, i, j) (INTEGER(x) [i-1 + \
        INTEGER(GET_DIM(x))[0] * (j-1)])

/* dangerous! */
#define RNG_INTEGER(low, high) ((long int)(unif_rand()*(high-low+1)+low))
#define RNG_REAL(low, high)    (unif_rand()*(high-low)+low)

/* -------------------------------------------------- */
/* BASIC GRAPH OPERATIONS                             */
/* -------------------------------------------------- */

SEXP REST_neighborlist_delete_vertices(SEXP neilist, SEXP vids);
SEXP REST_add_edges_adjacencylist(SEXP data, SEXP edges, 
				  SEXP pdirected);

/* -------------------------------------------------- */
/* Games                                              */
/* -------------------------------------------------- */

SEXP REST_ba_game(SEXP pn, SEXP pm, SEXP outseq, SEXP poutpref);

/* -------------------------------------------------- */
/* Structural properties                              */
/* -------------------------------------------------- */

SEXP REST_diameter(SEXP interface, SEXP graph, SEXP pdirected, SEXP punconn);
SEXP REST_closeness(SEXP interface, SEXP graph, SEXP nodes, SEXP pmode);
SEXP REST_clusters(SEXP interface, SEXP graph);
SEXP REST_betweenness (SEXP interface, SEXP graph, SEXP pdirected);
SEXP REST_edge_betweenness (SEXP interface, SEXP graph, SEXP pdirected);
SEXP REST_shortest_paths(SEXP interface, SEXP graph, SEXP from, SEXP pmode);
SEXP REST_get_shortest_paths(SEXP interface, SEXP graph, SEXP from, SEXP pmode);
SEXP REST_cocitation(SEXP interface, SEXP graph, SEXP mode);
SEXP REST_minimum_spanning_tree_unweighted(SEXP interface, SEXP graph);
SEXP REST_minimum_spanning_tree_prim(SEXP interface, SEXP graph);

/* -------------------------------------------------- */
/* Community Structure                                */
/* -------------------------------------------------- */

SEXP REST_eb_community(SEXP interface, SEXP graph, SEXP pdirected);

/* -------------------------------------------------- */
/* Read and write foreign formats                     */
/* -------------------------------------------------- */

SEXP REST_import_pajek(SEXP interface, SEXP lines, SEXP other,
		       SEXP pattributes);

/* -------------------------------------------------- */
/* Layouts, mostly from SNA                           */
/* -------------------------------------------------- */

SEXP REST_layout_kamadakawai(SEXP pn, SEXP pniter, 
			     SEXP pelen, SEXP pinitemp, SEXP pcoolexp, 
			     SEXP pkkconst, SEXP psigma, SEXP px, SEXP py);

/* -------------------------------------------------- */
/* Dynamics measurement                               */
/* -------------------------------------------------- */

SEXP REST_measure_dynamics_idage(SEXP interface, SEXP graph, SEXP st,
				 SEXP pagebins, SEXP pmaxind, SEXP psd);
SEXP REST_measure_dynamics_idage_st(SEXP interface, SEXP graph, SEXP akl,
				    SEXP pagebins, SEXP pmaxind);

/* -------------------------------------------------- */
/* Other, not graph related                           */
/* -------------------------------------------------- */

SEXP REST_running_mean(SEXP data, SEXP pbinwidth);

/* -------------------------------------------------- */
/* The C igraph interface                             */
/* -------------------------------------------------- */

typedef SEXP (*VCOUNT_t)(SEXP, SEXP);
typedef SEXP (*ECOUNT_t)(SEXP, SEXP);
typedef SEXP (*NEIGHBORS_t)(SEXP, SEXP, long int, SEXP);
typedef SEXP (*ADD_VERTICES_t)(SEXP, SEXP, long int);
typedef SEXP (*GRAPH_EMPTY_t)(SEXP, SEXP);
typedef SEXP (*ADD_EDGES_t)(SEXP, SEXP, SEXP);
typedef SEXP (*ADD_VERTEX_ATTRIBUTE_t)(SEXP, SEXP, const char*, SEXP);
typedef SEXP (*SET_VERTEX_ATTRIBUTE_t)(SEXP, SEXP, const char*, SEXP, SEXP);
typedef SEXP (*GET_EDGE_ATTRIBUTE_t)(SEXP, SEXP, SEXP, long int, long int);

typedef struct {
  VCOUNT_t vcount;
  ECOUNT_t ecount;
  NEIGHBORS_t neighbors;
  ADD_VERTICES_t add_vertices;
  GRAPH_EMPTY_t graph_empty;
  ADD_EDGES_t add_edges;
  ADD_VERTEX_ATTRIBUTE_t add_vertex_attribute;
  SET_VERTEX_ATTRIBUTE_t set_vertex_attribute;
  GET_EDGE_ATTRIBUTE_t get_edge_attribute;
} REST_i_ptrtable_t;

extern REST_i_ptrtable_t REST_i_table_default;
extern REST_i_ptrtable_t REST_i_table_adjacencylist;

REST_i_ptrtable_t REST_i_getptrtable(SEXP graph);

SEXP REST_i_default_vcount(SEXP, SEXP);
SEXP REST_i_default_ecount(SEXP, SEXP);
SEXP REST_i_default_neighbors(SEXP, SEXP, long int, SEXP);
SEXP REST_i_default_add_vertices(SEXP, SEXP, long int);
SEXP REST_i_default_graph_empty(SEXP, SEXP);
SEXP REST_i_default_add_edges(SEXP, SEXP, SEXP);
SEXP REST_i_default_add_vertex_attribute(SEXP, SEXP, const char*, SEXP);
SEXP REST_i_default_set_vertex_attribute(SEXP, SEXP, const char*, SEXP, SEXP);
SEXP REST_i_default_get_edge_attribute(SEXP, SEXP, SEXP, long int, long int);

SEXP REST_i_adjacencylist_vcount(SEXP, SEXP);
SEXP REST_i_adjacencylist_ecount(SEXP, SEXP);
SEXP REST_i_adjacencylist_neighbors(SEXP, SEXP, long int, SEXP);

/* -------------------------------------------------- */
/* INTERNALS                                          */
/* -------------------------------------------------- */

SEXP REST_i_get_list_element(SEXP list, const char *str);

/* -------------------------------------------------- */

/**
 */

typedef struct s_dqueue {
  long int *begin;
  long int *end;
  long int *stor_begin;
  long int *stor_end;
} dqueue_t;

int dqueue_init    (dqueue_t* q, long int size);
int dqueue_destroy (dqueue_t* q);
int dqueue_empty   (dqueue_t* q);
int dqueue_clear   (dqueue_t* q);
int dqueue_full    (dqueue_t* q);
long int dqueue_size    (dqueue_t* q);
long int dqueue_pop     (dqueue_t* q);
long int dqueue_pop_back(dqueue_t* q);
long int dqueue_head    (dqueue_t* q);
long int dqueue_back    (dqueue_t* q);
int dqueue_push    (dqueue_t* q, long int elem);

/**
 */

typedef struct s_clustset {
  long int allocated_size;
  long int no_of_nodes;
  long int no_of_clusters;
  long int *parent;
} clustset_t;

int clustset_init    (clustset_t* cs, long int size);
int clustset_destroy (clustset_t* cs);
long int clustset_in_which(clustset_t* cs, long int which);
long int clustset_addnode (clustset_t* cs);
long int clustset_rng     (clustset_t* cs);
long int clustset_merge   (clustset_t* cs, long int first,  long int second);
long int clustset_size    (clustset_t* cs);
long int clustset_nodes   (clustset_t* cs);
long int clustset_csize   (clustset_t* cs, long int which);
dqueue_t* clustset_csizes  (clustset_t* cs);

/**
 */

typedef struct s_vector {
  long int* stor_begin;
  long int* stor_end;
  long int* end;
} vector_t;

int vector_init      (vector_t* v, long int size);
int vector_destroy   (vector_t* v);
int vector_reserve   (vector_t* v, long int size);
int vector_empty     (vector_t* v);
long int vector_size      (vector_t* v);
int vector_clear     (vector_t* v);
int vector_null      (vector_t* v);
int vector_push_back (vector_t* v, long int e);
long int vector_e         (vector_t* v, long int pos);
int vector_set       (vector_t* v, long int pos, long int value);
int vector_add       (vector_t* v, long int pos, long int value);
int vector_replace_first(vector_t* v, long int old, long int new);
long int vector_pop_back(vector_t* v);
long int vector_find(vector_t* v, long int elem);
int vector_change(vector_t* v, long int pos1, long int pos2);

typedef struct s_stack {
  long int* stor_begin;
  long int* stor_end;
  long int* end;
} stack_t;

int stack_init       (stack_t* s, long int size);
int stack_destroy    (stack_t* s);
int stack_reserve    (stack_t* s, long int size);
int stack_empty      (stack_t* s);
long int stack_size       (stack_t* s);
int stack_clear      (stack_t* s);
int stack_push       (stack_t* s, long int elem);
long int stack_pop        (stack_t* s);

typedef struct s_multiset {
  long int* stor_begin;
  long int* stor_end;
  long int*end;
} multiset_t;

int multiset_init    (multiset_t* m, long int size);
int multiset_destroy (multiset_t* m);
int multiset_reserve (multiset_t* m, long int size);
int multiset_add     (multiset_t* m, long int elem);
int multiset_clear   (multiset_t* m);
long int multiset_choose  (multiset_t* m);
long int multiset_choose_random(multiset_t* m);
long int multiset_choose_remove (multiset_t* m);
long int multiset_choose_remove_random(multiset_t* m);
long int multiset_size(multiset_t* m);
int multiset_remove (multiset_t* m, long int elem);
int multiset_remove_all (multiset_t* m, long int elem);
long int* multiset_get_vector(multiset_t * m);
long int multiset_count(multiset_t* m, long int elem);
long int multiset_count_different(multiset_t* m, long int elem);
long int multiset_choose_random_different(multiset_t* m, long int elem);

typedef struct s_heap {
  long int* stor_begin;
  long int* stor_end;
  long int* end;
  int destroy;
} heap_t;

int heap_init           (heap_t* h, long int size);
int heap_init_array     (heap_t *t, long int* data, long int len);
int heap_destroy        (heap_t* h);
int heap_empty          (heap_t* h);
int heap_push           (heap_t* h, long int elem);
long int heap_max       (heap_t* h);
long int heap_delete_max(heap_t* h);
long int heap_size      (heap_t* h);
int heap_reserve        (heap_t* h, long int size);

int heap_i_build(long int* arr, long int size, long int head);
int heap_i_shift_up(long int* arr, long int size, long int elem);
int heap_i_sink(long int* arr, long int size, long int head);
int heap_i_switch(long int* arr, long int e1, long int e2);

typedef struct s_indheap {
  long int* stor_begin;
  long int* stor_end;
  long int* end;
  int destroy;
  long int* index_begin;
} indheap_t;

int indheap_init           (indheap_t* h, long int size);
int indheap_init_array     (indheap_t *t, long int* data, long int len);
int indheap_destroy        (indheap_t* h);
int indheap_empty          (indheap_t* h);
int indheap_push           (indheap_t* h, long int elem);
long int indheap_max       (indheap_t* h);
long int indheap_delete_max(indheap_t* h);
long int indheap_size      (indheap_t* h);
int indheap_reserve        (indheap_t* h, long int size);
long int indheap_max_index(indheap_t *h);

int indheap_i_build(indheap_t* h, long int head);
int indheap_i_shift_up(indheap_t* h, long int elem);
int indheap_i_sink(indheap_t* h, long int head);
int indheap_i_switch(indheap_t* h, long int e1, long int e2);

/* This is a heap containing double elements and 
   two indices, its intended usage is the storage of
   weighted edges.
*/

typedef struct s_indheap_d {
  double* stor_begin;
  double* stor_end;
  double* end;
  int destroy;
  long int* index_begin;
  long int* index2_begin;
} d_indheap_t;

int d_indheap_init           (d_indheap_t* h, long int size);
int d_indheap_destroy        (d_indheap_t* h);
int d_indheap_empty          (d_indheap_t* h);
int d_indheap_push           (d_indheap_t* h, double elem, 
			      long int idx, long int idx2);
double d_indheap_max       (d_indheap_t* h);
double d_indheap_delete_max(d_indheap_t* h);
long int d_indheap_size      (d_indheap_t* h);
int d_indheap_reserve        (d_indheap_t* h, long int size);
int d_indheap_max_index(d_indheap_t *h, long int *idx, long int *idx2);

int d_indheap_i_build(d_indheap_t* h, long int head);
int d_indheap_i_shift_up(d_indheap_t* h, long int elem);
int d_indheap_i_sink(d_indheap_t* h, long int head);
int d_indheap_i_switch(d_indheap_t* h, long int e1, long int e2);

#endif

