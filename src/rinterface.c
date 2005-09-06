/* -*- mode: C -*-  */
/* 
   IGraph library R interface.
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

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP R_matrix_to_SEXP(matrix_t *m) {

  SEXP result, dim; 
  long int nrow, ncol;
  
  PROTECT(result=NEW_NUMERIC(matrix_size(m)));
  matrix_copy_to(m, REAL(result));
  PROTECT(dim=NEW_INTEGER(2));
  INTEGER(dim)[0]=matrix_nrow(m);
  INTEGER(dim)[1]=matrix_ncol(m);
  SET_DIM(result, dim);

  UNPROTECT(2);
  return result;
}

SEXP R_igraph_to_SEXP(igraph_t *graph) {
  
  SEXP result;
  long int no_of_nodes=igraph_vcount(graph);
  long int no_of_edges=igraph_ecount(graph);
  
  PROTECT(result=NEW_LIST(8));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(1));
  SET_VECTOR_ELT(result, 1, NEW_LOGICAL(1));
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 3, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 4, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 5, NEW_NUMERIC(no_of_edges));
  SET_VECTOR_ELT(result, 6, NEW_NUMERIC(no_of_nodes+1));
  SET_VECTOR_ELT(result, 7, NEW_NUMERIC(no_of_nodes+1));
  
  REAL(VECTOR_ELT(result, 0))[0]=no_of_nodes;
  LOGICAL(VECTOR_ELT(result, 1))[0]=graph->directed;
  memcpy(REAL(VECTOR_ELT(result, 2)), graph->from.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 3)), graph->to.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 4)), graph->oi.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 5)), graph->ii.stor_begin, 
	 sizeof(real_t)*no_of_edges);
  memcpy(REAL(VECTOR_ELT(result, 6)), graph->os.stor_begin, 
	 sizeof(real_t)*(no_of_nodes+1));
  memcpy(REAL(VECTOR_ELT(result, 7)), graph->is.stor_begin, 
	 sizeof(real_t)*(no_of_nodes+1));
  
  SET_CLASS(result, ScalarString(CREATE_STRING_VECTOR("igraph")));
  
  UNPROTECT(1);
  return result;
}

int R_SEXP_to_matrix(SEXP pakl, matrix_t *akl) {
  akl->data=vector_as_vector(REAL(pakl), GET_LENGTH(pakl));
  akl->nrow=INTEGER(GET_DIM(pakl))[0];
  akl->ncol=INTEGER(GET_DIM(pakl))[1];

  return 0;
}
 
int R_SEXP_to_igraph(SEXP graph, igraph_t *res) {
  
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  res->from=vector_as_vector(REAL(VECTOR_ELT(graph, 2)), 
			    GET_LENGTH(VECTOR_ELT(graph, 2)));
  res->to=vector_as_vector(REAL(VECTOR_ELT(graph, 3)), 
			   GET_LENGTH(VECTOR_ELT(graph, 3)));
  res->oi=vector_as_vector(REAL(VECTOR_ELT(graph, 4)), 
			   GET_LENGTH(VECTOR_ELT(graph, 4)));
  res->ii=vector_as_vector(REAL(VECTOR_ELT(graph, 5)), 
			   GET_LENGTH(VECTOR_ELT(graph, 5)));
  res->os=vector_as_vector(REAL(VECTOR_ELT(graph, 6)), 
			   GET_LENGTH(VECTOR_ELT(graph, 6)));
  res->is=vector_as_vector(REAL(VECTOR_ELT(graph, 7)), 
			   GET_LENGTH(VECTOR_ELT(graph, 7)));
  return 0;
}

int R_SEXP_to_igraph_copy(SEXP graph, igraph_t *res) {
  
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  vector_init_copy(&res->from, REAL(VECTOR_ELT(graph, 2)), 
		   GET_LENGTH(VECTOR_ELT(graph, 2)));
  vector_init_copy(&res->to, REAL(VECTOR_ELT(graph, 3)), 
		   GET_LENGTH(VECTOR_ELT(graph, 3)));
  vector_init_copy(&res->oi, REAL(VECTOR_ELT(graph, 4)), 
		   GET_LENGTH(VECTOR_ELT(graph, 4)));
  vector_init_copy(&res->ii, REAL(VECTOR_ELT(graph, 5)), 
		   GET_LENGTH(VECTOR_ELT(graph, 5)));
  vector_init_copy(&res->os, REAL(VECTOR_ELT(graph, 6)), 
		   GET_LENGTH(VECTOR_ELT(graph, 6)));
  vector_init_copy(&res->is, REAL(VECTOR_ELT(graph, 7)),
		   GET_LENGTH(VECTOR_ELT(graph, 7)));
  return 0;
}

SEXP R_igraph_empty(SEXP pn, SEXP pdirected) {
  
  SEXP result;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  igraph_t g;

  igraph_empty(&g, n, directed);  
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_edges(SEXP graph, SEXP edges) {
  
  vector_t v;			/* do NOT destroy! */
  igraph_t g;
  SEXP result;

  v=vector_as_vector(REAL(edges), GET_LENGTH(edges));
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_edges(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_add_vertices(SEXP graph, SEXP pnv) {
  
  integer_t nv;
  igraph_t g;
  SEXP result;
  
  nv=REAL(pnv)[0];
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_add_vertices(&g, nv);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_vcount(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=igraph_vcount(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_ecount(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=igraph_ecount(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_neighbors(SEXP graph, SEXP pvid, SEXP pmode) {
  
  igraph_t g;
  vector_t neis;
  SEXP result;
  real_t vid;
  integer_t mode;
  
  vector_init(&neis, 0);
  vid=REAL(pvid)[0];
  mode=REAL(pmode)[0];
  R_SEXP_to_igraph(graph, &g);
  igraph_neighbors(&g, &neis, vid, mode);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&neis)));
  vector_copy_to(&neis, REAL(result));
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_delete_edges(SEXP graph, SEXP edges) {
  
  vector_t v;
  igraph_t g;
  SEXP result;
  
  v=vector_as_vector(REAL(edges), GET_LENGTH(edges));
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_delete_edges(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_delete_vertices(SEXP graph, SEXP vertices) {
  
  vector_t v;
  igraph_t g;
  SEXP result;
  
  v=vector_as_vector(REAL(vertices), GET_LENGTH(vertices));
  R_SEXP_to_igraph_copy(graph, &g);
  igraph_delete_vertices(&g, &v);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_is_directed(SEXP graph) {
  
  igraph_t g;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=igraph_is_directed(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_create(SEXP edges, SEXP pn, SEXP pdirected) {
  
  igraph_t g;
  vector_t v;
  integer_t n=REAL(pn)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  v=vector_as_vector(REAL(edges), GET_LENGTH(edges));
  igraph_create(&g, &v, n, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_degree(SEXP graph, SEXP vids, SEXP pmode, SEXP ploops) {
  
  igraph_t g;
  vector_t v;
  vector_t res;
  integer_t mode=REAL(pmode)[0];
  bool_t loops=LOGICAL(ploops)[0];
  SEXP result;

  v=vector_as_vector(REAL(vids), GET_LENGTH(vids));
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, vector_size(&v));
  igraph_degree(&g, &res, &v, mode, loops);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  
  vector_destroy(&res);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_clusters(SEXP graph, SEXP pmode) {

  igraph_t g;
  integer_t mode=REAL(pmode)[0];
  vector_t membership;
  vector_t csize;
  SEXP result, names;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&membership, 0);
  vector_init(&csize, 0);
  igraph_clusters(&g, &membership, &csize, mode);
  
  PROTECT(result=NEW_LIST(2));
  PROTECT(names=NEW_CHARACTER(2));
  SET_VECTOR_ELT(result, 0, NEW_NUMERIC(vector_size(&membership)));
  SET_VECTOR_ELT(result, 1, NEW_NUMERIC(vector_size(&csize)));
  SET_STRING_ELT(names, 0, CREATE_STRING_VECTOR("membership"));
  SET_STRING_ELT(names, 1, CREATE_STRING_VECTOR("csize"));
  SET_NAMES(result, names);
  vector_copy_to(&membership, REAL(VECTOR_ELT(result, 0)));
  vector_copy_to(&csize, REAL(VECTOR_ELT(result, 1)));
  
  vector_destroy(&membership);
  vector_destroy(&csize);
  UNPROTECT(2);
  return result;
}

SEXP R_igraph_is_connected(SEXP graph, SEXP pmode) {

  igraph_t g;
  integer_t mode=REAL(pmode)[0];
  bool_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  igraph_is_connected(&g, &res, mode);
  
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_diameter(SEXP graph, SEXP pdirected, SEXP punconnected) {
  
  igraph_t g;
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t unconnected=LOGICAL(punconnected)[0];
  integer_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  igraph_diameter(&g, &res, directed, unconnected);
  
  PROTECT(result=NEW_NUMERIC(1));
  REAL(result)[0]=res;
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_closeness(SEXP graph, SEXP pvids, SEXP pmode) {
  
  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  integer_t mode=REAL(pmode)[0];
  vector_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_closeness(&g, &res, &vids, mode);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_subcomponent(SEXP graph, SEXP pvertex, SEXP pmode) {
  
  igraph_t g;
  real_t vertex=REAL(pvertex)[0];
  integer_t mode=REAL(pmode)[0];
  vector_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_subcomponent(&g, &res, vertex, mode);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_betweenness(SEXP graph, SEXP pvids, SEXP pdirected) {
  
  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  bool_t directed=LOGICAL(pdirected)[0];
  vector_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_betweenness(&g, &res, &vids, directed);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_running_mean(SEXP pdata, SEXP pbinwidth) {
  
  vector_t data=vector_as_vector(REAL(pdata), GET_LENGTH(pdata));
  integer_t binwidth=REAL(pbinwidth)[0];
  vector_t res;
  SEXP result;
  
  vector_init(&res, 0);
  igraph_running_mean(&data, &res, binwidth);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  vector_destroy(&res);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_cocitation(SEXP graph, SEXP pvids) {

  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  matrix_t m;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&m, 0, 0);
  igraph_cocitation(&g, &m, &vids);
  
  PROTECT(result=R_matrix_to_SEXP(&m));
  matrix_destroy(&m);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_bibcoupling(SEXP graph, SEXP pvids) {

  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  matrix_t m;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&m, 0, 0);
  igraph_bibcoupling(&g, &m, &vids);
  
  PROTECT(result=R_matrix_to_SEXP(&m));
  matrix_destroy(&m);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_growing_random_game(SEXP pn, SEXP pm, SEXP pdirected,
				  SEXP pcitation) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t m=REAL(pm)[0];
  bool_t citation=LOGICAL(pcitation)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  igraph_growing_random_game(&g, n, m, directed, citation);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_shortest_paths(SEXP graph, SEXP pvids, SEXP pmode) {
  
  igraph_t g;
  vector_t vids=vector_as_vector(REAL(pvids), GET_LENGTH(pvids));
  integer_t mode=REAL(pmode)[0];
  matrix_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_shortest_paths(&g, &res, &vids, mode);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_lattice(SEXP pdimvector, SEXP pnei, SEXP pdirected,
		      SEXP pmutual, SEXP pcircular) {
  
  igraph_t g;
  vector_t dimvector=vector_as_vector(REAL(pdimvector), 
				      GET_LENGTH(pdimvector));
  integer_t nei=REAL(pnei)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  bool_t mutual=LOGICAL(pmutual)[0];
  bool_t circular=LOGICAL(pcircular)[0];  
  SEXP result;
  
  igraph_lattice(&g, &dimvector, nei, directed, mutual, circular);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;  
}

SEXP R_igraph_barabasi_game(SEXP pn, SEXP pm, SEXP poutseq,
			    SEXP poutpref, SEXP pdirected) {
  
  igraph_t g;
  integer_t n=REAL(pn)[0];
  integer_t m=REAL(pm)[0]; 
  vector_t outseq=vector_as_vector(REAL(poutseq), GET_LENGTH(poutseq));
  bool_t outpref=LOGICAL(poutpref)[0];
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;
  
  igraph_barabasi_game(&g, n, m, &outseq, outpref, directed);
  PROTECT(result=R_igraph_to_SEXP(&g));
  igraph_destroy(&g);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_layout_kamada_kawai(SEXP graph, SEXP pniter, SEXP pinitemp, 
				  SEXP pcoolexp, SEXP pkkconst, SEXP psigma) {
  
  igraph_t g;
  integer_t niter=REAL(pniter)[0];
  real_t initemp=REAL(pinitemp)[0];
  real_t coolexp=REAL(pcoolexp)[0];
  real_t kkconst=REAL(pkkconst)[0];
  real_t sigma=REAL(psigma)[0];
  matrix_t res;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&res, 0, 0);
  igraph_layout_kamada_kawai(&g, &res, niter, sigma, 
			     initemp, coolexp, kkconst);
  PROTECT(result=R_matrix_to_SEXP(&res));
  matrix_destroy(&res);
  
  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_minimum_spanning_tree_unweighted(SEXP graph) {
  
  igraph_t g;
  igraph_t mst;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_minimum_spanning_tree_unweighted(&g, &mst);
  PROTECT(result=R_igraph_to_SEXP(&mst));
  igraph_destroy(&mst);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_minimum_spanning_tree_prim(SEXP graph, SEXP pweights) {
  
  igraph_t g;
  igraph_t mst;
  vector_t weights=vector_as_vector(REAL(pweights), GET_LENGTH(pweights));
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  igraph_minimum_spanning_tree_prim(&g, &mst, &weights);
  PROTECT(result=R_igraph_to_SEXP(&mst));
  igraph_destroy(&mst);
  
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_edge_betweenness(SEXP graph, SEXP pdirected) {
  
  igraph_t g;
  vector_t res;
  bool_t directed=LOGICAL(pdirected)[0];
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  vector_init(&res, 0);
  igraph_edge_betweenness(&g, &res, directed);
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  
  vector_destroy(&res);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage(SEXP graph, SEXP pst, SEXP pagebins,
				     SEXP pmaxind, SEXP plsd) {
  
  igraph_t g;
  matrix_t akl, sd;
  vector_t st=vector_as_vector(REAL(pst), GET_LENGTH(pst));
  integer_t agebins=REAL(pagebins)[0];
  integer_t maxind=REAL(pmaxind)[0];
  bool_t lsd=LOGICAL(plsd)[0];
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&akl, 0, 0);
  matrix_init(&sd, 0, 0);
  igraph_measure_dynamics_idage(&g, &akl, &sd, &st, agebins, maxind, lsd);
  
  PROTECT(result=NEW_LIST(2));
  SET_VECTOR_ELT(result, 0, R_matrix_to_SEXP(&akl));
  matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_matrix_to_SEXP(&sd));
  matrix_destroy(&sd);

  UNPROTECT(1);
  return result;
}

SEXP R_igraph_measure_dynamics_idage_debug(SEXP graph, SEXP pst, 
					   SEXP pagebins, SEXP pmaxind, 
					   SEXP plsd, SEXP pest_ind, 
					   SEXP pest_age) {
  igraph_t g;
  matrix_t akl, sd;
  vector_t st=vector_as_vector(REAL(pst), GET_LENGTH(pst));
  integer_t agebins=REAL(pagebins)[0];
  integer_t maxind=REAL(pmaxind)[0];
  bool_t lsd=LOGICAL(plsd)[0];
  integer_t est_ind=REAL(pest_ind)[0];
  integer_t est_age=REAL(pest_age)[0];
  vector_t estimates;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  matrix_init(&akl, 0, 0);
  matrix_init(&sd, 0, 0);
  vector_init(&estimates, 0);
  igraph_measure_dynamics_idage_debug(&g, &akl, &sd, &st, agebins, maxind, lsd,
				      &estimates, est_ind, est_age);
  
  PROTECT(result=NEW_LIST(3));
  SET_VECTOR_ELT(result, 0, R_matrix_to_SEXP(&akl));
  matrix_destroy(&akl);
  SET_VECTOR_ELT(result, 1, R_matrix_to_SEXP(&sd));
  matrix_destroy(&sd);
  SET_VECTOR_ELT(result, 2, NEW_NUMERIC(vector_size(&estimates)));
  vector_copy_to(&estimates, REAL(VECTOR_ELT(result, 2)));
  vector_destroy(&estimates);

  UNPROTECT(1);
  return result;
}
  
SEXP R_igraph_measure_dynamics_idage_st(SEXP graph, SEXP pakl) {
  
  igraph_t g;
  matrix_t akl;
  vector_t res;
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  R_SEXP_to_matrix(pakl, &akl);
  vector_init(&res, 0);
  igraph_measure_dynamics_idage_st(&g, &res, &akl);
  
  PROTECT(result=NEW_NUMERIC(vector_size(&res)));
  vector_copy_to(&res, REAL(result));
  
  vector_destroy(&res);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_get_shortest_paths(SEXP graph, SEXP pfrom, SEXP pmode) {

  igraph_t g;
  integer_t from=REAL(pfrom)[0];
  integer_t mode=REAL(pmode)[0];
  vector_t *vects;
  long int no_of_nodes, i;
  SEXP result;
  
  R_SEXP_to_igraph(graph, &g);
  no_of_nodes=igraph_vcount(&g);
  vects=Calloc(no_of_nodes, vector_t);
  for (i=0; i<no_of_nodes; i++) {
    vector_init(&vects[i], 0);
  }
  igraph_get_shortest_paths(&g, vects, from, mode);
  PROTECT(result=NEW_LIST(no_of_nodes));
  for (i=0; i<no_of_nodes; i++) {
    SET_VECTOR_ELT(result, i, NEW_NUMERIC(vector_size(&vects[i])));
    vector_copy_to(&vects[i], REAL(VECTOR_ELT(result, i)));
    vector_destroy(&vects[i]);
  }
  
  Free(vects);
  UNPROTECT(1);
  return result;
}

SEXP R_igraph_are_connected(SEXP graph, SEXP pv1, SEXP pv2) {
  
  igraph_t g;
  integer_t v1=REAL(pv1)[0];
  integer_t v2=REAL(pv2)[0];
  SEXP result;

  R_SEXP_to_igraph(graph, &g);
  PROTECT(result=NEW_LOGICAL(1));
  LOGICAL(result)[0]=igraph_are_connected(&g, v1, v2);
  
  UNPROTECT(1);
  return result;
}
