/* -*- mode: C++ -*-  */
/* 
   IGraph library.
   Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
   Rue de l'Industrie 5, Lausanne 1005, Switzerland
   
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

#include "igraph_interface.h"
#include "igraph_community.h"
#include "igraph_memory.h"
#include "igraph_constructors.h"
#include "igraph_attributes.h"
#include "igraph_foreign.h"
#include "igraph_hrg.h"
#include "igraph_random.h"

#include "hrg_dendro.h"
#include "hrg_graph.h"
#include "hrg_graph_simp.h"

using namespace fitHRG;

namespace fitHRG {
  struct pblock { double L; int i; int j; };
}

int markovChainMonteCarlo(dendro *d, unsigned int period, 
			  igraph_hrg_t *hrg) {
  
  igraph_real_t bestL=d->getLikelihood();
  double  dL;
  bool    flag_taken;
  
  // Because moves in the dendrogram space are chosen (Monte
  // Carlo) so that we sample dendrograms with probability
  // proportional to their likelihood, a likelihood-proportional
  // sampling of the dendrogram models would be equivalent to a
  // uniform sampling of the walk itself. We would still have to
  // decide how often to sample the walk (at most once every n
  // steps is recommended) but for simplicity, the code here
  // simply runs the MCMC itself. To actually compute something 
  // over the set of sampled dendrogram models (in a Bayesian
  // model averaging sense), you'll need to code that yourself.
  
  // do 'period' MCMC moves before doing anything else
  for (unsigned int i=0; i<period; i++) {
    
    // make a MCMC move
    IGRAPH_CHECK(! d->monteCarloMove(dL, flag_taken, 1.0));
    
    // get likelihood of this D given G
    igraph_real_t cl= d->getLikelihood(); 
    if (cl > bestL) {
      // store the current best likelihood
      bestL = cl;
      // record the HRG structure
      d->recordDendrogramStructure(hrg);
    }
  }
  // corrects floating-point errors O(n)
  d->refreshLikelihood();	
  
  return 0;
}

int markovChainMonteCarlo2(dendro *d, int num_samples) {
  bool flag_taken;
  double dL, ptest = 1.0/(50.0*(double)(d->g->numNodes()));
  int sample_num=0, t=1, thresh = 200 * d->g->numNodes();
  
  // Since we're sampling uniformly at random over the equilibrium
  // walk, we just need to do a bunch of MCMC moves and let the
  // sampling happen on its own.
  while (sample_num < num_samples) {
    // Make a single MCMC move
    d->monteCarloMove(dL, flag_taken, 1.0);
		
    // We sample the dendrogram space once every n MCMC moves (on
    // average). Depending on the flags on the command line, we sample
    // different aspects of the dendrograph structure.
    if (t > thresh && RNG_UNIF01() < ptest) {
      sample_num++;
      d->sampleSplitLikelihoods(sample_num);
    }
		
    t++;

    // correct floating-point errors O(n)
    d->refreshLikelihood();	// TODO: less frequently    
  }
	
  return 0;
}

int MCMCEquilibrium_Find(dendro *d, igraph_hrg_t *hrg) {

  // We want to run the MCMC until we've found equilibrium; we
  // use the heuristic of the average log-likelihood (which is
  // exactly the entropy) over X steps being very close to the
  // average log-likelihood (entropy) over the X steps that
  // preceded those. In other words, we look for an apparent
  // local convergence of the entropy measure of the MCMC.

  bool flag_taken;
  igraph_real_t dL, Likeli;
  igraph_real_t oldMeanL;
  igraph_real_t newMeanL=-1e-49;
  
  while (1) {
    oldMeanL = newMeanL;
    newMeanL = 0.0;
    for (int i=0; i<65536; i++) {
      IGRAPH_CHECK(! d->monteCarloMove(dL, flag_taken, 1.0));
      Likeli = d->getLikelihood();
      newMeanL += Likeli;
    }
    // corrects floating-point errors O(n)
    d->refreshLikelihood();
    if (fabs(newMeanL-oldMeanL)/65536.0 < 1.0) { break; }	
  }
  
  // Record the result
  if (hrg) { d->recordDendrogramStructure(hrg); }

  return 0;
}

int igraph_i_hrg_getgraph(const igraph_t *igraph,
			  dendro *d) {
  
  int no_of_nodes = igraph_vcount(igraph);
  int no_of_edges = igraph_ecount(igraph);
  int i;
  
  // Create graph
  d->g=new graph(no_of_nodes);  

  // Add edges
  for (i=0; i<no_of_edges; i++) {
    int from=IGRAPH_FROM(igraph, i);
    int to=IGRAPH_TO(igraph, i);
    if (from==to) { continue; }
    if (!d->g->doesLinkExist(from, to)) { d->g->addLink(from, to); }
    if (!d->g->doesLinkExist(to, from)) { d->g->addLink(to, from); }
  }
  
  d->buildDendrogram();

  return 0;
}

int igraph_i_hrg_getsimplegraph(const igraph_t *igraph, 
				dendro *d, simpleGraph **sg, 
				int num_bins) {

  int no_of_nodes = igraph_vcount(igraph);
  int no_of_edges = igraph_ecount(igraph);
  int i;
  
  // Create graphs
  d->g = new graph(no_of_nodes, true);
  d->g->setAdjacencyHistograms(num_bins);
  (*sg) = new simpleGraph(no_of_nodes);
  
  for (i=0; i<no_of_edges; i++) {
    int from=IGRAPH_FROM(igraph, i);
    int to=IGRAPH_TO(igraph, i);
    if (from==to) { continue; }
    if (!d->g->doesLinkExist(from, to)) { d->g->addLink(from, to); }
    if (!d->g->doesLinkExist(to, from)) { d->g->addLink(to, from); }
    if (!(*sg)->doesLinkExist(from, to)) { (*sg)->addLink(from, to); }
    if (!(*sg)->doesLinkExist(to, from)) { (*sg)->addLink(to, from); }
  }

  d->buildDendrogram();
  
  return 0;
}

int igraph_hrg_init(igraph_hrg_t *hrg, int no_of_nodes) {
  IGRAPH_VECTOR_INIT_FINALLY(&hrg->left,      no_of_nodes-1);
  IGRAPH_VECTOR_INIT_FINALLY(&hrg->right,     no_of_nodes-1);
  IGRAPH_VECTOR_INIT_FINALLY(&hrg->prob,      no_of_nodes-1);  
  IGRAPH_VECTOR_INIT_FINALLY(&hrg->edges,     no_of_nodes-1);  
  IGRAPH_VECTOR_INIT_FINALLY(&hrg->vertices,  no_of_nodes-1);  
  IGRAPH_FINALLY_CLEAN(5);
  return 0;
}

void igraph_hrg_destroy(igraph_hrg_t *hrg) {
  igraph_vector_destroy(&hrg->left);
  igraph_vector_destroy(&hrg->right);
  igraph_vector_destroy(&hrg->prob);
  igraph_vector_destroy(&hrg->edges);
  igraph_vector_destroy(&hrg->vertices);
}

int igraph_hrg_size(const igraph_hrg_t *hrg) {
  return igraph_vector_size(&hrg->left)+1;
}

int igraph_hrg_resize(igraph_hrg_t *hrg, int newsize) {
  int origsize=igraph_hrg_size(hrg);
  int ret=0;
  igraph_error_handler_t *oldhandler =
    igraph_set_error_handler(igraph_error_handler_ignore);

  ret  = igraph_vector_resize(&hrg->left, newsize-1);
  ret |= igraph_vector_resize(&hrg->right, newsize-1);
  ret |= igraph_vector_resize(&hrg->prob, newsize-1);
  ret |= igraph_vector_resize(&hrg->edges, newsize-1);
  ret |= igraph_vector_resize(&hrg->vertices, newsize-1);

  igraph_set_error_handler(oldhandler);

  if (ret) { 
    igraph_vector_resize(&hrg->left, origsize);
    igraph_vector_resize(&hrg->right, origsize);
    igraph_vector_resize(&hrg->prob, origsize);
    igraph_vector_resize(&hrg->edges, origsize);
    igraph_vector_resize(&hrg->vertices, origsize);
    IGRAPH_ERROR("Cannot resize HRG", ret);
  }  

  return 0;
}

/**
 * \function igraph_hrg_fit
 * Fit a hierarchical random graph model to a network
 */

int igraph_hrg_fit(const igraph_t *graph, 
		   igraph_hrg_t *hrg,
		   igraph_bool_t start,
		   int steps) {
  
  int no_of_nodes=igraph_vcount(graph);
  dendro *d;

  IGRAPH_CHECK(igraph_hrg_resize(hrg, no_of_nodes));

  RNG_BEGIN();

  d = new dendro;  

  // Convert the igraph graph
  IGRAPH_CHECK(igraph_i_hrg_getgraph(graph, d));

  // Run fixed number of steps, or until convergence
  if (steps > 0) {
    IGRAPH_CHECK(markovChainMonteCarlo(d, steps, hrg));
  } else {
    IGRAPH_CHECK(MCMCEquilibrium_Find(d, hrg));
  }
  
  delete d;

  RNG_END();

  return 0;

}

/**
 * \function igraph_hrg_sample
 * Sample from a hierarchical random graph model
 */

int igraph_hrg_sample(const igraph_t *input_graph,
		      igraph_t *sample,
		      igraph_vector_ptr_t *samples,
		      int no_samples,
		      igraph_hrg_t *hrg,
		      igraph_bool_t start) {

  int i;
  dendro *d;

  if (no_samples < 0) {
    IGRAPH_ERROR("Number of samples must be non-negative", IGRAPH_EINVAL);
  }

  if (!sample && !samples) {
    IGRAPH_ERROR("Give at least one of `sample' and `samples'",
		 IGRAPH_EINVAL);
  }

  if (no_samples != 1 && sample) {
    IGRAPH_ERROR("Number of samples should be one if `sample' is given", 
		 IGRAPH_EINVAL);
  }

  if (no_samples > 1 && !samples) {
    IGRAPH_ERROR("`samples' must be non-null if number of samples "
		 "is larger than 1", IGRAPH_EINVAL);
  }

  if (!start && !input_graph) {
    IGRAPH_ERROR("Input graph must be given if initial HRG is not used", 
		 IGRAPH_EINVAL);
  }

  if (!start) {
    IGRAPH_CHECK(igraph_hrg_resize(hrg, igraph_vcount(input_graph)));
  }

  if (input_graph && igraph_hrg_size(hrg) != igraph_vcount(input_graph)) {
    IGRAPH_ERROR("Invalid HRG size, should match number of nodes", 
		 IGRAPH_EINVAL);
  }

  RNG_BEGIN();

  d = new dendro;

  // Need to find equilibrium first?
  if (start) {
    d->importDendrogramStructure(hrg);
  } else {
    IGRAPH_CHECK(MCMCEquilibrium_Find(d, hrg));
  }

  // TODO: free on error

  if (sample) {
    // A single graph
    d->makeRandomGraph();
    d->recordGraphStructure(sample);
    if (samples) {
      igraph_t *G=igraph_Calloc(1, igraph_t);
      if (!G) { IGRAPH_ERROR("Cannot sample HRG graphs", IGRAPH_ENOMEM); }
      d->recordGraphStructure(G);
      IGRAPH_CHECK(igraph_vector_ptr_resize(samples, 1));
      VECTOR(*samples)[0]=G;
    }
  } else {
    // Sample many     
    IGRAPH_CHECK(igraph_vector_ptr_resize(samples, no_samples));
    for (i=0; i<no_samples; i++) {
      igraph_t *G=igraph_Calloc(1, igraph_t);
      if (!G) { IGRAPH_ERROR("Cannot sample HRG graphs", IGRAPH_ENOMEM); }
      d->makeRandomGraph();
      d->recordGraphStructure(G);
      VECTOR(*samples)[i]=G;
    }
  }

  delete d;

  RNG_END();

  return 0;
}

/** 
 * \function igraph_hrg_game
 * Generate a hierarchical random graph
 */

int igraph_hrg_game(igraph_t *graph,
		    const igraph_hrg_t *hrg) {
  return igraph_hrg_sample(/* input_graph= */ 0, /* sample= */ graph, 
			   /* samples= */ 0, /* no_samples=*/ 1,
			   /* hrg= */ (igraph_hrg_t*) hrg, 
			   /* start= */ 1);
}

int igraph_hrg_dendrogram(igraph_t *graph,
			  const igraph_hrg_t *hrg) {
  
  int orig_nodes=igraph_hrg_size(hrg);
  int no_of_nodes=orig_nodes * 2 - 1;
  int no_of_edges=no_of_nodes-1;
  igraph_vector_t edges;
  int i, idx=0;
  igraph_vector_ptr_t vattrs;
  igraph_vector_t prob;
  igraph_attribute_record_t rec = { "probability", 
				    IGRAPH_ATTRIBUTE_NUMERIC,
				    &prob };
  
  // Probability labels, for leaf nodes they are IGRAPH_NAN
  IGRAPH_VECTOR_INIT_FINALLY(&prob, no_of_nodes);
  for (i=0; i<orig_nodes; i++) {
    VECTOR(prob)[i] = IGRAPH_NAN;
  }
  for (i=0; i<orig_nodes-1; i++) {
    VECTOR(prob)[orig_nodes+i] = VECTOR(hrg->prob)[i];
  }

  IGRAPH_VECTOR_INIT_FINALLY(&edges, no_of_edges * 2);
  IGRAPH_CHECK(igraph_vector_ptr_init(&vattrs, 1));
  IGRAPH_FINALLY(igraph_vector_ptr_destroy, &vattrs);
  VECTOR(vattrs)[0] = &rec;

  for (i=0; i<orig_nodes-1; i++) {
    int left=VECTOR(hrg->left)[i];
    int right=VECTOR(hrg->right)[i];

    VECTOR(edges)[idx++] = orig_nodes+i;
    VECTOR(edges)[idx++] = left < 0 ? orig_nodes-left-1 : left;
    VECTOR(edges)[idx++] = orig_nodes+i;
    VECTOR(edges)[idx++] = right < 0 ? orig_nodes-right-1 : right;
  }

  IGRAPH_CHECK(igraph_empty(graph, 0, IGRAPH_DIRECTED));
  IGRAPH_FINALLY(igraph_destroy, graph);
  IGRAPH_CHECK(igraph_add_vertices(graph, no_of_nodes, &vattrs));
  IGRAPH_CHECK(igraph_add_edges(graph, &edges, 0));

  igraph_vector_ptr_destroy(&vattrs);
  igraph_vector_destroy(&edges);
  igraph_vector_destroy(&prob);
  IGRAPH_FINALLY_CLEAN(4);	// + 1 for graph
  
  return 0;
}

int igraph_hrg_consensus(const igraph_t *graph,
			 igraph_vector_t *parents,
			 igraph_vector_t *weights,
			 igraph_hrg_t *hrg,
			 igraph_bool_t start, 
			 int num_samples) {

  dendro *d;

  if (start && !hrg) {
    IGRAPH_ERROR("`hrg' must be given is `start' is true", IGRAPH_EINVAL);
  }

  RNG_BEGIN();
  
  d = new dendro;
  
  IGRAPH_CHECK(igraph_i_hrg_getgraph(graph, d));

  if (start) {
    d->importDendrogramStructure(hrg);
  } else {
    if (hrg) { igraph_hrg_resize(hrg, igraph_vcount(graph)); }
    IGRAPH_CHECK(MCMCEquilibrium_Find(d, hrg));
  }

  IGRAPH_CHECK(markovChainMonteCarlo2(d, num_samples));

  d->recordConsensusTree(parents, weights);

  delete d;
  
  RNG_END();

  return 0;
}

int MCMCEquilibrium_Sample(dendro *d, int num_samples) {

  // Because moves in the dendrogram space are chosen (Monte
  // Carlo) so that we sample dendrograms with probability
  // proportional to their likelihood, a likelihood-proportional
  // sampling of the dendrogram models would be equivalent to a
  // uniform sampling of the walk itself. We would still have to
  // decide how often to sample the walk (at most once every n steps
  // is recommended) but for simplicity, the code here simply runs the
  // MCMC itself. To actually compute something over the set of
  // sampled dendrogram models (in a Bayesian model averaging sense),
  // you'll need to code that yourself.
  
  double dL;
  bool flag_taken;
  int sample_num=0;
  int t=1, thresh=100 * d->g->numNodes();
  double ptest=1.0/10.0/d->g->numNodes();

  while (sample_num < num_samples) {
    d->monteCarloMove(dL, flag_taken, 1.0);
    if (t > thresh && RNG_UNIF01() < ptest) {
      sample_num++;
      d->sampleAdjacencyLikelihoods();
    }
    d->refreshLikelihood();	// TODO: less frequently
    t++;
  }
  
  return 0;
}

int QsortPartition (pblock* array, int left, int right, int index) {
  pblock p_value, temp;
  p_value.L = array[index].L;
  p_value.i = array[index].i;
  p_value.j = array[index].j;
  
  // swap(array[p_value], array[right])
  temp.L = array[right].L;
  temp.i = array[right].i;
  temp.j = array[right].j;
  array[right].L = array[index].L;
  array[right].i = array[index].i;
  array[right].j = array[index].j;
  array[index].L = temp.L;
  array[index].i = temp.i;
  array[index].j = temp.j;
	
  int stored = left;
  for (int i=left; i<right; i++) {
    if (array[i].L <= p_value.L) {
      // swap(array[stored], array[i])
      temp.L = array[i].L;
      temp.i = array[i].i;
      temp.j = array[i].j;
      array[i].L = array[stored].L;
      array[i].i = array[stored].i;
      array[i].j = array[stored].j;
      array[stored].L = temp.L;
      array[stored].i = temp.i;
      array[stored].j = temp.j;
      stored++;
    }
  }
  // swap(array[right], array[stored])
  temp.L = array[stored].L;
  temp.i = array[stored].i;
  temp.j = array[stored].j;
  array[stored].L = array[right].L;
  array[stored].i = array[right].i;
  array[stored].j = array[right].j;
  array[right].L  = temp.L;
  array[right].i  = temp.i;
  array[right].j  = temp.j;
  
  return stored;
}

void QsortMain (pblock* array, int left, int right) {
  if (right > left) {
    int pivot = left;
    int part  = QsortPartition(array, left, right, pivot);
    QsortMain(array, left,   part-1);
    QsortMain(array, part+1, right  );
  }
  return;
}

int rankCandidatesByProbability(simpleGraph *sg, dendro *d, 
				pblock *br_list, int mk) {
  int mkk=0;
  int n=sg->getNumNodes();
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      if (sg->getAdjacency(i, j) < 0.5) {
	double temp=d->g->getAdjacencyAverage(i, j);
	br_list[mkk].L = temp * (1.0 + RNG_UNIF01()/1000.0);
	br_list[mkk].i = i;
	br_list[mkk].j = j;
	mkk++;
      }
    }
  }
  
  // Sort the candidates by their average probability
  QsortMain(br_list, 0, mk-1);  

  return 0;
}

int recordPredictions(dendro *d, pblock *br_list, igraph_vector_t *edges, 
		      igraph_vector_t *prob, int mk) {
    
  IGRAPH_CHECK(igraph_vector_resize(edges, mk*2));
  IGRAPH_CHECK(igraph_vector_resize(prob, mk));
  
  for (int i=mk-1, idx=0, idx2=0; i>=0; i--) {
    VECTOR(*edges)[idx++] = br_list[i].i;
    VECTOR(*edges)[idx++] = br_list[i].j;
    VECTOR(*prob)[idx2++] = br_list[i].L;
  }

  return 0;
}

int igraph_hrg_predict(const igraph_t *graph,
		       igraph_vector_t *edges,
		       igraph_vector_t *prob,
		       igraph_hrg_t *hrg,
		       igraph_bool_t start, 
		       int num_samples, 
		       int num_bins) {

  dendro *d;
  pblock *br_list;
  int mk;
  simpleGraph *sg;

  if (start && !hrg) {
    IGRAPH_ERROR("`hrg' must be given is `start' is true", IGRAPH_EINVAL);
  }

  RNG_BEGIN();

  d = new dendro;

  IGRAPH_CHECK(igraph_i_hrg_getsimplegraph(graph, d, &sg, num_bins));

  mk = sg->getNumNodes() * (sg->getNumNodes()-1) / 2 - sg->getNumLinks()/2;
  br_list = new pblock[mk];
  for (int i=0; i<mk; i++) { 
    br_list[i].L = 0.0; 
    br_list[i].i = -1; 
    br_list[i].j = -1; 
  }
  
  if (start) {
    d->importDendrogramStructure(hrg);
  } else {
    if (hrg) { igraph_hrg_resize(hrg, igraph_vcount(graph)); }
    IGRAPH_CHECK(MCMCEquilibrium_Find(d, hrg));
  }

  IGRAPH_CHECK(MCMCEquilibrium_Sample(d, num_samples));
  IGRAPH_CHECK(rankCandidatesByProbability(sg, d, br_list, mk));
  IGRAPH_CHECK(recordPredictions(d, br_list, edges, prob, mk));

  delete d;
  delete sg;
  delete [] br_list;

  RNG_END();

  return 0;
}

int igraph_hrg_create(igraph_hrg_t *hrg,
		      const igraph_t *graph, 
		      const igraph_vector_t *prob) {

  int no_of_nodes=igraph_vcount(graph);
  int no_of_internal=(no_of_nodes-1)/2;
  igraph_vector_t deg, idx;
  int root=0;
  int d0=0, d1=0, d2=0;
  int ii=0, il=0;
  igraph_vector_t neis;
  igraph_vector_t path;

  // --------------------------------------------------------
  // CHECKS
  // --------------------------------------------------------

  // At least three vertices are required
  if (no_of_nodes < 3) {
    IGRAPH_ERROR("HRG tree must have at least three vertices",
		 IGRAPH_EINVAL);
  }

  // Prob vector was given
  if (!prob) {
    IGRAPH_ERROR("Probability vector must be given for HRG", 
		 IGRAPH_EINVAL);
  }

  // Length of prob vector
  if (igraph_vector_size(prob) != no_of_nodes) {
    IGRAPH_ERROR("HRG probability vector of wrong size", IGRAPH_EINVAL);
  }

  // Must be a directed graph
  if (!igraph_is_directed(graph)) {
    IGRAPH_ERROR("HRG graph must be directed", IGRAPH_EINVAL);
  }

  // Number of nodes must be odd
  if (! no_of_nodes / 2) {
    IGRAPH_ERROR("Complete HRG graph must have odd number of vertices", 
		 IGRAPH_EINVAL);
  }

  IGRAPH_VECTOR_INIT_FINALLY(&deg, 0);

  // Every vertex, except for the root must have in-degree one.
  IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_IN,
			     IGRAPH_LOOPS));
  for (int i=0; i<no_of_nodes; i++) {
    int d=VECTOR(deg)[i];
    switch (d) {
    case 0: d0++; root=i; break;
    case 1: d1++; break;
    default:
      IGRAPH_ERROR("HRG nodes must have in-degree one, except for the "
		   "root vertex", IGRAPH_EINVAL);
    }
  }
  if (d1 != no_of_nodes-1 || d0 != 1) {
    IGRAPH_ERROR("HRG nodes must have in-degree one, except for the "
		 "root vertex", IGRAPH_EINVAL);
  }
  
  // Every internal vertex must have out-degree two,
  // leaves out-degree zero
  d0=d1=d2=0;
  IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_OUT, 
			     IGRAPH_LOOPS));
  for (int i=0; i<no_of_nodes; i++) {
    int d=VECTOR(deg)[i];
    switch (d) {
    case 0: d0++; break;
    case 2: d2++; break;
    default:
      IGRAPH_ERROR("HRG nodes must have out-degree 2 (internal nodes) or "
		   "degree 0 (leaves)", IGRAPH_EINVAL);
    }            
  }
  
  // Number of internal and external nodes is correct
  // This basically checks that the graph has one component
  if (d0 != d2+1) {
    IGRAPH_ERROR("HRG degrees are incorrect, maybe multiple components?",
		 IGRAPH_EINVAL);
  }
  
  // --------------------------------------------------------
  // Graph is good, do the conversion
  // --------------------------------------------------------

  // Create an index, that maps the root node as first, then
  // the internal nodes, then the leaf nodes
  IGRAPH_VECTOR_INIT_FINALLY(&idx, no_of_nodes);
  VECTOR(idx)[root] = - (ii++) - 1;
  for (int i=0; i<no_of_nodes; i++) {
    int d=VECTOR(deg)[i];
    if (i==root) { continue; }
    if (d==2) { VECTOR(idx)[i] = - (ii++) - 1; }
    if (d==0) { VECTOR(idx)[i] = (il++); }
  }

  igraph_hrg_resize(hrg, no_of_internal+1);
  IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
  for (int i=0; i<no_of_nodes; i++) {
    int ri=VECTOR(idx)[i];
    if (ri >= 0) { continue; }
    IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
    VECTOR(hrg->left )[-ri-1] = VECTOR(idx)[ (int) VECTOR(neis)[0] ];
    VECTOR(hrg->right)[-ri-1] = VECTOR(idx)[ (int) VECTOR(neis)[1] ];
    VECTOR(hrg->prob )[-ri-1] = VECTOR(*prob)[i];
  }

  // Calculate the number of vertices and edges in each subtree
  igraph_vector_null(&hrg->edges);
  igraph_vector_null(&hrg->vertices);
  IGRAPH_VECTOR_INIT_FINALLY(&path, 0);
  IGRAPH_CHECK(igraph_vector_push_back(&path, VECTOR(idx)[root]));
  while (!igraph_vector_empty(&path)) {
    int tail=igraph_vector_tail(&path);
    int ri=igraph_vector_tail(&path);
    int lc=VECTOR(hrg->left)[-ri-1];
    int rc=VECTOR(hrg->right)[-ri-1];
    if (lc < 0 && VECTOR(hrg->vertices)[-lc-1]==0) {
      // Go left
      IGRAPH_CHECK(igraph_vector_push_back(&path, lc));
    } else if (rc < 0 && VECTOR(hrg->vertices)[-rc-1]==0) {
      // Go right
      IGRAPH_CHECK(igraph_vector_push_back(&path, rc));
    } else {
      // Subtrees are done, update node and go up
      VECTOR(hrg->vertices)[-ri-1] += 
	lc < 0 ? VECTOR(hrg->vertices)[-lc-1] : 1;
      VECTOR(hrg->vertices)[-ri-1] += 
	rc < 0 ? VECTOR(hrg->vertices)[-rc-1] : 1;
      VECTOR(hrg->edges)[-ri-1] += lc < 0 ? VECTOR(hrg->edges)[-lc-1]+1 : 1;
      VECTOR(hrg->edges)[-ri-1] += rc < 0 ? VECTOR(hrg->edges)[-rc-1]+1 : 1;
      igraph_vector_pop_back(&path);
    }
  }

  igraph_vector_destroy(&path);
  igraph_vector_destroy(&neis);
  igraph_vector_destroy(&idx);
  igraph_vector_destroy(&deg);
  IGRAPH_FINALLY_CLEAN(4);

  return 0;
}
