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

#include "hrg_dendro.h"
#include "hrg_graph.h"

using namespace fitHRG;

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
  MTRand mtr;			// TODO: igraph RNG
  
  // Since we're sampling uniformly at random over the equilibrium
  // walk, we just need to do a bunch of MCMC moves and let the
  // sampling happen on its own.
  while (sample_num < num_samples) {
    // Make a single MCMC move
    d->monteCarloMove(dL, flag_taken, 1.0);
		
    // We sample the dendrogram space once every n MCMC moves (on
    // average). Depending on the flags on the command line, we sample
    // different aspects of the dendrograph structure.
    if (t > thresh and mtr.randExc() < ptest) {
      sample_num++;
      d->sampleSplitLikelihoods(sample_num);
    }
		
    t++;
  }
	
  // correct floating-point errors O(n)
  d->refreshLikelihood();
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
  igraph_real_t newMeanL=1e-49;
  igraph_real_t bestL=d->getLikelihood();
  
  while (1) {
    oldMeanL = newMeanL;
    newMeanL = 0.0;
    for (int i=0; i<65536; i++) {
      IGRAPH_CHECK(! d->monteCarloMove(dL, flag_taken, 1.0));
      Likeli = d->getLikelihood();
      if (Likeli > bestL) { bestL = Likeli; }
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
  igraph_real_t bestL;
  dendro *d;

  IGRAPH_CHECK(igraph_hrg_resize(hrg, no_of_nodes));

  d = new dendro;  

  // Convert the igraph graph
  IGRAPH_CHECK(igraph_i_hrg_getgraph(graph, d));

  bestL = d->getLikelihood();

  // Run fixed number of steps, or until convergence
  if (steps > 0) {
    IGRAPH_CHECK(markovChainMonteCarlo(d, steps, hrg));
  } else {
    IGRAPH_CHECK(MCMCEquilibrium_Find(d, hrg));
  }

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
    IGRAPH_ERROR("Give at least one of `sample' and `samples'",IGRAPH_EINVAL);
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
  igraph_real_t bestL;

  if (start && !hrg) {
    IGRAPH_ERROR("`hrg' must be given is `start' is true", IGRAPH_EINVAL);
  }

  d = new dendro;
  
  IGRAPH_CHECK(igraph_i_hrg_getgraph(graph, d));

  if (start) {
    d->importDendrogramStructure(hrg);
  } else {
    IGRAPH_CHECK(MCMCEquilibrium_Find(d, hrg));    
  }

  bestL=d->getLikelihood();

  IGRAPH_CHECK(markovChainMonteCarlo2(d, num_samples));
  
  d->recordConsensusTree(parents, weights);

  return 0;
}
