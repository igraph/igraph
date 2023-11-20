/* -*- mode: C++ -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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
#include "igraph_attributes.h"
#include "igraph_hrg.h"
#include "igraph_random.h"
#include "igraph_structural.h"

#include "hrg/dendro.h"
#include "hrg/graph.h"
#include "hrg/graph_simp.h"

#include "core/exceptions.h"

#include <memory>
#include <climits>

using namespace fitHRG;

/**
 * \section hrg_intro Introduction
 *
 * <para>A hierarchical random graph is an ensemble of undirected
 * graphs with \c n vertices. It is defined via a binary tree with \c
 * n leaf and \c n-1 internal vertices, where the
 * internal vertices are labeled with probabilities.
 * The probability that two vertices are connected in the random graph
 * is given by the probability label at their closest common
 * ancestor.
 * </para>
 *
 * <para>Please read the following two articles for more about
 * hierarchical random graphs: A. Clauset, C. Moore, and M.E.J. Newman.
 * Hierarchical structure and the prediction of missing links in networks.
 * Nature 453, 98 - 101 (2008); and A. Clauset, C. Moore, and M.E.J. Newman.
 * Structural Inference of Hierarchies in Networks. In E. M. Airoldi
 * et al. (Eds.): ICML 2006 Ws, Lecture Notes in Computer Science
 * 4503, 1-13. Springer-Verlag, Berlin Heidelberg (2007).
 * </para>
 *
 * <para>
 * igraph contains functions for fitting HRG models to a given network
 * (\ref igraph_hrg_fit), for generating networks from a given HRG
 * ensemble (\ref igraph_hrg_game, \ref igraph_hrg_sample), converting
 * an igraph graph to a HRG and back (\ref igraph_hrg_create, \ref
 * igraph_hrg_dendrogram), for calculating a consensus tree from a
 * set of sampled HRGs (\ref igraph_hrg_consensus) and for predicting
 * missing edges in a network based on its HRG models (\ref
 * igraph_hrg_predict).
 * </para>
 *
 * <para>The igraph HRG implementation is heavily based on the code
 * published by Aaron Clauset, at his website,
 * http://tuvalu.santafe.edu/~aaronc/hierarchy/
 * </para>
 */

namespace fitHRG {
struct pblock {
    double L;
    int i;
    int j;
};
}

static void markovChainMonteCarlo(dendro &d, const igraph_integer_t period,
                          igraph_hrg_t *hrg) {

    igraph_real_t bestL = d.getLikelihood();
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
    for (igraph_integer_t i = 0; i < period; i++) {

        // make a MCMC move
        d.monteCarloMove(dL, flag_taken, 1.0);

        // get likelihood of this D given G
        igraph_real_t cl = d.getLikelihood();
        if (cl > bestL) {
            // store the current best likelihood
            bestL = cl;
            // record the HRG structure
            d.recordDendrogramStructure(hrg);
        }
    }
    // corrects floating-point errors O(n)
    d.refreshLikelihood();
}

static void markovChainMonteCarlo2(dendro &d, const int num_samples) {
    bool flag_taken;
    double dL;
    const double ptest = 1.0 / (50.0 * static_cast<double>(d.getGraph()->numNodes()));
    igraph_integer_t sample_num = 0;
    int t = 1;
    const int thresh = 200 * d.getGraph()->numNodes();

    // Since we're sampling uniformly at random over the equilibrium
    // walk, we just need to do a bunch of MCMC moves and let the
    // sampling happen on its own.
    while (sample_num < num_samples) {
        // Make a single MCMC move
        d.monteCarloMove(dL, flag_taken, 1.0);

        // We sample the dendrogram space once every n MCMC moves (on
        // average). Depending on the flags on the command line, we sample
        // different aspects of the dendrograph structure.
        if (t > thresh && RNG_UNIF01() < ptest) {
            sample_num++;
            d.sampleSplitLikelihoods();
        }

        t++;

        // correct floating-point errors O(n)
        d.refreshLikelihood(); // TODO: less frequently
    }
}

static void MCMCEquilibrium_Find(dendro &d, igraph_hrg_t *hrg) {

    // We want to run the MCMC until we've found equilibrium; we
    // use the heuristic of the average log-likelihood (which is
    // exactly the entropy) over X steps being very close to the
    // average log-likelihood (entropy) over the X steps that
    // preceded those. In other words, we look for an apparent
    // local convergence of the entropy measure of the MCMC.

    bool flag_taken;
    igraph_real_t dL;
    igraph_real_t newMeanL = -1e-49;

    while (true) {
        const igraph_real_t oldMeanL = newMeanL;
        newMeanL = 0.0;
        for (int i = 0; i < 65536; i++) {
            d.monteCarloMove(dL, flag_taken, 1.0);
            const igraph_real_t Likeli = d.getLikelihood();
            newMeanL += Likeli;
        }
        // corrects floating-point errors O(n)
        d.refreshLikelihood();
        if (fabs(newMeanL - oldMeanL) / 65536.0 < 1.0) {
            break;
        }
    }

    // Record the result
    if (hrg) {
        d.recordDendrogramStructure(hrg);
    }
}

igraph_error_t dendro::setGraph(const igraph_t *igraph) {
    igraph_integer_t no_of_nodes = igraph_vcount(igraph);
    igraph_integer_t no_of_edges = igraph_ecount(igraph);

    if (no_of_nodes > INT_MAX) {
        IGRAPH_ERROR("Graph too large for the HRG module.", IGRAPH_EOVERFLOW);
    }

    // TODO: Can this be relaxed? buildDendrogram() creates a tree with n-2 internal edges,
    // i.e. zero internal edges for a 2-vertex graph. This is not handled at the moment.
    if (no_of_nodes < 3) {
        IGRAPH_ERROR("Graph must have at least 3 vertices for HRG, got only %" IGRAPH_PRId " vertices.", IGRAPH_EINVAL);
    }

    // Create graph
    g = new graph(no_of_nodes);

    // Add edges
    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        int from = IGRAPH_FROM(igraph, i);
        int to = IGRAPH_TO(igraph, i);
        if (from == to) {
            continue;
        }
        if (!g->doesLinkExist(from, to)) {
            g->addLink(from, to);
        }
        if (!g->doesLinkExist(to, from)) {
            g->addLink(to, from);
        }
    }

    buildDendrogram();

    return IGRAPH_SUCCESS;
}

static std::unique_ptr<simpleGraph> igraph_i_hrg_getsimplegraph(const igraph_t *igraph,
                                                                dendro &d,
                                                                igraph_integer_t num_bins) {

    const igraph_integer_t no_of_nodes = igraph_vcount(igraph);
    const igraph_integer_t no_of_edges = igraph_ecount(igraph);

    // TODO replace the following throw's with IGRAPH_ERROR

    if (no_of_nodes > INT_MAX) {
        throw std::runtime_error("Graph too large for the HRG module.");
    }

    // TODO: Can this be relaxed? buildDendrogram() creates a tree with n-2 internal edges,
    // i.e. zero internal edges for a 2-vertex graph. This is not handled at the moment.
    if (no_of_nodes < 3) {
        throw std::runtime_error("Graph must have at least 3 vertices for HRG.");
    }

    // Create graphs
    std::unique_ptr<graph> g(new graph(no_of_nodes, true));
    g->setAdjacencyHistograms(num_bins);

    std::unique_ptr<simpleGraph> sg(new simpleGraph(no_of_nodes));

    for (igraph_integer_t i = 0; i < no_of_edges; i++) {
        int from = (int) IGRAPH_FROM(igraph, i);
        int to = (int) IGRAPH_TO(igraph, i);
        if (from == to) {
            continue;
        }
        if (!g->doesLinkExist(from, to)) {
            g->addLink(from, to);
        }
        if (!g->doesLinkExist(to, from)) {
            g->addLink(to, from);
        }
        if (! sg->doesLinkExist(from, to)) {
            sg->addLink(from, to);
        }
        if (! sg->doesLinkExist(to, from)) {
            sg->addLink(to, from);
        }
    }

    d.setGraph(g.release());
    d.buildDendrogram();

    return sg;
}

/**
 * \function igraph_hrg_init
 * \brief Allocate memory for a HRG.
 *
 * This function must be called before passing an \ref igraph_hrg_t to
 * an igraph function.
 *
 * \param hrg Pointer to the HRG data structure to initialize.
 * \param n The number of vertices in the graph that is modeled by
 *    this HRG. It can be zero, if this is not yet known.
 * \return Error code.
 *
 * Time complexity: O(n), the number of vertices in the graph.
 */

igraph_error_t igraph_hrg_init(igraph_hrg_t *hrg, igraph_integer_t n) {
    if (n < 0) {
        IGRAPH_ERRORF("Number of vertices should not be negative, got %" IGRAPH_PRId ".", IGRAPH_EINVAL, n);
    }
    if (n == 0) {
        n = 1;
    }
    IGRAPH_VECTOR_INT_INIT_FINALLY(&hrg->left,      n - 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&hrg->right,     n - 1);
    IGRAPH_VECTOR_INIT_FINALLY    (&hrg->prob,      n - 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&hrg->edges,     n - 1);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&hrg->vertices,  n - 1);
    IGRAPH_FINALLY_CLEAN(5);
    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hrg_destroy
 * \brief Deallocate memory for an HRG.
 *
 * The HRG data structure can be reinitialized again with an \ref
 * igraph_hrg_destroy call.
 *
 * \param hrg Pointer to the HRG data structure to deallocate.
 *
 * Time complexity: operating system dependent.
 */

void igraph_hrg_destroy(igraph_hrg_t *hrg) {
    igraph_vector_int_destroy(&hrg->left);
    igraph_vector_int_destroy(&hrg->right);
    igraph_vector_destroy(&hrg->prob);
    igraph_vector_int_destroy(&hrg->edges);
    igraph_vector_int_destroy(&hrg->vertices);
}

/**
 * \function igraph_hrg_size
 * \brief Returns the size of the HRG, the number of leaf nodes.
 *
 * \param hrg Pointer to the HRG.
 * \return The number of leaf nodes in the HRG.
 *
 * Time complexity: O(1).
 */

igraph_integer_t igraph_hrg_size(const igraph_hrg_t *hrg) {
    return igraph_vector_int_size(&hrg->left) + 1;
}

/**
 * \function igraph_hrg_resize
 * \brief Resize a HRG.
 *
 * \param hrg Pointer to an initialized (see \ref igraph_hrg_init)
 *   HRG.
 * \param newsize The new size, i.e. the number of leaf nodes.
 * \return Error code.
 *
 * Time complexity: O(n), n is the new size.
 */

igraph_error_t igraph_hrg_resize(igraph_hrg_t *hrg, igraph_integer_t newsize) {
    igraph_integer_t origsize = igraph_hrg_size(hrg);

    /* The data structure must be left in a consistent state if resizing fails. */

#define CHECK_ERR(expr) \
    do { \
        igraph_error_t err = (expr); \
        if (err != IGRAPH_SUCCESS) { \
            igraph_vector_int_resize(&hrg->left, origsize); \
            igraph_vector_int_resize(&hrg->right, origsize); \
            igraph_vector_resize(&hrg->prob, origsize); \
            igraph_vector_int_resize(&hrg->edges, origsize); \
            igraph_vector_int_resize(&hrg->vertices, origsize); \
            IGRAPH_FINALLY_EXIT(); \
            IGRAPH_ERROR("Cannot resize HRG.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */ \
        } \
    } while (0)

    IGRAPH_FINALLY_ENTER();
    {
        CHECK_ERR(igraph_vector_int_resize(&hrg->left, newsize - 1));
        CHECK_ERR(igraph_vector_int_resize(&hrg->right, newsize - 1));
        CHECK_ERR(igraph_vector_resize(&hrg->prob, newsize - 1));
        CHECK_ERR(igraph_vector_int_resize(&hrg->edges, newsize - 1));
        CHECK_ERR(igraph_vector_int_resize(&hrg->vertices, newsize - 1));
    }
    IGRAPH_FINALLY_EXIT();

#undef CHECK_ERR

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hrg_fit
 * \brief Fit a hierarchical random graph model to a network.
 *
 * \param graph The igraph graph to fit the model to. Edge directions
 *   are ignored in directed graphs.
 * \param hrg Pointer to an initialized HRG, the result of the fitting
 *   is stored here. It can also be used to pass a HRG to the
 *   function, that can be used as the starting point of the Markov
 *   Chain Monte Carlo fitting, if the \p start argument is true.
 * \param start Logical, whether to start the fitting from the given
 *   HRG model.
 * \param steps Integer, the number of MCMC steps to take in the
 *   fitting procedure. If this is zero, then the fitting stops if a
 *   convergence criteria is fulfilled.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_hrg_fit(const igraph_t *graph,
                   igraph_hrg_t *hrg,
                   igraph_bool_t start,
                   igraph_integer_t steps) {

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN

    const igraph_integer_t no_of_nodes = igraph_vcount(graph);

    RNG_BEGIN();

    dendro d;

    // If we want to start from HRG
    if (start) {
        if (igraph_hrg_size(hrg) != no_of_nodes) {
            IGRAPH_ERROR("Invalid HRG to start from.", IGRAPH_EINVAL);
        }
        // Convert the igraph graph
        IGRAPH_CHECK(d.setGraph(graph));
        d.clearDendrograph();
        d.importDendrogramStructure(hrg);
    } else {
        // Convert the igraph graph
        IGRAPH_CHECK(d.setGraph(graph));
        IGRAPH_CHECK(igraph_hrg_resize(hrg, no_of_nodes));
    }

    // Run fixed number of steps, or until convergence
    if (steps > 0) {
        markovChainMonteCarlo(d, steps, hrg);
    } else {
        MCMCEquilibrium_Find(d, hrg);
    }

    RNG_END();

    return IGRAPH_SUCCESS;

    IGRAPH_HANDLE_EXCEPTIONS_END
}

/**
 * \function igraph_hrg_sample
 * \brief Sample from a hierarchical random graph model.
 *
 * This function draws a single sample from a hierarchical random graph model.
 *
 * \param hrg A HRG model to sample from
 * \param sample Pointer to an uninitialized graph; the sample is stored here.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_hrg_sample(const igraph_hrg_t *hrg, igraph_t *sample) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN
    dendro d;

    // TODO: error handling

    RNG_BEGIN();

    d.clearDendrograph();
    d.importDendrogramStructure(hrg);
    d.makeRandomGraph();
    IGRAPH_CHECK(d.recordGraphStructure(sample));

    RNG_END();

    return IGRAPH_SUCCESS;
    IGRAPH_HANDLE_EXCEPTIONS_END
}

/**
 * \function igraph_hrg_sample_many
 * \brief Draw multiple samples from a hierarchical random graph model.
 *
 * This function draws multiple samples from the hierarchical random graph
 * ensemble \p hrg.
 *
 * \param hrg A HRG model to sample from
 * \param samples An initialized graph list that will contain the sampled
 *   graphs. Note that existing graphs in the graph list are \em not removed
 *   so make sure you supply an empty list if you do not need the old contents
 *   of the list.
 * \param num_samples The number of samples to generate.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_hrg_sample_many(
    const igraph_hrg_t *hrg, igraph_graph_list_t *samples,
    igraph_integer_t num_samples
) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN
    igraph_t g;
    dendro d;

    if (num_samples < 0) {
        IGRAPH_ERROR("Number of samples must be non-negative.", IGRAPH_EINVAL);
    }

    if (num_samples == 0) {
        return IGRAPH_SUCCESS;
    }

    RNG_BEGIN();

    d.clearDendrograph();
    d.importDendrogramStructure(hrg);
    while (num_samples-- > 0) {
        d.makeRandomGraph();
        IGRAPH_CHECK(d.recordGraphStructure(&g));
        IGRAPH_FINALLY(igraph_destroy, &g);
        IGRAPH_CHECK(igraph_graph_list_push_back(samples, &g));
        IGRAPH_FINALLY_CLEAN(1);
    }

    RNG_END();

    return IGRAPH_SUCCESS;
    IGRAPH_HANDLE_EXCEPTIONS_END
}

/**
 * \function igraph_hrg_game
 * \brief Generate a hierarchical random graph.
 *
 * This function is a simple shortcut to \ref igraph_hrg_sample.
 * It creates a single graph from the given HRG.
 *
 * \param graph Pointer to an uninitialized graph, the new graph is
 *   created here.
 * \param hrg The hierarchical random graph model to sample from.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_hrg_game(igraph_t *graph,
                    const igraph_hrg_t *hrg) {
    return igraph_hrg_sample(hrg, graph);
}

/**
 * \function igraph_from_hrg_dendrogram
 * \brief Create a graph representation of the dendrogram of a hierarchical random graph model.
 *
 * Creates the igraph graph equivalent of the dendrogram encoded in an
 * \ref igraph_hrg_t data structure. The probabilities associated to the
 * nodes are returned in a vector so this function works without an
 * attribute handler.
 *
 * \param graph Pointer to an uninitialized graph, the result is
 *   stored here.
 * \param hrg The hierarchical random graph to convert.
 * \param prob Pointer to an \em initialized vector; the probabilities
 *   associated to the nodes of the dendrogram will be stored here. Leaf nodes
 *   will have an associated probability of \c IGRAPH_NAN .
 *   You may set this to \c NULL if you do not need the probabilities.
 * \return Error code.
 *
 * Time complexity: O(n), the number of vertices in the graph.
 */

igraph_error_t igraph_from_hrg_dendrogram(
    igraph_t *graph, const igraph_hrg_t *hrg, igraph_vector_t *prob
) {
    const igraph_integer_t orig_nodes = igraph_hrg_size(hrg);
    const igraph_integer_t no_of_nodes = orig_nodes * 2 - 1;
    const igraph_integer_t no_of_edges = no_of_nodes > 0 ? no_of_nodes - 1 : 0;
    igraph_vector_int_t edges;
    igraph_integer_t i, idx = 0;

    // Probability labels, for leaf nodes they are IGRAPH_NAN
    if (prob) {
        IGRAPH_CHECK(igraph_vector_resize(prob, no_of_nodes));
        for (i = 0; i < orig_nodes; i++) {
            VECTOR(*prob)[i] = IGRAPH_NAN;
        }
        for (i = 0; i < orig_nodes - 1; i++) {
            VECTOR(*prob)[orig_nodes + i] = VECTOR(hrg->prob)[i];
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);

    for (i = 0; i < orig_nodes - 1; i++) {
        igraph_integer_t left = VECTOR(hrg->left)[i];
        igraph_integer_t right = VECTOR(hrg->right)[i];

        VECTOR(edges)[idx++] = orig_nodes + i;
        VECTOR(edges)[idx++] = left < 0 ? orig_nodes - left - 1 : left;
        VECTOR(edges)[idx++] = orig_nodes + i;
        VECTOR(edges)[idx++] = right < 0 ? orig_nodes - right - 1 : right;
    }

    IGRAPH_CHECK(igraph_empty(graph, 0, IGRAPH_DIRECTED));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, no_of_nodes, NULL));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, NULL));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(2);  // + 1 for graph

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hrg_dendrogram
 * \brief Create a dendrogram from a hierarchical random graph.
 *
 * Creates the igraph graph equivalent of an \ref igraph_hrg_t data
 * structure.
 *
 * \param graph Pointer to an uninitialized graph, the result is
 *   stored here.
 * \param hrg The hierarchical random graph to convert.
 * \return Error code.
 *
 * Time complexity: O(n), the number of vertices in the graph.
 *
 * \deprecated-by igraph_from_hrg_dendrogram 0.10.5
 */
igraph_error_t igraph_hrg_dendrogram(igraph_t *graph, const igraph_hrg_t *hrg) {
    const igraph_integer_t orig_nodes = igraph_hrg_size(hrg);
    const igraph_integer_t no_of_nodes = orig_nodes * 2 - 1;
    const igraph_integer_t no_of_edges = no_of_nodes > 0 ? no_of_nodes - 1 : 0;
    igraph_vector_int_t edges;
    igraph_integer_t i, idx = 0;
    igraph_vector_ptr_t vattrs;
    igraph_vector_t prob;
    igraph_attribute_record_t rec = { "probability",
                                      IGRAPH_ATTRIBUTE_NUMERIC,
                                      &prob
                                    };

    // Probability labels, for leaf nodes they are IGRAPH_NAN
    IGRAPH_VECTOR_INIT_FINALLY(&prob, no_of_nodes);
    for (i = 0; i < orig_nodes; i++) {
        VECTOR(prob)[i] = IGRAPH_NAN;
    }
    for (i = 0; i < orig_nodes - 1; i++) {
        VECTOR(prob)[orig_nodes + i] = VECTOR(hrg->prob)[i];
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, no_of_edges * 2);
    IGRAPH_CHECK(igraph_vector_ptr_init(&vattrs, 1));
    IGRAPH_FINALLY(igraph_vector_ptr_destroy, &vattrs);
    VECTOR(vattrs)[0] = &rec;

    for (i = 0; i < orig_nodes - 1; i++) {
        igraph_integer_t left = VECTOR(hrg->left)[i];
        igraph_integer_t right = VECTOR(hrg->right)[i];

        VECTOR(edges)[idx++] = orig_nodes + i;
        VECTOR(edges)[idx++] = left < 0 ? orig_nodes - left - 1 : left;
        VECTOR(edges)[idx++] = orig_nodes + i;
        VECTOR(edges)[idx++] = right < 0 ? orig_nodes - right - 1 : right;
    }

    IGRAPH_CHECK(igraph_empty(graph, 0, IGRAPH_DIRECTED));
    IGRAPH_FINALLY(igraph_destroy, graph);
    IGRAPH_CHECK(igraph_add_vertices(graph, no_of_nodes, &vattrs));
    IGRAPH_CHECK(igraph_add_edges(graph, &edges, NULL));

    igraph_vector_ptr_destroy(&vattrs);
    igraph_vector_int_destroy(&edges);
    igraph_vector_destroy(&prob);
    IGRAPH_FINALLY_CLEAN(4);  // + 1 for graph

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hrg_consensus
 * \brief Calculate a consensus tree for a HRG.
 *
 * The calculation can be started from the given HRG (\p hrg), or (if
 * \p start is false), a HRG is first fitted to the given graph.
 *
 * \param graph The input graph.
 * \param parents An initialized vector, the results are stored
 *   here. For each vertex, the id of its parent vertex is stored, or
 *   -1, if the vertex is the root vertex in the tree. The first n
 *   vertex IDs (from 0) refer to the original vertices of the graph,
 *   the other IDs refer to vertex groups.
 * \param weights Numeric vector, counts the number of times a given
 *   tree split occured in the generated network samples, for each
 *   internal vertices. The order is the same as in \p parents.
 * \param hrg A hierarchical random graph. It is used as a starting
 *   point for the sampling, if the \p start argument is true. It is
 *   modified along the MCMC.
 * \param start Logical, whether to use the supplied HRG (in \p hrg)
 *   as a starting point for the MCMC.
 * \param num_samples The number of samples to generate for creating
 *   the consensus tree.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_hrg_consensus(const igraph_t *graph,
                         igraph_vector_int_t *parents,
                         igraph_vector_t *weights,
                         igraph_hrg_t *hrg,
                         igraph_bool_t start,
                         igraph_integer_t num_samples) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN

    if (start && !hrg) {
        IGRAPH_ERROR("`hrg' must be given if `start' is true.", IGRAPH_EINVAL);
    }

    RNG_BEGIN();

    dendro d;

    if (start) {
        IGRAPH_CHECK(d.setGraph(graph));
        d.clearDendrograph();
        d.importDendrogramStructure(hrg);
    } else {
        IGRAPH_CHECK(d.setGraph(graph));
        if (hrg) {
            igraph_hrg_resize(hrg, igraph_vcount(graph));
        }
        MCMCEquilibrium_Find(d, hrg);
    }

    markovChainMonteCarlo2(d, num_samples);

    d.recordConsensusTree(parents, weights);

    RNG_END();

    return IGRAPH_SUCCESS;

    IGRAPH_HANDLE_EXCEPTIONS_END
}

static void MCMCEquilibrium_Sample(dendro &d, igraph_integer_t num_samples) {

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
    igraph_integer_t sample_num = 0;
    igraph_integer_t t = 1, thresh = 100 * d.getGraph()->numNodes();
    double ptest = 1.0 / 10.0 / d.getGraph()->numNodes();

    while (sample_num < num_samples) {
        d.monteCarloMove(dL, flag_taken, 1.0);
        if (t > thresh && RNG_UNIF01() < ptest) {
            sample_num++;
            d.sampleAdjacencyLikelihoods();
        }
        d.refreshLikelihood(); // TODO: less frequently
        t++;
    }
}

static igraph_integer_t QsortPartition (pblock* array, igraph_integer_t left, igraph_integer_t right, igraph_integer_t index) {
    pblock p_value = array[index];

    std::swap(array[right], array[index]);

    igraph_integer_t stored = left;
    for (igraph_integer_t i = left; i < right; i++) {
        if (array[i].L <= p_value.L) {
            std::swap(array[i], array[stored]);
            stored++;
        }
    }
    std::swap(array[right], array[stored]);

    return stored;
}

static void QsortMain (pblock* array, igraph_integer_t left, igraph_integer_t right) {
    if (right > left) {
        igraph_integer_t pivot = left;
        igraph_integer_t part  = QsortPartition(array, left, right, pivot);
        QsortMain(array, left,   part - 1);
        QsortMain(array, part + 1, right  );
    }
}

static void rankCandidatesByProbability(const simpleGraph &sg, const dendro &d,
                                pblock *br_list, int mk) {
    int mkk = 0;
    int n = sg.getNumNodes();
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (sg.getAdjacency(i, j) < 0.5) {
                double temp = d.getGraph()->getAdjacencyAverage(i, j);
                br_list[mkk].L = temp * (1.0 + RNG_UNIF01() / 1000.0);
                br_list[mkk].i = i;
                br_list[mkk].j = j;
                mkk++;
            }
        }
    }

    // Sort the candidates by their average probability
    QsortMain(br_list, 0, mk - 1);
}

static igraph_error_t recordPredictions(const pblock *br_list, igraph_vector_int_t *edges,
                      igraph_vector_t *prob, int mk) {

    IGRAPH_CHECK(igraph_vector_int_resize(edges, mk * 2));
    IGRAPH_CHECK(igraph_vector_resize(prob, mk));

    for (int i = mk - 1, idx = 0, idx2 = 0; i >= 0; i--) {
        VECTOR(*edges)[idx++] = br_list[i].i;
        VECTOR(*edges)[idx++] = br_list[i].j;
        VECTOR(*prob)[idx2++] = br_list[i].L;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \function igraph_hrg_predict
 * \brief Predict missing edges in a graph, based on HRG models.
 *
 * Samples HRG models for a network, and estimated the probability
 * that an edge was falsely observed as non-existent in the network.
 *
 * \param graph The input graph.
 * \param edges The list of missing edges is stored here, the first
 *   two elements are the first edge, the next two the second edge,
 *   etc.
 * \param prob Vector of probabilies for the existence of missing
 *   edges, in the order corresponding to \c edges.
 * \param hrg A HRG, it is used as a starting point if \c start is
 *   true. It is also modified during the MCMC sampling.
 * \param start Logical, whether to start the MCMC from the given HRG.
 * \param num_samples The number of samples to generate.
 * \param num_bins Controls the resolution of the edge
 *   probabilities. Higher numbers result higher resolution.
 * \return Error code.
 *
 * Time complexity: TODO.
 */

igraph_error_t igraph_hrg_predict(const igraph_t *graph,
                       igraph_vector_int_t *edges,
                       igraph_vector_t *prob,
                       igraph_hrg_t *hrg,
                       igraph_bool_t start,
                       igraph_integer_t num_samples,
                       igraph_integer_t num_bins) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN

    if (start && !hrg) {
        IGRAPH_ERROR("`hrg' must be given when `start' is true", IGRAPH_EINVAL);
    }

    RNG_BEGIN();

    dendro d;

    std::unique_ptr<simpleGraph> sg = igraph_i_hrg_getsimplegraph(graph, d, num_bins);

    int mk = sg->getNumNodes() * (sg->getNumNodes() - 1) / 2 - sg->getNumLinks() / 2;
    std::unique_ptr<pblock []> br_list(new pblock[mk]);
    for (int i = 0; i < mk; i++) {
        br_list[i].L = 0.0;
        br_list[i].i = -1;
        br_list[i].j = -1;
    }

    if (start) {
        d.clearDendrograph();
        d.importDendrogramStructure(hrg);
    } else {
        if (hrg) {
            igraph_hrg_resize(hrg, igraph_vcount(graph));
        }
        MCMCEquilibrium_Find(d, hrg);
    }

    MCMCEquilibrium_Sample(d, num_samples);
    rankCandidatesByProbability(*sg, d, br_list.get(), mk);
    IGRAPH_CHECK(recordPredictions(br_list.get(), edges, prob, mk));

    RNG_END();

    return IGRAPH_SUCCESS;

    IGRAPH_HANDLE_EXCEPTIONS_END
}

/**
 * \function igraph_hrg_create
 * \brief Create a HRG from an igraph graph.
 *
 * \param hrg Pointer to an initialized \ref igraph_hrg_t. The result
 *    is stored here.
 * \param graph The igraph graph to convert. It must be a directed
 *    binary tree, with n-1 internal and n leaf vertices. The root
 *    vertex must have in-degree zero.
 * \param prob The vector of probabilities, this is used to label the
 *    internal nodes of the hierarchical random graph.
 * \return Error code.
 *
 * Time complexity: O(n), the number of vertices in the tree.
 */

igraph_error_t igraph_hrg_create(igraph_hrg_t *hrg,
                      const igraph_t *graph,
                      const igraph_vector_t *prob) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_internal = no_of_nodes > 0 ? (no_of_nodes - 1) / 2 : 0;
    igraph_vector_int_t deg, idx;
    igraph_integer_t root = 0;
    igraph_integer_t d0 = 0, d1 = 0, d2 = 0;
    igraph_integer_t ii = 0, il = 0;
    igraph_vector_int_t neis;
    igraph_vector_int_t path;
    igraph_bool_t simple;

    // --------------------------------------------------------
    // CHECKS
    // --------------------------------------------------------

    // At least three vertices are required
    if (no_of_nodes < 3) {
        IGRAPH_ERROR("HRG tree must have at least three vertices.",
                     IGRAPH_EINVAL);
    }

    // Prob vector was given
    if (!prob) {
        IGRAPH_ERROR("Probability vector must be given for HRG.",
                     IGRAPH_EINVAL);
    }

    // Length of prob vector
    if (igraph_vector_size(prob) != no_of_nodes / 2) {
        IGRAPH_ERRORF("HRG probability vector size (%" IGRAPH_PRId ") should be equal "
                "to the number of internal nodes (%" IGRAPH_PRId ").", IGRAPH_EINVAL,
                igraph_vector_size(prob), no_of_nodes / 2);
    }

    // Must be a directed graph
    if (!igraph_is_directed(graph)) {
        IGRAPH_ERROR("HRG graph must be directed.", IGRAPH_EINVAL);
    }

    // Number of nodes must be odd
    if (no_of_nodes % 2 == 0) {
        IGRAPH_ERROR("Complete HRG graph must have odd number of vertices.",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_is_simple(graph, &simple));
    if (!simple) {
        IGRAPH_ERROR("HRG graph must be a simple graph.", IGRAPH_EINVAL);
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&deg, 0);

    // Every vertex, except for the root must have in-degree one.
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_IN,
                               IGRAPH_LOOPS));
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t d = VECTOR(deg)[i];
        switch (d) {
        case 0: d0++; root = i; break;
        case 1: d1++; break;
        default:
            IGRAPH_ERROR("HRG nodes must have in-degree one, except for the "
                         "root vertex.", IGRAPH_EINVAL);
        }
    }
    if (d1 != no_of_nodes - 1 || d0 != 1) {
        IGRAPH_ERROR("HRG nodes must have in-degree one, except for the "
                     "root vertex.", IGRAPH_EINVAL);
    }

    // Every internal vertex must have out-degree two,
    // leaves out-degree zero
    d0 = d1 = d2 = 0;
    IGRAPH_CHECK(igraph_degree(graph, &deg, igraph_vss_all(), IGRAPH_OUT,
                               IGRAPH_LOOPS));
    for (int i = 0; i < no_of_nodes; i++) {
        igraph_integer_t d = VECTOR(deg)[i];
        switch (d) {
        case 0: d0++; break;
        case 2: d2++; break;
        default:
            IGRAPH_ERROR("HRG nodes must have out-degree 2 (internal nodes) or "
                         "degree 0 (leaves).", IGRAPH_EINVAL);
        }
    }

    // Number of internal and external nodes is correct
    // This basically checks that the graph has one component
    if (d0 != d2 + 1) {
        IGRAPH_ERROR("HRG degrees are incorrect, maybe multiple components?",
                     IGRAPH_EINVAL);
    }

    // --------------------------------------------------------
    // Graph is good, do the conversion
    // --------------------------------------------------------

    // Create an index, that maps the root node as first, then
    // the internal nodes, then the leaf nodes
    IGRAPH_VECTOR_INT_INIT_FINALLY(&idx, no_of_nodes);
    VECTOR(idx)[root] = - (ii++) - 1;
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t d = VECTOR(deg)[i];
        if (i == root) {
            continue;
        }
        if (d == 2) {
            VECTOR(idx)[i] = - (ii++) - 1;
        }
        if (d == 0) {
            VECTOR(idx)[i] = (il++);
        }
    }

    IGRAPH_CHECK(igraph_hrg_resize(hrg, no_of_internal + 1));
    IGRAPH_VECTOR_INT_INIT_FINALLY(&neis, 0);
    for (igraph_integer_t i = 0; i < no_of_nodes; i++) {
        igraph_integer_t ri = VECTOR(idx)[i];
        if (ri >= 0) {
            continue;
        }
        IGRAPH_CHECK(igraph_neighbors(graph, &neis, i, IGRAPH_OUT));
        VECTOR(hrg->left )[-ri - 1] = VECTOR(idx)[ VECTOR(neis)[0] ];
        VECTOR(hrg->right)[-ri - 1] = VECTOR(idx)[ VECTOR(neis)[1] ];
        VECTOR(hrg->prob )[-ri - 1] = VECTOR(*prob)[i];
    }

    // Calculate the number of vertices and edges in each subtree
    igraph_vector_int_null(&hrg->edges);
    igraph_vector_int_null(&hrg->vertices);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&path, 0);
    IGRAPH_CHECK(igraph_vector_int_push_back(&path, VECTOR(idx)[root]));
    while (!igraph_vector_int_empty(&path)) {
        igraph_integer_t ri = igraph_vector_int_tail(&path);
        igraph_integer_t lc = VECTOR(hrg->left)[-ri - 1];
        igraph_integer_t rc = VECTOR(hrg->right)[-ri - 1];
        if (lc < 0 && VECTOR(hrg->vertices)[-lc - 1] == 0) {
            // Go left
            IGRAPH_CHECK(igraph_vector_int_push_back(&path, lc));
        } else if (rc < 0 && VECTOR(hrg->vertices)[-rc - 1] == 0) {
            // Go right
            IGRAPH_CHECK(igraph_vector_int_push_back(&path, rc));
        } else {
            // Subtrees are done, update node and go up
            VECTOR(hrg->vertices)[-ri - 1] +=
                lc < 0 ? VECTOR(hrg->vertices)[-lc - 1] : 1;
            VECTOR(hrg->vertices)[-ri - 1] +=
                rc < 0 ? VECTOR(hrg->vertices)[-rc - 1] : 1;
            VECTOR(hrg->edges)[-ri - 1] += lc < 0 ? VECTOR(hrg->edges)[-lc - 1] + 1 : 1;
            VECTOR(hrg->edges)[-ri - 1] += rc < 0 ? VECTOR(hrg->edges)[-rc - 1] + 1 : 1;
            igraph_vector_int_pop_back(&path);
        }
    }

    igraph_vector_int_destroy(&path);
    igraph_vector_int_destroy(&neis);
    igraph_vector_int_destroy(&idx);
    igraph_vector_int_destroy(&deg);
    IGRAPH_FINALLY_CLEAN(4);

    return IGRAPH_SUCCESS;
}
