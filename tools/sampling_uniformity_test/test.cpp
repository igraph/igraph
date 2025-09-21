
// Written by Szabolcs Horvát <szhorvat@gmail.com>

// This program demonstrates how to verify that a stochastic graph generator
// samples *uniformly*. We generate random graphs with some method, convert them
// to their adjacency matrices, and count how many times each matrix occurs.
// The frequencies should be the same within expected statistical fluctuations.
// Note that this is a necessary but not sufficient condition for uniform
// sampling from a set of graphs. It must also be verified manually that
// *every* element of the set is generated.

// C++23 is required for the convenience of std::print()

#include <igraph.h>

#include <cmath>
#include <map>
#include <print>
#include <random>
#include <vector>

int main() {
    igraph_t g; // the graph
    igraph_matrix_t am; // adjacency matrix

    // Ensure a different random stream on each run.
    igraph_rng_seed(igraph_rng_default(), std::random_device()());

    igraph_matrix_init(&am, 0, 0);

    // N must be the number of vertices of the graph.
    // Optionally, define other constants used by the graph generator here (such as M below).
    const igraph_int_t N = 4;
    const igraph_int_t M = 4;

    // const igraph_int_t N1 = 2, N2 = 2; // bipartite
    // const igraph_int_t N = N1*N2;

    // Set whether the graph is directed,
    // and whether it can potentially have self-loops and/or multi-edges.
    const bool directed = true;
    const bool loops = true;

    igraph_edge_type_sw_t allowed_edge_types = IGRAPH_MULTI_SW;
    if (loops) allowed_edge_types |= IGRAPH_LOOPS_SW;

    // Use std::map, NOT std::unordered_map, so that the matrices are always
    // printed in the same order. This helps detect consistent patterns which
    // indicate non-uniform sampling.
    std::map<std::vector<int>, int> counts;

    // The adjacency matrix elements are stored contiguously in this vector.
    std::vector<int> v;

    const int samples = 1'000'000; // number of graphs to generate
    for (int iter = 0; iter < samples; iter++) {

        igraph_erdos_renyi_game_gnm(&g, N, M, directed, allowed_edge_types, IGRAPH_EDGE_UNLABELED); // unipartite
        // igraph_bipartite_game_gnm(&g, NULL, N1, N2, M, directed, IGRAPH_ALL, allowed_edge_types, IGRAPH_EDGE_UNLABELED); // bipartite

        igraph_get_adjacency(&g, &am, IGRAPH_GET_ADJACENCY_BOTH, NULL, IGRAPH_LOOPS_TWICE);

        v.clear();
        for (igraph_int_t j=0; j < N; j++) {
            for (igraph_int_t i=0; i < (directed ? N : j+1); i++) {
                if (!loops && i == j) continue;
                v.push_back(MATRIX(am, i, j));
            }
        }

        counts[v] += 1;

        igraph_destroy(&g);
    }

    // Prints the count of each graph/matrix along with the deviation from
    // the expected count, measured in standard deviations.
    {
        double p = 1.0 / counts.size();
        double expected = double(samples) * p;
        double var = double(samples) * p * (1 - p);
        double sd = std::sqrt(var);
        double chi2 = 0, entropy = 0;
        for (auto el : counts) {
            double deviation = expected - el.second;
            chi2 += deviation * deviation / el.second;
            double q = el.second / double(samples);
            if (q > 0) entropy += -q * std::log(q);
            std::println("{} {} {: .2g}", el.first, el.second, deviation / sd);
        }
        std::println();
        std::println("expected counts:    {:.1f}", expected);
        std::println("number of graphs:   {}", counts.size());

        std::println("normalized entropy: {}", entropy / std::log(counts.size()));
        // The normalized entropy should be close to 1.

        std::println("chi^2:              {:.3f}", chi2);
        // Use the chi^2 value in another system to get a p value.
        // A small p value (e.g. < 0.01) indicates that sampling is not uniform.
        // R: pchisq(chi2, ngraphs - 1, lower.tail=F)
        // Mathematica: SurvivalFunction[ChiSquareDistribution[ngraphs - 1], chi2]
        // Python: scipy.stats.chi2.sf(chi2, ngraphs - 1)
        // Rule of thumb: A chi^2 value notably larger than the number of graphs
        // is cause for concern.
    }

    igraph_matrix_destroy(&am);

    return 0;
}
