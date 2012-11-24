#include "prpack_preprocessed_gs_graph.h"
#include <algorithm>
using namespace prpack;
using namespace std;

void prpack_preprocessed_gs_graph::initialize() {
    heads = NULL;
    tails = NULL;
    vals = NULL;
    ii = NULL;
    d = NULL;
    num_outlinks = NULL;
}

void prpack_preprocessed_gs_graph::initialize_weighted(const prpack_base_graph* bg) {
    vals = new double[num_es];
    d = new double[num_vs];
    fill(d, d + num_vs, 1);
    for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
        tails[tails_i] = heads_i;
        ii[tails_i] = 0;
        const int start_j = bg->tails[tails_i];
        const int end_j = (tails_i + 1 != num_vs) ? bg->tails[tails_i + 1]: bg->num_es;
        for (int j = start_j; j < end_j; ++j) {
            if (tails_i == bg->heads[j])
                ii[tails_i] += bg->vals[j];
            else {
                heads[heads_i] = bg->heads[j];
                vals[heads_i] = bg->vals[j];
                ++heads_i;
            }
            d[bg->heads[j]] -= bg->vals[j];
        }
    }
}

void prpack_preprocessed_gs_graph::initialize_unweighted(const prpack_base_graph* bg) {
    num_outlinks = new double[num_vs];
    fill(num_outlinks, num_outlinks + num_vs, 0);
    for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
        tails[tails_i] = heads_i;
        ii[tails_i] = 0;
        const int start_j = bg->tails[tails_i];
        const int end_j = (tails_i + 1 != num_vs) ? bg->tails[tails_i + 1]: bg->num_es;
        for (int j = start_j; j < end_j; ++j) {
            if (tails_i == bg->heads[j])
                ++ii[tails_i];
            else
                heads[heads_i++] = bg->heads[j];
            ++num_outlinks[bg->heads[j]];
        }
    }
    for (int i = 0; i < num_vs; ++i) {
        if (num_outlinks[i] == 0)
            num_outlinks[i] = -1;
        ii[i] /= num_outlinks[i];
    }
}

prpack_preprocessed_gs_graph::prpack_preprocessed_gs_graph(const prpack_base_graph* bg) {
    initialize();
    num_vs = bg->num_vs;
    num_es = bg->num_es - bg->num_self_es;
    heads = new int[num_es];
    tails = new int[num_vs];
    ii = new double[num_vs];
    if (bg->vals != NULL)
        initialize_weighted(bg);
    else
        initialize_unweighted(bg);
}

prpack_preprocessed_gs_graph::~prpack_preprocessed_gs_graph() {
    delete[] heads;
    delete[] tails;
    delete[] vals;
    delete[] ii;
    delete[] d;
    delete[] num_outlinks;
}

