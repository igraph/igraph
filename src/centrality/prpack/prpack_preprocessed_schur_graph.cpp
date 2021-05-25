#include "prpack_preprocessed_schur_graph.h"
#include <algorithm>
#include <cstring>
using namespace prpack;
using namespace std;

void prpack_preprocessed_schur_graph::initialize() {
    heads = NULL;
    tails = NULL;
    vals = NULL;
    ii = NULL;
    d = NULL;
    num_outlinks = NULL;
    encoding = NULL;
    decoding = NULL;
}

void prpack_preprocessed_schur_graph::initialize_weighted(const prpack_base_graph* bg) {
    // permute d
    ii = d;
    d = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        d[encoding[i]] = ii[i];
    // convert bg to head/tail format
    for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
        ii[tails_i] = 0;
        tails[tails_i] = heads_i;
        const int decoded = decoding[tails_i];
        const int start_i = bg->tails[decoded];
        const int end_i = (decoded + 1 != num_vs) ? bg->tails[decoded + 1] : bg->num_es;
        for (int i = start_i; i < end_i; ++i) {
            if (decoded == bg->heads[i])
                ii[tails_i] += bg->vals[i];
            else {
                heads[heads_i] = encoding[bg->heads[i]];
                vals[heads_i] = bg->vals[i];
                ++heads_i;
            }
        }
    }
}

void prpack_preprocessed_schur_graph::initialize_unweighted(const prpack_base_graph* bg) {
    // permute num_outlinks
    ii = num_outlinks;
    num_outlinks = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        num_outlinks[encoding[i]] = (ii[i] == 0) ? -1 : ii[i];
    // convert bg to head/tail format
    for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
        ii[tails_i] = 0;
        tails[tails_i] = heads_i;
        const int decoded = decoding[tails_i];
        const int start_i = bg->tails[decoded];
        const int end_i = (decoded + 1 != num_vs) ? bg->tails[decoded + 1] : bg->num_es;
        for (int i = start_i; i < end_i; ++i) {
            if (decoded == bg->heads[i])
                ++ii[tails_i];
            else
                heads[heads_i++] = encoding[bg->heads[i]];
        }
        if (ii[tails_i] > 0)
            ii[tails_i] /= num_outlinks[tails_i];
    }
}

prpack_preprocessed_schur_graph::prpack_preprocessed_schur_graph(const prpack_base_graph* bg) {
    initialize();
    // initialize instance variables
    num_vs = bg->num_vs;
    num_es = bg->num_es - bg->num_self_es;
    tails = new int[num_vs];
    heads = new int[num_es];
    const bool weighted = bg->vals != NULL;
    if (weighted) {
        vals = new double[num_vs];
        d = new double[num_vs];
        fill(d, d + num_vs, 1);
        for (int i = 0; i < bg->num_es; ++i)
            d[bg->heads[i]] -= bg->vals[i];
    } else {
        num_outlinks = new double[num_vs];
        fill(num_outlinks, num_outlinks + num_vs, 0);
        for (int i = 0; i < bg->num_es; ++i)
            ++num_outlinks[bg->heads[i]];
    }
    // permute no-inlink vertices to the beginning, and no-outlink vertices to the end
    encoding = new int[num_vs];
    decoding = new int[num_vs];
    num_no_in_vs = num_no_out_vs = 0;
    for (int i = 0; i < num_vs; ++i) {
        if (bg->tails[i] == ((i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es)) {
            decoding[encoding[i] = num_no_in_vs] = i;
            ++num_no_in_vs;
        } else if ((weighted) ? (d[i] == 1) : (num_outlinks[i] == 0)) {
            decoding[encoding[i] = num_vs - 1 - num_no_out_vs] = i;
            ++num_no_out_vs;
        }
    }
    // permute everything else
    for (int i = 0, p = num_no_in_vs; i < num_vs; ++i)
        if (bg->tails[i] < ((i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es) && ((weighted) ? (d[i] < 1) : (num_outlinks[i] > 0)))
            decoding[encoding[i] = p++] = i;
    // continue initialization based off of weightedness
    if (weighted)
        initialize_weighted(bg);
    else
        initialize_unweighted(bg);
}

prpack_preprocessed_schur_graph::~prpack_preprocessed_schur_graph() {
    delete[] heads;
    delete[] tails;
    delete[] vals;
    delete[] ii;
    delete[] d;
    delete[] num_outlinks;
    delete[] encoding;
    delete[] decoding;
}
