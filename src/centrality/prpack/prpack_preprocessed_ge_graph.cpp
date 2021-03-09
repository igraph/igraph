#include "prpack_preprocessed_ge_graph.h"
#include <algorithm>
using namespace prpack;
using namespace std;

void prpack_preprocessed_ge_graph::initialize() {
    matrix = NULL;
    d = NULL;
}

void prpack_preprocessed_ge_graph::initialize_weighted(const prpack_base_graph* bg) {
    // initialize d
    fill(d, d + num_vs, 1);
    // fill in the matrix
    for (int i = 0, inum_vs = 0; i < num_vs; ++i, inum_vs += num_vs) {
        const int start_j = bg->tails[i];
        const int end_j = (i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es;
        for (int j = start_j; j < end_j; ++j) {
            matrix[inum_vs + bg->heads[j]] += bg->vals[j];
            d[bg->heads[j]] -= bg->vals[j];
        }
    }
}

void prpack_preprocessed_ge_graph::initialize_unweighted(const prpack_base_graph* bg) {
    // fill in the matrix
    for (int i = 0, inum_vs = 0; i < num_vs; ++i, inum_vs += num_vs) {
        const int start_j = bg->tails[i];
        const int end_j = (i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es;
        for (int j = start_j; j < end_j; ++j)
            ++matrix[inum_vs + bg->heads[j]];
    }
    // normalize the columns
    for (int j = 0; j < num_vs; ++j) {
        double sum = 0;
        for (int inum_vs = 0; inum_vs < num_vs*num_vs; inum_vs += num_vs)
            sum += matrix[inum_vs + j];
        if (sum > 0) {
            d[j] = 0;
            const double coeff = 1/sum;
            for (int inum_vs = 0; inum_vs < num_vs*num_vs; inum_vs += num_vs)
                matrix[inum_vs + j] *= coeff;
        } else {
            d[j] = 1;
        }
    }
}

prpack_preprocessed_ge_graph::prpack_preprocessed_ge_graph(const prpack_base_graph* bg) {
    initialize();
    num_vs = bg->num_vs;
    num_es = bg->num_es;
    matrix = new double[num_vs*num_vs];
    d = new double[num_vs];
    fill(matrix, matrix + num_vs*num_vs, 0);
    if (bg->vals != NULL)
        initialize_weighted(bg);
    else
        initialize_unweighted(bg);
}

prpack_preprocessed_ge_graph::~prpack_preprocessed_ge_graph() {
    delete[] matrix;
    delete[] d;
}

