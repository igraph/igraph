/*
 *
 * gengraph - generation of random simple connected graphs with prescribed
 *            degree sequence
 *
 * Copyright (C) 2006  Fabien Viger
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "gengraph_qsort.h"
#include "gengraph_hash.h"
#include "gengraph_degree_sequence.h"
#include "gengraph_graph_molloy_hash.h"

#include "igraph_constructors.h"
#include "igraph_error.h"
#include "igraph_progress.h"

#include <stdexcept>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

namespace gengraph {

//_________________________________________________________________________
void graph_molloy_hash::compute_neigh() {
    igraph_integer_t *p = links;
    for (igraph_integer_t i = 0; i < n; i++) {
        neigh[i] = p;
        p += HASH_SIZE(deg[i]);
    }
}

//_________________________________________________________________________
void graph_molloy_hash::compute_size() {
    size = 0;
    for (igraph_integer_t i = 0; i < n; i++) {
        size += HASH_SIZE(deg[i]);
    }
}

//_________________________________________________________________________
void graph_molloy_hash::init() {
    for (igraph_integer_t i = 0; i < size; i++) {
        links[i] = HASH_NONE;
    }
}

//_________________________________________________________________________
graph_molloy_hash::graph_molloy_hash(degree_sequence &degs) {
    alloc(degs);
}

//_________________________________________________________________________
igraph_integer_t graph_molloy_hash::alloc(degree_sequence &degs) {
    n = degs.size();
    a = degs.sum();
    assert(a % 2 == 0);

    deg = degs.seq();
    compute_size();
    deg = new igraph_integer_t[n + size];
    if (deg == NULL) {
        return 0;
    }
    igraph_integer_t i;
    for (i = 0; i < n; i++) {
        deg[i] = degs[i];
    }
    links = deg + n;
    init();
    neigh = new igraph_integer_t*[n];
    if (neigh == NULL) {
        return 0;
    }
    compute_neigh();
    return sizeof(igraph_integer_t *)*n + sizeof(igraph_integer_t) * (n + size);
}

//_________________________________________________________________________
graph_molloy_hash::~graph_molloy_hash() {
    if (deg != NULL) {
        delete[] deg;
    }
    if (neigh != NULL) {
        delete[] neigh;
    }
    deg = NULL;
    neigh = NULL;
}

//_________________________________________________________________________
graph_molloy_hash::graph_molloy_hash(igraph_integer_t *svg) {
    // Read n
    n = *(svg++);
    // Read a
    a = *(svg++);
    assert(a % 2 == 0);
    // Read degree sequence
    degree_sequence dd(n, svg);
    // Build neigh[] and alloc links[]
    alloc(dd);
    // Read links[]
    restore(svg + n);
}

//_________________________________________________________________________
igraph_integer_t *graph_molloy_hash::hard_copy() {
    igraph_integer_t *hc = new igraph_integer_t[2 + n + a / 2]; // to store n,a,deg[] and links[]
    hc[0] = n;
    hc[1] = a;
    memcpy(hc + 2, deg, sizeof(igraph_integer_t)*n);
    igraph_integer_t *p = hc + 2 + n;
    igraph_integer_t *l = links;
    for (igraph_integer_t i = 0; i < n; i++) for (igraph_integer_t j = HASH_SIZE(deg[i]); j--; l++) {
            igraph_integer_t d;
            if ((d = *l) != HASH_NONE && d >= i) {
                *(p++) = d;
            }
        }
    assert(p == hc + 2 + n + a / 2);
    return hc;
}

//_________________________________________________________________________
bool graph_molloy_hash::is_connected() {
    bool *visited = new bool[n];
    igraph_integer_t *buff = new igraph_integer_t[n];
    igraph_integer_t comp_size = depth_search(visited, buff);
    delete[] visited;
    delete[] buff;
    return (comp_size == n);
}

//_________________________________________________________________________
igraph_integer_t* graph_molloy_hash::backup() {
    igraph_integer_t *b = new igraph_integer_t[a / 2];
    igraph_integer_t *c = b;
    igraph_integer_t *p = links;
    for (igraph_integer_t i = 0; i < n; i++)
        for (igraph_integer_t d = HASH_SIZE(deg[i]); d--; p++) if (*p != HASH_NONE && *p > i) {
                *(c++) = *p;
            }
    assert(c == b + (a / 2));
    return b;
}

//_________________________________________________________________________
void graph_molloy_hash::restore(igraph_integer_t* b) {
    init();
    igraph_integer_t i;
    igraph_integer_t *dd = new igraph_integer_t[n];
    memcpy(dd, deg, sizeof(igraph_integer_t)*n);
    for (i = 0; i < n; i++) {
        deg[i] = 0;
    }
    for (i = 0; i < n - 1; i++) {
        while (deg[i] < dd[i]) {
            add_edge(i, *b, dd);
            b++;
        }
    }
    delete[] dd;
}

//_________________________________________________________________________
bool graph_molloy_hash::isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited) {
    if (K < 2) {
        return false;
    }
#ifdef OPT_ISOLATED
    if (K <= deg[v] + 1) {
        return false;
    }
#endif //OPT_ISOLATED
    igraph_integer_t *seen  = Kbuff;
    igraph_integer_t *known = Kbuff;
    igraph_integer_t *max   = Kbuff + K;
    *(known++) = v;
    visited[v] = true;
    bool is_isolated = true;

    while (known != seen) {
        v = *(seen++);
        igraph_integer_t *ww = neigh[v];
        igraph_integer_t w;
        for (igraph_integer_t d = HASH_SIZE(deg[v]); d--; ww++) if ((w = *ww) != HASH_NONE && !visited[w]) {
#ifdef OPT_ISOLATED
                if (K <= deg[w] + 1 || known == max) {
#else //OPT_ISOLATED
                if (known == max) {
#endif //OPT_ISOLATED
                    is_isolated = false;
                    goto end_isolated;
                }
                visited[w] = true;
                *(known++) = w;
            }
    }
end_isolated:
    // Undo the changes to visited[]...
    while (known != Kbuff) {
        visited[*(--known)] = false;
    }
    return is_isolated;
}

//_________________________________________________________________________
int graph_molloy_hash::random_edge_swap(igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited) {
    // Pick two random vertices a and c
    igraph_integer_t f1 = pick_random_vertex();
    igraph_integer_t f2 = pick_random_vertex();
    // Check that f1 != f2
    if (f1 == f2) {
        return 0;
    }
    // Get two random edges (f1,*f1t1) and (f2,*f2t2)
    igraph_integer_t *f1t1 = random_neighbour(f1);
    igraph_integer_t t1 = *f1t1;
    igraph_integer_t *f2t2 = random_neighbour(f2);
    igraph_integer_t t2 = *f2t2;
    // Check simplicity
    if (t1 == t2 || f1 == t2 || f2 == t1) {
        return 0;
    }
    if (is_edge(f1, t2) || is_edge(f2, t1)) {
        return 0;
    }
    // Swap
    igraph_integer_t *f1t2 = H_rpl(neigh[f1], deg[f1], f1t1, t2);
    igraph_integer_t *f2t1 = H_rpl(neigh[f2], deg[f2], f2t2, t1);
    igraph_integer_t *t1f2 = H_rpl(neigh[t1], deg[t1], f1, f2);
    igraph_integer_t *t2f1 = H_rpl(neigh[t2], deg[t2], f2, f1);
    // isolation test
    if (K <= 2) {
        return 1;
    }
    if ( !isolated(f1, K, Kbuff, visited) && !isolated(f2, K, Kbuff, visited) ) {
        return 1;
    }
    // undo swap
    H_rpl(neigh[f1], deg[f1], f1t2, t1);
    H_rpl(neigh[f2], deg[f2], f2t1, t2);
    H_rpl(neigh[t1], deg[t1], t1f2, f1);
    H_rpl(neigh[t2], deg[t2], t2f1, f2);
    return 0;
}

//_________________________________________________________________________
igraph_integer_t graph_molloy_hash::shuffle(igraph_integer_t times,
        igraph_integer_t maxtimes, int type) {
    igraph_progress("Shuffle", 0, 0);
    // assert(verify());
    // counters
    igraph_integer_t nb_swaps = 0;
    igraph_integer_t all_swaps = 0;
    unsigned long cost = 0;
    // window
    double T = double(((a < times) ? a : times) / 10);
    if (type == OPTIMAL_HEURISTICS) {
        T = double(optimal_window());
    }
    if (type == BRUTE_FORCE_HEURISTICS) {
        T = double(times * 2);
    }
    // isolation test parameter, and buffers
    double K = 2.4;
    igraph_integer_t *Kbuff = new igraph_integer_t[int(K) + 1];
    bool *visited = new bool[n];
    for (igraph_integer_t i = 0; i < n; i++) {
        visited[i] = false;
    }
    // Used for monitoring , active only if VERBOSE()
    igraph_integer_t failures = 0;
    igraph_integer_t successes = 0;
    double avg_K = 0;
    double avg_T = 0;
    unsigned long next = times;
    next = 0;

    // Shuffle: while #edge swap attempts validated by connectivity < times ...
    while (times > nb_swaps && maxtimes > all_swaps) {
        // Backup graph
        igraph_integer_t *save = backup();
        // Prepare counters, K, T
        igraph_integer_t swaps = 0;
        igraph_integer_t K_int = 0;
        if (type == FINAL_HEURISTICS || type == BRUTE_FORCE_HEURISTICS) {
            K_int = igraph_integer_t(K);
        }
        igraph_integer_t T_int = (igraph_integer_t)(floor(T));
        if (T_int < 1) {
            T_int = 1;
        }
        // compute cost
        cost += T_int;
        if (K_int > 2) {
            cost += K_int + T_int;
        }
        // Perform T edge swap attempts
        for (igraph_integer_t i = T_int; i > 0; i--) {
            // try one swap
            swaps += random_edge_swap(K_int, Kbuff, visited);
            all_swaps++;
            // Verbose
            if (nb_swaps + swaps > next) {
                next = (nb_swaps + swaps) + (times / 1000 > 100 ? times / 1000 : 100);
                int progress = int(double(nb_swaps + swaps) / double(times));
                igraph_progress("Shuffle",  progress, 0);
            }
        }
        // test connectivity
        cost += (a / 2);
        bool ok = is_connected();
        // performance monitor
        {
            avg_T += double(T_int); avg_K += double(K_int);
            if (ok) {
                successes++;
            } else {
                failures++;
            }
        }
        // restore graph if needed, and count validated swaps
        if (ok) {
            nb_swaps += swaps;
        } else {
            restore(save);
            next = nb_swaps;
        }
        delete[] save;
        // Adjust K and T following the heuristics.
        switch (type) {
            int steps;
        case GKAN_HEURISTICS:
            if (ok) {
                T += 1.0;
            } else {
                T *= 0.5;
            }
            break;
        case FAB_HEURISTICS:
            steps = 50 / (8 + failures + successes);
            if (steps < 1) {
                steps = 1;
            }
            while (steps--) if (ok) {
                    T *= 1.17182818;
                } else {
                    T *= 0.9;
                }
            if (T > double(5 * a)) {
                T = double(5 * a);
            }
            break;
        case FINAL_HEURISTICS:
            if (ok) {
                if ((K + 10.0)*T > 5.0 * double(a)) {
                    K /= 1.03;
                } else {
                    T *= 2;
                }
            } else {
                K *= 1.35;
                delete[] Kbuff;
                Kbuff = new igraph_integer_t[igraph_integer_t(K) + 1];
            }
            break;
        case OPTIMAL_HEURISTICS:
            if (ok) {
                T = double(optimal_window());
            }
            break;
        case BRUTE_FORCE_HEURISTICS:
            K *= 2; delete[] Kbuff; Kbuff = new igraph_integer_t[igraph_integer_t(K) + 1];
            break;
        default:
            throw std::invalid_argument("Error in graph_molloy_hash::shuffle(): Unknown heuristics type.");
        }
    }

    delete[] Kbuff;
    delete[] visited;

    if (maxtimes <= all_swaps) {
        IGRAPH_WARNING("Cannot shuffle graph, maybe it is the only realization of its degree sequence?");
    }

    return nb_swaps;
}

//_________________________________________________________________________

/*
void graph_molloy_hash::print(FILE *f) {
    igraph_integer_t i, j;
    for (i = 0; i < n; i++) {
        fprintf(f, "%" IGRAPH_PRId, i);
        for (j = 0; j < HASH_SIZE(deg[i]); j++) if (neigh[i][j] != HASH_NONE) {
                fprintf(f, " %" IGRAPH_PRId, neigh[i][j]);
            }
        fprintf(f, "\n");
    }
}
*/

igraph_error_t graph_molloy_hash::print(igraph_t *graph) {
    igraph_integer_t i, j;
    igraph_integer_t ptr = 0;
    igraph_vector_int_t edges;

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, a); // every edge is counted twice....

    for (i = 0; i < n; i++) {
        for (j = 0; j < HASH_SIZE(deg[i]); j++) {
            if (neigh[i][j] != HASH_NONE) {
                if (neigh[i][j] > i) {
                    VECTOR(edges)[ptr++] = i;
                    VECTOR(edges)[ptr++] = neigh[i][j];
                }
            }
        }
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, n, /*undirected=*/ 0));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

//_________________________________________________________________________
bool graph_molloy_hash::try_shuffle(igraph_integer_t T, igraph_integer_t K, igraph_integer_t *backup_graph) {
    // init all
    igraph_integer_t *Kbuff = NULL;
    bool *visited = NULL;
    if (K > 2) {
        Kbuff = new igraph_integer_t[K];
        visited = new bool[n];
        for (igraph_integer_t i = 0; i < n; i++) {
            visited[i] = false;
        }
    }
    igraph_integer_t *back = backup_graph;
    if (back == NULL) {
        back = backup();
    }
    // perform T edge swap attempts
    while (T--) {
        random_edge_swap(K, Kbuff, visited);
    }
    // clean
    if (visited != NULL) {
        delete[] visited;
    }
    if (Kbuff   != NULL) {
        delete[] Kbuff;
    }
    // check & restore
    bool yo = is_connected();
    restore(back);
    if (backup_graph == NULL) {
        delete[] back;
    }
    return yo;
}

//_________________________________________________________________________
#define _TRUST_BERNOULLI_LOWER 0.01

bool bernoulli_param_is_lower(int success, int trials, double param) {
    if (double(success) >= double(trials)*param) {
        return false;
    }
    double comb = 1.0;
    double fact = 1.0;
    for (int i = 0; i < success; i++) {
        comb *= double(trials - i);
        fact *= double(i + 1);
    }
    comb /= fact;
    comb *= pow(param, double(success)) * exp(double(trials - success) * log1p(-param));
    double sum = comb;
    while (success && sum < _TRUST_BERNOULLI_LOWER) {
        comb *= double(success) * (1.0 - param) / (double(trials - success) * param);
        sum += comb;
        success--;
    }
    // fprintf(stderr,"bernoulli test : %d/%d success against p=%f -> %s\n",success, trials, param, (sum < _TRUST_BERNOULLI_LOWER) ? "lower" : "can't say");
    return (sum < _TRUST_BERNOULLI_LOWER);
}

//_________________________________________________________________________
#define _MIN_SUCCESS_FOR_BERNOULLI_TRUST 100
double graph_molloy_hash::average_cost(igraph_integer_t T, igraph_integer_t *backup, double min_cost) {
    if (T < 1) {
        return 1e+99;
    }
    int successes = 0;
    int trials = 0;
    while (successes < _MIN_SUCCESS_FOR_BERNOULLI_TRUST &&
           !bernoulli_param_is_lower(successes, trials, 1.0 / min_cost)) {
        if (try_shuffle(T, 0, backup)) {
            successes++;
        }
        trials++;
    }
    if (successes >= _MIN_SUCCESS_FOR_BERNOULLI_TRUST) {
        return double(trials) / double(successes) * (1.0 + double(a / 2) / double(T));
    } else {
        return 2.0 * min_cost;
    }
}

//_________________________________________________________________________
igraph_integer_t graph_molloy_hash::optimal_window() {
    igraph_integer_t Tmax;
    igraph_integer_t optimal_T = 1;
    double min_cost = 1e+99;
    igraph_integer_t *back = backup();
    // on cherche une borne sup pour Tmax
    int been_greater = 0;
    for (Tmax = 1; Tmax <= 5 * a ; Tmax *= 2) {
        double c = average_cost(Tmax, back, min_cost);
        if (c > 1.5 * min_cost) {
            break;
        }
        if (c > 1.2 * min_cost && ++been_greater >= 3) {
            break;
        }
        if (c < min_cost) {
            min_cost = c;
            optimal_T = Tmax;
        }
    }
    // on cherche autour
    double span = 2.0;
    int try_again = 4;
    while (span > 1.05 && optimal_T <= 5 * a) {
        igraph_integer_t T_low  = igraph_integer_t(double(optimal_T) / span);
        igraph_integer_t T_high = igraph_integer_t(double(optimal_T) * span);
        double c_low  = average_cost(T_low, back, min_cost);
        double c_high = average_cost(T_high, back, min_cost);
        if (c_low < min_cost && c_high < min_cost) {
            if (try_again--) {
                continue;
            }
            delete[] back;
            return optimal_T;
        }
        if (c_low < min_cost) {
            optimal_T = T_low;
            min_cost = c_low;
        } else if (c_high < min_cost) {
            optimal_T = T_high;
            min_cost = c_high;
        };
        span = pow(span, 0.618);
    }
    delete[] back;
    return optimal_T;
}

//_________________________________________________________________________
double graph_molloy_hash::eval_K(int quality) {
    double K = 5.0;
    double avg_K = 1.0;
    for (int i = quality; i--; ) {
        int int_K = int(floor(K + 0.5));
        if (try_shuffle(a / (int_K + 1), int_K)) {
            K *= 0.8; /*fprintf(stderr,"+");*/
        } else {
            K *= 1.25; /*fprintf(stderr,"-");*/
        }
        if (i < quality / 2) {
            avg_K *= K;
        }
    }
    return pow(avg_K, 1.0 / double(quality / 2));
}

//_________________________________________________________________________
double graph_molloy_hash::effective_K(igraph_integer_t K, int quality) {
    if (K < 3) {
        return 0.0;
    }
    long sum_K = 0;
    igraph_integer_t *Kbuff = new igraph_integer_t[K];
    bool *visited = new bool[n];
    igraph_integer_t i;
    for (i = 0; i < n; i++) {
        visited[i] = false;
    }
    for (i = 0; i < quality; i++) {
        // assert(verify());
        igraph_integer_t f1, f2, t1, t2;
        igraph_integer_t *f1t1, *f2t2;
        do {
            // Pick two random vertices
            do {
                f1 = pick_random_vertex();
                f2 = pick_random_vertex();
            } while (f1 == f2);
            // Pick two random neighbours
            f1t1 = random_neighbour(f1);
            t1 = *f1t1;
            f2t2 = random_neighbour(f2);
            t2 = *f2t2;
            // test simplicity
        } while (t1 == t2 || f1 == t2 || f2 == t1 || is_edge(f1, t2) || is_edge(f2, t1));
        // swap
        swap_edges(f1, t2, f2, t1);
        // assert(verify());
        sum_K += effective_isolated(deg[f1] > deg[t2] ? f1 : t2, K, Kbuff, visited);
        // assert(verify());
        sum_K += effective_isolated(deg[f2] > deg[t1] ? f2 : t1, K, Kbuff, visited);
        // assert(verify());
        // undo swap
        swap_edges(f1, t2, f2, t1);
        // assert(verify());
    }
    delete[] Kbuff;
    delete[] visited;
    return double(sum_K) / double(2 * quality);
}

//_________________________________________________________________________
igraph_integer_t graph_molloy_hash::effective_isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited) {
    igraph_integer_t i;
    for (i = 0; i < K; i++) {
        Kbuff[i] = -1;
    }
    igraph_integer_t count = 0;
    igraph_integer_t left = K;
    igraph_integer_t *KB = Kbuff;
    //yapido = (my_random()%1000 == 0);
    depth_isolated(v, count, left, K, KB, visited);
    while (KB-- != Kbuff) {
        visited[*KB] = false;
    }
    //if(yapido) fprintf(stderr,"\n");
    return count;
}

//_________________________________________________________________________
void graph_molloy_hash::depth_isolated(igraph_integer_t v, igraph_integer_t &calls, igraph_integer_t &left_to_explore, igraph_integer_t dmax, igraph_integer_t * &Kbuff, bool *visited) {
    if (left_to_explore == 0) {
        return;
    }
//  if(yapido) fprintf(stderr,"%d ",deg[v]);
    if (--left_to_explore == 0) {
        return;
    }
    if (deg[v] + 1 >= dmax) {
        left_to_explore = 0;
        return;
    }
    *(Kbuff++) = v;
    visited[v] = true;
//  print();
//  fflush(stdout);
    calls++;
    igraph_integer_t *copy = NULL;
    igraph_integer_t *w = neigh[v];
    if (IS_HASH(deg[v])) {
        copy = new igraph_integer_t[deg[v]];
        H_copy(copy, w, deg[v]);
        w = copy;
    }
    qsort(deg, w, deg[v]);
    w += deg[v];
    for (igraph_integer_t i = deg[v]; i--; ) {
        if (visited[*--w]) {
            calls++;
        } else {
            depth_isolated(*w, calls, left_to_explore, dmax, Kbuff, visited);
        }
        if (left_to_explore == 0) {
            break;
        }
    }
    if (copy != NULL) {
        delete[] copy;
    }
}

//_________________________________________________________________________
igraph_integer_t graph_molloy_hash::depth_search(bool *visited, igraph_integer_t *buff, igraph_integer_t v0) {
    for (igraph_integer_t i = 0; i < n; i++) {
        visited[i] = false;
    }
    igraph_integer_t *to_visit = buff;
    igraph_integer_t nb_visited = 1;
    visited[v0] = true;
    *(to_visit++) = v0;
    while (to_visit != buff && nb_visited < n) {
        igraph_integer_t v = *(--to_visit);
        igraph_integer_t *ww = neigh[v];
        igraph_integer_t w;
        for (igraph_integer_t k = HASH_SIZE(deg[v]); k--; ww++) {
            if (HASH_NONE != (w = *ww) && !visited[w]) {
                visited[w] = true;
                nb_visited++;
                *(to_visit++) = w;
            }
        }
    }
    return nb_visited;
}

//_________________________________________________________________________
// bool graph_molloy_hash::verify() {
//   fprintf(stderr,"Warning: graph_molloy_hash::verify() called..\n");
//   fprintf(stderr,"   try to convert graph into graph_molloy_opt() instead\n");
//   return true;
// }

} // namespace gengraph
