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
#include "gengraph_definitions.h"
#include <stdexcept>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "gengraph_qsort.h"
#include "gengraph_hash.h"
#include "gengraph_degree_sequence.h"
#include "gengraph_graph_molloy_hash.h"

#include "config.h"
#include "core/math.h"
#include "igraph_constructors.h"
#include "igraph_error.h"
#include "igraph_statusbar.h"
#include "igraph_progress.h"

namespace gengraph {

//_________________________________________________________________________
void graph_molloy_hash::compute_neigh() {
    int *p = links;
    for (int i = 0; i < n; i++) {
        neigh[i] = p;
        p += HASH_SIZE(deg[i]);
    }
}

//_________________________________________________________________________
void graph_molloy_hash::compute_size() {
    size = 0;
    for (int i = 0; i < n; i++) {
        size += HASH_SIZE(deg[i]);
    }
}

//_________________________________________________________________________
void graph_molloy_hash::init() {
    for (int i = 0; i < size; i++) {
        links[i] = HASH_NONE;
    }
}

//_________________________________________________________________________
graph_molloy_hash::graph_molloy_hash(degree_sequence &degs) {
    igraph_status("Allocating memory for graph...", 0);
    int s = alloc(degs);
    igraph_statusf("%d bytes allocated successfully\n", 0, s);
}

//_________________________________________________________________________
int graph_molloy_hash::alloc(degree_sequence &degs) {
    n = degs.size();
    a = degs.sum();
    assert(a % 2 == 0);

    deg = degs.seq();
    compute_size();
    deg = new int[n + size];
    if (deg == NULL) {
        return 0;
    }
    int i;
    for (i = 0; i < n; i++) {
        deg[i] = degs[i];
    }
    links = deg + n;
    init();
    neigh = new int*[n];
    if (neigh == NULL) {
        return 0;
    }
    compute_neigh();
    return sizeof(int *)*n + sizeof(int) * (n + size);
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
graph_molloy_hash::graph_molloy_hash(int *svg) {
    // Read n
    n = *(svg++);
    // Read a
    a = *(svg++);
    assert(a % 2 == 0);
    // Read degree sequence
    degree_sequence dd(n, svg);
    // Build neigh[] and alloc links[]
    alloc(dd);
    dd.detach();
    // Read links[]
    restore(svg + n);
}

//_________________________________________________________________________
int *graph_molloy_hash::hard_copy() {
    int *hc = new int[2 + n + a / 2]; // to store n,a,deg[] and links[]
    hc[0] = n;
    hc[1] = a;
    memcpy(hc + 2, deg, sizeof(int)*n);
    int *p = hc + 2 + n;
    int *l = links;
    for (int i = 0; i < n; i++) for (int j = HASH_SIZE(deg[i]); j--; l++) {
            int d;
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
    int *buff = new int[n];
    int comp_size = depth_search(visited, buff);
    delete[] visited;
    delete[] buff;
    return (comp_size == n);
}

//_________________________________________________________________________
int* graph_molloy_hash::backup() {
    int *b = new int[a / 2];
    int *c = b;
    int *p = links;
    for (int i = 0; i < n; i++)
        for (int d = HASH_SIZE(deg[i]); d--; p++) if (*p != HASH_NONE && *p > i) {
                *(c++) = *p;
            }
    assert(c == b + (a / 2));
    return b;
}

//_________________________________________________________________________
void graph_molloy_hash::restore(int* b) {
    init();
    int i;
    int *dd = new int[n];
    memcpy(dd, deg, sizeof(int)*n);
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
bool graph_molloy_hash::isolated(int v, int K, int *Kbuff, bool *visited) {
    if (K < 2) {
        return false;
    }
#ifdef OPT_ISOLATED
    if (K <= deg[v] + 1) {
        return false;
    }
#endif //OPT_ISOLATED
    int *seen  = Kbuff;
    int *known = Kbuff;
    int *max   = Kbuff + K;
    *(known++) = v;
    visited[v] = true;
    bool is_isolated = true;

    while (known != seen) {
        v = *(seen++);
        int *ww = neigh[v];
        int w;
        for (int d = HASH_SIZE(deg[v]); d--; ww++) if ((w = *ww) != HASH_NONE && !visited[w]) {
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
int graph_molloy_hash::random_edge_swap(int K, int *Kbuff, bool *visited) {
    // Pick two random vertices a and c
    int f1 = pick_random_vertex();
    int f2 = pick_random_vertex();
    // Check that f1 != f2
    if (f1 == f2) {
        return 0;
    }
    // Get two random edges (f1,*f1t1) and (f2,*f2t2)
    int *f1t1 = random_neighbour(f1);
    int t1 = *f1t1;
    int *f2t2 = random_neighbour(f2);
    int t2 = *f2t2;
    // Check simplicity
    if (t1 == t2 || f1 == t2 || f2 == t1) {
        return 0;
    }
    if (is_edge(f1, t2) || is_edge(f2, t1)) {
        return 0;
    }
    // Swap
    int *f1t2 = H_rpl(neigh[f1], deg[f1], f1t1, t2);
    int *f2t1 = H_rpl(neigh[f2], deg[f2], f2t2, t1);
    int *t1f2 = H_rpl(neigh[t1], deg[t1], f1, f2);
    int *t2f1 = H_rpl(neigh[t2], deg[t2], f2, f1);
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
unsigned long graph_molloy_hash::shuffle(unsigned long times,
        unsigned long maxtimes, int type) {
    igraph_progress("Shuffle", 0, 0);
    // assert(verify());
    // counters
    unsigned long nb_swaps = 0;
    unsigned long all_swaps = 0;
    unsigned long cost = 0;
    // window
    double T = double(min((unsigned long)(a), times) / 10);
    if (type == OPTIMAL_HEURISTICS) {
        T = double(optimal_window());
    }
    if (type == BRUTE_FORCE_HEURISTICS) {
        T = double(times * 2);
    }
    // isolation test parameter, and buffers
    double K = 2.4;
    int *Kbuff = new int[int(K) + 1];
    bool *visited = new bool[n];
    for (int i = 0; i < n; i++) {
        visited[i] = false;
    }
    // Used for monitoring , active only if VERBOSE()
    int failures = 0;
    int successes = 0;
    double avg_K = 0;
    double avg_T = 0;
    unsigned long next = times;
    next = 0;

    // Shuffle: while #edge swap attempts validated by connectivity < times ...
    while (times > nb_swaps && maxtimes > all_swaps) {
        // Backup graph
        int *save = backup();
        // Prepare counters, K, T
        unsigned long swaps = 0;
        int K_int = 0;
        if (type == FINAL_HEURISTICS || type == BRUTE_FORCE_HEURISTICS) {
            K_int = int(K);
        }
        unsigned long T_int = (unsigned long)(floor(T));
        if (T_int < 1) {
            T_int = 1;
        }
        // compute cost
        cost += T_int;
        if (K_int > 2) {
            cost += (unsigned long)(K_int) * (unsigned long)(T_int);
        }
        // Perform T edge swap attempts
        for (int i = T_int; i > 0; i--) {
            // try one swap
            swaps += (unsigned long)(random_edge_swap(K_int, Kbuff, visited));
            all_swaps++;
            // Verbose
            if (nb_swaps + swaps > next) {
                next = (nb_swaps + swaps) + max((unsigned long)(100), (unsigned long)(times / 1000));
                int progress = int(double(nb_swaps + swaps) / double(times));
                igraph_progress("Shuffle",  progress, 0);
            }
        }
        // test connectivity
        cost += (unsigned long)(a / 2);
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
                Kbuff = new int[int(K) + 1];
            }
            break;
        case OPTIMAL_HEURISTICS:
            if (ok) {
                T = double(optimal_window());
            }
            break;
        case BRUTE_FORCE_HEURISTICS:
            K *= 2; delete[] Kbuff; Kbuff = new int[int(K) + 1];
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

    // Status report
    {
        igraph_status("*** Shuffle Monitor ***\n", 0);
        igraph_statusf(" - Average cost : %f / validated edge swap\n", 0,
                       double(cost) / double(nb_swaps));
        igraph_statusf(" - Connectivity tests : %d (%d successes, %d failures)\n",
                       0, successes + failures, successes, failures);
        igraph_statusf(" - Average window : %d\n", 0,
                       int(avg_T / double(successes + failures)));
        if (type == FINAL_HEURISTICS || type == BRUTE_FORCE_HEURISTICS)
            igraph_statusf(" - Average isolation test width : %f\n", 0,
                           avg_K / double(successes + failures));
    }
    return nb_swaps;
}

//_________________________________________________________________________
void graph_molloy_hash::print(FILE *f) {
    int i, j;
    for (i = 0; i < n; i++) {
        fprintf(f, "%d", i);
        for (j = 0; j < HASH_SIZE(deg[i]); j++) if (neigh[i][j] != HASH_NONE) {
                fprintf(f, " %d", neigh[i][j]);
            }
        fprintf(f, "\n");
    }
}

int graph_molloy_hash::print(igraph_t *graph) {
    int i, j;
    long int ptr = 0;
    igraph_vector_t edges;

    IGRAPH_VECTOR_INIT_FINALLY(&edges, a); // every edge is counted twice....

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
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

//_________________________________________________________________________
bool graph_molloy_hash::try_shuffle(int T, int K, int *backup_graph) {
    // init all
    int *Kbuff = NULL;
    bool *visited = NULL;
    if (K > 2) {
        Kbuff = new int[K];
        visited = new bool[n];
        for (int i = 0; i < n; i++) {
            visited[i] = false;
        }
    }
    int *back = backup_graph;
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
double graph_molloy_hash::average_cost(int T, int *backup, double min_cost) {
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
int graph_molloy_hash::optimal_window() {
    int Tmax;
    int optimal_T = 1;
    double min_cost = 1e+99;
    int *back = backup();
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
        igraph_statusf("Tmax = %d [%f]", 0, Tmax, min_cost);
    }
    // on cree Tmin
    int Tmin = int(0.5 * double(a) / (min_cost - 1.0));
    igraph_statusf("Optimal T is in [%d, %d]\n", 0, Tmin, Tmax);
    // on cherche autour
    double span = 2.0;
    int try_again = 4;
    while (span > 1.05 && optimal_T <= 5 * a) {
        igraph_statusf("Best T [cost]: %d [%f]", 0, optimal_T, min_cost);
        int T_low  = int(double(optimal_T) / span);
        int T_high = int(double(optimal_T) * span);
        double c_low  = average_cost(T_low, back, min_cost);
        double c_high = average_cost(T_high, back, min_cost);
        if (c_low < min_cost && c_high < min_cost) {
            if (try_again--) {
                continue;
            }
            {
                igraph_status("Warning: when looking for optimal T,\n", 0);
                igraph_statusf("Low: %d [%f]  Middle: %d [%f]  High: %d [%f]\n", 0,
                               T_low, c_low, optimal_T, min_cost, T_high, c_high);
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
double graph_molloy_hash::effective_K(int K, int quality) {
    if (K < 3) {
        return 0.0;
    }
    long sum_K = 0;
    int *Kbuff = new int[K];
    bool *visited = new bool[n];
    int i;
    for (i = 0; i < n; i++) {
        visited[i] = false;
    }
    for (int i = 0; i < quality; i++) {
        // assert(verify());
        int f1, f2, t1, t2;
        int *f1t1, *f2t2;
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
long graph_molloy_hash::effective_isolated(int v, int K, int *Kbuff, bool *visited) {
    int i;
    for (i = 0; i < K; i++) {
        Kbuff[i] = -1;
    }
    long count = 0;
    int left = K;
    int *KB = Kbuff;
    //yapido = (my_random()%1000 == 0);
    depth_isolated(v, count, left, K, KB, visited);
    while (KB-- != Kbuff) {
        visited[*KB] = false;
    }
    //if(yapido) fprintf(stderr,"\n");
    return count;
}

//_________________________________________________________________________
void graph_molloy_hash::depth_isolated(int v, long &calls, int &left_to_explore, int dmax, int * &Kbuff, bool *visited) {
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
    int *copy = NULL;
    int *w = neigh[v];
    if (IS_HASH(deg[v])) {
        copy = new int[deg[v]];
        H_copy(copy, w, deg[v]);
        w = copy;
    }
    qsort(deg, w, deg[v]);
    w += deg[v];
    for (int i = deg[v]; i--; ) {
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
int graph_molloy_hash::depth_search(bool *visited, int *buff, int v0) {
    for (int i = 0; i < n; i++) {
        visited[i] = false;
    }
    int *to_visit = buff;
    int nb_visited = 1;
    visited[v0] = true;
    *(to_visit++) = v0;
    while (to_visit != buff && nb_visited < n) {
        int v = *(--to_visit);
        int *ww = neigh[v];
        int w;
        for (int k = HASH_SIZE(deg[v]); k--; ww++) {
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


/*____________________________________________________________________________
  Not to use anymore : use graph_molloy_opt class instead

bool graph_molloy_hash::verify() {
int i;
  assert(neigh[0]==links);
  // verify edges count
  int sum = 0;
  for(i=0; i<n; i++) sum+=deg[i];
  assert(sum==a);
  // verify neigh[] and deg[] compatibility
  for(i=0; i<n-1; i++) assert(neigh[i]+HASH_SIZE(deg[i])==neigh[i+1]);
  // verify hash tables : do we see everyone ?
  for(i=0; i<n; i++) for(int j=HASH_SIZE(deg[i]); j--; )
    if(neigh[i][j]!=HASH_NONE) assert(H_is(neigh[i],deg[i],neigh[i][j]));
  degree_sequence dd(n,deg);
  graph_molloy_opt g(dd);
  dd.detach();
  int *bb = backup();
  g.restore(bb);
  delete[] bb;
  return g.verify();
}

graph_molloy_hash::graph_molloy_hash(FILE *f) {
  char *buff = new char[FBUFF_SIZE];
  // How many vertices ?
  if(VERBOSE()) fprintf(stderr,"Read file: #vertices=");
  int i;
  int n=0;
  while(fgets(buff,FBUFF_SIZE,f)) if(sscanf(buff,"%d",&i)==1 && i>n) n=i;
  n++;
  // degrees ?
  if(VERBOSE()) fprintf(stderr,"%d, #edges=",n);
  int *degs = new int[n];
  rewind(f);
  while(fgets(buff,FBUFF_SIZE,f)) {
    int d = 0;
    if(sscanf(buff,"%d",&i)==1) {
      char *b = buff;
      while(skip_int(b)) d++;
      degs[i]=d;
    }
  }
  // allocate memory
  degree_sequence dd(n,degs);
  if(VERBOSE()) fprintf(stderr,"%d\nAllocating memory...",dd.sum());
  alloc(dd);
  // add edges
  if(VERBOSE()) fprintf(stderr,"done\nCreating edges...");
  rewind(f);
  for(i=0; i<n; i++) deg[i]=0;
  int line=0;
  int j;
  while(fgets(buff,FBUFF_SIZE,f)) {
    line++;
    if(sscanf(buff,"%d",&i)==1) {
      char *b = buff;
      while(skip_int(b)) {
        if(sscanf(b,"%d",&j)!=1) {
          fprintf(stderr,"\nParse error at line %d, col %d : integer expected\n",line,int(b-buff));
          exit(6);
        }
        if(i<j) add_edge(i,j,dd.seq());
      }
    }
  }
  if(VERBOSE()) fprintf(stderr,"done\n");
  delete[] buff;
}


int graph_molloy_hash::max_degree() {
  int m=0;
  for(int k=0; k<n; k++) if(deg[k]>m) m=deg[k];
  return m;
}


bool graph_molloy_hash::havelhakimi() {

  int i;
  int dmax = max_degree()+1;
  // Sort vertices using basket-sort, in descending degrees
  int *nb = new int[dmax];
  int *sorted = new int[n];
  // init basket
  for(i=0; i<dmax; i++) nb[i]=0;
  // count basket
  for(i=0; i<n; i++) nb[deg[i]]++;
  // cumul
  int c = 0;
  for(i=dmax-1; i>=0; i--) {
    int t=nb[i];
    nb[i]=c;
    c+=t;
  }
  // sort
  for(i=0; i<n; i++) sorted[nb[deg[i]]++]=i;
  // Init edge count
  for(i=0; i<n; i++) deg[i] = 0;

// Binding process starts
  int first = 0;  // vertex with biggest residual degree
  int d = dmax-1; // maximum residual degree available

  for(c=a/2; c>0; ) {
    // pick a vertex. we could pick any, but here we pick the one with biggest degree
    int v = sorted[first];
    // look for current degree of v
    while(nb[d]<=first) d--;
    // store it in dv
    int dv = d;
    // bind it !
    c -= dv;
    int dc = d;         // residual degree of vertices we bind to
    int fc = ++first;   // position of the first vertex with degree dc

    while(dv>0 && dc>0) {
      int lc = nb[dc];
      if(lc!=fc) {
        while(dv>0 && lc>fc) {
          // binds v with sorted[--lc]
          dv--;
          int w = sorted[--lc];
          add_edge(v,w);
        }
        fc = nb[dc];
        nb[dc] = lc;
      }
      dc--;
    }
    if(dv != 0) { // We couldn't bind entirely v
      if(VERBOSE()) {
        fprintf(stderr,"Error in graph_molloy_hash::havelhakimi() :\n");
        fprintf(stderr,"Couldn't bind vertex %d entirely (%d edges remaining)\n",v,dv);
      }
      delete[] nb;
      delete[] sorted;
      return false;
    }
  }
  assert(c==0);
  delete[] nb;
  delete[] sorted;
  return true;
}


bool graph_molloy_hash::make_connected() {
  assert(verify());
  if(a/2 < n-1) {
    // fprintf(stderr,"\ngraph::make_connected() failed : #edges < #vertices-1\n");
    return false;
  }
  int i;

// Data struct for the visit :
// - buff[] contains vertices to visit
// - dist[V] is V's distance modulo 4 to the root of its comp, or -1 if it hasn't been visited yet
#define MC_BUFF_SIZE (n+2)
  int *buff = new int[MC_BUFF_SIZE];
  unsigned char * dist  = new unsigned char[n];
#define NOT_VISITED 255
#define FORBIDDEN   254
  for(i=n; i>0; dist[--i]=NOT_VISITED);

// Data struct to store components : either surplus trees or surplus edges are stored at buff[]'s end
// - A Tree is coded by one of its vertices
// - An edge (a,b) is coded by the TWO ints a and b
  int *ffub = buff+MC_BUFF_SIZE;
  edge *edges = (edge *) ffub;
  int *trees = ffub;
  int *min_ffub = buff+1+(MC_BUFF_SIZE%2 ? 0 : 1);

// There will be only one "fatty" component, and trees.
  edge fatty_edge;
  fatty_edge.from = -1;
  bool enough_edges = false;

  // start main loop
  for(int v0=0; v0<n; v0++) if(dist[v0]==NOT_VISITED) {
    // is v0 an isolated vertex?
    if(deg[v0]==0) {
#ifdef VERBOSE
      fprintf(stderr,"graph_molloy_opt::make_connected() returned FALSE : vertex %d has degree 0\n",v0);
#endif //VERBOSE
      delete[] dist;
      delete[] buff;
      return false;
    }
    dist[v0] = 0; // root
    int *to_visit = buff;
    int *current  = buff;
    *(to_visit++) = v0;

    // explore component connected to v0
    bool is_a_tree = true;
    while(current != to_visit) {
      int v = *(current++);
      unsigned char current_dist = dist[v];
      unsigned char next_dist = (current_dist+1) & 0x03;
      //unsigned char prev_dist = (current_dist-1) & 0x03;
      int* ww = neigh[v];
      int w;
      for(int k=HASH_SIZE(deg[v]); k--; ww++) if((w=*ww)!=HASH_NONE) {
        if(dist[w]==NOT_VISITED) {
          // we didn't visit w yet
          dist[w] = next_dist;
          *(to_visit++) = w;
          if(to_visit>min_ffub) min_ffub+=2; // update limit of ffub's storage
          //assert(verify());
        }
        else if(dist[w]==next_dist || (w!=HASH_NONE && w>v && dist[w]==current_dist)) {
          // we found a removable edge
          if(is_a_tree) {
            // we must first merge with the fatty component
            is_a_tree = false;
            if(fatty_edge.from < 0) {
              // we ARE the first component! fatty is us
              fatty_edge.from = v;
              fatty_edge.to   = w;
            }
            else {
              // we connect to fatty
              swap_edges(fatty_edge.from, fatty_edge.to, v, w);
              //assert(verify());
            }
          }
          else {
            // we have removable edges to give!
            if(trees!=ffub) {
              // some trees still.. Let's merge with them!
              assert(trees>=min_ffub);
              assert(edges==(edge *)ffub);
              swap_edges(v,w,*trees,neigh[*trees][0]);
              trees++;
              //assert(verify());
            }
            else if(!enough_edges) {
              // Store the removable edge for future use
              if(edges<=(edge *)min_ffub+1)
                enough_edges = true;
              else {
                edges--;
                edges->from = v;
                edges->to   = w;
              }
            }
          }
        }
      }
    }
    // Mark component
    while(to_visit!=buff) dist[*(--to_visit)] = FORBIDDEN;
    // Check if it is a tree
    if(is_a_tree ) {
      assert(deg[v0]!=0);
      if(edges!=(edge *)ffub) {
        // let's bind the tree we found with a removable edge in stock
        assert(trees == ffub);
        if(edges<(edge *)min_ffub) edges=(edge *)min_ffub;
        swap_edges(v0,neigh[v0][0],edges->from,edges->to);
        edges++;
        assert(verify());
    }
      else {
        // add the tree to the list of trees
        assert(trees>min_ffub);
        *(--trees) = v0;
        assert(verify());
      }
    }
  }
  delete[] buff;
  delete[] dist;
  return(trees == ffub);
}

int64_t graph_molloy_hash::slow_connected_shuffle(int64_t times) {
  assert(verify());
  int64_t nb_swaps = 0;
  int T = 1;

  while(times>nb_swaps) {
    // Backup graph
    int *save = backup();
    // Swaps
    int swaps = 0;
    for(int i=T; i>0; i--) {
      // Pick two random vertices a and c
      int f1 = pick_random_vertex();
      int f2 = pick_random_vertex();
      // Check that f1 != f2
      if(f1==f2) continue;
      // Get two random edges (f1,*f1t1) and (f2,*f2t2)
      int *f1t1 = random_neighbour(f1);
      int t1 = *f1t1;
      int *f2t2 = random_neighbour(f2);
      int t2 = *f2t2;
      // Check simplicity
      if(t1==t2 || f1==t2 || f2==t1) continue;
      if(is_edge(f1,t2) || is_edge(f2,t1)) continue;
      // Swap
      H_rpl(neigh[f1],deg[f1],f1t1,t2);
      H_rpl(neigh[f2],deg[f2],f2t2,t1);
      H_rpl(neigh[t1],deg[t1],f1,f2);
      H_rpl(neigh[t2],deg[t2],f2,f1);
      swaps++;
    }
    // test connectivity
    bool ok = is_connected();
    if(ok) {
      nb_swaps += swaps;
    }
    else {
      restore(save);
    }
    delete[] save;
  }
  return nb_swaps;
}


int graph_molloy_hash::width_search(unsigned char *dist, int *buff, int v0) {
  for(int i=0; i<n; i++) dist[i] = 0;
  int *to_visit = buff;
  int *to_add = buff;
  int nb_visited = 1;
  dist[v0]=1;
  *(to_add++)=v0;
  while(to_visit != to_add && nb_visited<n) {
    int v = *(to_visit++);
    int *ww = neigh[v];
    int w;
    unsigned char d = next_dist(dist[v]);
    for(int k=HASH_SIZE(deg[v]); k--; ww++) {
      if(HASH_NONE!=(w=*ww) && dist[w]==0) {
        dist[w]=d;
        nb_visited++;
        *(to_add++)=w;
      }
    }
  }
  return nb_visited;
}



int *graph_molloy_hash::vertex_betweenness_rsp(bool trivial_paths) {
  int i;
  unsigned char *dist = new unsigned char[n];
  int *buff = new int[n];
  int *b = new int[n];
  int *bb = new int[n];
  for(i=0; i<n; i++) b[i]=0;
  for(int v0 = 0; v0<n; v0++) {
    for(i=0; i<n; i++) bb[i]=0;
    int nb_vertices = width_search(dist, buff, v0);
    while(--nb_vertices) {
      int v=buff[nb_vertices];
      int d = prev_dist(dist[v]);
      int *adj = neigh[v];
      int adj_size = deg[v];
      int *ww;
      do ww=H_random(adj,adj_size); while(dist[*ww]!=d);
      if(trivial_paths || *ww!=v0) bb[*ww] += bb[v]+1;
      if(trivial_paths) bb[v]++;
    }
    for(i=0; i<n; i++) b[i]+=bb[i];
  }
  delete[] dist;
  delete[] buff;
  delete[] bb;
  return b;
}

double *graph_molloy_hash::vertex_betweenness_asp(bool trivial_paths) {
  int i;
  unsigned char *dist = new unsigned char[n];
  int *buff = new int[n];
  double *b = new double[n];
  double *bb = new double[n];
  for(i=0; i<n; i++) b[i]=0.0;
  for(int v0 = 0; v0<n; v0++) {
    for(i=0; i<n; i++) bb[i]=0.0;
    int nb_vertices = width_search(dist, buff, v0);
    if(!trivial_paths) dist[v0]=2;
    while(--nb_vertices) {
      int v=buff[nb_vertices];
      int d = prev_dist(dist[v]);
      int nb_father = 0;
      int *ww = neigh[v];
      int k;
      for(k=HASH_SIZE(deg[v]); k--; ww++) if(*ww != HASH_NONE && dist[*ww]==d) nb_father++;
      if(nb_father!=0) {
        double badd = (bb[v]+1.0)/double(nb_father);
        ww = neigh[v];
        for(k=HASH_SIZE(deg[v]); k--; ww++) if(*ww != HASH_NONE && dist[*ww]==d) bb[*ww]+=badd;
      }
      if(trivial_paths) bb[v]+=1.0;
    }
    for(i=0; i<n; i++) b[i]+=bb[i];
  }
  delete[] dist;
  delete[] buff;
  delete[] bb;
  return b;
}

//___________________________________________________________________________________
*/

} // namespace gengraph
