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
#include <cassert>
#include <cstdio>
#include <cmath>
#include <limits>

#include "gengraph_qsort.h"
#include "gengraph_box_list.h"
#include "gengraph_vertex_cover.h"
#include "gengraph_degree_sequence.h"
#include "gengraph_graph_molloy_optimized.h"

#include "igraph_error.h"
#include "igraph_statusbar.h"
#include "igraph_progress.h"


using namespace std;

namespace gengraph {

void graph_molloy_opt::breadth_search(int *dist, int v0, int *buff) {
    bool tmpbuff = (buff == NULL);
    if (tmpbuff) {
        buff = new int[n];
    }
    for (int i = 0; i < n; i++) {
        dist[i] = -1;
    }
    dist[v0] = 0;
    int *visited = buff;
    int *to_visit = buff;
    *to_visit++ = v0;
    while (visited != to_visit) {
        int v = *visited++;
        int *w = neigh[v];
        int dd = dist[v] + 1;
        for (int d = deg[v]; d--; w++) if (dist[*w] < 0) {
                dist[*w] = dd;
                *to_visit++ = *w;
            }
    }
    if (tmpbuff) {
        delete[] buff;
    }
}


int graph_molloy_opt::max_degree() {
    int m = 0;
    for (int k = 0; k < n; k++) if (deg[k] > m) {
            m = deg[k];
        }
    return m;
}

void graph_molloy_opt::compute_neigh() {
    int *p = links;
    for (int i = 0; i < n; i++) {
        neigh[i] = p;
        p += deg[i];
    }
}

void graph_molloy_opt::alloc(degree_sequence &degs) {
    n = degs.size();
    a = degs.sum();
    assert(a % 2 == 0);
    deg = new int[n + a];
    for (int i = 0; i < n; i++) {
        deg[i] = degs[i];
    }
    links = deg + n;
    neigh = new int*[n];
    compute_neigh();
}

graph_molloy_opt::graph_molloy_opt(degree_sequence &degs) {
    alloc(degs);
}

// graph_molloy_opt::graph_molloy_opt(FILE *f) {
//   char *buff = new char[FBUFF_SIZE];
//   // How many vertices ?
//   if(VERBOSE()) fprintf(stderr,"Read file: #vertices=");
//   int i;
//   int n=0;
//   while(fgets(buff,FBUFF_SIZE,f)) if(sscanf(buff,"%d",&i)==1 && i>n) n=i;
//   n++;
//   // degrees ?
//   if(VERBOSE()) fprintf(stderr,"%d, #edges=",n);
//   int *degs = new int[n];
//   for(i=0; i<n; i++) degs[i]=0;
//   rewind(f);
//   while(fgets(buff,FBUFF_SIZE,f)) {
//     int d = 0;
//     if(sscanf(buff,"%d",&i)==1) {
//       char *b = buff;
//       while(skip_int(b)) d++;
//       degs[i]=d;
//     }
//   }
//   // allocate memory
//   degree_sequence dd(n,degs);
//   a = dd.sum();
//   if(VERBOSE()) fprintf(stderr,"%d\nAllocating memory...",a);
//   alloc(dd);
//   // add edges
//   if(VERBOSE()) fprintf(stderr,"done\nCreating edges...");
//   rewind(f);
//   int line=0;
//   int j;
//   while(fgets(buff,FBUFF_SIZE,f)) {
//     line++;
//     if(sscanf(buff,"%d",&i)==1) {
//       char *b = buff;
//       while(skip_int(b)) {
//         if(sscanf(b,"%d",&j)!=1) {
//           fprintf(stderr,"\nParse error at line %d, col %d : integer expected\n",line,int(b-buff));
//           exit(6);
//         }
//         *(neigh[i]++) = j;
//       }
//     }
//   }
//   delete[] buff;
//   compute_neigh();
//   if(VERBOSE()) fprintf(stderr,"done\n");
// }

graph_molloy_opt::graph_molloy_opt(int *svg) {
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

void graph_molloy_opt::detach() {
    deg = NULL;
    neigh = NULL;
}

graph_molloy_opt::~graph_molloy_opt() {
    if (deg != NULL) {
        delete[] deg;
    }
    if (neigh != NULL) {
        delete[] neigh;
    }
    detach();
}

int* graph_molloy_opt::backup(int *b) {
    if (b == NULL) {
        b = new int[a / 2];
    }
    int *c = b;
    for (int i = 0; i < n; i++) {
        int *p = neigh[i];
        for (int d = deg[i]; d--; p++) {
            assert(*p != i);
            if (*p >= i) {
                *(c++) = *p;
            }
        }
    }
    assert(c == b + (a / 2));
    return b;
}

int *graph_molloy_opt::hard_copy() {
    int *hc = new int[2 + n + a / 2]; // to store n,a,deg[] and links[]
    hc[0] = n;
    hc[1] = a;
    memcpy(hc + 2, deg, sizeof(int)*n);
    int *c = hc + 2 + n;
    for (int i = 0; i < n; i++) {
        int *p = neigh[i];
        for (int d = deg[i]; d--; p++) {
            assert(*p != i);
            if (*p >= i) {
                *(c++) = *p;
            }
        }
    }
    assert(c == hc + 2 + n + a / 2);
    return hc;
}

void graph_molloy_opt::restore(int* b) {
    int i;
    for (i = 0; i < n; i++) {
        deg[i] = 0;
    }
    int *p = links;
    for (i = 0; i < n - 1; i++) {
        p += deg[i];
        deg[i] = int(neigh[i + 1] - neigh[i]);
        assert((neigh[i] + deg[i]) == neigh[i + 1]);
        while (p != neigh[i + 1]) {
            // b points to the current 'j'
            neigh[*b][deg[*b]++] = i;
            *(p++) = *(b++);
        }
    }
}

int* graph_molloy_opt::backup_degs(int *b) {
    if (b == NULL) {
        b = new int[n];
    }
    memcpy(b, deg, sizeof(int)*n);
    return b;
}

void graph_molloy_opt::restore_degs_only(int *b) {
    memcpy(deg, b, sizeof(int)*n);
    refresh_nbarcs();
}

void graph_molloy_opt::restore_degs_and_neigh(int *b) {
    restore_degs_only(b);
    compute_neigh();
}

void graph_molloy_opt::restore_degs(int last_degree) {
    a = last_degree;
    deg[n - 1] = last_degree;
    for (int i = n - 2; i >= 0; i--) {
        a += (deg[i] = int(neigh[i + 1] - neigh[i]));
    }
    refresh_nbarcs();
}

void graph_molloy_opt::clean() {
    int *b = hard_copy();
    replace(b);
    delete[] b;
}

void graph_molloy_opt::replace(int *_hardcopy) {
    delete[] deg;
    n = *(_hardcopy++);
    a = *(_hardcopy++);
    deg = new int[a + n];
    memcpy(deg, _hardcopy, sizeof(int)*n);
    links = deg + n;
    compute_neigh();
    restore(_hardcopy + n);
}

int* graph_molloy_opt::components(int *comp) {
    int i;
    // breadth-first search buffer
    int *buff = new int[n];
    // comp[i] will contain the index of the component that contains vertex i
    if (comp == NULL) {
        comp = new int[n];
    }
    memset(comp, 0, sizeof(int)*n);
    // current component index
    int curr_comp = 0;
    // loop over all non-visited vertices...
    for (int v0 = 0; v0 < n; v0++) if (comp[v0] == 0) {
            curr_comp++;
            // initiate breadth-first search
            int *to_visit = buff;
            int *visited = buff;
            *(to_visit++) = v0;
            comp[v0] = curr_comp;
            // breadth-first search
            while (visited != to_visit) {
                int v = *(visited++);
                int d = deg[v];
                for (int *w = neigh[v]; d--; w++) if (comp[*w] == 0) {
                        comp[*w] = curr_comp;
                        *(to_visit++) = *w;
                    }
            }
        }
    // compute component sizes and store them in buff[]
    int nb_comp = 0;
    memset(buff, 0, sizeof(int)*n);
    for (i = 0; i < n; i++)
        if (buff[comp[i] - 1]++ == 0 && comp[i] > nb_comp) {
            nb_comp = comp[i];
        }
    // box-sort sizes
    int offset = 0;
    int *box = pre_boxsort(buff, nb_comp, offset);
    for (i = nb_comp - 1; i >= 0; i--) {
        buff[i] = --box[buff[i] - offset];
    }
    delete[] box;
    // reassign component indexes
    for (int *c = comp + n; comp != c--; *c = buff[*c - 1]) { }
    // clean.. at last!
    delete[] buff;
    return comp;
}

void graph_molloy_opt::giant_comp() {
    int *comp = components();
    // Clear edges of all vertices that do not belong to comp 0
    for (int i = 0; i < n; i++) if (comp[i] != 0) {
            deg[i] = 0;
        }
    // Clean comp[]
    delete[] comp;
}

int graph_molloy_opt::nbvertices_comp() {
    int *comp = components();
    // Count all vertices that belong to comp 0
    int nb = 0;
    for (int i = 0; i < n; i++) if (comp[i] == 0) {
            nb++;
        }
    // Clean comp[]
    delete[] comp;
    return nb;
}

int graph_molloy_opt::nbarcs_comp() {
    int *comp = components();
    // Count all vertices that belong to comp 0
    int nb = 0;
    for (int i = 0; i < n; i++) if (comp[i] == 0) {
            nb += deg[i];
        }
    // Clean comp[]
    delete[] comp;
    return nb;
}

bool graph_molloy_opt::havelhakimi() {

    int i;
    int dmax = max_degree() + 1;
    // Sort vertices using basket-sort, in descending degrees
    int *nb = new int[dmax];
    int *sorted = new int[n];
    // init basket
    for (i = 0; i < dmax; i++) {
        nb[i] = 0;
    }
    // count basket
    for (i = 0; i < n; i++) {
        nb[deg[i]]++;
    }
    // cumul
    int c = 0;
    for (i = dmax - 1; i >= 0; i--) {
        c += nb[i];
        nb[i] = -nb[i] + c;
    }
    // sort
    for (i = 0; i < n; i++) {
        sorted[nb[deg[i]]++] = i;
    }

// Binding process starts
    int first = 0;  // vertex with biggest residual degree
    int d = dmax - 1; // maximum residual degree available

    for (c = a / 2; c > 0; ) {
        // pick a vertex. we could pick any, but here we pick the one with biggest degree
        int v = sorted[first];
        // look for current degree of v
        while (nb[d] <= first) {
            d--;
        }
        // store it in dv
        int dv = d;
        // bind it !
        c -= dv;
        int dc = d;         // residual degree of vertices we bind to
        int fc = ++first;   // position of the first vertex with degree dc

        while (dv > 0 && dc > 0) {
            int lc = nb[dc];
            if (lc != fc) {
                while (dv > 0 && lc > fc) {
                    // binds v with sorted[--lc]
                    dv--;
                    int w = sorted[--lc];
                    *(neigh[v]++) = w;
                    *(neigh[w]++) = v;
                }
                fc = nb[dc];
                nb[dc] = lc;
            }
            dc--;
        }
        if (dv != 0) { // We couldn't bind entirely v
            delete[] nb;
            delete[] sorted;
            compute_neigh();

            /* Cannot use IGRAPH_ERRORF() as this function does not return
             * an error code. This situation should only occur when the degree
             * sequence is not graphical, but that is already checked at the top
             * level. Therefore, we report EINTERNAL, as triggering this
             * indicates a bug. */
            igraph_errorf("Error in graph_molloy_opt::havelhakimi(): "
                          "Couldn't bind vertex %d entirely (%d edges remaining)",
                          IGRAPH_FILE_BASENAME, __LINE__,
                          IGRAPH_EINTERNAL, v, dv);
            return false;
        }
    }
    assert(c == 0);
    compute_neigh();
    delete[] nb;
    delete[] sorted;
    return true;
}

bool graph_molloy_opt::is_connected() {
    bool *visited = new bool[n];
    for (int i = n; i > 0; visited[--i] = false) { }
    int *to_visit = new int[n];
    int *stop = to_visit;
    int left = n - 1;
    *(to_visit++) = 0;
    visited[0] = true;
    while (left > 0 && to_visit != stop) {
        int v = *(--to_visit);
        int *w = neigh[v];
        for (int k = deg[v]; k--; w++) if (!visited[*w]) {
                visited[*w] = true;
                left--;
                *(to_visit++) = *w;
            }
    }
    delete[] visited;
    delete[] stop;
    assert(left >= 0);
    return (left == 0);
}


bool graph_molloy_opt::make_connected() {
    //assert(verify());
    if (a / 2 < n - 1) {
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
    for (i = n; i > 0; dist[--i] = NOT_VISITED) { }

// Data struct to store components : either surplus trees or surplus edges are stored at buff[]'s end
// - A Tree is coded by one of its vertices
// - An edge (a,b) is coded by the TWO ints a and b
    int *ffub = buff + MC_BUFF_SIZE;
    edge *edges = (edge *) ffub;
    int *trees = ffub;
    int *min_ffub = buff + 1 + (MC_BUFF_SIZE % 2 ? 0 : 1);

// There will be only one "fatty" component, and trees.
    edge fatty_edge = { -1, -1 };
    bool enough_edges = false;

    // start main loop
    for (int v0 = 0; v0 < n; v0++) if (dist[v0] == NOT_VISITED) {
            // is v0 an isolated vertex?
            if (deg[v0] == 0) {
                delete[] dist;
                delete[] buff;
                /* Cannot use IGRAPH_ERROR() as this function does not return an error code. */
                igraph_error("Cannot create connected graph from degree sequence: "
                              "vertex with degree 0 found.",
                              IGRAPH_FILE_BASENAME, __LINE__,
                              IGRAPH_EINVAL);
                return false;
            }
            dist[v0] = 0; // root
            int *to_visit = buff;
            int *current  = buff;
            *(to_visit++) = v0;

            // explore component connected to v0
            bool is_a_tree = true;
            while (current != to_visit) {
                int v = *(current++);
                unsigned char current_dist = dist[v];
                unsigned char next_dist = (current_dist + 1) & 0x03;
                //unsigned char prev_dist = (current_dist-1) & 0x03;
                int* ww = neigh[v];
                int w;
                for (int k = deg[v]; k--; ww++) {
                    if (dist[w = *ww] == NOT_VISITED) {
                        // we didn't visit *w yet
                        dist[w] = next_dist;
                        *(to_visit++) = w;
                        if (to_visit > min_ffub) {
                            min_ffub += 2;    // update limit of ffub's storage
                        }
                        //assert(verify());
                    } else if (dist[w] == next_dist || (w >= v && dist[w] == current_dist)) {
                        // we found a removable edge
                        if (trees != ffub) {
                            // some trees still.. Let's merge with them!
                            assert(trees >= min_ffub);
                            assert(edges == (edge *)ffub);
                            swap_edges(v, w, *trees, neigh[*trees][0]);
                            trees++;
                            //assert(verify());
                        } else if (is_a_tree) {
                            // we must merge with the fatty component
                            is_a_tree = false;
                            if (fatty_edge.from < 0) {
                                // we ARE the first component! fatty is us
                                fatty_edge.from = v;
                                fatty_edge.to   = w;
                            } else {
                                // we connect to fatty
                                swap_edges(fatty_edge.from, fatty_edge.to, v, w);
                                fatty_edge.to = w;
                                //assert(verify());
                            }
                        } else if (!enough_edges) {
                            // Store the removable edge for future use
                            if (edges <= (edge *)min_ffub + 1) {
                                enough_edges = true;
                            } else {
                                edges--;
                                edges->from = v;
                                edges->to   = w;
                            }
                        }
                    }
                }
            }
            // Mark component
            while (to_visit != buff) {
                dist[*(--to_visit)] = FORBIDDEN;
            }
            // Check if it is a tree
            if (is_a_tree ) {
                assert(deg[v0] != 0);
                if (edges != (edge *)ffub) {
                    // let's bind the tree we found with a removable edge in stock
                    assert(trees == ffub);
                    if (edges < (edge *)min_ffub) {
                        edges = (edge *)min_ffub;
                    }
                    swap_edges(v0, neigh[v0][0], edges->from, edges->to);
                    edges++;
                    assert(verify());
                } else if (fatty_edge.from >= 0) {
                    // if there is a fatty component, let's merge with it ! and discard fatty :-/
                    assert(trees == ffub);
                    swap_edges(v0, neigh[v0][0], fatty_edge.from, fatty_edge.to);
                    fatty_edge.from = -1;
                    fatty_edge.to = -1;
                    assert(verify());
                } else {
                    // add the tree to the list of trees
                    assert(trees > min_ffub);
                    *(--trees) = v0;
                    assert(verify());
                }
            }
        }
    delete[] buff;
    delete[] dist;
    // Should ALWAYS return true : either we have no tree left, or we are a unique, big tree
    return (trees == ffub || ((trees + 1) == ffub && fatty_edge.from < 0));
}

bool graph_molloy_opt::swap_edges_simple(int from1, int to1, int from2, int to2) {
    if (from1 == to1 || from1 == from2 || from1 == to2 || to1 == from2 || to1 == to2 || from2 == to2) {
        return false;
    }
    if (is_edge(from1, to2) || is_edge(from2, to1)) {
        return false;
    }
    swap_edges(from1, to1, from2, to2);
    return true;
}

long graph_molloy_opt::fab_connected_shuffle(long times) {
    //assert(verify());
    long nb_swaps = 0;
    double T = double(min(a, times)) / 10.0;
    double q1 = 1.131;
    double q2 = 0.9237;

    while (times > 0) {
        long iperiod = max(1, long(T));
        // Backup graph
        int *save = backup();
        //assert(verify());
        // Swaps
        long swaps = 0;
        for (long i = iperiod; i > 0; i--) {
            // Pick two random vertices
            int f1 = links[my_random() % a];
            int f2 = links[my_random() % a];
            if (f1 == f2) {
                continue;
            }
            // Pick two random neighbours
            int *f1t1 = neigh[f1] + my_random() % deg[f1];
            int *f2t2 = neigh[f2] + my_random() % deg[f2];
            int t1 = *f1t1;
            int t2 = *f2t2;
            // test simplicity
            if (t1 != t2 && f1 != t2 && f2 != t1 && is_edge(f1, t2) && !is_edge(f2, t1)) {
                // swap
                *f1t1 = t2;
                *f2t2 = t1;
                fast_rpl(neigh[t1], f1, f2);
                fast_rpl(neigh[t2], f2, f1);
                swaps++;
            }
        }
        //assert(verify());
        // test connectivity
        if (is_connected()) {
            nb_swaps += swaps;
            times -= iperiod;
            // adjust T
            T *= q1;
        } else {
            restore(save);
            //assert(verify());
            T *= q2;
        }
        delete[] save;
    }
    return nb_swaps;
}

long graph_molloy_opt::opt_fab_connected_shuffle(long times) {
    //assert(verify());
    long nb_swaps = 0;
    double T = double(min(a, times)) / 10.0;
    double q1 = 1.131;
    double q2 = 0.9237;

    while (times > 0) {
        long iperiod = max(1, long(T));
        // Backup graph
        int *save = backup();
        //assert(verify());
        // Swaps
        long swaps = 0;
        for (long i = iperiod; i > 0; i--) {
            // Pick two random vertices
            int f1 = links[my_random() % a];
            int f2 = links[my_random() % a];
            if (f1 == f2) {
                continue;
            }
            // Pick two random neighbours
            int *f1t1 = neigh[f1] + my_random() % deg[f1];
            int *f2t2 = neigh[f2] + my_random() % deg[f2];
            int t1 = *f1t1;
            int t2 = *f2t2;
            if (
                // test simplicity
                t1 != t2 && f1 != t2 && f2 != t1 && is_edge(f1, t2) && !is_edge(f2, t1) &&
                // test isolated pair
                (deg[f1] > 1 || deg[t2] > 1) && (deg[f2] > 1 || deg[t1] > 1)
            ) {
                // swap
                *f1t1 = t2;
                *f2t2 = t1;
                fast_rpl(neigh[t1], f1, f2);
                fast_rpl(neigh[t2], f2, f1);
                swaps++;
            }
        }
        //assert(verify());
        // test connectivity
        if (is_connected()) {
            nb_swaps += swaps;
            times -= iperiod;
            // adjust T
            T *= q1;
        } else {
            restore(save);
            //assert(verify());
            T *= q2;
        }
        delete[] save;
    }
    return nb_swaps;
}

long graph_molloy_opt::gkantsidis_connected_shuffle(long times) {
    //assert(verify());
    long nb_swaps = 0;
    long T = min(a, times) / 10;

    while (times > 0) {
        // Backup graph
        int *save = backup();
        //assert(verify());
        // Swaps
        long swaps = 0;
        for (int i = T; i > 0; i--) {
            // Pick two random vertices
            int f1 = links[my_random() % a];
            int f2 = links[my_random() % a];
            if (f1 == f2) {
                continue;
            }
            // Pick two random neighbours
            int *f1t1 = neigh[f1] + my_random() % deg[f1];
            int *f2t2 = neigh[f2] + my_random() % deg[f2];
            int t1 = *f1t1;
            int t2 = *f2t2;
            // test simplicity
            if (t1 != t2 && f1 != t2 && f2 != t1 && is_edge(f1, t2) && !is_edge(f2, t1)) {
                // swap
                *f1t1 = t2;
                *f2t2 = t1;
                fast_rpl(neigh[t1], f1, f2);
                fast_rpl(neigh[t2], f2, f1);
                swaps++;
            }
        }
        //assert(verify());
        // test connectivity
        if (is_connected()) {
            nb_swaps += swaps;
            times -= T;
            // adjust T
            T++;
        } else {
            restore(save);
            //assert(verify());
            T /= 2; if (T == 0) T = 1;
        }
        delete[] save;
    }
    return nb_swaps;
}

long graph_molloy_opt::slow_connected_shuffle(long times) {
    //assert(verify());
    long nb_swaps = 0;

    while (times--) {
        // Pick two random vertices
        int f1 = links[my_random() % a];
        int f2 = links[my_random() % a];
        if (f1 == f2) {
            continue;
        }
        // Pick two random neighbours
        int *f1t1 = neigh[f1] + my_random() % deg[f1];
        int *f2t2 = neigh[f2] + my_random() % deg[f2];
        int t1 = *f1t1;
        int t2 = *f2t2;
        // test simplicity
        if (t1 != t2 && f1 != t2 && f2 != t1 && is_edge(f1, t2) && !is_edge(f2, t1)) {
            // swap
            *f1t1 = t2;
            *f2t2 = t1;
            int *t1f1 = fast_rpl(neigh[t1], f1, f2);
            int *t2f2 = fast_rpl(neigh[t2], f2, f1);
            // test connectivity
            if (is_connected()) {
                nb_swaps++;
            } else {
                // undo swap
                *t1f1 = f1; *t2f2 = f2; *f1t1 = t1; *f2t2 = t2;
            }
        }
    }
    return nb_swaps;
}

void graph_molloy_opt::print(FILE *f, bool NOZERO) {
    int i, j;
    for (i = 0; i < n; i++) {
        if (!NOZERO || deg[i] > 0) {
            fprintf(f, "%d", i);
            for (j = 0; j < deg[i]; j++) {
                fprintf(f, " %d", neigh[i][j]);
            }
            fprintf(f, "\n");
        }
    }
}

long graph_molloy_opt::effective_isolated(int v, int K, int *Kbuff, bool *visited) {
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

void graph_molloy_opt::depth_isolated(int v, long &calls, int &left_to_explore, int dmax, int * &Kbuff, bool *visited) {
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
    calls++;
    int *w = neigh[v];
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
}

int graph_molloy_opt::depth_search(bool *visited, int *buff, int v0) {
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
        for (int k = deg[v]; k--; ww++) if (!visited[w = *ww]) {
                visited[w] = true;
                nb_visited++;
                *(to_visit++) = w;
            }
    }
    return nb_visited;
}

int graph_molloy_opt::width_search(unsigned char *dist, int *buff, int v0, int toclear) {
    if (toclear >= 0) for (int i = 0; i < toclear; i++) {
            dist[buff[i]] = 0;
        } else for (int i = 0; i < n; i++) {
            dist[i] = 0;
        }
    int *to_visit = buff;
    int *to_add = buff;
    int nb_visited = 1;
    dist[v0] = 1;
    *(to_add++) = v0;
    while (to_visit != to_add && nb_visited < n) {
        int v = *(to_visit++);
        int *ww = neigh[v];
        int w;
        unsigned char d = next_dist(dist[v]);
        for (int k = deg[v]; k--; ww++) if (dist[w = *ww] == 0) {
                dist[w] = d;
                nb_visited++;
                *(to_add++) = w;
            }
    }
    return nb_visited;
}

double graph_molloy_opt::avg_dist(unsigned char *dist, int *buff, int v0, int &nb_visited, int toclear) {
    nb_visited = width_search(dist, buff, v0, toclear);
    unsigned char curr_dist = 1;
    assert(curr_dist == dist[v0]);
    double total_dist = 0.0;
    int current_dist = 0;
    for (int p = 0; p < nb_visited; p++) {
        v0 = buff[p];
        if (dist[v0] != curr_dist) {
            current_dist++;
            curr_dist = dist[v0];
        }
        total_dist += double(current_dist);
    }
    nb_visited--;
    return total_dist / double(nb_visited);
}


void graph_molloy_opt::add_traceroute_edge(int v, int k, int *newdeg, double **edge_redudancy, double red) {
    int *ww = neigh[v] + k;
    int w = *ww;
    int k2 = 0;
    // Is neigh[v][k] a new edge ?
    if (k >= newdeg[v]) {
        int *p = neigh[v] + (newdeg[v]++);
        *ww = *p;
        *p = w;
        // Now, add the dual edge
        ww = neigh[w];
        p = ww + (newdeg[w]);
        while (ww != p && *ww != v) {
            ww++;
            k2++;
        }
        if (ww == p) {
            // dual edge was not discovered.. search it and add it.
            while (*ww != v) {
                ww++;
                k2++;
            }
            *ww = *p;
            *p = v;
            newdeg[w]++;
        }
    }
    // if edge redudancy is asked, look for dual edge
    else if (edge_redudancy != NULL)
        for (int *ww = neigh[w]; * (ww++) != v; k2++) { }
    // add edge redudancy
    if (edge_redudancy != NULL) {
        edge_redudancy[v][k]  += red;
        edge_redudancy[w][k2] += red;
    }
    assert(newdeg[v] <= deg[v]);
}

// dist[] MUST be full of zeros !!!!
int graph_molloy_opt::breadth_path_search(int src, int *buff, double *paths, unsigned char *dist) {
    unsigned char last_dist = 0;
    unsigned char curr_dist = 1;
    int *to_visit = buff;
    int *visited  = buff;
    *(to_visit++) = src;
    paths[src] = 1.0;
    dist[src]  = curr_dist;
    int nb_visited = 1;
    while (visited != to_visit) {
        int v = *(visited++);
        if (last_dist == (curr_dist = dist[v])) {
            break;
        }
        unsigned char nd = next_dist(curr_dist);
        int *ww = neigh[v];
        double p = paths[v];
        for (int k = deg[v]; k--;) {
            int w = *(ww++);
            unsigned char d = dist[w];
            if (d == 0) {
                // not visited yet !
                *(to_visit++) = w;
                dist[w] = nd;
                paths[w] = p;
                // is it the last one ?
                if (++nb_visited == n) {
                    last_dist = nd;
                }
            } else if (d == nd) {
                if ((paths[w] += p) == numeric_limits<double>::infinity()) {
                    IGRAPH_ERROR("Fatal error : too many (>MAX_DOUBLE) possible"
                                 " paths in graph", IGRAPH_EOVERFLOW);
                }
            }
        }
    }
    assert(to_visit == buff + nb_visited);
    return nb_visited;
}

// dist[] MUST be full of zeros !!!!
void graph_molloy_opt::explore_usp(double *target, int nb_vertices, int *buff, double *paths, unsigned char *dist, int *newdeg, double **edge_redudancy) {

    while (--nb_vertices) {
        int v = buff[nb_vertices];
        if (target[v] > 0.0) {
            unsigned char pd = prev_dist(dist[v]);
            int *ww = neigh[v];
            int k = 0;
            // pick ONE father at random
            double father_index = my_random01() * paths[v];
            double f = 0.0;
            int father = -1;
            while (f < father_index) {
                while (dist[father = ww[k++]] != pd) { }
                f += paths[father];
            }
            // increase target[] of father
            target[father] += target[v];
            // add edge, if necessary
            if (newdeg != NULL) {
                add_traceroute_edge(v, k - 1, newdeg, edge_redudancy, target[v]);
            }
        }
        // clear dist[]
        dist[v] = 0;
    }
    dist[buff[0]] = 0;
}

// dist[] MUST be full of zeros !!!!
void graph_molloy_opt::explore_asp(double *target, int nb_vertices, int *buff, double *paths, unsigned char *dist, int *newdeg, double **edge_redudancy) {

    while (--nb_vertices) {
        int v = buff[nb_vertices];
        if (target[v] > 0.0) {
            unsigned char pd = prev_dist(dist[v]);
            int *ww = neigh[v];
            int dv = deg[v];
            double f = target[v] / paths[v];
            // pick ALL fathers
            int father;
            for (int k = 0; k < dv; k++) if (dist[father = ww[k]] == pd) {
                    // increase target[] of father
                    target[father] += paths[father] * f;
                    // add edge, if necessary
                    if (newdeg != NULL) {
                        add_traceroute_edge(v, k, newdeg, edge_redudancy, target[v]);
                    }
                }
        }
        // clear dist[]
        dist[v] = 0;
    }
    dist[buff[0]] = 0;
}

// dist[] MUST be full of zeros !!!!
void graph_molloy_opt::explore_rsp(double *target, int nb_vertices, int *buff, double *paths, unsigned char *dist, int *newdeg, double** edge_redudancy) {

    while (--nb_vertices) {
        int v = buff[nb_vertices];
        if (target[v] > 0.0) {
            unsigned char pd = prev_dist(dist[v]);
            int *ww = neigh[v];
            // for all fathers : do we take it ?
            int paths_left = int(target[v]);
            double father_index = paths[v];
            int father;
            for (int k = 0; k < deg[v]; k++) if (dist[father = ww[k]] == pd) {
                    double pf = paths[father];
                    int to_add_to_father = my_binomial(pf / father_index, paths_left);
                    father_index -= pf;
                    if (to_add_to_father > 0) {
                        paths_left -= to_add_to_father;
                        // increase target[] of father
                        target[father] += to_add_to_father;
                        // add edge, if necessary
                        if (newdeg != NULL) {
                            add_traceroute_edge(v, k, newdeg, edge_redudancy, target[v]);
                        }
                    }
                }
        }
        // clear dist[]
        dist[v] = 0;
    }
    dist[buff[0]] = 0;
}

double *graph_molloy_opt::vertex_betweenness(int mode, bool trivial_paths) {
    char MODES[3] = {'U', 'A', 'R'};
    igraph_statusf("Computing vertex betweenness %cSP...", 0, MODES[mode]);

    // breadth-first search vertex fifo
    int *buff = new int[n];
    // breadth-first search path count
    double *paths = new double[n];
    // breadth-first search distance vector
    unsigned char *dist = new unsigned char[n];
    // global betweenness
    double *b = new double[n];
    // local betweenness (for one source)
    double *target = new double[n];
    // init all
    int progress = 0;
    memset(dist, 0, sizeof(unsigned char)*n);
    for (double *yo = target + n; (yo--) != target; *yo = 1.0) { }
    for (double *yo = b + n; (yo--) != b; *yo = 0.0) { }

    int progress_steps = max(1000, n / 10);
    // Main loop
    for (int v0 = 0; v0 < n; v0++) {
        // Verbose
        if (v0 > (progress * n) / progress_steps) {
            progress++;
            igraph_progressf("Computing vertex betweenness %cSP",
                             100.0 * double(progress) / double(progress_steps), 0,
                             MODES[mode]);
        }
        // Breadth-first search
        int nb_vertices = breadth_path_search(v0, buff, paths, dist);
        // initialize target[vertices in component] to 1
        //for(int *yo = buff+nb_vertices; (yo--)!=buff; target[*yo]=1.0);
        // backwards-cumulative exploration
        switch (mode) {
        case MODE_USP:
            explore_usp(target, nb_vertices, buff, paths, dist); break;
        case MODE_ASP:
            explore_asp(target, nb_vertices, buff, paths, dist); break;
        case MODE_RSP:
            explore_rsp(target, nb_vertices, buff, paths, dist); break;
        default:
            IGRAPH_WARNING("graph_molloy_opt::vertex_betweenness() "
                           "called with Invalid Mode");
        }
        // add targets[vertices in component] to global betweenness and reset targets[]
        if (nb_vertices == n) {
            // cache optimization if all vertices are in component
            double *bb = b;
            double *tt_end = target + n;
            if (trivial_paths) for (double *yo = target; yo != tt_end; * (bb++) += *(yo++)) {}
            else {
                for (double *yo = target; yo != tt_end; * (bb++) += (*(yo++) - 1.0)) { }
                b[*buff] -= (target[*buff] - 1.0);
            }
            for (double *yo = target; yo != tt_end; * (yo++) = 1.0) { }
        } else {
            if (trivial_paths)
                for (int *yo = buff + nb_vertices; (yo--) != buff; b[*yo] += target[*yo]) { }
            else
                for (int *yo = buff + nb_vertices; (--yo) != buff; b[*yo] += (target[*yo] - 1.0)) { }
            for (int *yo = buff + nb_vertices; (yo--) != buff; target[*yo] = 1.0) { }
        }
    }
    // Clean all & return
    delete[] target;
    delete[] dist;
    delete[] buff;
    delete[] paths;
    igraph_status("Done\n", 0);
    return b;
}

double graph_molloy_opt::traceroute_sample(int mode, int nb_src, int *src, int nb_dst, int* dst, double *redudancy, double **edge_redudancy) {
    // verify & verbose
    assert(verify());
    char MODES[3] = {'U', 'A', 'R'};
    igraph_statusf("traceroute %cSP on G(N=%d,M=%d) with %d src and %d dst...",
                   0, MODES[mode], nbvertices_real(), nbarcs(), nb_src, nb_dst);

    // create dst[] buffer if necessary
    bool newdist = dst == NULL;
    if (newdist) {
        dst = new int[n];
    }
    // breadth-first search vertex fifo
    int *buff = new int[n];
    // breadth-first search path count
    double *paths = new double[n];
    // breadth-first search distance vector
    unsigned char *dist = new unsigned char[n];
    // newdeg[] allows to tag discovered edges
    int *newdeg = new int[n];
    // target[v] is > 0 if v is a destination
    double *target = new double[n];

    // init all
    int i;
    memset(dist, 0, sizeof(unsigned char)*n);
    memset(newdeg, 0, sizeof(int)*n);
    for (double *yo = target + n; (yo--) != target; *yo = 0.0) { }
    if (redudancy != NULL)
        for (double *yo = redudancy + n; (yo--) != redudancy; *yo = 0.0) { }

    // src_0 counts the number of sources having degree 0
    int src_0 = 0;
    // nopath counts the number of pairs (src,dst) having no possible path
    int nopath = 0;
    // nb_paths & total_dist are for the average distance estimator
    int nb_paths = 0;
    double total_dist = 0;
    // s will be the current source
    int s;

    while (nb_src--) if (deg[s = *(src++)] == 0) {
            src_0++;
        } else {
            // breadth-first search
            int nb_vertices = breadth_path_search(s, buff, paths, dist);
            // do we have to pick new destinations ?
            if (newdist) {
                pick_random_dst(double(nb_dst), NULL, dst);
            }
            // mark reachable destinations as "targets"
            for (i = 0; i < nb_dst; i++) {
                if (dist[dst[i]] != 0) {
                    target[dst[i]] = 1.0;
                } else {
                    nopath++;
                }
            }
            // compute avg_dist estimator
            int current_dist = 0;
            unsigned char curr_dist = 1;
            for (int p = 1; p < nb_vertices; p++) {
                int v = buff[p];
                if (dist[v] != curr_dist) {
                    curr_dist = dist[v];
                    current_dist++;
                }
                if (target[v] > 0.0) {
                    total_dist += double(current_dist);
                    nb_paths++;
                }
            }
            // substract target[] to redudancy if needed
            if (redudancy != NULL) for (i = 1; i < nb_vertices; i++) {
                    redudancy[buff[i]] -= (target[buff[i]]);
                }
            // traceroute exploration
            switch (mode) {
            case MODE_USP:
                explore_usp(target, nb_vertices, buff, paths, dist, newdeg, edge_redudancy); break;
            case MODE_ASP:
                explore_asp(target, nb_vertices, buff, paths, dist, newdeg, edge_redudancy); break;
            case MODE_RSP:
                explore_rsp(target, nb_vertices, buff, paths, dist, newdeg, edge_redudancy); break;
            default:
                IGRAPH_WARNING("graph_molloy_opt::traceroute_sample() called "
                               "with Invalid Mode");
            }
            // add target[] to redudancy[] if needed
            if (redudancy != NULL) for (i = 1; i < nb_vertices; i++) {
                    redudancy[buff[i]] += (target[buff[i]]);
                }
            // clear target[]
            for (int *yo = buff + nb_vertices; yo-- != buff; target[*yo] = 0.0) { }
        }
    // update degrees
    for (i = 0; i < n; i++) {
        deg[i] = newdeg[i];
    }
    refresh_nbarcs();
    // clean all
    delete[] buff;
    delete[] paths;
    delete[] dist;
    delete[] newdeg;
    delete[] target;
    if (newdist) {
        delete[] dst;
    }
    {
        igraph_statusf("discovered %d vertices and %d edges\n", 0,
                       nbvertices_real(), nbarcs());
        if (src_0)  igraph_warningf("%d sources had degree 0\n", IGRAPH_FILE_BASENAME,
                                        __LINE__, -1, src_0);
        if (nopath) igraph_warningf("%d (src,dst) pairs had no possible path\n",
                                        IGRAPH_FILE_BASENAME, __LINE__, -1, nopath);
    }
    return total_dist / double(nb_paths);
}

int graph_molloy_opt::disconnecting_edges() {
    int removed = 0;
    while (is_connected()) {
        // replace random edge by loops
        int i;
        do {
            i = pick_random_vertex();
        } while (i < 0 || deg[i] < 1);
        int *p = neigh[i] + (my_random() % deg[i]);
        int j = *p; *p = i;
        fast_rpl(neigh[j], i, j);
        removed++;
    }
    return removed;
}

void graph_molloy_opt::vertex_covering() {
    vertex_cover(n, links, deg, neigh);
}


// optimisations a faire :
// 1/ arreter le breadth-first search qd on a vu toutes les dst
// 2/ faire une seule redescente pour toutes les dst.

double graph_molloy_opt::path_sampling(int *nb_dst, int *dst, double* redudancies, double **edge_redudancies) {
    assert(verify());
    // do we have to store the destinations (for one src) in a temp buffer?
    bool NOMEM = (dst == NULL);
    if (NOMEM) {
        dst = new int[n];
    }
    int i;
    int next_step = n + 1;
    {
        igraph_status("Sampling paths", 0);
        next_step = 0;
    }
    // breadth-first search buffers buff[] and dist[]
    int *buff = new int[n];
    unsigned char *dist = new unsigned char[n];
    for (i = 0; i < n; i++) {
        dist[i] = 0;
    }
    // nb_pos[] counts the number of possible paths to get to a vertex
    int *nb_pos = new int[n];
    for (i = 0; i < n; i++) {
        nb_pos[i] = 0;
    }
    // newdeg[i] is the number of edges of vertex i "seen" by traceroute
    int *newdeg = new int[n];
    for (i = 0; i < n; i++) {
        newdeg[i] = 0;
    }

    // src_0 counts the number of sources having degree 0
    int src_0 = 0;
    // nopath counts the number of pairs (src,dst) having no possible path
    int nopath = 0;
    // nb_paths & total_dist are for the average distance estimator
    int nb_paths = 0;
    unsigned int total_dist = 0;
    unsigned int total_dist64 = 0;

    // s is the source of the breadth-first search
    for (int s = 0; s < n; s++) if (nb_dst[s] > 0) {
            if (deg[s] == 0) {
                src_0++;
            } else {
                if (s > next_step) {
                    next_step = s + (n / 1000) + 1;
                    igraph_progress("Sampling paths", double(s) / double(n), 0);
                }
                int v;
                // breadth-first search
                int *to_visit = buff;
                int *visited = buff;
                *(to_visit++) = s;
                dist[s] = 1;
                nb_pos[s] = 1;
                while (visited != to_visit) {
                    v = *(visited++);
                    unsigned char n_dist = next_dist(dist[v]);
                    int *w0 = neigh[v];
                    for (int *w = w0 + deg[v]; w-- != w0; ) {
                        unsigned char d2 = dist[*w];
                        if (d2 == 0) {
                            dist[*w] = d2 = n_dist;
                            *(to_visit++) = *w;
                        }
                        if (d2 == n_dist) {
                            nb_pos[*w] += nb_pos[v];
                        }
                    }
                }

                // for every target, pick a random path.
                int t_index = nb_dst[s];
                // create dst[] if necessary
                if (NOMEM) {
                    pick_random_src(double(t_index), NULL, dst);
                }
                while (t_index--) if (dist[v = *(dst++)] == 0) {
                        nopath++;
                    } else {
#ifdef DEGSEQ_VL_DEBUG
                        igraph_statusf("Sampling path %d -> %d\n", 0, s, v);
#endif // DEGSEQ_VL_DEBUG
                        nb_paths++;
                        // while we haven't reached the source..
                        while (v != s) {
                            // pick a random father
                            int index = my_random() % nb_pos[v];
                            unsigned char p_dist = prev_dist(dist[v]);
                            int *w = neigh[v];
                            int k = 0;
                            int new_father;
                            while (dist[new_father = w[k]] != p_dist || (index -= nb_pos[new_father]) >= 0) {
                                k++;
                            }
                            // add edge
                            add_traceroute_edge(v, k, newdeg, edge_redudancies, 1.0);
                            if (redudancies != NULL && new_father != s) {
                                redudancies[new_father] += 1.0;
                            }
                            // step down to father
                            v = new_father;
                            // increase total distance
                            total_dist++;
                            if (total_dist == 0) {
                                total_dist64++;
                            }
                        }
                    }
                // reset (int *)dst if necessary
                if (NOMEM) {
                    dst -= nb_dst[s];
                }

                // clear breadth-first search buffers
                while (visited != buff) {
                    v = *(--visited);
                    dist[v] = 0;
                    nb_pos[v] = 0;
                }
            }
        }
    // update degrees
    for (i = 0; i < n; i++) {
        deg[i] = newdeg[i];
    }
    refresh_nbarcs();
    // clean
    delete[] newdeg;
    delete[] buff;
    delete[] dist;
    delete[] nb_pos;
    if (NOMEM) {
        delete[] dst;
    }
    if (VERBOSE()) {
        igraph_status("Sampling paths :  Done   \n", 0);
        if (src_0)  igraph_warningf("%d sources had degree 0", IGRAPH_FILE_BASENAME,
                                        __LINE__, -1, src_0);
        if (nopath) igraph_warningf("%d (src,dst) pairs had no possible path",
                                        IGRAPH_FILE_BASENAME, __LINE__, -1, nopath);
    }
    double tdist = double(total_dist64);
    if (total_dist64 > 0) {
        tdist *= 4294967296.0;
    }
    tdist += double(total_dist);
    return tdist / double(nb_paths);
}

int *graph_molloy_opt::vertices_real(int &nb_v) {
    int *yo;
    if (nb_v < 0) {
        nb_v = 0;
        for (yo = deg; yo != deg + n; ) if (*(yo++) > 0) {
                nb_v++;
            }
    }
    if (nb_v == 0) {
        IGRAPH_WARNING("graph is empty");
        return NULL;
    }
    int *buff = new int[nb_v];
    yo = buff;
    for (int i = 0; i < n; i++) if (deg[i] > 0) {
            *(yo++) = i;
        }
    if (yo != buff + nb_v) {
        igraph_warningf("wrong #vertices in graph_molloy_opt::vertices_real(%d)",
                        IGRAPH_FILE_BASENAME, __LINE__, -1, nb_v);
        delete[] buff;
        return NULL;
    } else {
        return buff;
    }
}

int *graph_molloy_opt::pick_random_vertices(int &k, int *output, int nb_v, int *among) {
    int i;
    bool CREATED_AMONG = false;
    if (among == NULL && k > 0) {
        among = vertices_real(nb_v);
        CREATED_AMONG = true;
    }
    if (k > nb_v) {
        igraph_warningf("Warning : tried to pick %d among %d vertices. "
                        "Picked only %d", IGRAPH_FILE_BASENAME, __LINE__, -1, k, nb_v, nb_v);
        k = nb_v;
    }
    if (k > 0) {
        if (output == NULL) {
            output = new int[k];
        }
        for (i = 0; i < k; i++) {
            int tmp = i + my_random() % (nb_v - i);
            output[i] = among[tmp];
            among[tmp] = among[i];
            among[i] = output[i];
        }
    }
    if (CREATED_AMONG) {
        delete[] among;
    }
    return output;
}

int *graph_molloy_opt::pick_random_src(double k, int *nb, int* buff, int nb_v, int* among) {
    bool AMONG_CREATED = false;
    if (among == NULL || nb_v < 0) {
        AMONG_CREATED = true;
        among = vertices_real(nb_v);
    }
    int kk = int(floor(0.5 + (k >= 1.0 ? k : k * double(nb_v))));
    if (kk == 0) {
        kk = 1;
    }
    int *yo = pick_random_vertices(kk, buff, nb_v, among);
    if (nb != NULL) {
        *nb = kk;
    }
    if (AMONG_CREATED) {
        delete[] among;
    }
    return yo;
}

int *graph_molloy_opt::pick_random_dst(double k, int *nb, int* buff, int nb_v, int* among) {
    bool AMONG_CREATED = false;
    if (among == NULL || nb_v < 0) {
        AMONG_CREATED = true;
        among = vertices_real(nb_v);
    }
    int kk = int(floor(0.5 + (k > 1.0 ? k : k * double(nb_v))));
    if (kk == 0) {
        kk = 1;
    }
    int *yo = pick_random_vertices(kk, buff, nb_v, among);
    if (nb != NULL) {
        *nb = kk;
    }
    if (AMONG_CREATED) {
        delete[] among;
    }
    return yo;
}

int graph_molloy_opt::core() {
    box_list b(n, deg);
    int v;
    int removed = 0;
    while ((v = b.get_one()) >= 0) {
        b.pop_vertex(v, neigh);
        deg[v] = 0;
        removed++;
    }
    refresh_nbarcs();
    return removed;
}

int graph_molloy_opt::try_disconnect(int K, int max_tries) {
    bool *visited = new bool[n];
    for (bool *p = visited + n; p != visited; * (--p) = false) { }
    int *Kbuff = new int[K];
    int tries = 0;
    int next_step = -1;
    if (VERBOSE()) {
        next_step = 0;
    }
    bool yo = true;
    while (yo && tries < max_tries) {
        if (tries == next_step) {
            igraph_statusf("Trying to disconnect the graph... "
                           "%d edges swaps done so far", 0, tries);
            next_step += 100;
        }
        int v1 = pick_random_vertex();
        int v2 = pick_random_vertex();
        int w1 = *(random_neighbour(v1));
        int w2 = *(random_neighbour(v2));
        if (swap_edges_simple(v1, w1, v2, w2)) {
            tries++;
            yo = (!isolated(v1, K, Kbuff, visited) && !isolated(v2, K, Kbuff, visited) && !is_connected());
            swap_edges(v1, w2, v2, w1);
        }
    }
    delete[] visited;
    delete[] Kbuff;
    return tries;
}

bool graph_molloy_opt::isolated(int v, int K, int *Kbuff, bool *visited) {
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
    int *max   = Kbuff + (K - 1);
    *(known++) = v;
    visited[v] = true;
    bool is_isolated = true;

    while (known != seen) {
        v = *(seen++);
        int *w = neigh[v];
        for (int d = deg[v]; d--; w++) if (!visited[*w]) {
#ifdef OPT_ISOLATED
                if (K <= deg[*w] + 1 || known == max) {
#else //OPT_ISOLATED
                if (known == max) {
#endif //OPT_ISOLATED
                    is_isolated = false;
                    goto end_isolated;
                }
                visited[*w] = true;
                *(known++) = *w;
            }
    }
end_isolated:
    // Undo the changes to visited[]...
    while (known != Kbuff) {
        visited[*(--known)] = false;
    }
    return is_isolated;
}

double graph_molloy_opt::rho(int mode, int nb_src, int *src, int nb_dst, int *dst) {
    assert(verify());

    // create dst[] buffer if necessary
    bool newdist = dst == NULL;
    if (newdist) {
        dst = new int[n];
    }
    // breadth-first search vertex fifo
    int *buff = new int[n];
    // breadth-first search path count
    double *paths = new double[n];
    // breadth-first search distance vector
    unsigned char *dist = new unsigned char[n];
    // target[v] is > 0 if v is a destination
    double *target = new double[n];
    // times_seen count the times we saw each vertex
    int *times_seen = new int[n];

    // init all
    int i;
    memset(dist, 0, sizeof(unsigned char)*n);
    memset(times_seen, 0, sizeof(int)*n);
    for (double *yo = target + n; (yo--) != target; *yo = 0.0) { }

    // src_0 counts the number of sources having degree 0
    int src_0 = 0;
    // nopath counts the number of pairs (src,dst) having no possible path
    int nopath = 0;
    // s will be the current source
    int s;

    for (int nsrc = 0; nsrc < nb_src; nsrc++) if (deg[s = *(src++)] == 0) {
            src_0++;
        } else {
            // breadth-first search
            int nb_vertices = breadth_path_search(s, buff, paths, dist);
            // do we have to pick new destinations ?
            if (newdist) {
                pick_random_dst(double(nb_dst), NULL, dst);
            }
            // mark reachable destinations as "targets" and substract one time_seen
            for (i = 0; i < nb_dst; i++) {
                if (dist[dst[i]] != 0) {
                    target[dst[i]] = 1.0;
                } else {
                    nopath++;
                }
            }
            // traceroute exploration
            switch (mode) {
            case MODE_USP:
                explore_usp(target, nb_vertices, buff, paths, dist); break;
            case MODE_ASP:
                explore_asp(target, nb_vertices, buff, paths, dist); break;
            case MODE_RSP:
                explore_rsp(target, nb_vertices, buff, paths, dist); break;
            default:
                IGRAPH_WARNING("graph_molloy_opt::rho() called with Invalid Mode");
            }
            // remove destinations that weren't discovered by a path coming through
            for (i = 0; i < nb_dst; i++) {
                int yo = dst[i];
                if (target[yo] == 1.0) {
                    target[yo] = 0.0;
                }
            }
            // add target[] to times_seen[]
            for (i = 1; i < nb_vertices; i++) {
                int yo = buff[i];
                if (target[yo] != 0.0) {
                    target[yo] = 0.0;
                    times_seen[yo]++;
                }
            }
            // also clear  the source
            target[buff[0]] = 0.0;
        }
    // clean all
    delete[] buff;
    delete[] paths;
    delete[] dist;
    delete[] target;
    if (newdist) {
        delete[] dst;
    }
    // compute rho
    double sum_nij = 0.0;
    double sum_ni = 0.0;
    for (i = 0; i < n; i++) {
        double d = double(times_seen[i]);
        sum_ni += d;
        sum_nij += d * d;
    }
    delete[] times_seen;
    {
        igraph_status("done\n", 0);
        if (src_0)  igraph_warningf("%d sources had degree 0", IGRAPH_FILE_BASENAME, __LINE__,
                                        -1, src_0);
        if (nopath) igraph_warningf("%d (src,dst) pairs had no possible path",
                                        IGRAPH_FILE_BASENAME, __LINE__, -1, nopath);
    }
    return (sum_nij - sum_ni) * double(n) * double(nb_src) / (sum_ni * sum_ni * double(nb_src - 1));
}

void graph_molloy_opt::sort() {
    for (int v = 0; v < n; v++) {
        qsort(neigh[v], deg[v]);
    }
}

int* graph_molloy_opt::sort_vertices(int *buff) {
    // pre-sort vertices by degrees
    buff = boxsort(deg, n, buff);
    // sort vertices having the same degrees
    int i = 0;
    while (i < n) {
        int d = deg[buff[i]];
        int j = i + 1;
        while (j < n && deg[buff[j]] == d) {
            j++;
        }
        lex_qsort(neigh, buff + i, j - i, d);
        i = j;
    }
    return buff;
}

int graph_molloy_opt::cycles(int v) {
    return v;
}

// void graph_molloy_opt::remove_vertex(int v) {
//   fprintf(stderr,"Warning : graph_molloy_opt::remove_vertex(%d) called",v);
// }

bool graph_molloy_opt::verify(int mode) {
    int i, j, k;
    assert(neigh[0] == links);
    // verify edges count
    if ((mode & VERIFY_NOARCS) == 0) {
        int sum = 0;
        for (i = 0; i < n; i++) {
            sum += deg[i];
        }
        assert(sum == a);
    }
    // verify neigh[] and deg[] compatibility
    if ((mode & VERIFY_NONEIGH) == 0)
        for (i = 0; i < n - 1; i++) {
            assert(neigh[i] + deg[i] == neigh[i + 1]);
        }
    // verify vertex range
    for (i = 0; i < a; i++) {
        assert(links[i] >= 0 && links[i] < n);
    }
    // verify simplicity
//  for(i=0; i<n; i++) for(j=0; j<deg[i]; j++) for(k=j+1; k<deg[i]; k++)
//    assert(neigh[i][j]!=neigh[i][k]);
    // verify symmetry
    for (i = 0; i < n; i++) for (j = 0; j < deg[i]; j++) {
            int v = neigh[i][j];
            int nb = 0;
            for (k = 0; k < deg[v]; k++) if (neigh[v][k] == i) {
                    nb++;
                }
            assert(nb > 0);
        }
    return true;
}

/*___________________________________________________________________________________
  Not to use anymore : use graph_molloy_hash class instead

void graph_molloy_opt::shuffle(long times) {
  while(times) {
    int f1 = links[my_random()%a];
    int f2 = links[my_random()%a];
    int t1 = neigh[f1][my_random()%deg[f1]];
    int t2 = neigh[f2][my_random()%deg[f2]];
    if(swap_edges_simple(f1,t1,f2,t2)) times--;
  }
}


long graph_molloy_opt::connected_shuffle(long times) {
  //assert(verify());
#ifdef PERFORMANCE_MONITOR
  long failures = 0;
  long successes = 0;
  double avg_K = 0.0;
  long avg_T = 0;
#endif //PERFORMANCE_MONITOR

  long nb_swaps = 0;
  long T = min(a,times)/10;
  double double_K = 1.0;
  int K = int(double_K);
  double Q1 = 1.35;
  double Q2 = 1.01;
  int *Kbuff = new int[K];
  bool *visited = new bool[n];
  for(int i=0; i<n; i++) visited[i] = false;

  while(times>nb_swaps) {
    // Backup graph
#ifdef PERFORMANCE_MONITOR
    avg_K+=double_K;
    avg_T+=T;
#endif //PERFORMANCE_MONITOR
    int *save = backup();
    //assert(verify());
    // Swaps
    long swaps = 0;
    for(int i=T; i>0; i--) {
      // Pick two random vertices
      int f1 = pick_random_vertex();
      int f2 = pick_random_vertex();
      if(f1==f2) continue;
      // Pick two random neighbours
      int *f1t1 = random_neighbour(f1);
      int t1 = *f1t1;
      int *f2t2 = random_neighbour(f2);
      int t2 = *f2t2;
      // test simplicity
      if(t1!=t2 && f1!=t2 && f2!=t1 && !is_edge(f1,t2) && !is_edge(f2,t1)) {
        // swap
        *f1t1 = t2;
        *f2t2 = t1;
        int *t1f1 = fast_rpl(neigh[t1],f1,f2);
        int *t2f2 = fast_rpl(neigh[t2],f2,f1);
        // isolation test
        if(isolated(f1, K, Kbuff, visited) || isolated(f2, K, Kbuff, visited)) {
          // undo swap
          *t1f1 = f1; *t2f2 = f2; *f1t1 = t1; *f2t2 = t2;
        }
        else swaps++;
      }
    }
    //assert(verify());
    // test connectivity
    bool ok = is_connected();
#ifdef PERFORMANCE_MONITOR
    if(ok) successes++; else failures++;
#endif //PERFORMANCE_MONITOR
    if(ok) {
      nb_swaps += swaps;
      // adjust K and T
      if((K+10)*T>5*a) {
        double_K/=Q2;
        K = int(double_K);
      }
      else T*=2;
    }
    else {
      restore(save);
      //assert(verify());
      double_K*=Q1;
      K = int(double_K);
      delete[] Kbuff;
      Kbuff = new int[K];
    }
    delete[] save;
  }
#ifdef PERFORMANCE_MONITOR
    fprintf(stderr,"\n*** Performance Monitor ***\n");
    fprintf(stderr," - Connectivity test successes : %ld\n",successes);
    fprintf(stderr," - Connectivity test failures  : %ld\n",failures);
    fprintf(stderr," - Average window : %ld\n",avg_T/long(successes+failures));
    fprintf(stderr," - Average isolation test width : %f\n",avg_K/double(successes+failures));
#endif //PERFORMANCE_MONITOR
  return nb_swaps;
}

bool graph_molloy_opt::try_shuffle(int T, int K) {
    int i;
    int *Kbuff = NULL;
    if(K>0) Kbuff = new int[K];
    bool *visited = new bool[n];
    for(i=0; i<n; i++) visited[i]=false;
    int *back=backup();
    for(i=T; i>0; i--) {
      // Pick two random vertices
      int f1 = pick_random_vertex();
      int f2 = pick_random_vertex();
      if(f1==f2) continue;
      // Pick two random neighbours
      int *f1t1 = random_neighbour(f1);
      int t1 = *f1t1;
      int *f2t2 = random_neighbour(f2);
      int t2 = *f2t2;
      // test simplicity
      if(t1!=t2 && f1!=t2 && f2!=t1 && is_edge(f1,t2) && !is_edge(f2,t1)) {
        // swap
        *f1t1 = t2;
        *f2t2 = t1;
        int *t1f1 = fast_rpl(neigh[t1],f1,f2);
        int *t2f2 = fast_rpl(neigh[t2],f2,f1);
        // isolation test
        if(isolated(f1, K, Kbuff, visited) || isolated(f2, K, Kbuff, visited)) {
          // undo swap
          *t1f1 = f1; *t2f2 = f2; *f1t1 = t1; *f2t2 = t2;
        }
      }
    }
    delete[] visited;
    if(Kbuff != NULL) delete[] Kbuff;
    bool yo = is_connected();
    restore(back);
    delete[] back;
    return yo;
}

double graph_molloy_opt::window(int K, double ratio) {
  int steps = 100;
  double T = double(a*10);
  double q2 = 0.1;
  double q1 = pow(q2,(ratio-1.0)/ratio);

  int failures = 0;
  int successes = 0;
  int *Kbuff = new int[K];
  bool *visited = new bool[n];

  while(successes<10*steps) {
    int *back=backup();
    for(int i=int(T); i>0; i--) {
      // Pick two random vertices
      int f1 = links[my_random()%a];
      int f2 = links[my_random()%a];
      if(f1==f2) continue;
      // Pick two random neighbours
      int *f1t1 = neigh[f1]+my_random()%deg[f1];
      int *f2t2 = neigh[f2]+my_random()%deg[f2];
      int t1 = *f1t1;
      int t2 = *f2t2;
      // test simplicity
      if(t1!=t2 && f1!=t2 && f2!=t1 && is_edge(f1,t2) && !is_edge(f2,t1)) {
        // swap
        *f1t1 = t2;
        *f2t2 = t1;
        int *t1f1 = fast_rpl(neigh[t1],f1,f2);
        int *t2f2 = fast_rpl(neigh[t2],f2,f1);
        // isolation test
        if(isolated(f1, K, Kbuff, visited) || isolated(f2, K, Kbuff, visited)) {
          // undo swap
          *t1f1 = f1; *t2f2 = f2; *f1t1 = t1; *f2t2 = t2;
        }
      }
    }
    if(is_connected()) {
      T *= q1;
      if(T>double(5*a)) T=double(5*a);
      successes++;
      if((successes%steps)==0) {
        q2 = sqrt(q2);
        q1 = sqrt(q1);
      }
    }
    else {
      T*=q2;
      failures++;
    }
    if(VERBOSE()) fprintf(stderr,".");
    restore(back);
    delete[] back;
  }
  delete[] Kbuff;
  delete[] visited;
  if(VERBOSE()) fprintf(stderr,"Failures:%d   Successes:%d\n",failures, successes);
  return T;
}


double graph_molloy_opt::eval_K(int quality) {
  double K = 5.0;
  double avg_K = 1.0;
  for(int i=quality; i--; ) {
    int int_K = int(floor(K+0.5));
    if(try_shuffle(a/(int_K+1),int_K)) {
      K*=0.8; fprintf(stderr,"+"); }
    else {
      K*=1.25; fprintf(stderr,"-"); }
    if(i<quality/2) avg_K *= K;
  }
  return pow(avg_K,1.0/double(quality/2));
}


double graph_molloy_opt::effective_K(int K, int quality) {
  if(K<3) return 0.0;
  long sum_K = 0;
  int *Kbuff = new int[K];
  bool *visited = new bool[n];
  int i;
  for(i=0; i<n; i++) visited[i] = false;
  for(int i=0; i<quality; i++) {
//    assert(verify());
    int f1,f2,t1,t2;
    int *f1t1, *f2t2;
    do {
      // Pick two random vertices
      do {
        f1 = pick_random_vertex();
        f2 = pick_random_vertex();
      } while(f1==f2);
      // Pick two random neighbours
      f1t1 = random_neighbour(f1);
      t1 = *f1t1;
      f2t2 = random_neighbour(f2);
      t2 = *f2t2;
      // test simplicity
    }
    while (t1==t2 || f1==t2 || f2==t1 || is_edge(f1,t2) || is_edge(f2,t1));
    // swap
    *f1t1 = t2;
    *f2t2 = t1;
    fast_rpl(neigh[t1],f1,f2);
    fast_rpl(neigh[t2],f2,f1);
    sum_K += effective_isolated(deg[f1]>deg[t2] ? f1 : t2, K, Kbuff, visited);
    sum_K += effective_isolated(deg[f2]>deg[t1] ? f2 : t1, K, Kbuff, visited);
    // undo swap
    swap_edges(f1,t2,f2,t1);
//    assert(verify());
  }
  delete[] Kbuff;
  delete[] visited;
  return double(sum_K)/double(2*quality);
}


//___________________________________________________________________________________
*/



/***** NOT USED ANYMORE (Modif 22/04/2005) ******

int64_t *graph_molloy_opt::vertex_betweenness_usp(bool trivial_paths) {
  if(VERBOSE()) fprintf(stderr,"Computing vertex betweenness USP...");
  int i;
  unsigned char *dist = new unsigned char[n];
  int *buff = new int[n];
  int64_t *b = new int64_t[n];
  int *bb = new int[n];
  int *dd = new int[max_degree()];
  for(i=0; i<n; i++) b[i]=0;
  int progress = 0;
  for(int v0 = 0; v0<n; v0++) {
    if(VERBOSE()==VERBOSE_LOTS && v0>(progress*n)/1000) {
      progress++;
      fprintf(stderr,"\rComputing vertex betweenness USP : %d.%d%% ",progress/10,progress%10);
    }
    int nb_vertices = width_search(dist, buff, v0);
    int nv = nb_vertices;
    for(i=0; i<nv; i++) bb[buff[i]]=0;
    while(--nv) {
      int v = buff[nv];
      unsigned char d = prev_dist(dist[v]);
      int n_father = 0;
      int *ww = neigh[v];
      for(int k=deg[v]; k--; ww++) if(dist[*ww]==d) dd[n_father++]=*ww;
      int w = dd[my_random()%n_father];
      if(trivial_paths || w!=v0) bb[w] += bb[v]+1;
      if(trivial_paths) bb[v]++;
    }
    for(i=0; i<nb_vertices; i++) b[buff[i]]+=(int64_t)(bb[buff[i]]);
  }
  delete[] dist;
  delete[] buff;
  delete[] bb;
  delete[] dd;
  return b;
}

int64_t *graph_molloy_opt::vertex_betweenness_rsp(bool trivial_paths) {
  if(VERBOSE()) fprintf(stderr,"Computing vertex betweenness RSP...");
  int i;
  unsigned char *dist = new unsigned char[n];
  int *buff = new int[n];
  int64_t *b = new int64_t[n];
  int *bb = new int[n];
  int *dd = new int[max_degree()];
  for(i=0; i<n; i++) b[i]=0;
  int progress = 0;
  for(int v0 = 0; v0<n; v0++) {
    if(VERBOSE()==VERBOSE_LOTS && v0>(progress*n)/1000) {
      progress++;
      fprintf(stderr,"\rComputing vertex betweenness RSP : %d.%d%% ",progress/10,progress%10);
    }
    int nb_vertices = width_search(dist, buff, v0);
    int nv = nb_vertices;
    for(i=0; i<nv; i++) bb[buff[i]]=0;
    while(--nv) {
      int v = buff[nv];
      unsigned char d = prev_dist(dist[v]);
      int n_father = 0;
      int *ww = neigh[v];
      for(int k=deg[v]; k--; ww++) if(dist[*ww]==d) dd[n_father++]=*ww;
      int to_give = bb[v]+1;
      if(dd[0]==v0) {
        if(trivial_paths) bb[v0]+= to_give;
      }
      else  {
        while(n_father>1 && to_give>2*n_father) {
          int o = rng.binomial(1.0/n_father,to_give);
          to_give -= o;
          bb[dd[--n_father]]+=o;
        }
        if(n_father==1) bb[dd[0]]+=to_give;
        else {
          while(to_give--) bb[dd[my_random()%n_father]]++;
        }
      }
      if(trivial_paths) bb[v]++;
    }
    for(i=0; i<nb_vertices; i++) b[buff[i]]+=(int64_t)(bb[buff[i]]);
  }
  delete[] dist;
  delete[] buff;
  delete[] bb;
  delete[] dd;
  return b;
}

double *graph_molloy_opt::vertex_betweenness_asp(bool trivial_paths) {
  if(VERBOSE()) fprintf(stderr,"Computing vertex betweenness ASP...");
  int i;
  unsigned char *dist = new unsigned char[n];
  int *buff = new int[n];
  double *b = new double[n];
  double *bb = new double[n];
  int *dd = new int[max_degree()];
  for(i=0; i<n; i++) b[i]=0.0;
  int progress = 0;
  for(int v0 = 0; v0<n; v0++) if(deg[v0]>0) {
    if(VERBOSE()==VERBOSE_LOTS && v0>(progress*n)/1000) {
      progress++;
      fprintf(stderr,"\rComputing vertex betweenness ASP : %d.%d%% ",progress/10,progress%10);
    }
    int nb_vertices = width_search(dist, buff, v0);
    if(!trivial_paths) dist[v0]=2;
    int nv = nb_vertices;
    for(i=0; i<nv; i++) bb[buff[i]]=0.0;
    while(--nv) {
      int v = buff[nv];
      unsigned char d = prev_dist(dist[v]);
      int n_father = 0;
      int *ww = neigh[v];
      for(int k=deg[v]; k--; ww++) if(dist[*ww]==d) dd[n_father++]=*ww;
      if(n_father!=0) {
        double badd = (bb[v]+1.0)/double(n_father);
        int *d2 = dd;
        while(n_father--) bb[*(d2++)]+=badd;
      }
      if(trivial_paths) bb[v]+=1.0;
    }
    for(i=0; i<nb_vertices; i++) b[buff[i]]+=bb[buff[i]];
  }
  delete[] dist;
  delete[] buff;
  delete[] bb;
  delete[] dd;
  if(VERBOSE()) fprintf(stderr,"done\n");
  return b;
}

*/

} // namespace gengraph
