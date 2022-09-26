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
#include <stdexcept>

#include "gengraph_qsort.h"
#include "gengraph_degree_sequence.h"
#include "gengraph_graph_molloy_optimized.h"

#include "igraph_error.h"
#include "igraph_progress.h"


using namespace std;

namespace gengraph {

igraph_integer_t graph_molloy_opt::max_degree() {
    igraph_integer_t m = 0;
    for (igraph_integer_t k = 0; k < n; k++) if (deg[k] > m) {
            m = deg[k];
        }
    return m;
}

void graph_molloy_opt::compute_neigh() {
    igraph_integer_t *p = links;
    for (igraph_integer_t i = 0; i < n; i++) {
        neigh[i] = p;
        p += deg[i];
    }
}

void graph_molloy_opt::alloc(degree_sequence &degs) {
    n = degs.size();
    a = degs.sum();
    assert(a % 2 == 0);
    deg = new igraph_integer_t[n + a];
    for (igraph_integer_t i = 0; i < n; i++) {
        deg[i] = degs[i];
    }
    links = deg + n;
    neigh = new igraph_integer_t*[n];
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

graph_molloy_opt::graph_molloy_opt(igraph_integer_t *svg) {
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

igraph_integer_t* graph_molloy_opt::backup(igraph_integer_t *b) {
    if (b == NULL) {
        b = new igraph_integer_t[a / 2];
    }
    igraph_integer_t *c = b;
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_integer_t *p = neigh[i];
        for (igraph_integer_t d = deg[i]; d--; p++) {
            assert(*p != i);
            if (*p >= i) {
                *(c++) = *p;
            }
        }
    }
    assert(c == b + (a / 2));
    return b;
}

igraph_integer_t *graph_molloy_opt::hard_copy() {
    igraph_integer_t *hc = new igraph_integer_t[2 + n + a / 2]; // to store n,a,deg[] and links[]
    hc[0] = n;
    hc[1] = a;
    memcpy(hc + 2, deg, sizeof(igraph_integer_t)*n);
    igraph_integer_t *c = hc + 2 + n;
    for (igraph_integer_t i = 0; i < n; i++) {
        igraph_integer_t *p = neigh[i];
        for (igraph_integer_t d = deg[i]; d--; p++) {
            assert(*p != i);
            if (*p >= i) {
                *(c++) = *p;
            }
        }
    }
    assert(c == hc + 2 + n + a / 2);
    return hc;
}

void graph_molloy_opt::restore(igraph_integer_t* b) {
    igraph_integer_t i;
    for (i = 0; i < n; i++) {
        deg[i] = 0;
    }
    igraph_integer_t *p = links;
    for (i = 0; i < n - 1; i++) {
        p += deg[i];
        deg[i] = igraph_integer_t(neigh[i + 1] - neigh[i]);
        assert((neigh[i] + deg[i]) == neigh[i + 1]);
        while (p != neigh[i + 1]) {
            // b points to the current 'j'
            neigh[*b][deg[*b]++] = i;
            *(p++) = *(b++);
        }
    }
}

igraph_integer_t* graph_molloy_opt::backup_degs(igraph_integer_t *b) {
    if (b == NULL) {
        b = new igraph_integer_t[n];
    }
    memcpy(b, deg, sizeof(igraph_integer_t)*n);
    return b;
}

void graph_molloy_opt::restore_degs_only(igraph_integer_t *b) {
    memcpy(deg, b, sizeof(igraph_integer_t)*n);
    refresh_nbarcs();
}

void graph_molloy_opt::restore_degs_and_neigh(igraph_integer_t *b) {
    restore_degs_only(b);
    compute_neigh();
}

void graph_molloy_opt::restore_degs(igraph_integer_t last_degree) {
    a = last_degree;
    deg[n - 1] = last_degree;
    for (igraph_integer_t i = n - 2; i >= 0; i--) {
        a += (deg[i] = igraph_integer_t(neigh[i + 1] - neigh[i]));
    }
    refresh_nbarcs();
}

void graph_molloy_opt::clean() {
    igraph_integer_t *b = hard_copy();
    replace(b);
    delete[] b;
}

void graph_molloy_opt::replace(igraph_integer_t *_hardcopy) {
    delete[] deg;
    n = *(_hardcopy++);
    a = *(_hardcopy++);
    deg = new igraph_integer_t[a + n];
    memcpy(deg, _hardcopy, sizeof(igraph_integer_t)*n);
    links = deg + n;
    compute_neigh();
    restore(_hardcopy + n);
}

igraph_integer_t* graph_molloy_opt::components(igraph_integer_t *comp) {
    igraph_integer_t i;
    // breadth-first search buffer
    igraph_integer_t *buff = new igraph_integer_t[n];
    // comp[i] will contain the index of the component that contains vertex i
    if (comp == NULL) {
        comp = new igraph_integer_t[n];
    }
    memset(comp, 0, sizeof(igraph_integer_t)*n);
    // current component index
    igraph_integer_t curr_comp = 0;
    // loop over all non-visited vertices...
    for (igraph_integer_t v0 = 0; v0 < n; v0++) if (comp[v0] == 0) {
            curr_comp++;
            // initiate breadth-first search
            igraph_integer_t *to_visit = buff;
            igraph_integer_t *visited = buff;
            *(to_visit++) = v0;
            comp[v0] = curr_comp;
            // breadth-first search
            while (visited != to_visit) {
                igraph_integer_t v = *(visited++);
                igraph_integer_t d = deg[v];
                for (igraph_integer_t *w = neigh[v]; d--; w++) if (comp[*w] == 0) {
                        comp[*w] = curr_comp;
                        *(to_visit++) = *w;
                    }
            }
        }
    // compute component sizes and store them in buff[]
    igraph_integer_t nb_comp = 0;
    memset(buff, 0, sizeof(igraph_integer_t)*n);
    for (i = 0; i < n; i++)
        if (buff[comp[i] - 1]++ == 0 && comp[i] > nb_comp) {
            nb_comp = comp[i];
        }
    // box-sort sizes
    igraph_integer_t offset = 0;
    igraph_integer_t *box = pre_boxsort(buff, nb_comp, offset);
    for (i = nb_comp - 1; i >= 0; i--) {
        buff[i] = --box[buff[i] - offset];
    }
    delete[] box;
    // reassign component indexes
    for (igraph_integer_t *c = comp + n; comp != c--; *c = buff[*c - 1]) { }
    // clean.. at last!
    delete[] buff;
    return comp;
}

bool graph_molloy_opt::havelhakimi() {

    igraph_integer_t i;
    igraph_integer_t dmax = max_degree() + 1;
    // Sort vertices using basket-sort, in descending degrees
    igraph_integer_t *nb = new igraph_integer_t[dmax];
    igraph_integer_t *sorted = new igraph_integer_t[n];
    // init basket
    for (i = 0; i < dmax; i++) {
        nb[i] = 0;
    }
    // count basket
    for (i = 0; i < n; i++) {
        nb[deg[i]]++;
    }
    // cumul
    igraph_integer_t c = 0;
    for (i = dmax - 1; i >= 0; i--) {
        c += nb[i];
        nb[i] = -nb[i] + c;
    }
    // sort
    for (i = 0; i < n; i++) {
        sorted[nb[deg[i]]++] = i;
    }

// Binding process starts
    igraph_integer_t first = 0;  // vertex with biggest residual degree
    igraph_integer_t d = dmax - 1; // maximum residual degree available

    for (c = a / 2; c > 0; ) {
        // pick a vertex. we could pick any, but here we pick the one with biggest degree
        igraph_integer_t v = sorted[first];
        // look for current degree of v
        while (nb[d] <= first) {
            d--;
        }
        // store it in dv
        igraph_integer_t dv = d;
        // bind it !
        c -= dv;
        igraph_integer_t dc = d;         // residual degree of vertices we bind to
        igraph_integer_t fc = ++first;   // position of the first vertex with degree dc

        while (dv > 0 && dc > 0) {
            igraph_integer_t lc = nb[dc];
            if (lc != fc) {
                while (dv > 0 && lc > fc) {
                    // binds v with sorted[--lc]
                    dv--;
                    igraph_integer_t w = sorted[--lc];
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
             * level. Therefore, we use IGRAPH_FATAL(), as triggering this
             * indicates a bug. */
            IGRAPH_FATALF("Error in graph_molloy_opt::havelhakimi(): "
                          "Couldn't bind vertex %" IGRAPH_PRId " entirely (%" IGRAPH_PRId " edges remaining)",
                          v, dv);
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
    for (igraph_integer_t i = n; i > 0; visited[--i] = false) { }
    igraph_integer_t *to_visit = new igraph_integer_t[n];
    igraph_integer_t *stop = to_visit;
    igraph_integer_t left = n - 1;
    *(to_visit++) = 0;
    visited[0] = true;
    while (left > 0 && to_visit != stop) {
        igraph_integer_t v = *(--to_visit);
        igraph_integer_t *w = neigh[v];
        for (igraph_integer_t k = deg[v]; k--; w++) {
            if (!visited[*w]) {
                visited[*w] = true;
                left--;
                *(to_visit++) = *w;
            }
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
    igraph_integer_t i;

// Data struct for the visit :
// - buff[] contains vertices to visit
// - dist[V] is V's distance modulo 4 to the root of its comp, or -1 if it hasn't been visited yet
#define MC_BUFF_SIZE (n+2)
    igraph_integer_t *buff = new igraph_integer_t[MC_BUFF_SIZE];
    unsigned char * dist  = new unsigned char[n];
#define NOT_VISITED 255
#define FORBIDDEN   254
    for (i = n; i > 0; dist[--i] = NOT_VISITED) { }

// Data struct to store components : either surplus trees or surplus edges are stored at buff[]'s end
// - A Tree is coded by one of its vertices
// - An edge (a,b) is coded by the TWO ints a and b
    igraph_integer_t *ffub = buff + MC_BUFF_SIZE;
    edge *edges = (edge *) ffub;
    igraph_integer_t *trees = ffub;
    igraph_integer_t *min_ffub = buff + 1 + (MC_BUFF_SIZE % 2 ? 0 : 1);

// There will be only one "fatty" component, and trees.
    edge fatty_edge = { -1, -1 };
    bool enough_edges = false;

    // start main loop
    for (igraph_integer_t v0 = 0; v0 < n; v0++) if (dist[v0] == NOT_VISITED) {
            // is v0 an isolated vertex?
            if (deg[v0] == 0) {
                delete[] dist;
                delete[] buff;
                // 0-degree vertex found, cannot create connected graph
                return false;
            }
            dist[v0] = 0; // root
            igraph_integer_t *to_visit = buff;
            igraph_integer_t *current  = buff;
            *(to_visit++) = v0;

            // explore component connected to v0
            bool is_a_tree = true;
            while (current != to_visit) {
                igraph_integer_t v = *(current++);
                unsigned char current_dist = dist[v];
                unsigned char next_dist = (current_dist + 1) & 0x03;
                //unsigned char prev_dist = (current_dist-1) & 0x03;
                igraph_integer_t* ww = neigh[v];
                igraph_integer_t w;
                for (igraph_integer_t k = deg[v]; k--; ww++) {
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

bool graph_molloy_opt::swap_edges_simple(igraph_integer_t from1, igraph_integer_t to1, igraph_integer_t from2, igraph_integer_t to2) {
    if (from1 == to1 || from1 == from2 || from1 == to2 || to1 == from2 || to1 == to2 || from2 == to2) {
        return false;
    }
    if (is_edge(from1, to2) || is_edge(from2, to1)) {
        return false;
    }
    swap_edges(from1, to1, from2, to2);
    return true;
}

void graph_molloy_opt::print(FILE *f, bool NOZERO) {
    igraph_integer_t i, j;
    for (i = 0; i < n; i++) {
        if (!NOZERO || deg[i] > 0) {
            fprintf(f, "%" IGRAPH_PRId, i);
            for (j = 0; j < deg[i]; j++) {
                fprintf(f, " %" IGRAPH_PRId, neigh[i][j]);
            }
            fprintf(f, "\n");
        }
    }
}

igraph_integer_t graph_molloy_opt::effective_isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited) {
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

void graph_molloy_opt::depth_isolated(igraph_integer_t v, igraph_integer_t &calls, igraph_integer_t &left_to_explore, igraph_integer_t dmax, igraph_integer_t * &Kbuff, bool *visited) {
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
    igraph_integer_t *w = neigh[v];
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
}

igraph_integer_t graph_molloy_opt::depth_search(bool *visited, igraph_integer_t *buff, igraph_integer_t v0) {
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
        for (igraph_integer_t k = deg[v]; k--; ww++) if (!visited[w = *ww]) {
                visited[w] = true;
                nb_visited++;
                *(to_visit++) = w;
            }
    }
    return nb_visited;
}

igraph_integer_t graph_molloy_opt::width_search(unsigned char *dist, igraph_integer_t *buff, igraph_integer_t v0, igraph_integer_t toclear) {
    if (toclear >= 0) for (igraph_integer_t i = 0; i < toclear; i++) {
            dist[buff[i]] = 0;
        } else for (igraph_integer_t i = 0; i < n; i++) {
            dist[i] = 0;
        }
    igraph_integer_t *to_visit = buff;
    igraph_integer_t *to_add = buff;
    igraph_integer_t nb_visited = 1;
    dist[v0] = 1;
    *(to_add++) = v0;
    while (to_visit != to_add && nb_visited < n) {
        igraph_integer_t v = *(to_visit++);
        igraph_integer_t *ww = neigh[v];
        igraph_integer_t w;
        unsigned char d = next_dist(dist[v]);
        for (igraph_integer_t k = deg[v]; k--; ww++) if (dist[w = *ww] == 0) {
                dist[w] = d;
                nb_visited++;
                *(to_add++) = w;
            }
    }
    return nb_visited;
}

// dist[] MUST be full of zeros !!!!
igraph_integer_t graph_molloy_opt::breadth_path_search(igraph_integer_t src, igraph_integer_t *buff, double *paths, unsigned char *dist) {
    unsigned char last_dist = 0;
    unsigned char curr_dist = 1;
    igraph_integer_t *to_visit = buff;
    igraph_integer_t *visited  = buff;
    *(to_visit++) = src;
    paths[src] = 1.0;
    dist[src]  = curr_dist;
    igraph_integer_t nb_visited = 1;
    while (visited != to_visit) {
        igraph_integer_t v = *(visited++);
        if (last_dist == (curr_dist = dist[v])) {
            break;
        }
        unsigned char nd = next_dist(curr_dist);
        igraph_integer_t *ww = neigh[v];
        double p = paths[v];
        for (igraph_integer_t k = deg[v]; k--;) {
            igraph_integer_t w = *(ww++);
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
                    throw std::runtime_error("Fatal error: too many (>MAX_DOUBLE) possible paths in graph.");
                }
            }
        }
    }
    assert(to_visit == buff + nb_visited);
    return nb_visited;
}

igraph_integer_t *graph_molloy_opt::vertices_real(igraph_integer_t &nb_v) {
    igraph_integer_t *yo;
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
    igraph_integer_t *buff = new igraph_integer_t[nb_v];
    yo = buff;
    for (igraph_integer_t i = 0; i < n; i++) if (deg[i] > 0) {
            *(yo++) = i;
        }
    if (yo != buff + nb_v) {
        IGRAPH_WARNINGF("wrong #vertices in graph_molloy_opt::vertices_real(%" IGRAPH_PRId ")", nb_v);
        delete[] buff;
        return NULL;
    } else {
        return buff;
    }
}

bool graph_molloy_opt::isolated(igraph_integer_t v, igraph_integer_t K, igraph_integer_t *Kbuff, bool *visited) {
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
    igraph_integer_t *max   = Kbuff + (K - 1);
    *(known++) = v;
    visited[v] = true;
    bool is_isolated = true;

    while (known != seen) {
        v = *(seen++);
        igraph_integer_t *w = neigh[v];
        for (igraph_integer_t d = deg[v]; d--; w++) if (!visited[*w]) {
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

void graph_molloy_opt::sort() {
    for (int v = 0; v < n; v++) {
        qsort(neigh[v], deg[v]);
    }
}

// void graph_molloy_opt::remove_vertex(int v) {
//   fprintf(stderr,"Warning : graph_molloy_opt::remove_vertex(%d) called",v);
// }

bool graph_molloy_opt::verify(int mode) {
    IGRAPH_UNUSED(mode);
#ifndef NDEBUG
    igraph_integer_t i, j, k;
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
            igraph_integer_t v = neigh[i][j];
            igraph_integer_t nb = 0;
            for (k = 0; k < deg[v]; k++) if (neigh[v][k] == i) {
                    nb++;
                }
            assert(nb > 0);
        }
#endif
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
