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
#include "gengraph_random.h"
#include "gengraph_powerlaw.h"
#include "gengraph_degree_sequence.h"
#include "gengraph_hash.h"

#include "igraph_statusbar.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include <stdexcept>

// using namespace __gnu_cxx;
using namespace std;

namespace gengraph {

// shuffle an int[] randomly
void random_permute(int *a, int n);

// sort an array of positive integers in time & place O(n + max)
void cumul_sort(int *q, int n);


void degree_sequence::detach() {
    deg = NULL;
}

degree_sequence::~degree_sequence() {
    if (deg != NULL) {
        delete[] deg;
    }
    deg = NULL;
}

void degree_sequence::make_even(int mini, int maxi) {
    if (total % 2 == 0) {
        return;
    }
    if (maxi < 0) {
        maxi = 0x7FFFFFFF;
    }
    int i;
    for (i = 0; i < n; i++) {
        if (deg[i] > mini) {
            deg[i]--;
            total--;
            break;
        } else if (deg[i] < maxi) {
            deg[i]++;
            total++;
            break;
        }
    }
    if (i == n) {
        IGRAPH_WARNING("Warning: degree_sequence::make_even() forced one "
                       "degree to go over degmax");
        deg[0]++;
        total++;
    }
}

void degree_sequence::shuffle() {
    random_permute(deg, n);
}

void degree_sequence::sort() {
    cumul_sort(deg, n);
}

void degree_sequence::compute_total() {
    total = 0;
    for (int i = 0; i < n; i++) {
        total += deg[i];
    }
}

degree_sequence::
degree_sequence(int n0, int *degs) {
    deg = degs;
    n = n0;
    compute_total();
}

degree_sequence::
degree_sequence(const igraph_vector_t *out_seq) {
    n = igraph_vector_size(out_seq);
    deg = new int[n];
    for (long int i = 0; i < n; i++) {
        deg[i] = VECTOR(*out_seq)[i];
    }
    compute_total();
}

#ifndef FBUFF_SIZE
    #define FBUFF_SIZE 999
#endif //FBUFF_SIZE

// degree_sequence::degree_sequence(FILE *f, bool DISTRIB) {
//   n = 0;
//   total = 0;
//   char *buff = new char[FBUFF_SIZE];
//   char *c;
//   vector<int> degree;
//   if(!DISTRIB) {
//     // Input is a 'raw' degree sequence d0 d1 d2 d3 ...
//     while(fgets(buff, FBUFF_SIZE, f)) {
//       int d = strtol(buff, &c, 10);
//       if(c == buff) continue;
//       degree.push_back(d);
//       total += d;
//     }
//     n = int(degree.size());
//     deg = new int[n];
//     int *yo = deg;
//     vector<int>::iterator end = degree.end();
//     for(vector<int>::iterator it=degree.begin(); it!=end; *(yo++) = *(it++));
//   }
//   else {
//     // Input is a degree distribution : d0 #(degree=d0), d1 #(degree=d1), ...
//     vector<int> n_with_degree;
//     int line = 0;
//     int syntax  = 0;
//     int ignored = 0;
//     int first_syntax  = 0;
//     int first_ignored = 0;
//     while(fgets(buff, FBUFF_SIZE, f)) {
//       line++;
//       int d = strtol(buff, &c, 10);
//       if(c == buff) { ignored++; first_ignored = line; continue; }
//       char *cc;
//       int i = strtol(c, &cc, 10);
//       if(cc == c) { syntax++; first_syntax = line; continue; }
//       n += i;
//       total += i*d;
//       degree.push_back(d);
//       n_with_degree.push_back(i);
//       if( cc != c) {  syntax++; first_syntax = line; }
//     }
//     if(VERBOSE()) {
//       if(ignored > 0) fprintf(stderr,"Ignored %d lines (first was line #%d)\n", ignored, first_ignored);
//       if(syntax > 0) fprintf(stderr,"Found %d probable syntax errors (first was line #%d)\n", syntax, first_syntax);
//     }
//     deg = new int[n];
//     int *yo = deg;
//     vector<int>::iterator it_n = n_with_degree.begin();
//     for(vector<int>::iterator it = degree.begin(); it != degree.end(); it++)
//       for(int k = *(it_n++); k--; *yo++ = *it);
//   }
//   if(VERBOSE()) {
//     if(total % 2 != 0) fprintf(stderr,"Warning: degree sequence is odd\n");
//     fprintf(stderr,"Degree sequence created. N=%d, 2M=%d\n", n, total);
//   }
// }

// n vertices, exponent, min degree, max degree, average degree (optional, default is -1)
degree_sequence::
degree_sequence(int _n, double exp, int degmin, int degmax, double z) {

    n = _n;
    if (exp == 0.0) {
        // Binomial distribution
        if (z < 0) {
            throw std::invalid_argument(
                        "Fatal error in degree_sequence constructor: "
                        "positive average degree must be specified.");
        }
        if (degmax < 0) {
            degmax = n - 1;
        }
        total = int(floor(double(n) * z + 0.5));
        deg = new int[n];
        KW_RNG::RNG myrand;
        double p = (z - double(degmin)) / double(n);
        total = 0;
        for (int i = 0; i < n; i++) {
            do {
                deg[i] = 1 + myrand.binomial(p, n);
            } while (deg[i] > degmax);
            total += deg[i];
        }
    } else {
        // Power-law distribution
        igraph_status("Creating powerlaw sampler...", 0);
        powerlaw pw(exp, degmin, degmax);
        if (z == -1.0) {
            pw.init();
            igraph_statusf("done. Mean=%f\n", 0, pw.mean());
        } else {
            double offset = pw.init_to_mean(z);
            igraph_statusf("done. Offset=%f, Mean=%f\n", 0, offset, pw.mean());
        }

        deg = new int[n];
        total = 0;
        int i;

        igraph_statusf("Sampling %d random numbers...", 0, n);
        for (i = 0; i < n; i++) {
            deg[i] = pw.sample();
            total += deg[i];
        }

        igraph_status("done\nSimple statistics on degrees...", 0);
        int wanted_total = int(floor(z * n + 0.5));
        sort();
        igraph_statusf("done : Max=%d, Total=%d.\n", 0, deg[0], total);
        if (z != -1.0)  {
            igraph_statusf("Adjusting total to %d...", 0, wanted_total);
            int iterations = 0;

            while (total != wanted_total) {
                sort();
                for (i = 0; i < n && total > wanted_total; i++) {
                    total -= deg[i];
                    if (total + degmin <= wanted_total) {
                        deg[i] = wanted_total - total;
                    } else {
                        deg[i] = pw.sample();
                    }
                    total += deg[i];
                }
                iterations += i;
                for (i = n - 1; i > 0 && total < wanted_total; i--) {
                    total -= deg[i];
                    if (total + (deg[0] >> 1) >= wanted_total) {
                        deg[i] = wanted_total - total;
                    } else {
                        deg[i] = pw.sample();
                    }
                    total += deg[i];
                }
                iterations += n - 1 - i;
            }
            igraph_statusf("done(%d iterations).", 0, iterations);
            igraph_statusf("  Now, degmax = %d\n", 0, dmax());
        }

        shuffle();
    }
}

// void degree_sequence::print() {
//   for(int i=0; i<n; i++) printf("%d\n",deg[i]);
// }

// void degree_sequence::print_cumul() {
//   if(n==0) return;
//   int dmax = deg[0];
//   int dmin = deg[0];
//   int i;
//   for(i=1; i<n; i++) if(dmax<deg[i]) dmax=deg[i];
//   for(i=1; i<n; i++) if(dmin>deg[i]) dmin=deg[i];
//   int *dd = new int[dmax-dmin+1];
//   for(i=dmin; i<=dmax; i++) dd[i-dmin]=0;
//   if(VERBOSE()) fprintf(stderr,"Computing cumulative distribution...");
//   for(i=0; i<n; i++) dd[deg[i]-dmin]++;
//   if(VERBOSE()) fprintf(stderr,"done\n");
//   for(i=dmin; i<=dmax; i++) if(dd[i-dmin]>0) printf("%d %d\n",i,dd[i-dmin]);
//   delete[] dd;
// }

bool degree_sequence::havelhakimi() {

    int i;
    int dm = dmax() + 1;
    // Sort vertices using basket-sort, in descending degrees
    int *nb = new int[dm];
    int *sorted = new int[n];
    // init basket
    for (i = 0; i < dm; i++) {
        nb[i] = 0;
    }
    // count basket
    for (i = 0; i < n; i++) {
        nb[deg[i]]++;
    }
    // cumul
    int c = 0;
    for (i = dm - 1; i >= 0; i--) {
        int t = nb[i];
        nb[i] = c;
        c += t;
    }
    // sort
    for (i = 0; i < n; i++) {
        sorted[nb[deg[i]]++] = i;
    }

// Binding process starts
    int first = 0;  // vertex with biggest residual degree
    int d = dm - 1; // maximum residual degree available

    for (c = total / 2; c > 0; ) {
        // We design by 'v' the vertex of highest degree (indexed by first)
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
                    lc--;
                }
                fc = nb[dc];
                nb[dc] = lc;
            }
            dc--;
        }
        if (dv != 0) { // We couldn't bind entirely v
            delete[] nb;
            delete[] sorted;
            return false;
        }
    }
    delete[] nb;
    delete[] sorted;
    return true;
}

//*************************
// Subroutines definitions
//*************************

inline int int_adjust(double x) {
    return (int(floor(x + random_float())));
}

void random_permute(int *a, int n) {
    int j, tmp;
    for (int i = 0; i < n - 1; i++) {
        j = i + my_random() % (n - i);
        tmp = a[i];
        a[i] = a[j];
        a[j] = tmp;
    }
}

void cumul_sort(int *q, int n) {
    // looks for the maximum q[i] and minimum
    if (n == 0) {
        return;
    }
    int qmax = q[0];
    int qmin = q[0];
    int i;
    for (i = 0; i < n; i++) if (q[i] > qmax) {
            qmax = q[i];
        }
    for (i = 0; i < n; i++) if (q[i] < qmin) {
            qmin = q[i];
        }

    // counts #q[i] with given q
    int *nb = new int[qmax - qmin + 1];
    for (int *onk = nb + (qmax - qmin + 1); onk != nb; * (--onk) = 0) { }
    for (i = 0; i < n; i++) {
        nb[q[i] - qmin]++;
    }

    // counts cumulative distribution
    for (i = qmax - qmin; i > 0; i--) {
        nb[i - 1] += nb[i];
    }

    // sort by q[i]
    int last_q;
    int tmp;
    int modifier = qmax - qmin + 1;
    for (int current = 0; current < n; current++) {
        tmp = q[current];
        if (tmp >= qmin && tmp <= qmax) {
            last_q = qmin;
            do {
                q[current] = last_q + modifier;
                last_q = tmp;
                current = --nb[last_q - qmin];
            } while ((tmp = q[current]) >= qmin && tmp <= qmax);
            q[current] = last_q + modifier;
        }
    }
    delete[] nb;
    for (i = 0; i < n; i++) {
        q[i] = q[i] - modifier;
    }
}

} // namespace gengraph
