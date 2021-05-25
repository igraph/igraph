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
#include "gengraph_box_list.h"
#include <cassert>

namespace gengraph {

void box_list::insert(int v) {
    int d = deg[v];
    if (d < 1) {
        return;
    }
    if (d > dmax) {
        dmax = d;
    }
    int yo = list[d - 1];
    list[d - 1] = v;
    prev[v] = -1;
    next[v] = yo;
    if (yo >= 0) {
        prev[yo] = v;
    }
}

void box_list::pop(int v) {
    int p = prev[v];
    int nxt = next[v];
    if (p < 0) {
        int d = deg[v];
        assert(list[d - 1] == v);
        list[d - 1] = nxt;
        if (d == dmax && nxt < 0) {
            do {
                dmax--;
            } while (dmax > 0 && list[dmax - 1] < 0);
        }
    } else {
        next[p] = nxt;
    }
    if (nxt >= 0) {
        prev[nxt] = p;
    }
}

box_list::box_list(int n0, int *deg0) : n(n0), deg(deg0) {
    next = new int[n];
    prev = new int[n];
    dmax = -1;
    int i;
    for (i = 0; i < n; i++) if (deg[i] > dmax) {
            dmax = deg[i];
        }
    list = new int[dmax];
    for (i = 0; i < dmax; i++) {
        list[i] = -1;
    }
    for (i = 0; i < n; i++) {
        insert(i);
    }
}

box_list::~box_list() {
    delete[] prev;
    delete[] next;
    delete[] list;
}

void box_list::pop_vertex(int v, int **neigh) {
    int k = deg[v];
    if (k < 1) {
        return;
    }
    pop(v);
    int *w = neigh[v];
    while (k--) {
        int v2 = *(w++);
        int *w2 = neigh[v2];
        while (*w2 != v) {
            w2++;
        }
        int *w3 = neigh[v2] + (deg[v2] - 1);
        assert(w2 <= w3);
        int tmp = *w3;
        *w3 = *w2;
        *w2 = tmp;
        pop(v2);
        deg[v2]--;
        insert(v2);
    }
}

} // namespace gengraph
