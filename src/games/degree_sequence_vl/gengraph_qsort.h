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
#ifndef QSORT_H
#define QSORT_H

#include "igraph_types.h"

#include <assert.h>
#include <stdio.h>

namespace gengraph {

//___________________________________________________________________________
// check if every element is zero
inline bool check_zero(igraph_integer_t *mem, igraph_integer_t n) {
    for (igraph_integer_t *v = mem + n; v != mem; ) {
        if (*(--v) != 0) {
            return false;
        }
    }
    return true;
}

//___________________________________________________________________________
//  Sort simple integer arrays in ASCENDING order
//___________________________________________________________________________
inline igraph_integer_t med3(igraph_integer_t a, igraph_integer_t b, igraph_integer_t c) {
    if (a < b) {
        if (c < b) {
            return (a < c) ? c : a;
        } else {
            return b;
        }
    } else {
        if (c < a) {
            return (b < c) ? c : b;
        } else {
            return a;
        }
    }
}

inline void isort(igraph_integer_t *v, igraph_integer_t t) {
    if (t < 2) {
        return;
    }
    for (igraph_integer_t i = 1; i < t; i++) {
        igraph_integer_t *w = v + i;
        igraph_integer_t tmp = *w;
        while (w != v && *(w - 1) > tmp) {
            *w = *(w - 1);
            w--;
        }
        *w = tmp;
    }
}

inline igraph_integer_t partitionne(igraph_integer_t *v, igraph_integer_t t, igraph_integer_t p) {
    igraph_integer_t i = 0;
    igraph_integer_t j = t - 1;
    while (i < j) {
        while (i <= j && v[i] < p) {
            i++;
        }
        while (i <= j && v[j] > p) {
            j--;
        }
        if (i < j) {
            igraph_integer_t tmp = v[i];
            v[i++] = v[j];
            v[j--] = tmp;
        }
    }
    if (i == j && v[i] < p) {
        i++;
    }
    assert(i != 0 && i != t);
    return i;
}

inline void qsort(igraph_integer_t *v, igraph_integer_t t) {
    if (t < 15) {
        isort(v, t);
    } else {
        igraph_integer_t x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
        qsort(v, x);
        qsort(v + x, t - x);
    }
}

inline igraph_integer_t qsort_median(igraph_integer_t *v, igraph_integer_t t, igraph_integer_t pos) {
    if (t < 10) {
        isort(v, t);
        return v[pos];
    }
    igraph_integer_t x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
    if (pos < x) {
        return qsort_median(v, x, pos);
    } else {
        return qsort_median(v + x, t - x, pos - x);
    }
}

inline igraph_integer_t qsort_median(igraph_integer_t *v, igraph_integer_t t) {
    return qsort_median(v, t, t / 2);
}

//___________________________________________________________________________
//  Sort simple double arrays in ASCENDING order
//___________________________________________________________________________
inline double med3(double a, double b, double c) {
    if (a < b) {
        if (c < b) {
            return (a < c) ? c : a;
        } else {
            return b;
        }
    } else {
        if (c < a) {
            return (b < c) ? c : b;
        } else {
            return a;
        }
    }
}

inline void isort(double *v, igraph_integer_t t) {
    if (t < 2) {
        return;
    }
    for (igraph_integer_t i = 1; i < t; i++) {
        double *w = v + i;
        double tmp = *w;
        while (w != v && *(w - 1) > tmp) {
            *w = *(w - 1);
            w--;
        }
        *w = tmp;
    }
}

inline igraph_integer_t partitionne(double *v, igraph_integer_t t, double p) {
    igraph_integer_t i = 0;
    igraph_integer_t j = t - 1;
    while (i < j) {
        while (i <= j && v[i] < p) {
            i++;
        }
        while (i <= j && v[j] > p) {
            j--;
        }
        if (i < j) {
            double tmp = v[i];
            v[i++] = v[j];
            v[j--] = tmp;
        }
    }
    if (i == j && v[i] < p) {
        i++;
    }
    assert(i != 0 && i != t);
    return i;
}

inline void qsort(double *v, igraph_integer_t t) {
    if (t < 15) {
        isort(v, t);
    } else {
        igraph_integer_t x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
        qsort(v, x);
        qsort(v + x, t - x);
    }
}

inline double qsort_median(double *v, igraph_integer_t t, igraph_integer_t pos) {
    if (t < 10) {
        isort(v, t);
        return v[pos];
    }
    igraph_integer_t x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
    if (pos < x) {
        return qsort_median(v, x, pos);
    } else {
        return qsort_median(v + x, t - x, pos - x);
    }
}

inline double qsort_median(double *v, igraph_integer_t t) {
    return qsort_median(v, t, t / 2);
}

//___________________________________________________________________________
// Sort integer arrays according to value stored in mem[], in ASCENDING order
inline void isort(igraph_integer_t *mem, igraph_integer_t *v, igraph_integer_t t) {
    if (t < 2) {
        return;
    }
    for (igraph_integer_t i = 1; i < t; i++) {
        igraph_integer_t vtmp = v[i];
        igraph_integer_t tmp = mem[vtmp];
        igraph_integer_t j;
        for (j = i; j > 0 && tmp < mem[v[j - 1]]; j--) {
            v[j] = v[j - 1];
        }
        v[j] = vtmp;
    }
}

inline void qsort(igraph_integer_t *mem, igraph_integer_t *v, igraph_integer_t t) {
    if (t < 15) {
        isort(mem, v, t);
    } else {
        igraph_integer_t p = med3(mem[v[t >> 1]], mem[v[(t >> 2) + 3]], mem[v[t - (t >> 1) - 3]]);
        igraph_integer_t i = 0;
        igraph_integer_t j = t - 1;
        while (i < j) {
            while (i <= j && mem[v[i]] < p) {
                i++;
            }
            while (i <= j && mem[v[j]] > p) {
                j--;
            }
            if (i < j) {
                igraph_integer_t tmp = v[i];
                v[i++] = v[j];
                v[j--] = tmp;
            }
        }
        if (i == j && mem[v[i]] < p) {
            i++;
        }
        assert(i != 0 && i != t);
        qsort(mem, v, i);
        qsort(mem, v + i, t - i);
    }
}

//Box-Sort 1..n according to value stored in mem[], in DESCENDING order.
inline igraph_integer_t *pre_boxsort(igraph_integer_t *mem, igraph_integer_t n, igraph_integer_t &offset) {
    igraph_integer_t *yo;
    // maximum and minimum
    igraph_integer_t mx = mem[0];
    igraph_integer_t mn = mem[0];
    for (yo = mem + n - 1; yo != mem; yo--) {
        igraph_integer_t x = *yo;
        if (x > mx) {
            mx = x;
        }
        if (x < mn) {
            mn = x;
        }
    }
    // box
    igraph_integer_t c = mx - mn + 1;
    igraph_integer_t *box = new igraph_integer_t[c];
    for (yo = box + c; yo != box; * (--yo) = 0) { }
    for (yo = mem + n; yo != mem; box[*(--yo) - mn]++) { }
    // cumul sum
    igraph_integer_t sum = 0;
    for (yo = box + c; yo != box; ) {
        sum += *(--yo);
        *yo = sum;
    }
    offset = mn;
    return box;
}

inline igraph_integer_t *boxsort(igraph_integer_t *mem, igraph_integer_t n, igraph_integer_t *buff = NULL) {
    igraph_integer_t i;
    if (n <= 0) {
        return buff;
    }
    igraph_integer_t offset = 0;
    igraph_integer_t *box = pre_boxsort(mem, n, offset);
    // sort
    if (buff == NULL) {
        buff = new igraph_integer_t[n];
    }
    for (i = 0; i < n; i++) {
        buff[--box[mem[i] - offset]] = i;
    }
    // clean
    delete[] box;
    return buff;
}

} // namespace gengraph

#endif //QSORT_H
