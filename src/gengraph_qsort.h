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

#include <assert.h>
#include <stdio.h>

namespace gengraph {

//___________________________________________________________________________
// check if every element is zero
inline bool check_zero(int *mem, int n) {
    for (int *v = mem + n; v != mem; ) if (*(--v) != 0) {
            return false;
        }
    return true;
}

//___________________________________________________________________________
//  Sort simple integer arrays in ASCENDING order
//___________________________________________________________________________
inline int med3(int a, int b, int c) {
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

inline void isort(int *v, int t) {
    if (t < 2) {
        return;
    }
    for (int i = 1; i < t; i++) {
        int *w = v + i;
        int tmp = *w;
        while (w != v && *(w - 1) > tmp) {
            *w = *(w - 1);
            w--;
        }
        *w = tmp;
    }
}

inline int partitionne(int *v, int t, int p) {
    int i = 0;
    int j = t - 1;
    while (i < j) {
        while (i <= j && v[i] < p) {
            i++;
        }
        while (i <= j && v[j] > p) {
            j--;
        }
        if (i < j) {
            int tmp = v[i];
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

inline void qsort(int *v, int t) {
    if (t < 15) {
        isort(v, t);
    } else {
        int x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
        qsort(v, x);
        qsort(v + x, t - x);
    }
}

inline int qsort_median(int *v, int t, int pos) {
    if (t < 10) {
        isort(v, t);
        return v[pos];
    }
    int x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
    if (pos < x) {
        return qsort_median(v, x, pos);
    } else {
        return qsort_median(v + x, t - x, pos - x);
    }
}

inline int qsort_median(int *v, int t) {
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

inline void isort(double *v, int t) {
    if (t < 2) {
        return;
    }
    for (int i = 1; i < t; i++) {
        double *w = v + i;
        double tmp = *w;
        while (w != v && *(w - 1) > tmp) {
            *w = *(w - 1);
            w--;
        }
        *w = tmp;
    }
}

inline int partitionne(double *v, int t, double p) {
    int i = 0;
    int j = t - 1;
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

inline void qsort(double *v, int t) {
    if (t < 15) {
        isort(v, t);
    } else {
        int x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
        qsort(v, x);
        qsort(v + x, t - x);
    }
}

inline double qsort_median(double *v, int t, int pos) {
    if (t < 10) {
        isort(v, t);
        return v[pos];
    }
    int x = partitionne(v, t, med3(v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2]));
    if (pos < x) {
        return qsort_median(v, x, pos);
    } else {
        return qsort_median(v + x, t - x, pos - x);
    }
}

inline double qsort_median(double *v, int t) {
    return qsort_median(v, t, t / 2);
}

//___________________________________________________________________________
// Sort integer arrays according to value stored in mem[], in ASCENDING order
inline void isort(int *mem, int *v, int t) {
    if (t < 2) {
        return;
    }
    for (int i = 1; i < t; i++) {
        int vtmp = v[i];
        int tmp = mem[vtmp];
        int j;
        for (j = i; j > 0 && tmp < mem[v[j - 1]]; j--) {
            v[j] = v[j - 1];
        }
        v[j] = vtmp;
    }
}

inline void qsort(int *mem, int *v, int t) {
    if (t < 15) {
        isort(mem, v, t);
    } else {
        int p = med3(mem[v[t >> 1]], mem[v[(t >> 2) + 3]], mem[v[t - (t >> 1) - 3]]);
        int i = 0;
        int j = t - 1;
        while (i < j) {
            while (i <= j && mem[v[i]] < p) {
                i++;
            }
            while (i <= j && mem[v[j]] > p) {
                j--;
            }
            if (i < j) {
                int tmp = v[i];
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
inline int *pre_boxsort(int *mem, int n, int &offset) {
    int *yo;
    // maximum and minimum
    int mx = mem[0];
    int mn = mem[0];
    for (yo = mem + n - 1; yo != mem; yo--) {
        int x = *yo;
        if (x > mx) {
            mx = x;
        }
        if (x < mn) {
            mn = x;
        }
    }
    // box
    int c = mx - mn + 1;
    int *box = new int[c];
    for (yo = box + c; yo != box; * (--yo) = 0) { }
    for (yo = mem + n; yo != mem; box[*(--yo) - mn]++) { }
    // cumul sum
    int sum = 0;
    for (yo = box + c; yo != box; ) {
        sum += *(--yo);
        *yo = sum;
    }
    offset = mn;
    return box;
}

inline int *boxsort(int *mem, int n, int *buff = NULL) {
    int i;
    if (n <= 0) {
        return buff;
    }
    int offset = 0;
    int *box = pre_boxsort(mem, n, offset);
    // sort
    if (buff == NULL) {
        buff = new int[n];
    }
    for (i = 0; i < n; i++) {
        buff[--box[mem[i] - offset]] = i;
    }
    // clean
    delete[] box;
    return buff;
}

// merge two sorted arays in their intersection. Store the result in first array, and return length
inline int intersect(int *a, int a_len, int *b, int b_len) {
    if (a_len == 0 || b_len == 0) {
        return 0;
    }
    int *asup = a + a_len;
    int *bsup = b + b_len;
    int len = 0;
    int *p = a;
    do {
        if (*a == *b) {
            p[len++] = *a;
        }
        do if (++a == asup) {
                return len;
            } while (*a < *b);
        if (*a == *b) {
            p[len++] = *a;
        }
        do if (++b == bsup) {
                return len;
            } while (*b < *a);
    } while (true);
}

// merge two sorted arays in their union, store result in m
inline int unify(int *m, int *a, int a_len, int *b, int b_len) {
    int *asup = a + a_len;
    int *bsup = b + b_len;
    int len = 0;
    while (a != asup && b != bsup) {
        if (*a < *b) {
            m[len++] = *(a++);
        } else {
            if (*a == *b) {
                a++;
            }
            m[len++] = *(b++);
        }
    }
    while (a != asup) {
        m[len++] = *(a++);
    }
    while (b != asup) {
        m[len++] = *(b++);
    }
    return len;
}

// lexicographic compare
inline int lex_comp(int *v1, int *v2, int n) {
    int *stop = v1 + n;
    while (v1 != stop && *v1 == *v2) {
        v1++;
        v2++;
    };
    if (v1 == stop) {
        return 0;
    } else if (*v1 < *v2) {
        return -1;
    } else {
        return 1;
    }
}
// lexicographic median of three
inline int *lex_med3(int *a, int *b, int *c, int s) {
    int ab = lex_comp(a, b, s);
    if (ab == 0) {
        return a;
    } else {
        int cb = lex_comp(c, b, s);
        if (cb == 0) {
            return b;
        }
        int ca = lex_comp(c, a, s);
        if (ab < 0) {
            if (cb > 0) {
                return b;
            } else {
                return (ca > 0) ? c : a;
            }
        } else     {
            if (cb < 0) {
                return b;
            } else {
                return (ca < 0) ? c : a;
            }
        }
    }
}

// Lexicographic sort
inline void lex_isort(int **l, int *v, int t, int s) {
    if (t < 2) {
        return;
    }
    for (int i = 1; i < t; i++) {
        int *w = v + i;
        int tmp = *w;
        while (w != v && lex_comp(l[tmp], l[*(w - 1)], s) < 0) {
            *w = *(w - 1);
            w--;
        }
        *w = tmp;
    }
}

#ifdef _STABLE_SORT_ONLY
    #define _CRITICAL_SIZE_QSORT 0x7FFFFFFF
    #warning "lex_qsort will be replaced by lex_isort"
#else
    #define _CRITICAL_SIZE_QSORT 15
#endif

inline void lex_qsort(int **l, int *v, int t, int s) {

    if (t < _CRITICAL_SIZE_QSORT) {
        lex_isort(l, v, t, s);
    } else {
        int *p = lex_med3(l[v[t >> 1]], l[v[(t >> 2) + 2]], l[v[t - (t >> 1) - 2]], s);
        int i = 0;
        int j = t - 1;
//    printf("pivot = %d\n",p);
        while (i < j) {
//      for(int k=0; k<t; k++) printf("%d ",v[k]);
            while (i <= j && lex_comp(l[v[i]], p, s) < 0) {
                i++;
            }
            while (i <= j && lex_comp(l[v[j]], p, s) > 0) {
                j--;
            }
            if (i < j) {
//        printf("  swap %d[%d] with %d[%d]\n",i,v[i],j,v[j]);
                int tmp = v[i];
                v[i++] = v[j];
                v[j--] = tmp;
            }
        }
        if (i == j && lex_comp(l[v[i]], p, s) < 0) {
            i++;
        }
        assert(i != 0 && i != t);
        lex_qsort(l, v, i, s);
        lex_qsort(l, v + i, t - i, s);
    }
}

// lexicographic indirect compare
inline int lex_comp_indirect(int *key, int *v1, int *v2, int n) {
    int *stop = v1 + n;
    while (v1 != stop && key[*v1] == key[*v2]) {
        v1++;
        v2++;
    };
    if (v1 == stop) {
        return 0;
    } else if (key[*v1] < key[*v2]) {
        return -1;
    } else {
        return 1;
    }
}

inline int qsort_min(const int a, const int b) {
    return a <= b ? a : b;
}

// mix indirect compare
inline int mix_comp_indirect(int *key, int a, int b, int **neigh, int *degs) {
    if (key[a] < key[b]) {
        return -1;
    } else if (key[a] > key[b]) {
        return 1;
    } else {
        int cmp = lex_comp_indirect(key, neigh[a], neigh[b], qsort_min(degs[a], degs[b]));
        if (cmp == 0) {
            if (degs[a] > degs[b]) {
                return -1;
            }
            if (degs[a] < degs[b]) {
                return 1;
            }
        }
        return cmp;
    }
}
// lexicographic indirect median of three
inline int mix_med3_indirect(int *key, int a, int b, int c, int **neigh, int *degs) {
    int ab = mix_comp_indirect(key, a, b, neigh, degs);
    if (ab == 0) {
        return a;
    } else {
        int cb = mix_comp_indirect(key, c, b, neigh, degs);
        if (cb == 0) {
            return b;
        }
        int ca = mix_comp_indirect(key, c, a, neigh, degs);
        if (ab < 0) {
            if (cb > 0) {
                return b;
            } else {
                return (ca > 0) ? c : a;
            }
        } else     {
            if (cb < 0) {
                return b;
            } else {
                return (ca < 0) ? c : a;
            }
        }
    }
}

// Sort integer arrays in ASCENDING order
inline void mix_isort_indirect(int *key, int *v, int t, int **neigh, int *degs) {
    if (t < 2) {
        return;
    }
    for (int i = 1; i < t; i++) {
        int *w = v + i;
        int tmp = *w;
        while (w != v && mix_comp_indirect(key, tmp, *(w - 1), neigh, degs) < 0) {
            *w = *(w - 1);
            w--;
        }
        *w = tmp;
    }
}

inline void mix_qsort_indirect(int *key, int *v, int t, int **neigh, int *degs) {
    if (t < 15) {
        mix_isort_indirect(key, v, t, neigh, degs);
    } else {
        int p = mix_med3_indirect(key, v[t >> 1], v[(t >> 2) + 2], v[t - (t >> 1) - 2], neigh, degs);
        int i = 0;
        int j = t - 1;
//    printf("pivot = %d\n",p);
        while (i < j) {
//      for(int k=0; k<t; k++) printf("%d ",v[k]);
            while (i <= j && mix_comp_indirect(key, v[i], p, neigh, degs) < 0) {
                i++;
            }
            while (i <= j && mix_comp_indirect(key, v[j], p, neigh, degs) > 0) {
                j--;
            }
            if (i < j) {
//        printf("  swap %d[%d] with %d[%d]\n",i,v[i],j,v[j]);
                int tmp = v[i];
                v[i++] = v[j];
                v[j--] = tmp;
            }
        }
        if (i == j && mix_comp_indirect(key, v[i], p, neigh, degs) < 0) {
            i++;
        }
        assert(i != 0 && i != t);
        mix_qsort_indirect(key, v, i, neigh, degs);
        mix_qsort_indirect(key, v + i, t - i, neigh, degs);
    }
}

} // namespace gengraph

#endif //QSORT_H
