/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2020  The igraph development team

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_error.h"

#include "core/indheap.h"

#include <string.h>         /* memcpy & co. */
#include <stdlib.h>

/* -------------------------------------------------- */
/* Indexed heap                                       */
/* -------------------------------------------------- */

#define PARENT(x)     (((x)+1)/2-1)
#define LEFTCHILD(x)  (((x)+1)*2-1)
#define RIGHTCHILD(x) (((x)+1)*2)

static void igraph_indheap_i_build(igraph_indheap_t* h, long int head);
static void igraph_indheap_i_shift_up(igraph_indheap_t* h, long int elem);
static void igraph_indheap_i_sink(igraph_indheap_t* h, long int head);
static void igraph_indheap_i_switch(igraph_indheap_t* h, long int e1, long int e2);

/**
 * \ingroup indheap
 * \brief Initializes an indexed heap (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_indheap_init(igraph_indheap_t* h, long int alloc_size) {
    if (alloc_size <= 0 ) {
        alloc_size = 1;
    }
    h->stor_begin = IGRAPH_CALLOC(alloc_size, igraph_real_t);
    if (h->stor_begin == 0) {
        h->index_begin = 0;
        IGRAPH_ERROR("indheap init failed", IGRAPH_ENOMEM);
    }
    h->index_begin = IGRAPH_CALLOC(alloc_size, long int);
    if (h->index_begin == 0) {
        IGRAPH_FREE(h->stor_begin);
        h->stor_begin = 0;
        IGRAPH_ERROR("indheap init failed", IGRAPH_ENOMEM);
    }

    h->stor_end = h->stor_begin + alloc_size;
    h->end = h->stor_begin;
    h->destroy = 1;

    return 0;
}

int igraph_indheap_clear(igraph_indheap_t *h) {
    h->end = h->stor_begin;
    return 0;
}

/**
 * \ingroup indheap
 * \brief Initializes and build an indexed heap from a C array (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_indheap_init_array     (igraph_indheap_t *h, igraph_real_t* data, long int len) {
    long int i;

    h->stor_begin = IGRAPH_CALLOC(len, igraph_real_t);
    if (h->stor_begin == 0) {
        h->index_begin = 0;
        IGRAPH_ERROR("indheap init from array failed", IGRAPH_ENOMEM);
    }
    h->index_begin = IGRAPH_CALLOC(len, long int);
    if (h->index_begin == 0) {
        IGRAPH_FREE(h->stor_begin);
        h->stor_begin = 0;
        IGRAPH_ERROR("indheap init from array failed", IGRAPH_ENOMEM);
    }
    h->stor_end = h->stor_begin + len;
    h->end = h->stor_end;
    h->destroy = 1;

    memcpy(h->stor_begin, data, (size_t) len * sizeof(igraph_real_t));
    for (i = 0; i < len; i++) {
        h->index_begin[i] = i + 1;
    }

    igraph_indheap_i_build (h, 0);

    return 0;
}

/**
 * \ingroup indheap
 * \brief Destroys an initialized indexed heap.
 */

void igraph_indheap_destroy        (igraph_indheap_t* h) {
    IGRAPH_ASSERT(h != 0);
    if (h->destroy) {
        if (h->stor_begin != 0) {
            IGRAPH_FREE(h->stor_begin);
            h->stor_begin = 0;
        }
        if (h->index_begin != 0) {
            IGRAPH_FREE(h->index_begin);
            h->index_begin = 0;
        }
    }
}

/**
 * \ingroup indheap
 * \brief Checks whether a heap is empty.
 */

igraph_bool_t igraph_indheap_empty          (igraph_indheap_t* h) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);
    return h->stor_begin == h->end;
}

/**
 * \ingroup indheap
 * \brief Adds an element to an indexed heap.
 */

int igraph_indheap_push           (igraph_indheap_t* h, igraph_real_t elem) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);

    /* full, allocate more storage */
    if (h->stor_end == h->end) {
        long int new_size = igraph_indheap_size(h) * 2;
        if (new_size == 0) {
            new_size = 1;
        }
        IGRAPH_CHECK(igraph_indheap_reserve(h, new_size));
    }

    *(h->end) = elem;
    h->end += 1;
    *(h->index_begin + igraph_indheap_size(h) - 1) = igraph_indheap_size(h) - 1;

    /* maintain indheap */
    igraph_indheap_i_shift_up(h, igraph_indheap_size(h) - 1);

    return 0;
}

/**
 * \ingroup indheap
 * \brief Adds an element to an indexed heap with a given index.
 */

int igraph_indheap_push_with_index(igraph_indheap_t* h, long int idx, igraph_real_t elem) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);

    /* full, allocate more storage */
    if (h->stor_end == h->end) {
        long int new_size = igraph_indheap_size(h) * 2;
        if (new_size == 0) {
            new_size = 1;
        }
        IGRAPH_CHECK(igraph_indheap_reserve(h, new_size));
    }

    *(h->end) = elem;
    h->end += 1;
    *(h->index_begin + igraph_indheap_size(h) - 1) = idx;

    /* maintain indheap */
    igraph_indheap_i_shift_up(h, igraph_indheap_size(h) - 1);

    return 0;
}

/**
 * \ingroup indheap
 * \brief Modifies an element in an indexed heap.
 */

int igraph_indheap_modify(igraph_indheap_t* h, long int idx, igraph_real_t elem) {
    long int i, n;

    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);

    n = igraph_indheap_size(h);
    for (i = 0; i < n; i++)
        if (h->index_begin[i] == idx) {
            h->stor_begin[i] = elem;
            break;
        }

    if (i == n) {
        return 0;
    }

    /* maintain indheap */
    igraph_indheap_i_build(h, 0);

    return 0;
}

/**
 * \ingroup indheap
 * \brief Returns the largest element in an indexed heap.
 */

igraph_real_t igraph_indheap_max       (igraph_indheap_t* h) {
    IGRAPH_ASSERT(h != NULL);
    IGRAPH_ASSERT(h->stor_begin != NULL);
    IGRAPH_ASSERT(h->stor_begin != h->end);

    return h->stor_begin[0];
}

/**
 * \ingroup indheap
 * \brief Removes the largest element from an indexed heap.
 */

igraph_real_t igraph_indheap_delete_max(igraph_indheap_t* h) {
    igraph_real_t tmp;

    IGRAPH_ASSERT(h != NULL);
    IGRAPH_ASSERT(h->stor_begin != NULL);

    tmp = h->stor_begin[0];
    igraph_indheap_i_switch(h, 0, igraph_indheap_size(h) - 1);
    h->end -= 1;
    igraph_indheap_i_sink(h, 0);

    return tmp;
}

/**
 * \ingroup indheap
 * \brief Gives the number of elements in an indexed heap.
 */

long int igraph_indheap_size      (igraph_indheap_t* h) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);
    return h->end - h->stor_begin;
}

/**
 * \ingroup indheap
 * \brief Reserves more memory for an indexed heap.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_indheap_reserve        (igraph_indheap_t* h, long int size) {
    long int actual_size = igraph_indheap_size(h);
    igraph_real_t *tmp1;
    long int *tmp2;
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);

    if (size <= actual_size) {
        return 0;
    }

    tmp1 = IGRAPH_CALLOC(size, igraph_real_t);
    if (tmp1 == 0) {
        IGRAPH_ERROR("indheap reserve failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmp1);
    tmp2 = IGRAPH_CALLOC(size, long int);
    if (tmp2 == 0) {
        IGRAPH_ERROR("indheap reserve failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmp2);
    memcpy(tmp1, h->stor_begin, (size_t) actual_size * sizeof(igraph_real_t));
    memcpy(tmp2, h->index_begin, (size_t) actual_size * sizeof(long int));
    IGRAPH_FREE(h->stor_begin);
    IGRAPH_FREE(h->index_begin);

    h->stor_begin = tmp1;
    h->index_begin = tmp2;
    h->stor_end = h->stor_begin + size;
    h->end = h->stor_begin + actual_size;

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

/**
 * \ingroup indheap
 * \brief Returns the index of the largest element in an indexed heap.
 */

long int igraph_indheap_max_index(igraph_indheap_t *h) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);
    return h->index_begin[0];
}

/**
 * \ingroup indheap
 * \brief Builds an indexed heap, this function should not be called
 * directly.
 */

static void igraph_indheap_i_build(igraph_indheap_t* h, long int head) {

    long int size = igraph_indheap_size(h);
    if (RIGHTCHILD(head) < size) {
        /* both subtrees */
        igraph_indheap_i_build(h, LEFTCHILD(head) );
        igraph_indheap_i_build(h, RIGHTCHILD(head));
        igraph_indheap_i_sink(h, head);
    } else if (LEFTCHILD(head) < size) {
        /* only left */
        igraph_indheap_i_build(h, LEFTCHILD(head));
        igraph_indheap_i_sink(h, head);
    } else {
        /* none */
    }
}

/**
 * \ingroup indheap
 * \brief Moves an element up in the heap, don't call this function
 * directly.
 */

static void igraph_indheap_i_shift_up(igraph_indheap_t *h, long int elem) {

    if (elem == 0 || h->stor_begin[elem] < h->stor_begin[PARENT(elem)]) {
        /* at the top */
    } else {
        igraph_indheap_i_switch(h, elem, PARENT(elem));
        igraph_indheap_i_shift_up(h, PARENT(elem));
    }
}

/**
 * \ingroup indheap
 * \brief Moves an element down in the heap, don't call this function
 * directly.
 */

static void igraph_indheap_i_sink(igraph_indheap_t* h, long int head) {

    long int size = igraph_indheap_size(h);
    if (LEFTCHILD(head) >= size) {
        /* no subtrees */
    } else if (RIGHTCHILD(head) == size ||
               h->stor_begin[LEFTCHILD(head)] >= h->stor_begin[RIGHTCHILD(head)]) {
        /* sink to the left if needed */
        if (h->stor_begin[head] < h->stor_begin[LEFTCHILD(head)]) {
            igraph_indheap_i_switch(h, head, LEFTCHILD(head));
            igraph_indheap_i_sink(h, LEFTCHILD(head));
        }
    } else {
        /* sink to the right */
        if (h->stor_begin[head] < h->stor_begin[RIGHTCHILD(head)]) {
            igraph_indheap_i_switch(h, head, RIGHTCHILD(head));
            igraph_indheap_i_sink(h, RIGHTCHILD(head));
        }
    }
}

/**
 * \ingroup indheap
 * \brief Switches two elements in a heap, don't call this function
 * directly.
 */

static void igraph_indheap_i_switch(igraph_indheap_t* h, long int e1, long int e2) {
    if (e1 != e2) {
        igraph_real_t tmp = h->stor_begin[e1];
        h->stor_begin[e1] = h->stor_begin[e2];
        h->stor_begin[e2] = tmp;

        tmp = h->index_begin[e1];
        h->index_begin[e1] = h->index_begin[e2];
        h->index_begin[e2] = (long int) tmp;
    }
}

/*************************************************/

/* -------------------------------------------------- */
/* Doubly indexed heap                                */
/* -------------------------------------------------- */

/* static void igraph_d_indheap_i_build(igraph_d_indheap_t* h, long int head); */ /* Unused function */
static void igraph_d_indheap_i_shift_up(igraph_d_indheap_t* h, long int elem);
static void igraph_d_indheap_i_sink(igraph_d_indheap_t* h, long int head);
static void igraph_d_indheap_i_switch(igraph_d_indheap_t* h, long int e1, long int e2);

/**
 * \ingroup doubleindheap
 * \brief Initializes an empty doubly indexed heap object (constructor).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_d_indheap_init           (igraph_d_indheap_t* h, long int alloc_size) {
    if (alloc_size <= 0 ) {
        alloc_size = 1;
    }
    h->stor_begin = IGRAPH_CALLOC(alloc_size, igraph_real_t);
    if (h->stor_begin == 0) {
        h->index_begin = 0;
        h->index2_begin = 0;
        IGRAPH_ERROR("d_indheap init failed", IGRAPH_ENOMEM);
    }
    h->stor_end = h->stor_begin + alloc_size;
    h->end = h->stor_begin;
    h->destroy = 1;
    h->index_begin = IGRAPH_CALLOC(alloc_size, long int);
    if (h->index_begin == 0) {
        IGRAPH_FREE(h->stor_begin);
        h->stor_begin = 0;
        h->index2_begin = 0;
        IGRAPH_ERROR("d_indheap init failed", IGRAPH_ENOMEM);
    }
    h->index2_begin = IGRAPH_CALLOC(alloc_size, long int);
    if (h->index2_begin == 0) {
        IGRAPH_FREE(h->stor_begin);
        IGRAPH_FREE(h->index_begin);
        h->stor_begin = 0;
        h->index_begin = 0;
        IGRAPH_ERROR("d_indheap init failed", IGRAPH_ENOMEM);
    }

    return 0;
}

/**
 * \ingroup doubleindheap
 * \brief Destroys an initialized doubly indexed heap object.
 */

void igraph_d_indheap_destroy        (igraph_d_indheap_t* h) {
    IGRAPH_ASSERT(h != 0);
    if (h->destroy) {
        if (h->stor_begin != 0) {
            IGRAPH_FREE(h->stor_begin);
            h->stor_begin = 0;
        }
        if (h->index_begin != 0) {
            IGRAPH_FREE(h->index_begin);
            h->index_begin = 0;
        }
        if (h->index2_begin != 0) {
            IGRAPH_FREE(h->index2_begin);
            h->index2_begin = 0;
        }
    }
}

/**
 * \ingroup doubleindheap
 * \brief Decides whether a heap is empty.
 */

igraph_bool_t igraph_d_indheap_empty          (igraph_d_indheap_t* h) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);
    return h->stor_begin == h->end;
}

/**
 * \ingroup doubleindheap
 * \brief Adds an element to the heap.
 */

int igraph_d_indheap_push           (igraph_d_indheap_t* h, igraph_real_t elem,
                                     long int idx, long int idx2) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);

    /* full, allocate more storage */
    if (h->stor_end == h->end) {
        long int new_size = igraph_d_indheap_size(h) * 2;
        if (new_size == 0) {
            new_size = 1;
        }
        IGRAPH_CHECK(igraph_d_indheap_reserve(h, new_size));
    }

    *(h->end) = elem;
    h->end += 1;
    *(h->index_begin + igraph_d_indheap_size(h) - 1) = idx ;
    *(h->index2_begin + igraph_d_indheap_size(h) - 1) = idx2 ;

    /* maintain d_indheap */
    igraph_d_indheap_i_shift_up(h, igraph_d_indheap_size(h) - 1);

    return 0;
}

/**
 * \ingroup doubleindheap
 * \brief Returns the largest element in the heap.
 */

igraph_real_t igraph_d_indheap_max       (igraph_d_indheap_t* h) {
    IGRAPH_ASSERT(h != NULL);
    IGRAPH_ASSERT(h->stor_begin != NULL);
    IGRAPH_ASSERT(h->stor_begin != h->end);

    return h->stor_begin[0];
}

/**
 * \ingroup doubleindheap
 * \brief Removes the largest element from the heap.
 */

igraph_real_t igraph_d_indheap_delete_max(igraph_d_indheap_t* h) {
    igraph_real_t tmp;

    IGRAPH_ASSERT(h != NULL);
    IGRAPH_ASSERT(h->stor_begin != NULL);

    tmp = h->stor_begin[0];
    igraph_d_indheap_i_switch(h, 0, igraph_d_indheap_size(h) - 1);
    h->end -= 1;
    igraph_d_indheap_i_sink(h, 0);

    return tmp;
}

/**
 * \ingroup doubleindheap
 * \brief Gives the number of elements in the heap.
 */

long int igraph_d_indheap_size      (igraph_d_indheap_t* h) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);
    return h->end - h->stor_begin;
}

/**
 * \ingroup doubleindheap
 * \brief Allocates memory for a heap.
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_d_indheap_reserve        (igraph_d_indheap_t* h, long int size) {
    long int actual_size = igraph_d_indheap_size(h);
    igraph_real_t *tmp1;
    long int *tmp2, *tmp3;
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);

    if (size <= actual_size) {
        return 0;
    }

    tmp1 = IGRAPH_CALLOC(size, igraph_real_t);
    if (tmp1 == 0) {
        IGRAPH_ERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmp1);
    tmp2 = IGRAPH_CALLOC(size, long int);
    if (tmp2 == 0) {
        IGRAPH_ERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmp2);
    tmp3 = IGRAPH_CALLOC(size, long int);
    if (tmp3 == 0) {
        IGRAPH_ERROR("d_indheap reserve failed", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, tmp3);

    memcpy(tmp1, h->stor_begin, (size_t) actual_size * sizeof(igraph_real_t));
    memcpy(tmp2, h->index_begin, (size_t) actual_size * sizeof(long int));
    memcpy(tmp3, h->index2_begin, (size_t) actual_size * sizeof(long int));
    IGRAPH_FREE(h->stor_begin);
    IGRAPH_FREE(h->index_begin);
    IGRAPH_FREE(h->index2_begin);

    h->stor_begin = tmp1;
    h->stor_end = h->stor_begin + size;
    h->end = h->stor_begin + actual_size;
    h->index_begin = tmp2;
    h->index2_begin = tmp3;

    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

/**
 * \ingroup doubleindheap
 * \brief Gives the indices of the maximal element in the heap.
 */

void igraph_d_indheap_max_index(igraph_d_indheap_t *h, long int *idx, long int *idx2) {
    IGRAPH_ASSERT(h != 0);
    IGRAPH_ASSERT(h->stor_begin != 0);
    (*idx) = h->index_begin[0];
    (*idx2) = h->index2_begin[0];
}

/**
 * \ingroup doubleindheap
 * \brief Builds the heap, don't call it directly.
 */

/* Unused function, temporarily disabled */
#if 0
static void igraph_d_indheap_i_build(igraph_d_indheap_t* h, long int head) {

    long int size = igraph_d_indheap_size(h);
    if (RIGHTCHILD(head) < size) {
        /* both subtrees */
        igraph_d_indheap_i_build(h, LEFTCHILD(head) );
        igraph_d_indheap_i_build(h, RIGHTCHILD(head));
        igraph_d_indheap_i_sink(h, head);
    } else if (LEFTCHILD(head) < size) {
        /* only left */
        igraph_d_indheap_i_build(h, LEFTCHILD(head));
        igraph_d_indheap_i_sink(h, head);
    } else {
        /* none */
    }
}
#endif

/**
 * \ingroup doubleindheap
 * \brief Moves an element up in the heap, don't call it directly.
 */

static void igraph_d_indheap_i_shift_up(igraph_d_indheap_t *h, long int elem) {

    if (elem == 0 || h->stor_begin[elem] < h->stor_begin[PARENT(elem)]) {
        /* at the top */
    } else {
        igraph_d_indheap_i_switch(h, elem, PARENT(elem));
        igraph_d_indheap_i_shift_up(h, PARENT(elem));
    }
}

/**
 * \ingroup doubleindheap
 * \brief Moves an element down in the heap, don't call it directly.
 */

static void igraph_d_indheap_i_sink(igraph_d_indheap_t* h, long int head) {

    long int size = igraph_d_indheap_size(h);
    if (LEFTCHILD(head) >= size) {
        /* no subtrees */
    } else if (RIGHTCHILD(head) == size ||
               h->stor_begin[LEFTCHILD(head)] >= h->stor_begin[RIGHTCHILD(head)]) {
        /* sink to the left if needed */
        if (h->stor_begin[head] < h->stor_begin[LEFTCHILD(head)]) {
            igraph_d_indheap_i_switch(h, head, LEFTCHILD(head));
            igraph_d_indheap_i_sink(h, LEFTCHILD(head));
        }
    } else {
        /* sink to the right */
        if (h->stor_begin[head] < h->stor_begin[RIGHTCHILD(head)]) {
            igraph_d_indheap_i_switch(h, head, RIGHTCHILD(head));
            igraph_d_indheap_i_sink(h, RIGHTCHILD(head));
        }
    }
}

/**
 * \ingroup doubleindheap
 * \brief Switches two elements in the heap, don't call it directly.
 */

static void igraph_d_indheap_i_switch(igraph_d_indheap_t* h, long int e1, long int e2) {
    if (e1 != e2) {
        long int tmpi;
        igraph_real_t tmp = h->stor_begin[e1];
        h->stor_begin[e1] = h->stor_begin[e2];
        h->stor_begin[e2] = tmp;

        tmpi = h->index_begin[e1];
        h->index_begin[e1] = h->index_begin[e2];
        h->index_begin[e2] = tmpi;

        tmpi = h->index2_begin[e1];
        h->index2_begin[e1] = h->index2_begin[e2];
        h->index2_begin[e2] = tmpi;
    }
}

/*************************************************/

/* -------------------------------------------------- */
/* Two-way indexed heap                               */
/* -------------------------------------------------- */

#undef PARENT
#undef LEFTCHILD
#undef RIGHTCHILD
#define PARENT(x)     (((x)+1)/2-1)
#define LEFTCHILD(x)  (((x)+1)*2-1)
#define RIGHTCHILD(x) (((x)+1)*2)

/* This is a smart indexed heap. In addition to the "normal" indexed heap
   it allows to access every element through its index in O(1) time.
   In other words, for this heap the indexing operation is O(1), the
   normal heap does this in O(n) time.... */

static void igraph_i_2wheap_switch(igraph_2wheap_t *h,
                                   long int e1, long int e2) {
    if (e1 != e2) {
        long int tmp1, tmp2;
        igraph_real_t tmp3 = VECTOR(h->data)[e1];
        VECTOR(h->data)[e1] = VECTOR(h->data)[e2];
        VECTOR(h->data)[e2] = tmp3;

        tmp1 = VECTOR(h->index)[e1];
        tmp2 = VECTOR(h->index)[e2];

        VECTOR(h->index2)[tmp1] = e2 + 2;
        VECTOR(h->index2)[tmp2] = e1 + 2;

        VECTOR(h->index)[e1] = tmp2;
        VECTOR(h->index)[e2] = tmp1;
    }
}

static void igraph_i_2wheap_shift_up(igraph_2wheap_t *h,
                                     long int elem) {
    if (elem == 0 || VECTOR(h->data)[elem] < VECTOR(h->data)[PARENT(elem)]) {
        /* at the top */
    } else {
        igraph_i_2wheap_switch(h, elem, PARENT(elem));
        igraph_i_2wheap_shift_up(h, PARENT(elem));
    }
}

static void igraph_i_2wheap_sink(igraph_2wheap_t *h,
                                 long int head) {
    long int size = igraph_2wheap_size(h);
    if (LEFTCHILD(head) >= size) {
        /* no subtrees */
    } else if (RIGHTCHILD(head) == size ||
               VECTOR(h->data)[LEFTCHILD(head)] >= VECTOR(h->data)[RIGHTCHILD(head)]) {
        /* sink to the left if needed */
        if (VECTOR(h->data)[head] < VECTOR(h->data)[LEFTCHILD(head)]) {
            igraph_i_2wheap_switch(h, head, LEFTCHILD(head));
            igraph_i_2wheap_sink(h, LEFTCHILD(head));
        }
    } else {
        /* sink to the right */
        if (VECTOR(h->data)[head] < VECTOR(h->data)[RIGHTCHILD(head)]) {
            igraph_i_2wheap_switch(h, head, RIGHTCHILD(head));
            igraph_i_2wheap_sink(h, RIGHTCHILD(head));
        }
    }
}

/* ------------------ */
/* These are public   */
/* ------------------ */

int igraph_2wheap_init(igraph_2wheap_t *h, long int size) {
    h->size = size;
    /* We start with the biggest */
    IGRAPH_CHECK(igraph_vector_long_init(&h->index2, size));
    IGRAPH_FINALLY(igraph_vector_long_destroy, &h->index2);
    IGRAPH_VECTOR_INIT_FINALLY(&h->data, 0);
    IGRAPH_CHECK(igraph_vector_long_init(&h->index, 0));
    /* IGRAPH_FINALLY(igraph_vector_long_destroy, &h->index); */

    IGRAPH_FINALLY_CLEAN(2);
    return 0;
}

void igraph_2wheap_destroy(igraph_2wheap_t *h) {
    igraph_vector_destroy(&h->data);
    igraph_vector_long_destroy(&h->index);
    igraph_vector_long_destroy(&h->index2);
}

int igraph_2wheap_clear(igraph_2wheap_t *h) {
    igraph_vector_clear(&h->data);
    igraph_vector_long_clear(&h->index);
    igraph_vector_long_null(&h->index2);
    return 0;
}

igraph_bool_t igraph_2wheap_empty(const igraph_2wheap_t *h) {
    return igraph_vector_empty(&h->data);
}

int igraph_2wheap_push_with_index(igraph_2wheap_t *h,
                                  long int idx, igraph_real_t elem) {

    /*   printf("-> %.2g [%li]\n", elem, idx); */

    long int size = igraph_vector_size(&h->data);
    IGRAPH_CHECK(igraph_vector_push_back(&h->data, elem));
    IGRAPH_CHECK(igraph_vector_long_push_back(&h->index, idx));
    VECTOR(h->index2)[idx] = size + 2;

    /* maintain heap */
    igraph_i_2wheap_shift_up(h, size);
    return 0;
}

long int igraph_2wheap_size(const igraph_2wheap_t *h) {
    return igraph_vector_size(&h->data);
}

long int igraph_2wheap_max_size(const igraph_2wheap_t *h) {
    return h->size;
}

igraph_real_t igraph_2wheap_max(const igraph_2wheap_t *h) {
    return VECTOR(h->data)[0];
}

long int igraph_2wheap_max_index(const igraph_2wheap_t *h) {
    return VECTOR(h->index)[0];
}

igraph_bool_t igraph_2wheap_has_elem(const igraph_2wheap_t *h, long int idx) {
    return VECTOR(h->index2)[idx] != 0;
}

igraph_bool_t igraph_2wheap_has_active(const igraph_2wheap_t *h, long int idx) {
    return VECTOR(h->index2)[idx] > 1;
}

igraph_real_t igraph_2wheap_get(const igraph_2wheap_t *h, long int idx) {
    long int i = VECTOR(h->index2)[idx] - 2;
    return VECTOR(h->data)[i];
}

igraph_real_t igraph_2wheap_delete_max(igraph_2wheap_t *h) {

    igraph_real_t tmp = VECTOR(h->data)[0];
    long int tmpidx = VECTOR(h->index)[0];
    igraph_i_2wheap_switch(h, 0, igraph_2wheap_size(h) - 1);
    igraph_vector_pop_back(&h->data);
    igraph_vector_long_pop_back(&h->index);
    VECTOR(h->index2)[tmpidx] = 0;
    igraph_i_2wheap_sink(h, 0);

    /*   printf("<-max %.2g\n", tmp); */

    return tmp;
}

igraph_real_t igraph_2wheap_deactivate_max(igraph_2wheap_t *h) {

    igraph_real_t tmp = VECTOR(h->data)[0];
    long int tmpidx = VECTOR(h->index)[0];
    igraph_i_2wheap_switch(h, 0, igraph_2wheap_size(h) - 1);
    igraph_vector_pop_back(&h->data);
    igraph_vector_long_pop_back(&h->index);
    VECTOR(h->index2)[tmpidx] = 1;
    igraph_i_2wheap_sink(h, 0);

    return tmp;
}

igraph_real_t igraph_2wheap_delete_max_index(igraph_2wheap_t *h, long int *idx) {

    igraph_real_t tmp = VECTOR(h->data)[0];
    long int tmpidx = VECTOR(h->index)[0];
    igraph_i_2wheap_switch(h, 0, igraph_2wheap_size(h) - 1);
    igraph_vector_pop_back(&h->data);
    igraph_vector_long_pop_back(&h->index);
    VECTOR(h->index2)[tmpidx] = 0;
    igraph_i_2wheap_sink(h, 0);

    if (idx) {
        *idx = tmpidx;
    }
    return tmp;
}

int igraph_2wheap_modify(igraph_2wheap_t *h, long int idx, igraph_real_t elem) {

    long int pos = VECTOR(h->index2)[idx] - 2;

    /*   printf("-- %.2g -> %.2g\n", VECTOR(h->data)[pos], elem); */

    VECTOR(h->data)[pos] = elem;
    igraph_i_2wheap_sink(h, pos);
    igraph_i_2wheap_shift_up(h, pos);

    return 0;
}

/* Check that the heap is in a consistent state */

int igraph_2wheap_check(igraph_2wheap_t *h) {
    long int size = igraph_2wheap_size(h);
    long int i;
    igraph_bool_t error = 0;

    /* Check the heap property */
    for (i = 0; i < size; i++) {
        if (LEFTCHILD(i) >= size) {
            break;
        }
        if (VECTOR(h->data)[LEFTCHILD(i)] > VECTOR(h->data)[i]) {
            error = 1; break;
        }
        if (RIGHTCHILD(i) >= size) {
            break;
        }
        if (VECTOR(h->data)[RIGHTCHILD(i)] > VECTOR(h->data)[i]) {
            error = 1; break;
        }
    }

    if (error) {
        IGRAPH_ERROR("Inconsistent heap", IGRAPH_EINTERNAL);
    }

    return 0;
}
