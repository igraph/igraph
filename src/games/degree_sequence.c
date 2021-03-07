/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph library.
   Copyright (C) 2003-2021 The igraph development team

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

#include "igraph_games.h"

#include "igraph_adjlist.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_graphicality.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_vector_ptr.h"

#include "core/interruption.h"
#include "core/set.h"

static int igraph_i_degree_sequence_game_simple(igraph_t *graph,
                                       const igraph_vector_t *out_seq,
                                       const igraph_vector_t *in_seq) {

    long int outsum = 0, insum = 0;
    igraph_bool_t directed = (in_seq != 0 && igraph_vector_size(in_seq) != 0);
    igraph_bool_t degseq_ok;
    long int no_of_nodes, no_of_edges;
    long int *bag1 = 0, *bag2 = 0;
    long int bagp1 = 0, bagp2 = 0;
    igraph_vector_t edges = IGRAPH_VECTOR_NULL;
    long int i, j;

    IGRAPH_CHECK(igraph_is_graphical(out_seq, in_seq, IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR(in_seq ? "No directed graph can realize the given degree sequences" :
                     "No undirected graph can realize the given degree sequence", IGRAPH_EINVAL);
    }

    outsum = (long int) igraph_vector_sum(out_seq);
    if (directed) {
        insum = (long int) igraph_vector_sum(in_seq);
    }

    no_of_nodes = igraph_vector_size(out_seq);
    no_of_edges = directed ? outsum : outsum / 2;

    bag1 = IGRAPH_CALLOC(outsum, long int);
    if (bag1 == 0) {
        IGRAPH_ERROR("degree sequence game (simple)", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, bag1);

    for (i = 0; i < no_of_nodes; i++) {
        for (j = 0; j < VECTOR(*out_seq)[i]; j++) {
            bag1[bagp1++] = i;
        }
    }
    if (directed) {
        bag2 = IGRAPH_CALLOC(insum, long int);
        if (bag2 == 0) {
            IGRAPH_ERROR("degree sequence game (simple)", IGRAPH_ENOMEM);
        }
        IGRAPH_FINALLY(igraph_free, bag2);
        for (i = 0; i < no_of_nodes; i++) {
            for (j = 0; j < VECTOR(*in_seq)[i]; j++) {
                bag2[bagp2++] = i;
            }
        }
    }

    IGRAPH_VECTOR_INIT_FINALLY(&edges, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&edges, no_of_edges * 2));

    RNG_BEGIN();

    if (directed) {
        for (i = 0; i < no_of_edges; i++) {
            long int from = RNG_INTEGER(0, bagp1 - 1);
            long int to = RNG_INTEGER(0, bagp2 - 1);
            igraph_vector_push_back(&edges, bag1[from]); /* safe, already reserved */
            igraph_vector_push_back(&edges, bag2[to]);   /* ditto */
            bag1[from] = bag1[bagp1 - 1];
            bag2[to] = bag2[bagp2 - 1];
            bagp1--; bagp2--;
        }
    } else {
        for (i = 0; i < no_of_edges; i++) {
            long int from = RNG_INTEGER(0, bagp1 - 1);
            long int to;
            igraph_vector_push_back(&edges, bag1[from]); /* safe, already reserved */
            bag1[from] = bag1[bagp1 - 1];
            bagp1--;
            to = RNG_INTEGER(0, bagp1 - 1);
            igraph_vector_push_back(&edges, bag1[to]);   /* ditto */
            bag1[to] = bag1[bagp1 - 1];
            bagp1--;
        }
    }

    RNG_END();

    IGRAPH_FREE(bag1);
    IGRAPH_FINALLY_CLEAN(1);
    if (directed) {
        IGRAPH_FREE(bag2);
        IGRAPH_FINALLY_CLEAN(1);
    }

    IGRAPH_CHECK(igraph_create(graph, &edges, (igraph_integer_t) no_of_nodes,
                               directed));
    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_degree_sequence_game_no_multiple_undirected(
    igraph_t *graph, const igraph_vector_t *seq) {

    igraph_vector_t stubs = IGRAPH_VECTOR_NULL;
    igraph_vector_int_t *neis;
    igraph_vector_t residual_degrees = IGRAPH_VECTOR_NULL;
    igraph_set_t incomplete_vertices;
    igraph_adjlist_t al;
    igraph_bool_t finished, failed;
    igraph_integer_t from, to, dummy;
    long int i, j, k;
    long int no_of_nodes, outsum = 0;
    igraph_bool_t degseq_ok;

    IGRAPH_CHECK(igraph_is_graphical(seq, 0, IGRAPH_SIMPLE_SW, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR("No simple undirected graph can realize the given degree sequence",
                     IGRAPH_EINVAL);
    }

    outsum = (long int) igraph_vector_sum(seq);
    no_of_nodes = igraph_vector_size(seq);

    /* Allocate required data structures */
    IGRAPH_CHECK(igraph_adjlist_init_empty(&al, (igraph_integer_t) no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    IGRAPH_VECTOR_INIT_FINALLY(&stubs, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&stubs, outsum));
    IGRAPH_VECTOR_INIT_FINALLY(&residual_degrees, no_of_nodes);
    IGRAPH_CHECK(igraph_set_init(&incomplete_vertices, 0));
    IGRAPH_FINALLY(igraph_set_destroy, &incomplete_vertices);

    /* Start the RNG */
    RNG_BEGIN();

    /* Outer loop; this will try to construct a graph several times from scratch
     * until it finally succeeds. */
    finished = 0;
    while (!finished) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* Be optimistic :) */
        failed = 0;

        /* Clear the adjacency list to get rid of the previous attempt (if any) */
        igraph_adjlist_clear(&al);

        /* Initialize the residual degrees from the degree sequence */
        IGRAPH_CHECK(igraph_vector_update(&residual_degrees, seq));

        /* While there are some unconnected stubs left... */
        while (!finished && !failed) {
            /* Construct the initial stub vector */
            igraph_vector_clear(&stubs);
            for (i = 0; i < no_of_nodes; i++) {
                for (j = 0; j < VECTOR(residual_degrees)[i]; j++) {
                    igraph_vector_push_back(&stubs, i);
                }
            }

            /* Clear the skipped stub counters and the set of incomplete vertices */
            igraph_vector_null(&residual_degrees);
            igraph_set_clear(&incomplete_vertices);

            /* Shuffle the stubs in-place */
            igraph_vector_shuffle(&stubs);

            /* Connect the stubs where possible */
            k = igraph_vector_size(&stubs);
            for (i = 0; i < k; ) {
                from = (igraph_integer_t) VECTOR(stubs)[i++];
                to = (igraph_integer_t) VECTOR(stubs)[i++];

                if (from > to) {
                    dummy = from; from = to; to = dummy;
                }

                neis = igraph_adjlist_get(&al, from);
                if (from == to || igraph_vector_int_binsearch(neis, to, &j)) {
                    /* Edge exists already */
                    VECTOR(residual_degrees)[from]++;
                    VECTOR(residual_degrees)[to]++;
                    IGRAPH_CHECK(igraph_set_add(&incomplete_vertices, from));
                    IGRAPH_CHECK(igraph_set_add(&incomplete_vertices, to));
                } else {
                    /* Insert the edge */
                    IGRAPH_CHECK(igraph_vector_int_insert(neis, j, to));
                }
            }

            finished = igraph_set_empty(&incomplete_vertices);

            if (!finished) {
                /* We are not done yet; check if the remaining stubs are feasible. This
                 * is done by enumerating all possible pairs and checking whether at
                 * least one feasible pair is found. */
                i = 0;
                failed = 1;
                while (failed && igraph_set_iterate(&incomplete_vertices, &i, &from)) {
                    j = 0;
                    while (igraph_set_iterate(&incomplete_vertices, &j, &to)) {
                        if (from == to) {
                            /* This is used to ensure that each pair is checked once only */
                            break;
                        }
                        if (from > to) {
                            dummy = from; from = to; to = dummy;
                        }
                        neis = igraph_adjlist_get(&al, from);
                        if (!igraph_vector_int_binsearch(neis, to, 0)) {
                            /* Found a suitable pair, so we can continue */
                            failed = 0;
                            break;
                        }
                    }
                }
            }
        }
    }

    /* Finish the RNG */
    RNG_END();

    /* Clean up */
    igraph_set_destroy(&incomplete_vertices);
    igraph_vector_destroy(&residual_degrees);
    igraph_vector_destroy(&stubs);
    IGRAPH_FINALLY_CLEAN(3);

    /* Create the graph. We cannot use IGRAPH_ALL here for undirected graphs
     * because we did not add edges in both directions in the adjacency list.
     * We will use igraph_to_undirected in an extra step. */
    IGRAPH_CHECK(igraph_adjlist(graph, &al, IGRAPH_OUT, 1));
    IGRAPH_CHECK(igraph_to_undirected(graph, IGRAPH_TO_UNDIRECTED_EACH, 0));

    /* Clear the adjacency list */
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

static int igraph_i_degree_sequence_game_no_multiple_directed(igraph_t *graph,
        const igraph_vector_t *out_seq, const igraph_vector_t *in_seq) {
    igraph_adjlist_t al;
    igraph_bool_t deg_seq_ok, failed, finished;
    igraph_vector_t in_stubs = IGRAPH_VECTOR_NULL;
    igraph_vector_t out_stubs = IGRAPH_VECTOR_NULL;
    igraph_vector_int_t *neis;
    igraph_vector_t residual_in_degrees = IGRAPH_VECTOR_NULL;
    igraph_vector_t residual_out_degrees = IGRAPH_VECTOR_NULL;
    igraph_set_t incomplete_in_vertices;
    igraph_set_t incomplete_out_vertices;
    igraph_integer_t from, to;
    long int i, j, k;
    long int no_of_nodes, outsum;

    IGRAPH_CHECK(igraph_is_graphical(out_seq, in_seq, IGRAPH_SIMPLE_SW, &deg_seq_ok));
    if (!deg_seq_ok) {
        IGRAPH_ERROR("No simple directed graph can realize the given degree sequence",
                     IGRAPH_EINVAL);
    }

    outsum = (long int) igraph_vector_sum(out_seq);
    no_of_nodes = igraph_vector_size(out_seq);

    /* Allocate required data structures */
    IGRAPH_CHECK(igraph_adjlist_init_empty(&al, (igraph_integer_t) no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &al);
    IGRAPH_VECTOR_INIT_FINALLY(&out_stubs, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&out_stubs, outsum));
    IGRAPH_VECTOR_INIT_FINALLY(&in_stubs, 0);
    IGRAPH_CHECK(igraph_vector_reserve(&in_stubs, outsum));
    IGRAPH_VECTOR_INIT_FINALLY(&residual_out_degrees, no_of_nodes);
    IGRAPH_VECTOR_INIT_FINALLY(&residual_in_degrees, no_of_nodes);
    IGRAPH_CHECK(igraph_set_init(&incomplete_out_vertices, 0));
    IGRAPH_FINALLY(igraph_set_destroy, &incomplete_out_vertices);
    IGRAPH_CHECK(igraph_set_init(&incomplete_in_vertices, 0));
    IGRAPH_FINALLY(igraph_set_destroy, &incomplete_in_vertices);

    /* Start the RNG */
    RNG_BEGIN();

    /* Outer loop; this will try to construct a graph several times from scratch
     * until it finally succeeds. */
    finished = 0;
    while (!finished) {
        IGRAPH_ALLOW_INTERRUPTION();

        /* Be optimistic :) */
        failed = 0;

        /* Clear the adjacency list to get rid of the previous attempt (if any) */
        igraph_adjlist_clear(&al);

        /* Initialize the residual degrees from the degree sequences */
        IGRAPH_CHECK(igraph_vector_update(&residual_out_degrees, out_seq));
        IGRAPH_CHECK(igraph_vector_update(&residual_in_degrees, in_seq));

        /* While there are some unconnected stubs left... */
        while (!finished && !failed) {
            /* Construct the initial stub vectors */
            igraph_vector_clear(&out_stubs);
            igraph_vector_clear(&in_stubs);
            for (i = 0; i < no_of_nodes; i++) {
                for (j = 0; j < VECTOR(residual_out_degrees)[i]; j++) {
                    igraph_vector_push_back(&out_stubs, i);
                }
                for (j = 0; j < VECTOR(residual_in_degrees)[i]; j++) {
                    igraph_vector_push_back(&in_stubs, i);
                }
            }

            /* Clear the skipped stub counters and the set of incomplete vertices */
            igraph_vector_null(&residual_out_degrees);
            igraph_vector_null(&residual_in_degrees);
            igraph_set_clear(&incomplete_out_vertices);
            igraph_set_clear(&incomplete_in_vertices);

            /* Shuffle the out-stubs in-place */
            igraph_vector_shuffle(&out_stubs);

            /* Connect the stubs where possible */
            k = igraph_vector_size(&out_stubs);
            for (i = 0; i < k; i++) {
                from = (igraph_integer_t) VECTOR(out_stubs)[i];
                to = (igraph_integer_t) VECTOR(in_stubs)[i];

                neis = igraph_adjlist_get(&al, from);
                if (from == to || igraph_vector_int_binsearch(neis, to, &j)) {
                    /* Edge exists already */
                    VECTOR(residual_out_degrees)[from]++;
                    VECTOR(residual_in_degrees)[to]++;
                    IGRAPH_CHECK(igraph_set_add(&incomplete_out_vertices, from));
                    IGRAPH_CHECK(igraph_set_add(&incomplete_in_vertices, to));
                } else {
                    /* Insert the edge */
                    IGRAPH_CHECK(igraph_vector_int_insert(neis, j, to));
                }
            }

            /* Are we finished? */
            finished = igraph_set_empty(&incomplete_out_vertices);

            if (!finished) {
                /* We are not done yet; check if the remaining stubs are feasible. This
                 * is done by enumerating all possible pairs and checking whether at
                 * least one feasible pair is found. */
                i = 0;
                failed = 1;
                while (failed && igraph_set_iterate(&incomplete_out_vertices, &i, &from)) {
                    j = 0;
                    while (igraph_set_iterate(&incomplete_in_vertices, &j, &to)) {
                        neis = igraph_adjlist_get(&al, from);
                        if (from != to && !igraph_vector_int_binsearch(neis, to, 0)) {
                            /* Found a suitable pair, so we can continue */
                            failed = 0;
                            break;
                        }
                    }
                }
            }
        }
    }

    /* Finish the RNG */
    RNG_END();

    /* Clean up */
    igraph_set_destroy(&incomplete_in_vertices);
    igraph_set_destroy(&incomplete_out_vertices);
    igraph_vector_destroy(&residual_in_degrees);
    igraph_vector_destroy(&residual_out_degrees);
    igraph_vector_destroy(&in_stubs);
    igraph_vector_destroy(&out_stubs);
    IGRAPH_FINALLY_CLEAN(6);

    /* Create the graph */
    IGRAPH_CHECK(igraph_adjlist(graph, &al, IGRAPH_OUT, 1));

    /* Clear the adjacency list */
    igraph_adjlist_destroy(&al);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/* swap two elements of a vector_int */
#define SWAP_INT_ELEM(vec, i, j) \
    { \
        igraph_integer_t temp; \
        temp = VECTOR(vec)[i]; \
        VECTOR(vec)[i] = VECTOR(vec)[j]; \
        VECTOR(vec)[j] = temp; \
    }

static int igraph_i_degree_sequence_game_no_multiple_undirected_uniform(igraph_t *graph, const igraph_vector_t *degseq) {
    igraph_vector_int_t stubs;
    igraph_vector_t edges;
    igraph_bool_t degseq_ok;
    igraph_vector_ptr_t adjlist;
    long i, j;
    long vcount, ecount, stub_count;

    IGRAPH_CHECK(igraph_is_graphical(degseq, NULL, IGRAPH_SIMPLE_SW, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR("No simple undirected graph can realize the given degree sequence", IGRAPH_EINVAL);
    }

    stub_count = (long) igraph_vector_sum(degseq);
    ecount = stub_count / 2;
    vcount = igraph_vector_size(degseq);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&stubs, stub_count);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, stub_count);

    /* Fill stubs vector. */
    {
        long k = 0;
        for (i = 0; i < vcount; ++i) {
            long deg = (long) VECTOR(*degseq)[i];
            for (j = 0; j < deg; ++j) {
                VECTOR(stubs)[k++] = i;
            }
        }
    }

    /* Build an adjacency list in terms of sets; used to check for multi-edges. */
    IGRAPH_CHECK(igraph_vector_ptr_init(&adjlist, vcount));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&adjlist, igraph_set_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &adjlist);
    for (i = 0; i < vcount; ++i) {
        igraph_set_t *set = IGRAPH_CALLOC(1, igraph_set_t);
        if (! set) {
            IGRAPH_ERROR("Out of memory", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_set_init(set, 0));
        VECTOR(adjlist)[i] = set;
        IGRAPH_CHECK(igraph_set_reserve(set, (long) VECTOR(*degseq)[i]));
    }

    RNG_BEGIN();

    for (;;) {
        igraph_bool_t success = 1;

        /* Shuffle stubs vector with Fisher-Yates and check for self-loops and multi-edges as we go. */
        for (i = 0; i < ecount; ++i) {
            long k, from, to;

            k = RNG_INTEGER(2*i, stub_count-1);
            SWAP_INT_ELEM(stubs, 2*i, k);

            k = RNG_INTEGER(2*i+1, stub_count-1);
            SWAP_INT_ELEM(stubs, 2*i+1, k);

            from = VECTOR(stubs)[2*i];
            to   = VECTOR(stubs)[2*i+1];

            /* self-loop, fail */
            if (from == to) {
                success = 0;
                break;
            }

            /* multi-edge, fail */
            if (igraph_set_contains((igraph_set_t *) VECTOR(adjlist)[to], from)) {
                success = 0;
                break;
            }

            /* sets are already reserved */
            igraph_set_add((igraph_set_t *) VECTOR(adjlist)[to], from);
            igraph_set_add((igraph_set_t *) VECTOR(adjlist)[from], to);

            /* register edge */
            VECTOR(edges)[2 * i]   = from;
            VECTOR(edges)[2 * i + 1] = to;
        }

        if (success) {
            break;
        }

        /* Clear adjacency list. */
        for (j = 0; j < vcount; ++j) {
            igraph_set_clear((igraph_set_t *) VECTOR(adjlist)[j]);
        }

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    igraph_vector_ptr_destroy_all(&adjlist);
    igraph_vector_int_destroy(&stubs);
    IGRAPH_FINALLY_CLEAN(2);

    IGRAPH_CHECK(igraph_create(graph, &edges, vcount, /* directed = */ 0));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


static int igraph_i_degree_sequence_game_no_multiple_directed_uniform(
    igraph_t *graph, const igraph_vector_t *out_deg, const igraph_vector_t *in_deg) {
    igraph_vector_int_t out_stubs, in_stubs;
    igraph_vector_t edges;
    igraph_bool_t degseq_ok;
    igraph_vector_ptr_t adjlist;
    long i, j;
    long vcount, ecount;

    IGRAPH_CHECK(igraph_is_graphical(out_deg, in_deg, IGRAPH_SIMPLE_SW, &degseq_ok));
    if (!degseq_ok) {
        IGRAPH_ERROR("No simple directed graph can realize the given degree sequence", IGRAPH_EINVAL);
    }

    ecount = (long) igraph_vector_sum(out_deg);
    vcount = igraph_vector_size(out_deg);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&out_stubs, ecount);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&in_stubs, ecount);
    IGRAPH_VECTOR_INIT_FINALLY(&edges, 2 * ecount);

    /* Fill in- and out-stubs vectors. */
    {
        long k = 0, l = 0;
        for (i = 0; i < vcount; ++i) {
            long dout, din;

            dout = (long) VECTOR(*out_deg)[i];
            for (j = 0; j < dout; ++j) {
                VECTOR(out_stubs)[k++] = i;
            }

            din  = (long) VECTOR(*in_deg)[i];
            for (j = 0; j < din; ++j) {
                VECTOR(in_stubs)[l++] = i;
            }
        }
    }

    /* Build an adjacency list in terms of sets; used to check for multi-edges. */
    IGRAPH_CHECK(igraph_vector_ptr_init(&adjlist, vcount));
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&adjlist, igraph_set_destroy);
    IGRAPH_FINALLY(igraph_vector_ptr_destroy_all, &adjlist);
    for (i = 0; i < vcount; ++i) {
        igraph_set_t *set = IGRAPH_CALLOC(1, igraph_set_t);
        if (! set) {
            IGRAPH_ERROR("Out of memory", IGRAPH_ENOMEM);
        }
        IGRAPH_CHECK(igraph_set_init(set, 0));
        VECTOR(adjlist)[i] = set;
        IGRAPH_CHECK(igraph_set_reserve(set, (long) VECTOR(*out_deg)[i]));
    }

    RNG_BEGIN();

    for (;;) {
        igraph_bool_t success = 1;

        /* Shuffle out-stubs vector with Fisher-Yates and check for self-loops and multi-edges as we go. */
        for (i = 0; i < ecount; ++i) {
            long k, from, to;
            igraph_set_t *set;

            k = RNG_INTEGER(i, ecount-1);
            SWAP_INT_ELEM(out_stubs, i, k);

            from = VECTOR(out_stubs)[i];
            to   = VECTOR(in_stubs)[i];

            /* self-loop, fail */
            if (to == from) {
                success = 0;
                break;
            }

            /* multi-edge, fail */
            set = (igraph_set_t *) VECTOR(adjlist)[from];
            if (igraph_set_contains(set, to)) {
                success = 0;
                break;
            }

            /* sets are already reserved */
            igraph_set_add(set, to);

            /* register edge */
            VECTOR(edges)[2 * i]   = from;
            VECTOR(edges)[2 * i + 1] = to;
        }

        if (success) {
            break;
        }

        /* Clear adjacency list. */
        for (j = 0; j < vcount; ++j) {
            igraph_set_clear((igraph_set_t *) VECTOR(adjlist)[j]);
        }

        IGRAPH_ALLOW_INTERRUPTION();
    }

    RNG_END();

    igraph_vector_ptr_destroy_all(&adjlist);
    igraph_vector_int_destroy(&out_stubs);
    igraph_vector_int_destroy(&in_stubs);
    IGRAPH_FINALLY_CLEAN(3);

    IGRAPH_CHECK(igraph_create(graph, &edges, vcount, /* directed = */ 1));

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

#undef SWAP_INT_ELEM


/* This is in gengraph_mr-connected.cpp */

int igraph_degree_sequence_game_vl(igraph_t *graph,
                                   const igraph_vector_t *out_seq,
                                   const igraph_vector_t *in_seq);

/**
 * \ingroup generators
 * \function igraph_degree_sequence_game
 * \brief Generates a random graph with a given degree sequence.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param out_deg The degree sequence for an undirected graph (if
 *        \p in_seq is \c NULL or of length zero), or the out-degree
 *        sequence of a directed graph (if \p in_deq is not
 *        of length zero).
 * \param in_deg It is either a zero-length vector or
 *        \c NULL (if an undirected
 *        graph is generated), or the in-degree sequence.
 * \param method The method to generate the graph. Possible values:
 *        \clist
 *          \cli IGRAPH_DEGSEQ_SIMPLE
 *          This method implements the configuration model.
 *          For undirected graphs, it puts all vertex IDs in a bag
 *          such that the multiplicity of a vertex in the bag is the same as
 *          its degree. Then it draws pairs from the bag until the bag becomes
 *          empty. This method may generate both loop (self) edges and multiple
 *          edges. For directed graphs, the algorithm is basically the same,
 *          but two separate bags are used for the in- and out-degrees.
 *          Undirected graphs are generated with probability proportional to
 *          <code>(\prod_{i&lt;j} A_{ij} ! \prod_i A_{ii} !!)^{-1}</code>,
 *          where \c A denotes the adjacency matrix and <code>!!</code> denotes
 *          the double factorial. Here \c A is assumed to have twice the number of
 *          self-loops on its diagonal.
 *          The corresponding  expression for directed graphs is
 *          <code>(\prod_{i,j} A_{ij}!)^{-1}</code>.
 *          Thus the probability of all simple graphs (which only have 0s and 1s
 *          in the adjacency matrix) is the same, while that of
 *          non-simple ones depends on their edge and self-loop multiplicities.
 *          \cli IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE
 *          This method generates simple graphs.
 *          It is similar to \c IGRAPH_DEGSEQ_SIMPLE
 *          but tries to avoid multiple and loop edges and restarts the
 *          generation from scratch if it gets stuck. It can generate all simple
 *          realizations of a degree sequence, but it is not guaranteed
 *          to sample them uniformly. This method is relatively fast and it will
 *          eventually succeed if the provided degree sequence is graphical,
 *          but there is no upper bound on the number of iterations.
 *          \cli IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM
 *          This method is identical to \c IGRAPH_DEGSEQ_SIMPLE, but if the
 *          generated graph is not simple, it rejects it and re-starts the
 *          generation. It generates all simple graphs with the same probability.
 *          \cli IGRAPH_DEGSEQ_VL
 *          This method samples undirected \em connected graphs approximately
 *          uniformly. It is a Monte Carlo method based on degree-preserving
 *          edge swaps.
 *          This generator should be favoured if undirected and connected
 *          graphs are to be generated and execution time is not a concern.
 *          igraph uses the original implementation of Fabien Viger; for the algorithm,
 *          see https://www-complexnetworks.lip6.fr/~latapy/FV/generation.html
 *          and the paper https://arxiv.org/abs/cs/0502085
 *        \endclist
 * \return Error code:
 *          \c IGRAPH_ENOMEM: there is not enough
 *           memory to perform the operation.
 *          \c IGRAPH_EINVAL: invalid method parameter, or
 *           invalid in- and/or out-degree vectors. The degree vectors
 *           should be non-negative, \p out_deg should sum
 *           up to an even integer for undirected graphs; the length
 *           and sum of \p out_deg and
 *           \p in_deg
 *           should match for directed graphs.
 *
 * Time complexity: O(|V|+|E|), the number of vertices plus the number of edges
 *                  for \c IGRAPH_DEGSEQ_SIMPLE. The time complexity of the
 *                  other modes is not known.
 *
 * \sa \ref igraph_barabasi_game(), \ref igraph_erdos_renyi_game(),
 *     \ref igraph_is_graphical()
 *
 * \example examples/simple/igraph_degree_sequence_game.c
 */

int igraph_degree_sequence_game(igraph_t *graph, const igraph_vector_t *out_deg,
                                const igraph_vector_t *in_deg,
                                igraph_degseq_t method) {
    if (in_deg && igraph_vector_empty(in_deg) && !igraph_vector_empty(out_deg)) {
        in_deg = 0;
    }

    switch (method) {
    case IGRAPH_DEGSEQ_SIMPLE:
        return igraph_i_degree_sequence_game_simple(graph, out_deg, in_deg);

    case IGRAPH_DEGSEQ_VL:
        return igraph_degree_sequence_game_vl(graph, out_deg, in_deg);

    case IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE:
        if (in_deg == 0) {
            return igraph_i_degree_sequence_game_no_multiple_undirected(graph, out_deg);
        } else {
            return igraph_i_degree_sequence_game_no_multiple_directed(graph, out_deg, in_deg);
        }

    case IGRAPH_DEGSEQ_SIMPLE_NO_MULTIPLE_UNIFORM:
        if (in_deg == 0) {
            return igraph_i_degree_sequence_game_no_multiple_undirected_uniform(graph, out_deg);
        } else {
            return igraph_i_degree_sequence_game_no_multiple_directed_uniform(graph, out_deg, in_deg);
        }

    default:
        IGRAPH_ERROR("Invalid degree sequence game method", IGRAPH_EINVAL);
    }
}
