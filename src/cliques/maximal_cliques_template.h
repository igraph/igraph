/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

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

#ifdef IGRAPH_MC_ORIG
#define RESTYPE igraph_vector_int_list_t *res
#define RESNAME res
#define SUFFIX
#define RECORD do {                         \
        IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(res, R));     \
    } while (0)
#define PREPARE do {                    \
        igraph_vector_int_list_clear(res);           \
    } while (0)
#define CLEANUP
#define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++)
#define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_COUNT
    #define RESTYPE igraph_integer_t *res
    #define RESNAME res
    #define SUFFIX _count
    #define RECORD (*res)++
    #define PREPARE *res=0;
    #define CLEANUP
    #define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++)
    #define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_FILE
    #define RESTYPE FILE *res
    #define RESNAME res
    #define SUFFIX _file
    #define RECORD igraph_vector_int_fprint(R, res)
    #define PREPARE
    #define CLEANUP
    #define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++)
    #define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_FULL
#define RESTYPE                 \
    const igraph_vector_int_t *subset,            \
    igraph_vector_int_list_t *res,           \
    igraph_integer_t *no,           \
    FILE *outfile
#define RESNAME subset, res, no, outfile
#define SUFFIX _subset
#define RECORD do {                         \
        if (res) {                                \
            IGRAPH_CHECK(igraph_vector_int_list_push_back_copy(res, R));      \
        }                                 \
        if (no) { (*no)++; }                              \
        if (outfile) { igraph_vector_int_fprint(R, outfile); }        \
    } while (0)
#define PREPARE do {                        \
        if (res) {                                 \
            igraph_vector_int_list_clear(res);     \
        }                             \
        if (no) { *no=0; }                        \
    } while (0)
#define CLEANUP
#define FOR_LOOP_OVER_VERTICES                  \
    nn= subset ? igraph_vector_int_size(subset) : no_of_nodes;    \
    for (ii=0; ii<nn; ii++)
#define FOR_LOOP_OVER_VERTICES_PREPARE do {  \
        i= subset ? VECTOR(*subset)[ii] : ii;    \
    } while (0)
#endif

#ifdef IGRAPH_MC_CALLBACK
#define RESTYPE \
    igraph_clique_handler_t *cliquehandler_fn, \
    void *arg
#define RESNAME cliquehandler_fn, arg
#define SUFFIX _callback
#define RECORD do { \
        igraph_error_t cliquehandler_retval; \
        cliquehandler_retval = cliquehandler_fn(R, arg); \
        if (cliquehandler_retval == IGRAPH_STOP) { \
            return IGRAPH_STOP; \
        } else if (cliquehandler_retval) { \
            IGRAPH_ERROR("Cannot list maximal cliques", cliquehandler_retval); \
        } \
    } while (0)
#define PREPARE
#define CLEANUP
#define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++)
#define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_HIST
#define RESTYPE igraph_vector_t *hist
#define RESNAME hist
#define SUFFIX _hist
#define RECORD do { \
        igraph_integer_t hsize = igraph_vector_size(hist); \
        if (clsize > hsize) { \
            igraph_integer_t hcapacity = igraph_vector_capacity(hist); \
            igraph_integer_t j; \
            igraph_error_t err; \
            if (hcapacity < clsize && clsize < 2*hcapacity) \
                err = igraph_vector_reserve(hist, 2*hcapacity); \
            err = igraph_vector_resize(hist, clsize); \
            if (err != IGRAPH_SUCCESS) \
                IGRAPH_ERROR("Cannot count maximal cliques", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */ \
            for (j=hsize; j < clsize; j++) \
                VECTOR(*hist)[j] = 0; \
        } \
        VECTOR(*hist)[clsize-1] += 1; \
    } while (0)
#define PREPARE \
    igraph_vector_clear(hist); \
    IGRAPH_CHECK(igraph_vector_reserve(hist, 50)); /* initially reserve space for 50 elements */
#define CLEANUP
#define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++)
#define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

static igraph_error_t FUNCTION(igraph_i_maximal_cliques_bk, SUFFIX)(
    igraph_vector_int_t *PX, igraph_integer_t PS, igraph_integer_t PE,
    igraph_integer_t XS, igraph_integer_t XE, igraph_integer_t oldPS, igraph_integer_t oldXE,
    igraph_vector_int_t *R,
    igraph_vector_int_t *pos,
    igraph_adjlist_t *adjlist,
    RESTYPE,
    igraph_vector_int_t *nextv,
    igraph_vector_int_t *H,
    igraph_integer_t min_size, igraph_integer_t max_size) {

    igraph_error_t err;

    IGRAPH_CHECK(igraph_vector_int_push_back(H, -1)); /* boundary */

    if (PS > PE && XS > XE) {
        /* Found a maximum clique, report it */
        igraph_integer_t clsize = igraph_vector_int_size(R);
        if (min_size <= clsize && (clsize <= max_size || max_size <= 0)) {
            RECORD;
        }
    } else if (PS <= PE) {
        /* Select a pivot element */
        igraph_integer_t pivot, mynextv;
        IGRAPH_CHECK(igraph_i_maximal_cliques_select_pivot(
            PX, PS, PE, XS, XE, pos, adjlist, &pivot, nextv, oldPS, oldXE
        ));
        while ((mynextv = igraph_vector_int_pop_back(nextv)) != -1) {
            igraph_integer_t newPS, newXE;

            /* Going down, prepare */
            IGRAPH_CHECK(igraph_i_maximal_cliques_down(
                PX, PS, PE, XS, XE, pos, adjlist, mynextv, R, &newPS, &newXE
            ));
            /* Recursive call */
            err = FUNCTION(igraph_i_maximal_cliques_bk, SUFFIX)(
                      PX, newPS, PE, XS, newXE, PS, XE, R,
                      pos, adjlist, RESNAME, nextv, H,
                      min_size, max_size);

            if (err == IGRAPH_STOP) {
                return err;
            } else {
                IGRAPH_CHECK(err);
            }
            /* Putting v from P to X */
            if (igraph_vector_int_tail(nextv) != -1) {
                IGRAPH_CHECK(igraph_i_maximal_cliques_PX(
                    PX, PS, &PE, &XS, XE, pos, adjlist, mynextv, H
                ));
            }
        }
    }

    /* Putting back vertices from X to P, see notes in H */
    IGRAPH_CHECK(igraph_i_maximal_cliques_up(PX, PS, PE, XS, XE, pos, adjlist, R, H));

    return IGRAPH_SUCCESS;
}

igraph_error_t FUNCTION(igraph_maximal_cliques, SUFFIX)(
    const igraph_t *graph,
    RESTYPE,
    igraph_integer_t min_size,
    igraph_integer_t max_size) {

    /* Implementation details. TODO */

    igraph_vector_int_t PX, R, H, pos, nextv;
    igraph_vector_int_t coreness;
    igraph_vector_int_t order;
    igraph_vector_int_t rank; /* TODO: this is not needed */
    igraph_integer_t i, ii, nn, no_of_nodes = igraph_vcount(graph);
    igraph_adjlist_t adjlist, fulladjlist;
    igraph_real_t pgreset = round(no_of_nodes / 100.0), pg = pgreset, pgc = 0;
    igraph_error_t err;
    IGRAPH_UNUSED(nn);

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored for maximal clique "
                       "calculation");
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&order, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&rank, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&coreness, no_of_nodes);
    IGRAPH_CHECK(igraph_coreness(graph, &coreness, /*mode=*/ IGRAPH_ALL));
    IGRAPH_CHECK(igraph_vector_int_qsort_ind(&coreness, &order, IGRAPH_ASCENDING));
    for (ii = 0; ii < no_of_nodes; ii++) {
        igraph_integer_t v = VECTOR(order)[ii];
        VECTOR(rank)[v] = ii;
    }

    igraph_vector_int_destroy(&coreness);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);

    IGRAPH_CHECK(igraph_adjlist_init(graph, &fulladjlist, IGRAPH_ALL, IGRAPH_NO_LOOPS, IGRAPH_NO_MULTIPLE));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &fulladjlist);

    IGRAPH_VECTOR_INT_INIT_FINALLY(&PX, 20);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&R, 20);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&H, 100);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&pos, no_of_nodes);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&nextv, 100);

    PREPARE;

    FOR_LOOP_OVER_VERTICES {
        igraph_integer_t v;
        igraph_integer_t vrank;
        igraph_vector_int_t *vneis;
        igraph_integer_t vdeg;
        igraph_integer_t Pptr, Xptr, PS, PE, XS, XE;
        igraph_integer_t j;

        FOR_LOOP_OVER_VERTICES_PREPARE;

        v = VECTOR(order)[i];
        vrank = VECTOR(rank)[v];
        vneis = igraph_adjlist_get(&fulladjlist, v);
        vdeg = igraph_vector_int_size(vneis);
        Pptr = 0; Xptr = vdeg - 1; PS = 0; XE = vdeg - 1;

        pg--;
        if (pg <= 0) {
            IGRAPH_PROGRESS("Maximal cliques: ", pgc++, NULL);
            pg = pgreset;
        }

        IGRAPH_ALLOW_INTERRUPTION();

        IGRAPH_CHECK(igraph_vector_int_resize(&PX, vdeg));
        IGRAPH_CHECK(igraph_vector_int_resize(&R, 1));
        IGRAPH_CHECK(igraph_vector_int_resize(&H, 1));
        igraph_vector_int_null(&pos); /* TODO: makes it quadratic? */
        IGRAPH_CHECK(igraph_vector_int_resize(&nextv, 1));

        VECTOR(H)[0] = -1;      /* marks the end of the recursion */
        VECTOR(nextv)[0] = -1;

        /* ================================================================*/
        /* P <- G(v[i]) intersect { v[i+1], ..., v[n-1] }
           X <- G(v[i]) intersect { v[0], ..., v[i-1] } */

        VECTOR(R)[0] = v;
        for (j = 0; j < vdeg; j++) {
            igraph_integer_t vx = VECTOR(*vneis)[j];
            if (VECTOR(rank)[vx] > vrank) {
                VECTOR(PX)[Pptr] = vx;
                VECTOR(pos)[vx] = Pptr + 1;
                Pptr++;
            } else if (VECTOR(rank)[vx] < vrank) {
                VECTOR(PX)[Xptr] = vx;
                VECTOR(pos)[vx] = Xptr + 1;
                Xptr--;
            }
        }

        PE = Pptr - 1; XS = Xptr + 1; /* end of P, start of X in PX */

        /* Create an adjacency list that is specific to the
           v vertex. It only contains 'v' and its neighbors. Moreover, we
           only deal with the vertices in P and X (and R). */
        IGRAPH_CHECK(igraph_vector_int_update(
            igraph_adjlist_get(&adjlist, v),
            igraph_adjlist_get(&fulladjlist, v)
        ));
        for (j = 0; j <= vdeg - 1; j++) {
            igraph_integer_t vv = VECTOR(PX)[j];
            igraph_vector_int_t *fadj = igraph_adjlist_get(&fulladjlist, vv);
            igraph_vector_int_t *radj = igraph_adjlist_get(&adjlist, vv);
            igraph_integer_t k, fn = igraph_vector_int_size(fadj);
            igraph_vector_int_clear(radj);
            for (k = 0; k < fn; k++) {
                igraph_integer_t nei = VECTOR(*fadj)[k];
                igraph_integer_t neipos = VECTOR(pos)[nei] - 1;
                if (neipos >= PS && neipos <= XE) {
                    IGRAPH_CHECK(igraph_vector_int_push_back(radj, nei));
                }
            }
        }

        /* Reorder the adjacency lists, according to P and X. */
        IGRAPH_CHECK(igraph_i_maximal_cliques_reorder_adjlists(
            &PX, PS, PE, XS, XE, &pos, &adjlist
        ));

        err = FUNCTION(igraph_i_maximal_cliques_bk, SUFFIX)(
                &PX, PS, PE, XS, XE, PS, XE, &R, &pos,
                &adjlist, RESNAME, &nextv, &H, min_size,
                max_size);
        if (err == IGRAPH_STOP) {
            break;
        } else {
            IGRAPH_CHECK(err);
        }
    }

    IGRAPH_PROGRESS("Maximal cliques: ", 100.0, NULL);

    CLEANUP;

    igraph_vector_int_destroy(&nextv);
    igraph_vector_int_destroy(&pos);
    igraph_vector_int_destroy(&H);
    igraph_vector_int_destroy(&R);
    igraph_vector_int_destroy(&PX);
    igraph_adjlist_destroy(&fulladjlist);
    igraph_adjlist_destroy(&adjlist);
    igraph_vector_int_destroy(&rank);
    igraph_vector_int_destroy(&order);
    IGRAPH_FINALLY_CLEAN(9);

    return IGRAPH_SUCCESS;
}

#undef RESTYPE
#undef RESNAME
#undef SUFFIX
#undef RECORD
#undef PREPARE
#undef CLEANUP
#undef FOR_LOOP_OVER_VERTICES
#undef FOR_LOOP_OVER_VERTICES_PREPARE
