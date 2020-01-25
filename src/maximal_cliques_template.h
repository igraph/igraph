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
#define RESTYPE igraph_vector_ptr_t *res
#define RESNAME res
#define SUFFIX
#define RECORD do {                         \
        igraph_vector_t *cl=igraph_Calloc(1, igraph_vector_t);      \
        int j;                              \
        if (!cl) {                              \
            IGRAPH_ERROR("Cannot list maximal cliques", IGRAPH_ENOMEM);   \
        }                                   \
        IGRAPH_CHECK(igraph_vector_ptr_push_back(res, cl));             \
        IGRAPH_CHECK(igraph_vector_init(cl, clsize));                   \
        for (j=0; j<clsize; j++) { VECTOR(*cl)[j] = VECTOR(*R)[j]; }    \
    } while (0)
#define FINALLY do {                    \
        igraph_vector_ptr_clear(res);           \
        IGRAPH_FINALLY(igraph_i_maximal_cliques_free, res); \
    } while (0)
#define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++) {
#define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_COUNT
    #define RESTYPE igraph_integer_t *res
    #define RESNAME res
    #define SUFFIX _count
    #define RECORD (*res)++
    #define FINALLY *res=0;
    #define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++) {
    #define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_FILE
    #define RESTYPE FILE *res
    #define RESNAME res
    #define SUFFIX _file
    #define RECORD igraph_vector_int_fprint(R, res)
    #define FINALLY
    #define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++) {
    #define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_FULL
#define RESTYPE                 \
    igraph_vector_int_t *subset,            \
    igraph_vector_ptr_t *res,           \
    igraph_integer_t *no,           \
    FILE *outfile
#define RESNAME subset, res, no, outfile
#define SUFFIX _subset
#define RECORD do {                         \
        if (res) {                                \
            igraph_vector_t *cl=igraph_Calloc(1, igraph_vector_t);      \
            int j;                              \
            if (!cl) {                              \
                IGRAPH_ERROR("Cannot list maximal cliques", IGRAPH_ENOMEM);   \
            }                                   \
            IGRAPH_CHECK(igraph_vector_ptr_push_back(res, cl));             \
            IGRAPH_CHECK(igraph_vector_init(cl, clsize));                   \
            for (j=0; j<clsize; j++) { VECTOR(*cl)[j] = VECTOR(*R)[j]; }    \
        }                                 \
        if (no) { (*no)++; }                              \
        if (outfile) { igraph_vector_int_fprint(R, outfile); }        \
    } while (0)
#define FINALLY do {                        \
        if (res) {                            \
            igraph_vector_ptr_clear(res);               \
            IGRAPH_FINALLY(igraph_i_maximal_cliques_free_full, res);    \
        }                             \
        if (no) { *no=0; }                        \
    } while (0)
#define FOR_LOOP_OVER_VERTICES                  \
    nn= subset ? igraph_vector_int_size(subset) : no_of_nodes;    \
    for (ii=0; ii<nn; ii++) {
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
        igraph_vector_t *cl=igraph_Calloc(1, igraph_vector_t); \
        long j; \
        if (!cl) { \
            IGRAPH_ERROR("Cannot list maximal cliques", IGRAPH_ENOMEM); \
        } \
        IGRAPH_CHECK(igraph_vector_init(cl, clsize)); \
        for (j=0; j<clsize; j++) { VECTOR(*cl)[j] = VECTOR(*R)[j]; } \
        if (!cliquehandler_fn(cl, arg)) \
            return IGRAPH_STOP; \
    } while (0)
#define FINALLY
#define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++) {
#define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_HIST
#define RESTYPE igraph_vector_t *hist
#define RESNAME hist
#define SUFFIX _hist
#define RECORD do { \
        long hsize = igraph_vector_size(hist); \
        if (clsize > hsize) { \
            long hcapacity = igraph_vector_capacity(hist); \
            long j; \
            int err; \
            if (hcapacity < clsize && clsize < 2*hcapacity) \
                err = igraph_vector_reserve(hist, 2*hcapacity); \
            err = igraph_vector_resize(hist, clsize); \
            if (err != IGRAPH_SUCCESS) \
                IGRAPH_ERROR("Cannot count maximal cliques", IGRAPH_ENOMEM); \
            for (j=hsize; j < clsize; j++) \
                VECTOR(*hist)[j] = 0; \
        } \
        VECTOR(*hist)[clsize-1] += 1; \
    } while (0)
#define FINALLY \
    igraph_vector_clear(hist); \
    igraph_vector_reserve(hist, 50); /* initially reserve space for 50 elements */
#define FOR_LOOP_OVER_VERTICES for (i=0; i<no_of_nodes; i++) {
#define FOR_LOOP_OVER_VERTICES_PREPARE
#endif

#ifdef IGRAPH_MC_ORIG
void igraph_i_maximal_cliques_free(void *ptr) {
    igraph_vector_ptr_t *res = (igraph_vector_ptr_t*) ptr;
    int i, n = igraph_vector_ptr_size(res);
    for (i = 0; i < n; i++) {
        igraph_vector_t *v = VECTOR(*res)[i];
        if (v) {
            igraph_Free(v);
            igraph_vector_destroy(v);
        }
    }
    igraph_vector_ptr_clear(res);
}
#endif

#ifdef IGRAPH_MC_FULL
void igraph_i_maximal_cliques_free_full(void *ptr) {
    if (ptr) {
        igraph_vector_ptr_t *res = (igraph_vector_ptr_t*) ptr;
        int i, n = igraph_vector_ptr_size(res);
        for (i = 0; i < n; i++) {
            igraph_vector_t *v = VECTOR(*res)[i];
            if (v) {
                igraph_Free(v);
                igraph_vector_destroy(v);
            }
        }
        igraph_vector_ptr_clear(res);
    }
}
#endif

int FUNCTION(igraph_i_maximal_cliques_bk, SUFFIX)(
    igraph_vector_int_t *PX, int PS, int PE,
    int XS, int XE, int oldPS, int oldXE,
    igraph_vector_int_t *R,
    igraph_vector_int_t *pos,
    igraph_adjlist_t *adjlist,
    RESTYPE,
    igraph_vector_int_t *nextv,
    igraph_vector_int_t *H,
    int min_size, int max_size) {

    int err;

    igraph_vector_int_push_back(H, -1); /* boundary */

    if (PS > PE && XS > XE) {
        /* Found a maximum clique, report it */
        int clsize = igraph_vector_int_size(R);
        if (min_size <= clsize && (clsize <= max_size || max_size <= 0)) {
            RECORD;
        }
    } else if (PS <= PE) {
        /* Select a pivot element */
        int pivot, mynextv;
        igraph_i_maximal_cliques_select_pivot(PX, PS, PE, XS, XE, pos,
                                              adjlist, &pivot, nextv,
                                              oldPS, oldXE);
        while ((mynextv = igraph_vector_int_pop_back(nextv)) != -1) {
            int newPS, newXE;

            /* Going down, prepare */
            igraph_i_maximal_cliques_down(PX, PS, PE, XS, XE, pos, adjlist,
                                          mynextv, R, &newPS, &newXE);
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
                igraph_i_maximal_cliques_PX(PX, PS, &PE, &XS, XE, pos, adjlist,
                                            mynextv, H);
            }
        }
    }

    /* Putting back vertices from X to P, see notes in H */
    igraph_i_maximal_cliques_up(PX, PS, PE, XS, XE, pos, adjlist, R, H);

    return 0;
}

int FUNCTION(igraph_maximal_cliques, SUFFIX)(
    const igraph_t *graph,
    RESTYPE,
    igraph_integer_t min_size,
    igraph_integer_t max_size) {

    /* Implementation details. TODO */

    igraph_vector_int_t PX, R, H, pos, nextv;
    igraph_vector_t coreness, order;
    igraph_vector_int_t rank; /* TODO: this is not needed */
    int i, ii, nn, no_of_nodes = igraph_vcount(graph);
    igraph_adjlist_t adjlist, fulladjlist;
    igraph_real_t pgreset = round(no_of_nodes / 100.0), pg = pgreset, pgc = 0;
    int err;
    IGRAPH_UNUSED(nn);

    if (igraph_is_directed(graph)) {
        IGRAPH_WARNING("Edge directions are ignored for maximal clique "
                       "calculation");
    }

    igraph_vector_init(&order, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_destroy, &order);
    igraph_vector_int_init(&rank, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &rank);
    igraph_vector_init(&coreness, no_of_nodes);
    igraph_coreness(graph, &coreness, /*mode=*/ IGRAPH_ALL);
    IGRAPH_FINALLY(igraph_vector_destroy, &coreness);
    igraph_vector_qsort_ind(&coreness, &order, /*descending=*/ 0);
    for (ii = 0; ii < no_of_nodes; ii++) {
        int v = VECTOR(order)[ii];
        VECTOR(rank)[v] = ii;
    }

    igraph_vector_destroy(&coreness);
    IGRAPH_FINALLY_CLEAN(1);

    igraph_adjlist_init(graph, &adjlist, IGRAPH_ALL);

    igraph_adjlist_simplify(&adjlist);
    igraph_adjlist_init(graph, &fulladjlist, IGRAPH_ALL);
    IGRAPH_FINALLY(igraph_adjlist_destroy, &fulladjlist);
    igraph_adjlist_simplify(&fulladjlist);
    igraph_vector_int_init(&PX, 20);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &PX);
    igraph_vector_int_init(&R,  20);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &R);
    igraph_vector_int_init(&H, 100);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &H);
    igraph_vector_int_init(&pos, no_of_nodes);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &pos);
    igraph_vector_int_init(&nextv, 100);
    IGRAPH_FINALLY(igraph_vector_int_destroy, &nextv);

    FINALLY;

    FOR_LOOP_OVER_VERTICES
    int v;
    int vrank;
    igraph_vector_int_t *vneis;
    int vdeg;
    int Pptr, Xptr, PS, PE, XS, XE;
    int j;

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

    igraph_vector_int_resize(&PX, vdeg);
    igraph_vector_int_resize(&R, 1);
    igraph_vector_int_resize(&H, 1);
    igraph_vector_int_null(&pos); /* TODO: makes it quadratic? */
    igraph_vector_int_resize(&nextv, 1);

    VECTOR(H)[0] = -1;      /* marks the end of the recursion */
    VECTOR(nextv)[0] = -1;

    /* ================================================================*/
    /* P <- G(v[i]) intersect { v[i+1], ..., v[n-1] }
       X <- G(v[i]) intersect { v[0], ..., v[i-1] } */

    VECTOR(R)[0] = v;
    for (j = 0; j < vdeg; j++) {
        int vx = VECTOR(*vneis)[j];
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
    igraph_vector_int_update(igraph_adjlist_get(&adjlist, v),
                             igraph_adjlist_get(&fulladjlist, v));
    for (j = 0; j <= vdeg - 1; j++) {
        int vv = VECTOR(PX)[j];
        igraph_vector_int_t *fadj = igraph_adjlist_get(&fulladjlist, vv);
        igraph_vector_int_t *radj = igraph_adjlist_get(&adjlist, vv);
        int k, fn = igraph_vector_int_size(fadj);
        igraph_vector_int_clear(radj);
        for (k = 0; k < fn; k++) {
            int nei = VECTOR(*fadj)[k];
            int neipos = VECTOR(pos)[nei] - 1;
            if (neipos >= PS && neipos <= XE) {
                igraph_vector_int_push_back(radj, nei);
            }
        }
    }

    /* Reorder the adjacency lists, according to P and X. */
    igraph_i_maximal_cliques_reorder_adjlists(&PX, PS, PE, XS, XE, &pos,
            &adjlist);

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

igraph_vector_int_destroy(&nextv);
igraph_vector_int_destroy(&pos);
igraph_vector_int_destroy(&H);
igraph_vector_int_destroy(&R);
igraph_vector_int_destroy(&PX);
igraph_adjlist_destroy(&fulladjlist);
igraph_adjlist_destroy(&adjlist);
igraph_vector_int_destroy(&rank);
igraph_vector_destroy(&order);
IGRAPH_FINALLY_CLEAN(10); /* + res */

return 0;
}

#undef RESTYPE
#undef RESNAME
#undef SUFFIX
#undef RECORD
#undef FINALLY
#undef FOR_LOOP_OVER_VERTICES
#undef FOR_LOOP_OVER_VERTICES_PREPARE
