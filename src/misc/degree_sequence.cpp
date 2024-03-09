/*
  IGraph library.
  Constructing realizations of degree sequences and bi-degree sequences.
  Copyright (C) 2018-2024  The igraph development team <igraph@igraph.org>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "igraph_constructors.h"

#include "core/exceptions.h"
#include "math/safe_intop.h"

#include <vector>
#include <list>
#include <algorithm>
#include <utility>

#define IGRAPH_I_MULTI_EDGES_SW 0x02 /* 010, more than one edge allowed between distinct vertices */
#define IGRAPH_I_MULTI_LOOPS_SW 0x04 /* 100, more than one self-loop allowed on the same vertex   */

/******************************/
/***** Helper constructs ******/
/******************************/

// (vertex, degree) pair
struct vd_pair {
    igraph_integer_t vertex;
    igraph_integer_t degree;

    vd_pair(igraph_integer_t vertex, igraph_integer_t degree) : vertex(vertex), degree(degree) {}
};

// (indegree, outdegree)
typedef std::pair<igraph_integer_t, igraph_integer_t> bidegree;

// (vertex, bidegree) pair
struct vbd_pair {
    igraph_integer_t vertex;
    bidegree degree;

    vbd_pair(igraph_integer_t vertex, bidegree degree) : vertex(vertex), degree(degree) {}
};

// Comparison function for vertex-degree pairs.
// Also used for lexicographic sorting of bi-degrees.
template<typename T> inline bool degree_greater(const T &a, const T &b) {
    return a.degree > b.degree;
}

template<typename T> inline bool degree_less(const T &a, const T &b) {
    return a.degree < b.degree;
}


/*************************************/
/***** Undirected simple graphs ******/
/*************************************/

// Generate simple undirected realization as edge-list.
// If largest=true, always choose the vertex with the largest remaining degree to connect up next.
// Otherwise, always choose the one with the smallest remaining degree.
static igraph_error_t igraph_i_havel_hakimi(const igraph_vector_int_t *deg, igraph_vector_int_t *edges, bool largest) {
    igraph_integer_t n = igraph_vector_int_size(deg);

    igraph_integer_t ec = 0; // number of edges added so far

    std::vector<vd_pair> vertices;
    vertices.reserve(n);
    for (igraph_integer_t i = 0; i < n; ++i) {
        vertices.push_back(vd_pair(i, VECTOR(*deg)[i]));
    }

    while (! vertices.empty()) {
        if (largest) {
            std::stable_sort(vertices.begin(), vertices.end(), degree_less<vd_pair>);
        } else {
            std::stable_sort(vertices.begin(), vertices.end(), degree_greater<vd_pair>);
        }

        // take the next vertex to be connected up
        vd_pair vd = vertices.back();
        vertices.pop_back();

        if (vd.degree == 0) {
            continue;
        }

        if (vertices.size() < size_t(vd.degree)) {
            goto fail;
        }

        if (largest) {
            for (igraph_integer_t i = 0; i < vd.degree; ++i) {
                if (--(vertices[vertices.size() - 1 - i].degree) < 0) {
                    goto fail;
                }

                VECTOR(*edges)[2 * (ec + i)] = vd.vertex;
                VECTOR(*edges)[2 * (ec + i) + 1] = vertices[vertices.size() - 1 - i].vertex;
            }
        } else {
            // this loop can only be reached if all zero-degree nodes have already been removed
            // therefore decrementing remaining degrees is safe
            for (igraph_integer_t i = 0; i < vd.degree; ++i) {
                vertices[i].degree--;

                VECTOR(*edges)[2 * (ec + i)] = vd.vertex;
                VECTOR(*edges)[2 * (ec + i) + 1] = vertices[i].vertex;
            }
        }

        ec += vd.degree;
    }

    return IGRAPH_SUCCESS;

fail:
    IGRAPH_ERROR("The given degree sequence cannot be realized as a simple graph.", IGRAPH_EINVAL);
}


// Choose vertices in the order of their IDs.
static igraph_error_t igraph_i_havel_hakimi_index(const igraph_vector_int_t *deg, igraph_vector_int_t *edges) {
    igraph_integer_t n = igraph_vector_int_size(deg);

    igraph_integer_t ec = 0; // number of edges added so far

    typedef std::list<vd_pair> vlist;
    vlist vertices;
    for (igraph_integer_t i = 0; i < n; ++i) {
        vertices.push_back(vd_pair(i, VECTOR(*deg)[i]));
    }

    std::vector<vlist::iterator> pointers;
    pointers.reserve(n);
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        pointers.push_back(it);
    }

    for (const auto &pt : pointers) {
        vertices.sort(degree_greater<vd_pair>);

        vd_pair vd = *pt;
        vertices.erase(pt);

        if (vd.degree == 0) {
            continue;
        }

        igraph_integer_t k;
        vlist::iterator it;
        for (it = vertices.begin(), k = 0;
             k != vd.degree && it != vertices.end();
             ++it, ++k) {
            if (--(it->degree) < 0) {
                goto fail;
            }

            VECTOR(*edges)[2 * (ec + k)] = vd.vertex;
            VECTOR(*edges)[2 * (ec + k) + 1] = it->vertex;
        }
        if (it == vertices.end() && k < vd.degree) {
            goto fail;
        }

        ec += vd.degree;
    }

    return IGRAPH_SUCCESS;

fail:
    IGRAPH_ERROR("The given degree sequence cannot be realized as a simple graph.", IGRAPH_EINVAL);
}


/***********************************/
/***** Undirected multigraphs ******/
/***********************************/

// Given a sequence that is sorted, except for its first element,
// move the first element to the correct position fully sort the sequence.
template<typename It, typename Compare>
static void bubble_up(It first, It last, Compare comp) {
    if (first == last)
        return;
    It it = first;
    it++;
    while (it != last) {
        if (comp(*first, *it)) {
            break;
        } else {
            std::swap(*first, *it);
        }
        first = it;
        it++;
    }
}

// In each step, choose a vertex (the largest degree one if largest=true,
// the smallest degree one otherwise) and connect it to the largest remaining degree vertex.
// This will create a connected loopless multigraph, if one exists.
// If loops=true, and a loopless multigraph does not exist, complete the procedure
// by adding loops on the last vertex.
// If largest=false, and the degree sequence was potentially connected, the resulting
// graph will be connected.
static igraph_error_t igraph_i_realize_undirected_multi(const igraph_vector_int_t *deg, igraph_vector_int_t *edges, bool loops, bool largest) {
    igraph_integer_t vcount = igraph_vector_int_size(deg);

    if (vcount == 0)
        return IGRAPH_SUCCESS;

    std::vector<vd_pair> vertices;
    vertices.reserve(vcount);
    for (igraph_integer_t i = 0; i < vcount; ++i) {
        igraph_integer_t d = VECTOR(*deg)[i];
        vertices.push_back(vd_pair(i, d));
    }

    // Initial sort in non-increasing order.
    std::stable_sort(vertices.begin(), vertices.end(), degree_greater<vd_pair>);

    igraph_integer_t ec = 0;
    while (! vertices.empty()) {
        // Remove any zero degrees, and error on negative ones.

        vd_pair &w = vertices.back();

        if (w.degree == 0) {
            vertices.pop_back();
            continue;
        }

        // If only one vertex remains, then the degree sequence cannot be realized as
        // a loopless multigraph. We either complete the graph by adding loops on this vertex
        // or throw an error, depending on the 'loops' setting.
        if (vertices.size() == 1) {
            if (loops) {
                for (igraph_integer_t i=0; i < w.degree/2; ++i) {
                    VECTOR(*edges)[2*ec]   = w.vertex;
                    VECTOR(*edges)[2*ec+1] = w.vertex;
                    ec++;
                }
                break;
            } else {
                IGRAPH_ERROR("The given degree sequence cannot be realized as a loopless multigraph.", IGRAPH_EINVAL);
            }
        }

        // At this point we are guaranteed to have at least two remaining vertices.

        vd_pair *u, *v;
        if (largest) {
            u = &vertices[0];
            v = &vertices[1];
        } else {
            u = &vertices.front();
            v = &vertices.back();
        }

        u->degree -= 1;
        v->degree -= 1;

        VECTOR(*edges)[2*ec]   = u->vertex;
        VECTOR(*edges)[2*ec+1] = v->vertex;
        ec++;

        // Now the first element may be out of order.
        // If largest=true, the first two elements may be out of order.
        // Restore the sorted order using a single step of bubble sort.
        if (largest) {
            bubble_up(vertices.begin()+1, vertices.end(), degree_greater<vd_pair>);
        }
        bubble_up(vertices.begin(), vertices.end(), degree_greater<vd_pair>);
    }

    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_i_realize_undirected_multi_index(const igraph_vector_int_t *deg, igraph_vector_int_t *edges, bool loops) {
    igraph_integer_t vcount = igraph_vector_int_size(deg);

    if (vcount == 0)
        return IGRAPH_SUCCESS;

    typedef std::list<vd_pair> vlist;
    vlist vertices;
    for (igraph_integer_t i = 0; i < vcount; ++i) {
        vertices.push_back(vd_pair(i, VECTOR(*deg)[i]));
    }

    std::vector<vlist::iterator> pointers;
    pointers.reserve(vcount);
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        pointers.push_back(it);
    }

    // Initial sort
    vertices.sort(degree_greater<vd_pair>);

    igraph_integer_t ec = 0;
    for (const auto &pt : pointers) {
        vd_pair vd = *pt;
        vertices.erase(pt);

        while (vd.degree > 0) {
            auto uit = vertices.begin();

            if (vertices.empty() || uit->degree == 0) {
                // We are out of non-zero degree vertices to connect to.
                if (loops) {
                    for (igraph_integer_t i=0; i < vd.degree/2; ++i) {
                        VECTOR(*edges)[2*ec]   = vd.vertex;
                        VECTOR(*edges)[2*ec+1] = vd.vertex;
                        ec++;
                    }
                    return IGRAPH_SUCCESS;
                } else {
                    IGRAPH_ERROR("The given degree sequence cannot be realized as a loopless multigraph.", IGRAPH_EINVAL);
                }
            }

            vd.degree   -= 1;
            uit->degree -= 1;

            VECTOR(*edges)[2*ec]   = vd.vertex;
            VECTOR(*edges)[2*ec+1] = uit->vertex;
            ec++;

            // If there are at least two elements, and the first two are not in order,
            // re-sort the list. A possible optimization would be a version of
            // bubble_up() that can exchange list nodes instead of swapping their values.
            if (vertices.size() > 1) {
                auto wit = uit;
                ++wit;

                if (wit->degree > uit->degree) {
                    vertices.sort(degree_greater<vd_pair>);
                }
            }
        }
    }

    return IGRAPH_SUCCESS;
}


/***********************************/
/***** Directed simple graphs ******/
/***********************************/

inline bool is_nonzero_outdeg(const vbd_pair &vd) {
    return (vd.degree.second != 0);
}


// The below implementations of the Kleitman-Wang algorithm follow the description in https://arxiv.org/abs/0905.4913

// Realize bi-degree sequence as edge list
// If smallest=true, always choose the vertex with "smallest" bi-degree for connecting up next,
// otherwise choose the "largest" (based on lexicographic bi-degree ordering).
static igraph_error_t igraph_i_kleitman_wang(const igraph_vector_int_t *outdeg, const igraph_vector_int_t *indeg, igraph_vector_int_t *edges, bool smallest) {
    igraph_integer_t n = igraph_vector_int_size(indeg); // number of vertices

    igraph_integer_t ec = 0; // number of edges added so far

    std::vector<vbd_pair> vertices;
    vertices.reserve(n);
    for (igraph_integer_t i = 0; i < n; ++i) {
        vertices.push_back(vbd_pair(i, bidegree(VECTOR(*indeg)[i], VECTOR(*outdeg)[i])));
    }

    while (true) {
        // sort vertices by (in, out) degree pairs in decreasing order
        std::stable_sort(vertices.begin(), vertices.end(), degree_greater<vbd_pair>);

        // remove (0,0)-degree vertices
        while (!vertices.empty() && vertices.back().degree == bidegree(0, 0)) {
            vertices.pop_back();
        }

        // if no vertices remain, stop
        if (vertices.empty()) {
            break;
        }

        // choose a vertex the out-stubs of which will be connected
        // note: a vertex with non-zero out-degree is guaranteed to exist
        // because there are _some_ non-zero degrees and the sum of in- and out-degrees
        // is the same
        vbd_pair *vdp;
        if (smallest) {
            vdp = &*std::find_if(vertices.rbegin(), vertices.rend(), is_nonzero_outdeg);
        } else {
            vdp = &*std::find_if(vertices.begin(), vertices.end(), is_nonzero_outdeg);
        }

        // are there a sufficient number of other vertices to connect to?
        if (static_cast<igraph_integer_t>(vertices.size()) - 1 < vdp->degree.second) {
            goto fail;
        }

        // create the connections
        igraph_integer_t k = 0;
        for (auto it = vertices.begin();
             k < vdp->degree.second;
             ++it) {
            if (it->vertex == vdp->vertex) {
                continue;    // do not create a self-loop
            }
            if (--(it->degree.first) < 0) {
                goto fail;
            }

            VECTOR(*edges)[2 * (ec + k)] = vdp->vertex;
            VECTOR(*edges)[2 * (ec + k) + 1] = it->vertex;

            k++;
        }

        ec += vdp->degree.second;
        vdp->degree.second = 0;
    }

    return IGRAPH_SUCCESS;

fail:
    IGRAPH_ERROR("The given directed degree sequences cannot be realized as a simple graph.", IGRAPH_EINVAL);
}


// Choose vertices in the order of their IDs.
static igraph_error_t igraph_i_kleitman_wang_index(const igraph_vector_int_t *outdeg, const igraph_vector_int_t *indeg, igraph_vector_int_t *edges) {
    igraph_integer_t n = igraph_vector_int_size(indeg); // number of vertices

    igraph_integer_t ec = 0; // number of edges added so far

    typedef std::list<vbd_pair> vlist;
    vlist vertices;
    for (igraph_integer_t i = 0; i < n; ++i) {
        vertices.push_back(vbd_pair(i, bidegree(VECTOR(*indeg)[i], VECTOR(*outdeg)[i])));
    }

    std::vector<vlist::iterator> pointers;
    pointers.reserve(n);
    for (auto it = vertices.begin(); it != vertices.end(); ++it) {
        pointers.push_back(it);
    }

    for (const auto &pt : pointers) {
        // sort vertices by (in, out) degree pairs in decreasing order
        // note: std::list::sort does a stable sort
        vertices.sort(degree_greater<vbd_pair>);

        // choose a vertex the out-stubs of which will be connected
        vbd_pair &vd = *pt;

        if (vd.degree.second == 0) {
            continue;
        }

        igraph_integer_t k = 0;
        vlist::iterator it;
        for (it = vertices.begin();
             k != vd.degree.second && it != vertices.end();
             ++it) {
            if (it->vertex == vd.vertex) {
                continue;
            }

            if (--(it->degree.first) < 0) {
                goto fail;
            }

            VECTOR(*edges)[2 * (ec + k)] = vd.vertex;
            VECTOR(*edges)[2 * (ec + k) + 1] = it->vertex;

            ++k;
        }
        if (it == vertices.end() && k < vd.degree.second) {
            goto fail;
        }

        ec += vd.degree.second;
        vd.degree.second = 0;
    }

    return IGRAPH_SUCCESS;

fail:
    IGRAPH_ERROR("The given directed degree sequences cannot be realized as a simple graph.", IGRAPH_EINVAL);
}


/**************************/
/***** Main functions *****/
/**************************/

static igraph_error_t igraph_i_realize_undirected_degree_sequence(
        igraph_t *graph,
        const igraph_vector_int_t *deg,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_realize_degseq_t method)
{
    igraph_integer_t node_count = igraph_vector_int_size(deg);
    igraph_integer_t deg_sum;

    IGRAPH_CHECK(igraph_i_safe_vector_int_sum(deg, &deg_sum));

    if (deg_sum % 2 != 0) {
        IGRAPH_ERROR("The sum of degrees must be even for an undirected graph.", IGRAPH_EINVAL);
    }

    if (node_count > 0 && igraph_vector_int_min(deg) < 0) {
        IGRAPH_ERROR("Vertex degrees must be non-negative.", IGRAPH_EINVAL);
    }

    igraph_vector_int_t edges;
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, deg_sum);

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;
    if ( (allowed_edge_types & IGRAPH_LOOPS_SW) && (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) && (allowed_edge_types & IGRAPH_I_MULTI_LOOPS_SW ) )
    {
        switch (method) {
        case IGRAPH_REALIZE_DEGSEQ_SMALLEST:
            IGRAPH_CHECK(igraph_i_realize_undirected_multi(deg, &edges, true, false));
            break;
        case IGRAPH_REALIZE_DEGSEQ_LARGEST:
            IGRAPH_CHECK(igraph_i_realize_undirected_multi(deg, &edges, true, true));
            break;
        case IGRAPH_REALIZE_DEGSEQ_INDEX:
            IGRAPH_CHECK(igraph_i_realize_undirected_multi_index(deg, &edges, true));
            break;
        default:
            IGRAPH_ERROR("Invalid degree sequence realization method.", IGRAPH_EINVAL);
        }
    }
    else if ( ! (allowed_edge_types & IGRAPH_LOOPS_SW) && (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) )
    {
        switch (method) {
        case IGRAPH_REALIZE_DEGSEQ_SMALLEST:
            IGRAPH_CHECK(igraph_i_realize_undirected_multi(deg, &edges, false, false));
            break;
        case IGRAPH_REALIZE_DEGSEQ_LARGEST:
            IGRAPH_CHECK(igraph_i_realize_undirected_multi(deg, &edges, false, true));
            break;
        case IGRAPH_REALIZE_DEGSEQ_INDEX:
            IGRAPH_CHECK(igraph_i_realize_undirected_multi_index(deg, &edges, false));
            break;
        default:
            IGRAPH_ERROR("Invalid degree sequence realization method.", IGRAPH_EINVAL);
        }
    }
    else if ( (allowed_edge_types & IGRAPH_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) )
    {
        IGRAPH_ERROR("Graph realization with at most one self-loop per vertex is not implemented.", IGRAPH_UNIMPLEMENTED);
    }
    else if ( ! (allowed_edge_types & IGRAPH_LOOPS_SW) && ! (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) )
    {
        switch (method) {
        case IGRAPH_REALIZE_DEGSEQ_SMALLEST:
            IGRAPH_CHECK(igraph_i_havel_hakimi(deg, &edges, false));
            break;
        case IGRAPH_REALIZE_DEGSEQ_LARGEST:
            IGRAPH_CHECK(igraph_i_havel_hakimi(deg, &edges, true));
            break;
        case IGRAPH_REALIZE_DEGSEQ_INDEX:
            IGRAPH_CHECK(igraph_i_havel_hakimi_index(deg, &edges));
            break;
        default:
            IGRAPH_ERROR("Invalid degree sequence realization method.", IGRAPH_EINVAL);
        }
    }
    else
    {
        /* Remainig cases:
         *  - At most one self-loop per vertex but multi-edges between distinct vertices allowed.
         *  - At most one edge between distinct vertices but multi-self-loops allowed.
         * These cases cannot currently be requested through the documented API,
         * so no explanatory error message for now. */
        return IGRAPH_UNIMPLEMENTED;
    }
    IGRAPH_HANDLE_EXCEPTIONS_END;

    IGRAPH_CHECK(igraph_create(graph, &edges, node_count, false));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


static igraph_error_t igraph_i_realize_directed_degree_sequence(
        igraph_t *graph,
        const igraph_vector_int_t *outdeg,
        const igraph_vector_int_t *indeg,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_realize_degseq_t method)
{
    igraph_integer_t node_count = igraph_vector_int_size(outdeg);
    igraph_integer_t edge_count, edge_count2, indeg_sum;

    IGRAPH_CHECK(igraph_i_safe_vector_int_sum(outdeg, &edge_count));

    if (igraph_vector_int_size(indeg) != node_count) {
        IGRAPH_ERROR("In- and out-degree sequences must have the same length.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_safe_vector_int_sum(indeg, &indeg_sum));
    if (indeg_sum != edge_count) {
        IGRAPH_ERROR("In- and out-degree sequences do not sum to the same value.", IGRAPH_EINVAL);
    }

    if (node_count > 0 && (igraph_vector_int_min(outdeg) < 0 || igraph_vector_int_min(indeg) < 0)) {
        IGRAPH_ERROR("Vertex degrees must be non-negative.", IGRAPH_EINVAL);
    }

    /* TODO implement loopless and loopy multigraph case */
    if (allowed_edge_types != IGRAPH_SIMPLE_SW) {
        IGRAPH_ERROR("Realizing directed degree sequences as non-simple graphs is not implemented.", IGRAPH_UNIMPLEMENTED);
    }

    igraph_vector_int_t edges;
    IGRAPH_SAFE_MULT(edge_count, 2, &edge_count2);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, edge_count2);

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;
    switch (method) {
    case IGRAPH_REALIZE_DEGSEQ_SMALLEST:
        IGRAPH_CHECK(igraph_i_kleitman_wang(outdeg, indeg, &edges, true));
        break;
    case IGRAPH_REALIZE_DEGSEQ_LARGEST:
        IGRAPH_CHECK(igraph_i_kleitman_wang(outdeg, indeg, &edges, false));
        break;
    case IGRAPH_REALIZE_DEGSEQ_INDEX:
        IGRAPH_CHECK(igraph_i_kleitman_wang_index(outdeg, indeg, &edges));
        break;
    default:
        IGRAPH_ERROR("Invalid directed degree sequence realization method.", IGRAPH_EINVAL);
    }
    IGRAPH_HANDLE_EXCEPTIONS_END;

    IGRAPH_CHECK(igraph_create(graph, &edges, node_count, true));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup generators
 * \function igraph_realize_degree_sequence
 * \brief Generates a graph with the given degree sequence.
 *
 * This function generates an undirected graph that realizes a given degree sequence,
 * or a directed graph that realized a given pair of out- and in-degree sequences.
 *
 * </para><para>
 * Simple undirected graphs are constructed using the Havel-Hakimi algorithm
 * (undirected case), or the analogous Kleitman-Wang algorithm (directed case).
 * These algorithms work by choosing an arbitrary vertex and connecting all its stubs
 * to other vertices of highest degree.  In the directed case, the "highest" (in, out) degree
 * pairs are determined based on lexicographic ordering. This step is repeated until all degrees
 * have been connected up.
 *
 * </para><para>
 * Loopless multigraphs are generated using an analogous algorithm: an arbitrary vertex is chosen,
 * and it is connected with a single connection to a highest remaining degee vertex. If self-loops
 * are also allowed, the same algorithm is used, but if a non-zero vertex remains at the end of the
 * procedure, the graph is completed by adding self-loops to it. Thus, the result will contain at most
 * one vertex with self-loops.
 *
 * </para><para>
 * The \c method parameter controls the order in which the vertices to be connected are chosen.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * V. Havel,
 * Poznámka o existenci konečných grafů (A remark on the existence of finite graphs),
 * Časopis pro pěstování matematiky 80, 477-480 (1955).
 * http://eudml.org/doc/19050
 *
 * </para><para>
 * S. L. Hakimi,
 * On Realizability of a Set of Integers as Degrees of the Vertices of a Linear Graph,
 * Journal of the SIAM 10, 3 (1962).
 * https://www.jstor.org/stable/2098770
 *
 * </para><para>
 * D. J. Kleitman and D. L. Wang,
 * Algorithms for Constructing Graphs and Digraphs with Given Valences and Factors,
 * Discrete Mathematics 6, 1 (1973).
 * https://doi.org/10.1016/0012-365X%2873%2990037-X
 *
 * </para><para>
 * Sz. Horvát and C. D. Modes,
 * Connectedness matters: construction and exact random sampling of connected networks (2021).
 * https://doi.org/10.1088/2632-072X/abced5
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param outdeg The degree sequence of an undirected graph
 *        (if \p indeg is NULL), or the out-degree sequence of
 *        a directed graph (if \p indeg is given).
 * \param indeg The in-degree sequence of a directed graph.
 *        Pass \c NULL to generate an undirected graph.
 * \param allowed_edge_types The types of edges to allow in the graph. For directed graphs,
 *        only \c IGRAPH_SIMPLE_SW is implemented at this moment. For undirected
 *        graphs, the following values are valid:
 *        \clist
 *          \cli IGRAPH_SIMPLE_SW
 *          simple graphs (i.e. no self-loops or multi-edges allowed).
 *          \cli IGRAPH_LOOPS_SW
 *          single self-loops are allowed, but not multi-edges; currently not implemented.
 *          \cli IGRAPH_MULTI_SW
 *          multi-edges are allowed, but not self-loops.
 *          \cli IGRAPH_LOOPS_SW | IGRAPH_MULTI_SW
 *          both self-loops and multi-edges are allowed.
 *        \endclist
 * \param method The method to generate the graph. Possible values:
 *        \clist
 *          \cli IGRAPH_REALIZE_DEGSEQ_SMALLEST
 *          The vertex with smallest remaining degree is selected first. The result is usually
 *          a graph with high negative degree assortativity. In the undirected case, this method
 *          is guaranteed to generate a connected graph, regardless of whether multi-edges are allowed,
 *          provided that a connected realization exists (see Horvát and Modes, 2021, as well as
 *          http://szhorvat.net/pelican/hh-connected-graphs.html).
 *          In the directed case it tends to generate weakly connected graphs, but this is not
 *          guaranteed.
 *          \cli IGRAPH_REALIZE_DEGSEQ_LARGEST
 *          The vertex with the largest remaining degree is selected first. The result
 *          is usually a graph with high positive degree assortativity, and is often disconnected.
 *          \cli IGRAPH_REALIZE_DEGSEQ_INDEX
 *          The vertices are selected in order of their index (i.e. their position in the degree vector).
 *          Note that sorting the degree vector and using the \c INDEX method is not equivalent
 *          to the \c SMALLEST method above, as \c SMALLEST uses the smallest \em remaining
 *          degree for selecting vertices, not the smallest \em initial degree.
 *         \endclist
 * \return Error code:
 *          \clist
 *          \cli IGRAPH_UNIMPLEMENTED
 *           The requested method is not implemented.
 *          \cli IGRAPH_ENOMEM
 *           There is not enough memory to perform the operation.
 *          \cli IGRAPH_EINVAL
 *           Invalid method parameter, or invalid in- and/or out-degree vectors.
 *           The degree vectors should be non-negative, the length
 *           and sum of \p outdeg and \p indeg should match for directed graphs.
 *          \endclist
 *
 * \sa  \ref igraph_is_graphical() to test graphicality without generating a graph;
 *      \ref igraph_degree_sequence_game() to generate random graphs with a given degree sequence;
 *      \ref igraph_k_regular_game() to generate random regular graphs;
 *      \ref igraph_rewire() to randomly rewire the edges of a graph while preserving its degree sequence.
 *
 * \example examples/simple/igraph_realize_degree_sequence.c
 */

igraph_error_t igraph_realize_degree_sequence(
        igraph_t *graph,
        const igraph_vector_int_t *outdeg, const igraph_vector_int_t *indeg,
        igraph_edge_type_sw_t allowed_edge_types,
        igraph_realize_degseq_t method)
{
    bool directed = indeg != NULL;

    if (directed) {
        return igraph_i_realize_directed_degree_sequence(graph, outdeg, indeg, allowed_edge_types, method);
    } else {
        return igraph_i_realize_undirected_degree_sequence(graph, outdeg, allowed_edge_types, method);
    }
}


// Uses index order to construct an undirected bipartite graph.
// degree1 is considered to range from index [0, len(degree1)[,
// so for this implementation degree1 is always the source degree
// sequence and degree2 is always the dest degree sequence.
static igraph_error_t igraph_i_realize_undirected_bipartite_index(
    igraph_t *graph,
    const igraph_vector_int_t *degree1, const igraph_vector_int_t *degree2,
    igraph_bool_t multiedges
) {
    igraph_integer_t ec = 0; // The number of edges added so far
    igraph_integer_t n1 = igraph_vector_int_size(degree1);
    igraph_integer_t n2 = igraph_vector_int_size(degree2);
    igraph_vector_int_t edges;
    igraph_integer_t ds1_sum;
    igraph_integer_t ds2_sum;

    std::vector<vd_pair> vertices1;
    std::vector<vd_pair> vertices2;
    std::vector<vd_pair> *src_vs = &vertices1;
    std::vector<vd_pair> *dest_vs = &vertices2;

    IGRAPH_CHECK(igraph_i_safe_vector_int_sum(degree1, &ds1_sum));
    IGRAPH_CHECK(igraph_i_safe_vector_int_sum(degree2, &ds2_sum));

    if (ds1_sum != ds2_sum) {
        goto fail;
    }

    // If both degree sequences are empty, it's bigraphical
    if (!(n1 == 0 && n2 == 0)) {
        if (igraph_vector_int_min(degree1) < 0 || igraph_vector_int_min(degree2) < 0){
            goto fail;
        }
    }

    vertices1.reserve(n1);
    vertices2.reserve(n2);

    for (igraph_integer_t i = 0; i < n1; i++) {
        vertices1.push_back(vd_pair(i, VECTOR(*degree1)[i]));
    }
    for (igraph_integer_t i=0; i < n2; i++) {
        vertices2.push_back(vd_pair(i+n1, VECTOR(*degree2)[i]));
    }

	IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ds1_sum + ds2_sum);

    while (!vertices1.empty() && !vertices2.empty()) {
        // Go by index, so we start in ds1, so ds2 needs to be sorted.
        std::stable_sort(vertices2.begin(), vertices2.end(), degree_greater<vd_pair>);
        // No sorting of ds1 needed for index case
        vd_pair vd_src = vertices1.front();
        // No multiedges - Take the first vertex, connect to the largest delta in opposite partition
        if (!multiedges) {
            // Remove the source degrees
            src_vs->erase(src_vs->begin());

            if (vd_src.degree == 0) {
                continue;
            }

            if (dest_vs->size() < size_t(vd_src.degree)) {
                goto fail;
            }

            for (igraph_integer_t i=0;i<vd_src.degree;i++) {
                IGRAPH_ASSERT((*dest_vs)[i].degree - 1 >= 0);

                (*dest_vs)[i].degree--;

                VECTOR(edges)[2*(ec + i)] = vd_src.vertex;
                VECTOR(edges)[2*(ec + i) + 1] = (*dest_vs)[i].vertex;
            }
            ec += vd_src.degree;
        }
        // If multiedges are allowed
        else {
            // If this is the last edge to be created from this vertex, we remove it.
            if (src_vs->front().degree <= 1) {
                src_vs->erase(src_vs->begin());
            } else {
                src_vs->front().degree--;
            }

            if (vd_src.degree == 0) {
                continue;
            }

            if (dest_vs->size() < size_t(1)) {
                goto fail;
            }
            // We should never decrement below zero, but check just in case.
            IGRAPH_ASSERT((*dest_vs)[0].degree - 1 >= 0);

            // Connect to the opposite partition
            (*dest_vs)[0].degree--;

            VECTOR(edges)[2 * ec] = vd_src.vertex;
            VECTOR(edges)[2 * ec + 1] = (*dest_vs)[0].vertex;
            ec++;
        }
    }
    IGRAPH_CHECK(igraph_create(graph, &edges, n1+n2, false));

    igraph_vector_int_destroy(&edges);
	IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;

fail:
    IGRAPH_ERRORF("The given bidegree sequence cannot be realized as a bipartite %sgraph.",
                  IGRAPH_EINVAL, multiedges ? "multi" : "simple ");
}

/**
 * \function igraph_realize_bipartite_degree_sequence
 * \brief Generates a bipartite graph with the given bidegree sequence.
 *
 * \experimental
 *
 * This function generates a bipartite graph with the given bidegree sequence,
 * using a Havel-Hakimi-like construction algorithm. The order in which vertices
 * are connected up is controlled by the \p method parameter. When using the
 * \c IGRAPH_REALIZE_DEGSEQ_SMALLEST method, it is ensured that the graph will be
 * connected if and only if the given bidegree sequence is potentially connected.
 *
 * </para><para>
 * The vertices of the graph will be ordered so that those having \p degrees1
 * come first, followed by \p degrees2.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param degrees1 The degree sequence of the first partition.
 * \param degrees2 The degree sequence of the second partition.
 * \param allowed_edge_types The types of edges to allow in the graph.
 *        \clist
 *          \cli IGRAPH_SIMPLE_SW
 *          simple graph (i.e. no multi-edges allowed).
 *          \cli IGRAPH_MULTI_SW
 *          multi-edges are allowed
 *        \endclist
 * \param method Controls the order in which vertices are selected for connection.
 *        Possible values:
 *        \clist
 *          \cli IGRAPH_REALIZE_DEGSEQ_SMALLEST
 *          The vertex with smallest remaining degree is selected first, from either
 *          partition. The result is usually a graph with high negative degree
 *          assortativity. This method is guaranteed to generate a connected graph,
 *          if one exists.
 *          \cli IGRAPH_REALIZE_DEGSEQ_LARGEST
 *          The vertex with the largest remaining degree is selected first, from
 *          either parition. The result is usually a graph with high positive degree
 *          assortativity, and is often disconnected.
 *          \cli IGRAPH_REALIZE_DEGSEQ_INDEX
 *          The vertices are selected in order of their index.
 *         \endclist
 * \return Error code.
 * \sa \ref igraph_is_bigraphical() to test bigraphicality without generating a graph.
 */

igraph_error_t igraph_realize_bipartite_degree_sequence(
    igraph_t *graph,
    const igraph_vector_int_t *degrees1, const igraph_vector_int_t *degrees2,
    const igraph_edge_type_sw_t allowed_edge_types, const igraph_realize_degseq_t method
) {
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_integer_t ec = 0; // The number of edges added so far
    igraph_integer_t n1 = igraph_vector_int_size(degrees1);
    igraph_integer_t n2 = igraph_vector_int_size(degrees2);
    igraph_vector_int_t edges;
    igraph_integer_t ds1_sum;
    igraph_integer_t ds2_sum;
    igraph_bool_t multiedges;
    igraph_bool_t largest;
    std::vector<vd_pair> vertices1;
    std::vector<vd_pair> vertices2;

    // Bipartite graphs can't have self loops, so we ignore those.
    if (allowed_edge_types & IGRAPH_I_MULTI_EDGES_SW) {
        // Multiedges allowed
        multiedges = true;
    } else {
        // No multiedges
        multiedges = false;
    }

    switch (method) {
    case IGRAPH_REALIZE_DEGSEQ_SMALLEST:
        largest = false;
        break;
    case IGRAPH_REALIZE_DEGSEQ_LARGEST:
        largest = true;
        break;
    case IGRAPH_REALIZE_DEGSEQ_INDEX:
        return igraph_i_realize_undirected_bipartite_index(graph, degrees1, degrees2, multiedges);
    default:
        IGRAPH_ERROR("Invalid bipartite degree sequence realization method.", IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_i_safe_vector_int_sum(degrees1, &ds1_sum));
    IGRAPH_CHECK(igraph_i_safe_vector_int_sum(degrees2, &ds2_sum));

    // Degree sequences of the two partitions must sum to the same value
    if (ds1_sum != ds2_sum) {
        goto fail;
    }

    // If both degree sequences are empty, it's bigraphical
    if (!(n1 == 0 && n2 == 0)) {
        if (igraph_vector_int_min(degrees1) < 0 || igraph_vector_int_min(degrees2) < 0){
            goto fail;
        }
    }

    vertices1.reserve(n1);
    vertices2.reserve(n2);

    for (igraph_integer_t i = 0; i < n1; i++) {
        vertices1.push_back(vd_pair(i, VECTOR(*degrees1)[i]));
    }
    for (igraph_integer_t i=0; i < n2; i++) {
        vertices2.push_back(vd_pair(i+n1, VECTOR(*degrees2)[i]));
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, ds1_sum+ds2_sum);


    std::vector<vd_pair> *src_vs;
    std::vector<vd_pair> *dest_vs;

    while (!vertices1.empty() && !vertices2.empty()) {
        // Sort in non-increasing order.
        // Note: for the smallest method, we can skip sorting the smaller ds, minor optimization.
        // (i.e., we only need to sort the dest partition, since we always just remove the back of the min partition)
        std::stable_sort(vertices1.begin(), vertices1.end(), degree_greater<vd_pair>);
        std::stable_sort(vertices2.begin(), vertices2.end(), degree_greater<vd_pair>);

        vd_pair vd_src(-1, -1);

        if (!largest) {
            vd_pair min1 = vertices1.back();
            vd_pair min2 = vertices2.back();
            if (min1.degree <= min2.degree) {
                src_vs = &vertices1;
                dest_vs = &vertices2;
            } else {
                src_vs = &vertices2;
                dest_vs = &vertices1;
            }

            vd_src = src_vs->back();

        } else {
            vd_pair max1 = vertices1.front();
            vd_pair max2 = vertices2.front();

            if (max1.degree >= max2.degree) {
                src_vs = &vertices1;
                dest_vs = &vertices2;
            } else {
                src_vs = &vertices2;
                dest_vs = &vertices1;
            }

            vd_src = src_vs->front();
        }

        IGRAPH_ASSERT(vd_src.degree != -1);

        if (!multiedges) {
            // Remove the smallest element
            if (!largest) {
                src_vs->pop_back();
            } else {
                // Remove the largest element.
                src_vs->erase(src_vs->begin());
            }

            if (vd_src.degree == 0) {
                continue;
            }
            if (dest_vs->size() < size_t(vd_src.degree)) {
                goto fail;
            }
            for (igraph_integer_t i=0;i < vd_src.degree; i++) {
                // decrement the degree of the delta largest vertices in the opposite partition

                // We should never decrement below zero, but check just in case.
                IGRAPH_ASSERT((*dest_vs)[i].degree - 1 >= 0);

                (*dest_vs)[i].degree--;

                VECTOR(edges)[2*(ec + i)] = vd_src.vertex;
                VECTOR(edges)[2*(ec + i) + 1] = (*dest_vs)[i].vertex;
            }
            ec += vd_src.degree;
        }
        // If multiedges are allowed
        else {
            // The smallest degree is in the back, and we know it is in vertices1
            // If this is the last edge to be created from this vertex, we remove it.
            if (!largest) {
                if (src_vs->back().degree <= 1) {
                    src_vs->pop_back();
                } else {
                    // Otherwise we decrement its degrees by 1 for the edge we are about to create.
                    src_vs->back().degree--;
                }
            } else {
                if (src_vs->front().degree <= 1) {
                    src_vs->erase(src_vs->begin());
                } else {
                    src_vs->front().degree--;
                }
            }

            if (vd_src.degree == 0) {
                continue;
            }

            if (dest_vs->size() < size_t(1)) {
                goto fail;
            }
            // We should never decrement below zero, but check just in case.
            IGRAPH_ASSERT((*dest_vs)[0].degree - 1 >= 0);

            // Connect to the opposite partition
            (*dest_vs)[0].degree--;

            VECTOR(edges)[2 * ec] = vd_src.vertex;
            VECTOR(edges)[2 * ec + 1] = (*dest_vs)[0].vertex;
            ec++;
        }
    }
    IGRAPH_CHECK(igraph_create(graph, &edges, n1+n2, IGRAPH_UNDIRECTED));

    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;

fail:
    IGRAPH_ERRORF("The given bidegree sequence cannot be realized as a bipartite %sgraph.",
                 IGRAPH_EINVAL, multiedges ? "multi" : "simple ");

    IGRAPH_HANDLE_EXCEPTIONS_END;
}
