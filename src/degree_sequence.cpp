/*
  Constructing realizations of degree sequences and bi-degree sequences.
  Copyright (C) 2018 Szabolcs Horvat <szhorvat@gmail.com>

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

#include "igraph_constructors.h"
#include "igraph_interface.h"

#include <vector>
#include <list>
#include <algorithm>
#include <utility>


// (vertex, degree) pair
struct vd_pair {
    long vertex;
    igraph_integer_t degree;

    vd_pair(long vertex, igraph_integer_t degree) : vertex(vertex), degree(degree) {}
};

// (indegree, outdegree)
typedef std::pair<igraph_integer_t, igraph_integer_t> bidegree;

// (vertex, bidegree) pair
struct vbd_pair {
    long vertex;
    bidegree degree;

    vbd_pair(long vertex, bidegree degree) : vertex(vertex), degree(degree) {}
};

// Comparison function for vertex-degree pairs.
// Also used for lexicographic sorting of bi-degrees.
template<typename T> inline bool degree_greater(const T &a, const T &b) {
    return a.degree > b.degree;
}

template<typename T> inline bool degree_less(const T &a, const T &b) {
    return a.degree < b.degree;
}


// Generate undirected realization as edge-list.
// If largest=true, always choose the vertex with the largest remaining degree to connect up next.
// Otherwise, always choose the one with the smallest remaining degree.
static int igraph_i_havel_hakimi(const igraph_vector_t *deg, igraph_vector_t *edges, bool largest) {
    long n = igraph_vector_size(deg);

    long ec = 0; // number of edges added so far

    std::vector<vd_pair> vertices;
    vertices.reserve(n);
    for (int i = 0; i < n; ++i) {
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

        if (vd.degree < 0) {
            IGRAPH_ERROR("Vertex degrees must be positive", IGRAPH_EINVAL);
        }

        if (vd.degree == 0) {
            continue;
        }

        if (vertices.size() < size_t(vd.degree)) {
            goto fail;
        }

        if (largest) {
            for (int i = 0; i < vd.degree; ++i) {
                if (--(vertices[vertices.size() - 1 - i].degree) < 0) {
                    goto fail;
                }

                VECTOR(*edges)[2 * (ec + i)] = vd.vertex;
                VECTOR(*edges)[2 * (ec + i) + 1] = vertices[vertices.size() - 1 - i].vertex;
            }
        } else {
            // this loop can only be reached if all zero-degree nodes have already been removed
            // therefore decrementing remaining degrees is safe
            for (int i = 0; i < vd.degree; ++i) {
                vertices[i].degree--;

                VECTOR(*edges)[2 * (ec + i)] = vd.vertex;
                VECTOR(*edges)[2 * (ec + i) + 1] = vertices[i].vertex;
            }
        }

        ec += vd.degree;
    }

    return IGRAPH_SUCCESS;

fail:
    IGRAPH_ERROR("The given degree sequence is not realizable", IGRAPH_EINVAL);
}


// Choose vertices in the order of their IDs.
static int igraph_i_havel_hakimi_index(const igraph_vector_t *deg, igraph_vector_t *edges) {
    long n = igraph_vector_size(deg);

    long ec = 0; // number of edges added so far

    typedef std::list<vd_pair> vlist;
    vlist vertices;
    for (int i = 0; i < n; ++i) {
        vertices.push_back(vd_pair(i, VECTOR(*deg)[i]));
    }

    std::vector<vlist::iterator> pointers;
    pointers.reserve(n);
    for (vlist::iterator it = vertices.begin(); it != vertices.end(); ++it) {
        pointers.push_back(it);
    }

    for (std::vector<vlist::iterator>::iterator pt = pointers.begin(); pt != pointers.end(); ++pt) {
        vertices.sort(degree_greater<vd_pair>);

        vd_pair vd = **pt;
        vertices.erase(*pt);

        if (vd.degree < 0) {
            IGRAPH_ERROR("Vertex degrees must be positive", IGRAPH_EINVAL);
        }

        if (vd.degree == 0) {
            continue;
        }

        int k;
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
    IGRAPH_ERROR("The given degree sequence is not realizable", IGRAPH_EINVAL);
}


inline bool is_nonzero_outdeg(const vbd_pair &vd) {
    return (vd.degree.second != 0);
}


// The below implementations of the Kleitman-Wang algorithm follow the description in https://arxiv.org/abs/0905.4913

// Realize bi-degree sequence as edge list
// If smallest=true, always choose the vertex with "smallest" bi-degree for connecting up next,
// otherwise choose the "largest" (based on lexicographic bi-degree ordering).
static int igraph_i_kleitman_wang(const igraph_vector_t *outdeg, const igraph_vector_t *indeg, igraph_vector_t *edges, bool smallest) {
    long n = igraph_vector_size(indeg); // number of vertices

    long ec = 0; // number of edges added so far

    std::vector<vbd_pair> vertices;
    vertices.reserve(n);
    for (int i = 0; i < n; ++i) {
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
        vbd_pair *vdp;
        if (smallest) {
            vdp = &*std::find_if(vertices.rbegin(), vertices.rend(), is_nonzero_outdeg);
        } else {
            vdp = &*std::find_if(vertices.begin(), vertices.end(), is_nonzero_outdeg);
        }


        if (vdp->degree.first < 0 || vdp->degree.second < 0) {
            IGRAPH_ERROR("Vertex degrees must be positive", IGRAPH_EINVAL);
        }

        // are there a sufficient number of other vertices to connect to?
        if (vertices.size() < vdp->degree.second - 1) {
            goto fail;
        }

        // create the connections
        int k = 0;
        for (std::vector<vbd_pair>::iterator it = vertices.begin();
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
    IGRAPH_ERROR("The given directed degree sequence is not realizable", IGRAPH_EINVAL);
}


// Choose vertices in the order of their IDs.
static int igraph_i_kleitman_wang_index(const igraph_vector_t *outdeg, const igraph_vector_t *indeg, igraph_vector_t *edges) {
    long n = igraph_vector_size(indeg); // number of vertices

    long ec = 0; // number of edges added so far

    typedef std::list<vbd_pair> vlist;
    vlist vertices;
    for (int i = 0; i < n; ++i) {
        vertices.push_back(vbd_pair(i, bidegree(VECTOR(*indeg)[i], VECTOR(*outdeg)[i])));
    }

    std::vector<vlist::iterator> pointers;
    pointers.reserve(n);
    for (vlist::iterator it = vertices.begin(); it != vertices.end(); ++it) {
        pointers.push_back(it);
    }

    for (std::vector<vlist::iterator>::iterator pt = pointers.begin(); pt != pointers.end(); ++pt) {
        // sort vertices by (in, out) degree pairs in decreasing order
        // note: std::list::sort does a stable sort
        vertices.sort(degree_greater<vbd_pair>);

        // choose a vertex the out-stubs of which will be connected
        vbd_pair &vd = **pt;

        if (vd.degree.second == 0) {
            continue;
        }

        if (vd.degree.first < 0 || vd.degree.second < 0) {
            IGRAPH_ERROR("Vertex degrees must be positive", IGRAPH_EINVAL);
        }

        int k = 0;
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
    IGRAPH_ERROR("The given directed degree sequence is not realizable", IGRAPH_EINVAL);
}


static int igraph_i_realize_undirected_degree_sequence(
    igraph_t *graph,
    const igraph_vector_t *deg,
    igraph_realize_degseq_t method) {
    long node_count = igraph_vector_size(deg);
    long deg_sum = long(igraph_vector_sum(deg));

    if (deg_sum % 2 != 0) {
        IGRAPH_ERROR("The sum of degrees must be even for an undirected graph", IGRAPH_EINVAL);
    }

    igraph_vector_t edges;
    IGRAPH_CHECK(igraph_vector_init(&edges, deg_sum));
    IGRAPH_FINALLY(igraph_vector_destroy, &edges);

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
        IGRAPH_ERROR("Invalid degree sequence realization method", IGRAPH_EINVAL);
    }

    igraph_create(graph, &edges, igraph_integer_t(node_count), false);

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


static int igraph_i_realize_directed_degree_sequence(
    igraph_t *graph,
    const igraph_vector_t *outdeg,
    const igraph_vector_t *indeg,
    igraph_realize_degseq_t method) {
    long node_count = igraph_vector_size(outdeg);
    long edge_count = long(igraph_vector_sum(outdeg));

    if (igraph_vector_size(indeg) != node_count) {
        IGRAPH_ERROR("In- and out-degree sequences must have the same length", IGRAPH_EINVAL);
    }
    if (igraph_vector_sum(indeg) != edge_count) {
        IGRAPH_ERROR("In- and out-degree sequences do not sum to the same value", IGRAPH_EINVAL);
    }

    igraph_vector_t edges;
    IGRAPH_CHECK(igraph_vector_init(&edges, 2 * edge_count));
    IGRAPH_FINALLY(igraph_vector_destroy, &edges);

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
        IGRAPH_ERROR("Invalid bi-degree sequence realization method", IGRAPH_EINVAL);
    }

    igraph_create(graph, &edges, igraph_integer_t(node_count), true);

    igraph_vector_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}


/**
 * \ingroup generators
 * \function igraph_realize_degree_sequence
 * \brief Generates a graph with the given degree sequence
 *
 * This function constructs a simple graph that realizes the given degree sequence
 * using the Havel-Hakimi algorithm, or the given (directed) out- and in-degree
 * sequences using the related Kleitman-Wang algorithm.
 *
 * The algorithms work by choosing an arbitrary vertex and connecting all its stubs
 * to other vertices of highest degree.  In the directed case, the "highest" (in, out) degree
 * pairs are determined based on lexicographic ordering.
 *
 * The \c method parameter controls the order in which the vertices to be connected are chosen.
 *
 * \param graph Pointer to an uninitialized graph object.
 * \param outdeg The degree sequence for a simple undirected graph
 *        (if \p indeg is NULL or of length zero), or the out-degree sequence of
 *        a directed graph (if \p indeg is of nonzero size).
 * \param indeg It is either a zero-length vector or \c NULL (if an undirected graph
 *        is generated), or the in-degree sequence.
 * \param method The method to generate the graph. Possible values:
 *        \clist
 *          \cli IGRAPH_REALIZE_DEGSEQ_SMALLEST
 *          The vertex with smallest remaining degree is selected first. The result is usually
 *          a graph with high negative degree assortativity. In the undirected case, this method
 *          is guaranteed to generate a connected graph, provided that a connected realization exists.
 *          See http://szhorvat.net/pelican/hh-connected-graphs.html for a proof.
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
 *          \cli IGRAPH_ENOMEM
 *           There is not enough memory to perform the operation.
 *          \cli IGRAPH_EINVAL
 *           Invalid method parameter, or invalid in- and/or out-degree vectors.
 *           The degree vectors should be non-negative, the length
 *           and sum of \p outdeg and \p indeg should match for directed graphs.
 *          \endclist
 *
 * \sa  \ref igraph_is_graphical_degree_sequence()
 *      \ref igraph_degree_sequence_game()
 *      \ref igraph_k_regular_game()
 *      \ref igraph_rewire()
 *
 */

int igraph_realize_degree_sequence(
    igraph_t *graph,
    const igraph_vector_t *outdeg, const igraph_vector_t *indeg,
    igraph_realize_degseq_t method) {
    long n = igraph_vector_size(outdeg);
    if (n != igraph_integer_t(n)) { // does the vector size fit into an igraph_integer_t ?
        IGRAPH_ERROR("Degree sequence vector too long", IGRAPH_EINVAL);
    }

    bool directed = bool(indeg) && igraph_vector_size(indeg) != 0;

    try {
        if (directed) {
            return igraph_i_realize_directed_degree_sequence(graph, outdeg, indeg, method);
        } else {
            return igraph_i_realize_undirected_degree_sequence(graph, outdeg, method);
        }
    } catch (const std::bad_alloc &) {
        IGRAPH_ERROR("Cannot realize degree sequence due to insufficient memory", IGRAPH_ENOMEM);
    }
}
