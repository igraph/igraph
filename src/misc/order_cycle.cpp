/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#include "misc/order_cycle.h"

#include "igraph_interface.h"

#include "core/exceptions.h"

#include <map>
#include <utility>

// Initialized to {-1, -1}
struct eid_pair_t : public std::pair<igraph_integer_t, igraph_integer_t> {
    eid_pair_t() : std::pair<igraph_integer_t, igraph_integer_t>(-1, -1) { }
};

/**
 * \function igraph_i_order_cycle
 * \brief Reorders edges of a cycle in cycle order
 *
 * This function takes \p cycle, a vector of arbitrarily ordered edge IDs,
 * representing a graph cycle. It produces a vector \p res containing the
 * same IDs in cycle order. \p res must be initialized when calling this function.
 */
igraph_error_t igraph_i_order_cycle(
        const igraph_t *graph,
        const igraph_vector_int_t *cycle,
        igraph_vector_int_t *res) {

    IGRAPH_HANDLE_EXCEPTIONS_BEGIN;

    igraph_integer_t n = igraph_vector_int_size(cycle);
    IGRAPH_ASSERT(n > 0);

    std::map<igraph_integer_t, eid_pair_t> inclist;
    for (igraph_integer_t i=0; i < n; ++i) {
        igraph_integer_t eid = VECTOR(*cycle)[i];

        {
            igraph_integer_t from = IGRAPH_FROM(graph, eid);
            auto &p = inclist[from];
            if (p.first < 0) {
                p.first = eid;
            } else {
                IGRAPH_ASSERT(p.second < 0);
                p.second = eid;
            }
        }

        {
            igraph_integer_t to = IGRAPH_TO(graph, eid);
            auto &p = inclist[to];
            if (p.first < 0) {
                p.first = eid;
            } else {
                IGRAPH_ASSERT(p.second < 0);
                p.second = eid;
            }
        }
    }

    igraph_vector_int_clear(res);
    IGRAPH_CHECK(igraph_vector_int_reserve(res, igraph_vector_int_size(cycle)));
    igraph_integer_t current_e = VECTOR(*cycle)[0];
    igraph_integer_t current_v = IGRAPH_FROM(graph, current_e);
    for (igraph_integer_t i=0; i < n; ++i) {
        const auto &p = inclist.at(current_v);
        igraph_vector_int_push_back(res, current_e); /* reserved */
        igraph_integer_t next_e = p.first;
        if (next_e == current_e) {
            next_e = p.second;
        }
        current_e = next_e;
        igraph_integer_t next_v = IGRAPH_FROM(graph, current_e);
        if (next_v == current_v) {
            next_v = IGRAPH_TO(graph, current_e);
        }
        current_v = next_v;
    }

    IGRAPH_HANDLE_EXCEPTIONS_END;

    return IGRAPH_SUCCESS;
}
