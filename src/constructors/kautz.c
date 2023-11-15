/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2005-2021 The igraph development team

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

#include "core/interruption.h"
#include "math/safe_intop.h"

/**
 * \function igraph_kautz
 * \brief Generate a Kautz graph.
 *
 * A Kautz graph is a labeled graph, vertices are labeled by strings
 * of length \c n+1 above an alphabet with \c m+1 letters, with
 * the restriction that every two consecutive letters in the string
 * must be different. There is a directed edge from a vertex \c v to
 * another vertex \c w if it is possible to transform the string of
 * \c v into the string of \c w by removing the first letter and
 * appending a letter to it. For string length 1 the new letter
 * cannot equal the old letter, so there are no loops.
 *
 * </para><para>
 * Kautz graphs have some interesting properties, see e.g. Wikipedia
 * for details.
 *
 * </para><para>
 * Vincent Matossian wrote the first version of this function in R,
 * thanks.
 * \param graph Pointer to an uninitialized graph object, the result
 * will be stored here.
 * \param m Integer, \c m+1 is the number of letters in the alphabet.
 * \param n Integer, \c n+1 is the length of the strings.
 * \return Error code.
 *
 * \sa \ref igraph_de_bruijn().
 *
 * Time complexity: O(|V|* [(m+1)/m]^n +|E|), in practice it is more
 * like O(|V|+|E|). |V| is the number of vertices, |E| is the number
 * of edges and \c m and \c n are the corresponding arguments.
 */
igraph_error_t igraph_kautz(igraph_t *graph, igraph_integer_t m, igraph_integer_t n) {

    /* m+1 - number of symbols */
    /* n+1 - length of strings */

    igraph_integer_t no_of_nodes, no_of_edges;
    igraph_integer_t allstrings;
    igraph_integer_t i, j, idx = 0;
    igraph_vector_int_t edges;
    igraph_vector_int_t digits, table;
    igraph_vector_int_t index1, index2;
    igraph_integer_t actb = 0;
    igraph_integer_t actvalue = 0;
    int iter = 0;

    if (m < 0 || n < 0) {
        IGRAPH_ERROR("`m' and `n' should be non-negative in a Kautz graph",
                     IGRAPH_EINVAL);
    }

    if (n == 0) {
        return igraph_full(graph, m + 1, IGRAPH_DIRECTED, IGRAPH_NO_LOOPS);
    }
    if (m == 0) {
        return igraph_empty(graph, 0, IGRAPH_DIRECTED);
    }

    /* no_of_nodes = ((m + 1) * pow(m, n)) */
    {
        igraph_real_t m_to_pow_n_real = pow(m, n);
        igraph_integer_t m_to_pow_n = m_to_pow_n_real;
        if (m_to_pow_n != m_to_pow_n_real) {
            IGRAPH_ERRORF("Parameters (%" IGRAPH_PRId ", %" IGRAPH_PRId ") too large for Kautz graph.", IGRAPH_EINVAL,
                          m, n);
        }
        IGRAPH_SAFE_MULT(m+1, m_to_pow_n, &no_of_nodes);
    }
    /* no_of_edges = m * no_of_nodes */
    IGRAPH_SAFE_MULT(no_of_nodes, m, &no_of_edges);

    {
        igraph_real_t allstrings_real = pow(m + 1, n + 1);
        allstrings = allstrings_real;
        if (allstrings != allstrings_real) {
            IGRAPH_ERRORF("Parameters (%" IGRAPH_PRId ", %" IGRAPH_PRId ") too large for Kautz graph.", IGRAPH_EINVAL,
                          m, n);
        }
    }

    IGRAPH_VECTOR_INT_INIT_FINALLY(&edges, 0);

    IGRAPH_CHECK(igraph_vector_int_init(&table, n + 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &table);
    j = 1;
    for (i = n; i >= 0; i--) {
        VECTOR(table)[i] = j;
        j *= (m + 1);
    }

    IGRAPH_CHECK(igraph_vector_int_init(&digits, n + 1));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &digits);
    IGRAPH_CHECK(igraph_vector_int_init(&index1, pow(m + 1, n + 1)));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &index1);
    IGRAPH_CHECK(igraph_vector_int_init(&index2, no_of_nodes));
    IGRAPH_FINALLY(igraph_vector_int_destroy, &index2);

    /* Fill the index tables*/
    while (true) {
        /* at the beginning of the loop, 0:actb contain the valid prefix */
        /* we might need to fill it to get a valid string */
        igraph_integer_t z = 0;
        if (VECTOR(digits)[actb] == 0) {
            z = 1;
        }
        for (actb++; actb <= n; actb++) {
            VECTOR(digits)[actb] = z;
            actvalue += z * VECTOR(table)[actb];
            z = 1 - z;
        }
        actb = n;

        /* ok, we have a valid string now */
        VECTOR(index1)[actvalue] = idx + 1;
        VECTOR(index2)[idx] = actvalue;
        idx++;

        /* finished? */
        if (idx >= no_of_nodes) {
            break;
        }

        /* not yet, we need a valid prefix now */
        while (true) {
            /* try to increase digits at position actb */
            igraph_integer_t next = VECTOR(digits)[actb] + 1;
            if (actb != 0 && VECTOR(digits)[actb - 1] == next) {
                next++;
            }
            if (next <= m) {
                /* ok, no problem */
                actvalue += (next - VECTOR(digits)[actb]) * VECTOR(table)[actb];
                VECTOR(digits)[actb] = next;
                break;
            } else {
                /* bad luck, try the previous digit */
                actvalue -= VECTOR(digits)[actb] * VECTOR(table)[actb];
                actb--;
            }
        }
    }

    {
        igraph_integer_t no_of_edges2;
        IGRAPH_SAFE_MULT(no_of_edges, 2, &no_of_edges2);
        IGRAPH_CHECK(igraph_vector_int_reserve(&edges, no_of_edges2));
    }

    /* Now come the edges at last */
    for (i = 0; i < no_of_nodes; i++) {
        igraph_integer_t fromvalue = VECTOR(index2)[i];
        igraph_integer_t lastdigit = fromvalue % (m + 1);
        igraph_integer_t basis = (fromvalue * (m + 1)) % allstrings;
        for (j = 0; j <= m; j++) {
            igraph_integer_t tovalue, to;
            if (j == lastdigit) {
                continue;
            }
            tovalue = basis + j;
            to = VECTOR(index1)[tovalue] - 1;
            if (to < 0) {
                continue;
            }
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, i));
            IGRAPH_CHECK(igraph_vector_int_push_back(&edges, to));
        }
        IGRAPH_ALLOW_INTERRUPTION_LIMITED(iter, 1 << 10);
    }

    igraph_vector_int_destroy(&index2);
    igraph_vector_int_destroy(&index1);
    igraph_vector_int_destroy(&digits);
    igraph_vector_int_destroy(&table);
    IGRAPH_FINALLY_CLEAN(4);

    IGRAPH_CHECK(igraph_create(graph, &edges, no_of_nodes, IGRAPH_DIRECTED));
    igraph_vector_int_destroy(&edges);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}
