/*
   igraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

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

#include <igraph.h>

#include "test_utilities.h"

int main(void) {
    igraph_vector_t outdeg, indeg;
    igraph_bool_t simple, multi;
    igraph_t g;

    igraph_rng_seed(igraph_rng_default(), 42);

    /* Zero-legth input */

    igraph_vector_init(&outdeg, 0);

    igraph_chung_lu_game(&g, &outdeg, NULL, true, IGRAPH_CHUNG_LU_ORIGINAL);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, &outdeg, true, IGRAPH_CHUNG_LU_ORIGINAL);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_destroy(&g);

    igraph_vector_destroy(&outdeg);

    /* ORIGINAL */

    /* Must use floating point literals here! 2.0 instead of 2 */
    igraph_vector_init_real_end(&outdeg, -1.0,
                                1.0, 0.0, 2.5, 2.0, 3.0, 2.0, 1.5,
                                -1.0);
    igraph_vector_init_real_end(&indeg, -1.0,
                                2.0, 2.0, 2.0, 2.0, 0.0, 2.0, 2.0,
                                -1.0);

    igraph_chung_lu_game(&g, &outdeg, NULL, false, IGRAPH_CHUNG_LU_ORIGINAL);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, NULL, true, IGRAPH_CHUNG_LU_ORIGINAL);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_has_multiple(&g, &multi);
    IGRAPH_ASSERT(!multi);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, &indeg, false, IGRAPH_CHUNG_LU_ORIGINAL);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, &indeg, true, IGRAPH_CHUNG_LU_ORIGINAL);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_has_multiple(&g, &multi);
    IGRAPH_ASSERT(!multi);
    igraph_destroy(&g);

    igraph_vector_destroy(&indeg);
    igraph_vector_destroy(&outdeg);

    /* GRG */

    /* Must use floating point literals here! 2.0 instead of 2 */
    igraph_vector_init_real_end(&outdeg, -1.0,
                                189.0, 0.0, 2.5, 12.0, 3.0, 2.0, 1.5,
                                -1.0);
    igraph_vector_init_real_end(&indeg, -1.0,
                                2.0, 2.0, 2.0, 2.0, 0.0, 200.0, 2.0,
                                -1.0);

    igraph_chung_lu_game(&g, &outdeg, NULL, false, IGRAPH_CHUNG_LU_MAXENT);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, NULL, true, IGRAPH_CHUNG_LU_MAXENT);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_has_multiple(&g, &multi);
    IGRAPH_ASSERT(!multi);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, &indeg, false, IGRAPH_CHUNG_LU_MAXENT);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, &indeg, true, IGRAPH_CHUNG_LU_MAXENT);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_has_multiple(&g, &multi);
    IGRAPH_ASSERT(!multi);
    igraph_destroy(&g);

    /* NR */

    igraph_chung_lu_game(&g, &outdeg, NULL, false, IGRAPH_CHUNG_LU_NR);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, NULL, true, IGRAPH_CHUNG_LU_NR);
    IGRAPH_ASSERT(! igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_has_multiple(&g, &multi);
    IGRAPH_ASSERT(!multi);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, &indeg, false, IGRAPH_CHUNG_LU_NR);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_is_simple(&g, &simple, IGRAPH_DIRECTED);
    IGRAPH_ASSERT(simple);
    igraph_destroy(&g);

    igraph_chung_lu_game(&g, &outdeg, &indeg, true, IGRAPH_CHUNG_LU_NR);
    IGRAPH_ASSERT(igraph_is_directed(&g));
    IGRAPH_ASSERT(igraph_vcount(&g) == igraph_vector_size(&outdeg));
    igraph_has_multiple(&g, &multi);
    IGRAPH_ASSERT(!multi);
    igraph_destroy(&g);

    /* Invalid input */

    /* Bad variant */
    CHECK_ERROR(igraph_chung_lu_game(&g, &outdeg, &indeg, true, (igraph_chung_lu_t) -1), IGRAPH_EINVAL);

    /* Inconsistent sum */
    VECTOR(outdeg)[0] = 0;
    CHECK_ERROR(igraph_chung_lu_game(&g, &outdeg, &indeg, true, IGRAPH_CHUNG_LU_ORIGINAL), IGRAPH_EINVAL);

    /* Negative weight */
    VECTOR(outdeg)[0] = -1;
    CHECK_ERROR(igraph_chung_lu_game(&g, &outdeg, NULL, true, IGRAPH_CHUNG_LU_ORIGINAL), IGRAPH_EINVAL);

    /* Non-finite weight */
    VECTOR(outdeg)[0] = IGRAPH_INFINITY;
    CHECK_ERROR(igraph_chung_lu_game(&g, &outdeg, NULL, true, IGRAPH_CHUNG_LU_ORIGINAL), IGRAPH_EINVAL);

    igraph_vector_destroy(&indeg);
    igraph_vector_destroy(&outdeg);

    VERIFY_FINALLY_STACK();

    return 0;
}
