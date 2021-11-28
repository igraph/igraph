/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

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
#include "igraph_umap.h"
#include "test_utilities.inc"

int main() {
	igraph_matrix_t data;
	igraph_matrix_t layout;
	igraph_real_t data_arr[] = {0, 1, 1, 0};

	matrix_init_real_row_major(&data, 2, 2, data_arr);
	igraph_matrix_init(&layout, 0, 0);

    IGRAPH_ASSERT(igraph_umap(&data, &layout) == IGRAPH_SUCCESS);

    VERIFY_FINALLY_STACK();
    return 0;
}
