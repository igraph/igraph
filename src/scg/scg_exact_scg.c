/*
 *  SCGlib : A C library for the spectral coarse graining of matrices
 *  as described in the paper: Shrinking Matrices while preserving their
 *  eigenpairs with Application to the Spectral Coarse Graining of Graphs.
 *  Preprint available at <http://people.epfl.ch/david.morton>
 *
 *  Copyright (C) 2008 David Morton de Lachapelle <david.morton@a3.epfl.ch>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 *  02110-1301 USA
 *
 *  DESCRIPTION
 *  -----------
 *    The exact_coarse_graining function labels all the objects whose
 *    components in 'v' are equal. The result is stored in 'gr'. Labels
 *    are positive consecutive integers starting from 0.
 *    See also Section 5.4.1 (last paragraph) of the above reference.
 */

#include "scg_headers.h"

#include "igraph_memory.h"
#include "igraph_qsort.h"

#include <math.h>

int igraph_i_exact_coarse_graining(const igraph_real_t *v,
                                   int *gr, int n) {
    int i, gr_nb;
    igraph_i_scg_indval_t *w = IGRAPH_CALLOC(n, igraph_i_scg_indval_t);

    if (!w) {
        IGRAPH_ERROR("SCG error", IGRAPH_ENOMEM);
    }
    IGRAPH_FINALLY(igraph_free, w);

    for (i = 0; i < n; i++) {
        w[i].val = v[i];
        w[i].ind = i;
    }

    igraph_qsort(w, (size_t) n, sizeof(igraph_i_scg_indval_t), igraph_i_compare_ind_val);

    gr_nb = 0;
    gr[w[0].ind] = gr_nb;
    for (i = 1; i < n; i++) {
        if ( fabs(w[i].val - w[i - 1].val) > 1e-14 ) {
            gr_nb++;
        }
        gr[w[i].ind] = gr_nb;
    }

    IGRAPH_FREE(w);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}
