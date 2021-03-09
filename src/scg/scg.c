/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-12  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge, MA, 02138 USA

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
 *    The grouping function takes as argument 'nev' eigenvectors and
 *    and tries to minimize the eigenpair shifts induced by the coarse
 *    graining (Section 5 of the above reference). The eigenvectors are
 *    stored in a 'nev'x'n' matrix 'v'.
 *    The 'algo' parameter can take the following values
 *      1  ->  Optimal method (sec. 5.3.1)
 *      2  ->  Intervals+k-means (sec. 5.3.3)
 *      3  ->  Intervals (sec. 5.3.2)
 *      4  ->  Exact SCG (sec. 5.4.1--last paragraph)
 *    'nt' is a vector of length 'nev' giving either the size of the
 *    partitions (if algo = 1) or the number of intervals to cut the
 *    eigenvectors if algo = 2 or algo = 3. When algo = 4 this parameter
 *    is ignored. 'maxiter' fixes the maximum number of iterations of
 *    the k-means algorithm, and is only considered when algo = 2.
 *    All the algorithms try to find a minimizing partition of
 *    ||v_i-Pv_i|| where P is a problem-specific projector and v_i denotes
 *    the eigenvectors stored in v. The final partition is worked out
 *    as decribed in Method 1 of Section 5.4.2.
 *    'matrix' provides the type of SCG (i.e. the form of P). So far,
 *    the options are those described in section 6, that is:
 *      1  ->  Symmetric (sec. 6.1)
 *      2  ->  Laplacian (sec. 6.2)
 *      3  ->  Stochastic (sec. 6.3)
 *    In the stochastic case, a valid distribution probability 'p' must be
 *    provided. In all other cases, 'p' is ignored and can be set to NULL.
 *    The group labels in the final partition are given in 'gr' as positive
 *    consecutive integers starting from 0.
 */

#include "igraph_scg.h"

#include "igraph_eigen.h"
#include "igraph_interface.h"
#include "igraph_structural.h"
#include "igraph_community.h"
#include "igraph_constructors.h"
#include "igraph_conversion.h"
#include "igraph_memory.h"
#include "igraph_qsort.h"

#include "misc/conversion_internal.h"

#include "scg_headers.h"

#include "math.h"

/**
 * \section about_scg
 *
 * <para>
 * The SCG functions provide a framework, called Spectral Coarse Graining
 * (SCG), for reducing large graphs while preserving their
 * <emphasis>spectral-related features</emphasis>, that is features
 * closely related with the eigenvalues and eigenvectors of a graph
 * matrix (which for now can be the adjacency, the stochastic, or the
 * Laplacian matrix).
 * </para>
 *
 * <para>
 * Common examples of such features comprise the first-passage-time of
 * random walkers on Markovian graphs, thermodynamic properties of
 * lattice models in statistical physics (e.g. Ising model), and the
 * epidemic threshold of epidemic network models (SIR and SIS models).
 * </para>
 *
 * <para>
 * SCG differs from traditional clustering schemes by producing a
 * <emphasis>coarse-grained graph</emphasis> (not just a partition of
 * the vertices), representative of the original one. As shown in [1],
 * Principal Component Analysis can be viewed as a particular SCG,
 * called <emphasis>exact SCG</emphasis>, where the matrix to be
 * coarse-grained is the covariance matrix of some data set.
 * </para>
 *
 * <para>
 * SCG should be of interest to practitioners of various
 * fields dealing with problems where matrix eigenpairs play an important
 * role, as for instance is the case of dynamical processes on networks.
 * </para>
 *
 * <section id="scg-in-brief"><title>SCG in brief</title>
 * <para>
 * The main idea of SCG is to operate on a matrix a shrinkage operation
 * specifically designed to preserve some of the matrix eigenpairs while
 * not altering other important matrix features (such as its structure).
 * Mathematically, this idea was expressed as follows. Consider a
 * (complex) n x n matrix M and form the product
 * <blockquote><para><phrase role="math">
 *   M'=LMR*,
 * </phrase></para></blockquote>
 * where n' &lt; n and L, R are from C[n'xn]} and are such
 * that LR*=I[n'] (R* denotes the conjugate transpose of R). Under
 * these assumptions, it can be shown that P=R*L is an n'-rank
 * projector and that, if (lambda, v) is a (right)
 * eigenpair of M (i.e. Mv=lambda v} and P is orthogonal, there exists
 * an eigenvalue lambda' of M' such that
 * <blockquote><para><phrase role="math">
 *   |lambda-lambda'| &lt;= const ||e[P](v)||
 *   [1+O(||e[P](v)||<superscript>2</superscript>)],
 * </phrase></para></blockquote>
 * where ||e[P](v)||=||v-Pv||. Hence, if P (or equivalently
 * L, R) is chosen so as to make ||e[P](v)|| as small as possible, one
 * can preserve to any desired level the original eigenvalue
 * lambda in the coarse-grained matrix M';
 * under extra assumptions on M, this result can be generalized to
 * eigenvectors [1]. This leads to the following generic definition of a
 * SCG problem.
 * </para>
 *
 * <para>
 * Given M (C[nxn]) and (lambda, v), a (right) eigenpair of M to be
 * preserved by the coarse graining, the problem is to find a projector
 * P' solving
 * <blockquote><para><phrase role="math">
 *   min(||e[P](v)||, p in Omega),
 * </phrase></para></blockquote>
 * where Omega is a set of projectors in C[nxn] described by some
 * ad hoc constraints c[1], ..., c[r]
 * (e.g. c[1]: P in R[nxn], c[2]: P=t(P), c[3]: P[i,j] >= 0}, etc).
 * </para>
 *
 * <para>
 * Choosing pertinent constraints to solve the SCG problem is of great
 * importance in applications. For instance, in the absence of
 * constraints the SCG problem is solved trivially by
 * P'=vv* (v is assumed normalized). We have designed a particular
 * constraint, called <emphasis>homogeneous mixing</emphasis>, which
 * ensures that vertices belonging to the same group are merged
 * consistently from a physical point of view (see [1] for
 * details). Under this constraint the SCG problem reduces to finding
 * the partition of 1, ..., n (labeling the original vertices)
 * minimizing
 * <blockquote><para><phrase role="math">
 *   ||e[P](v)||<superscript>2</superscript> =
 *   sum([v(i)-(Pv)(i)]<superscript>2</superscript>;
 *   alpha=1,...,n', i in alpha),
 * </phrase></para></blockquote>
 * where alpha denotes a group (i.e. a block) in a partition of
 * {1, ..., n}, and |alpha| is the number of elements in alpha.
 * </para>
 *
 * <para>
 * If M is symmetric or stochastic, for instance, then it may be
 * desirable (or mandatory) to choose L, R so that M' is symmetric or
 * stochastic as well. This <emphasis>structural constraint</emphasis>
 * has led to the construction of particular semi-projectors for
 * symmetric [1], stochastic [3] and Laplacian [2] matrices, that are
 * made available.
 * </para>
 *
 * <para>
 * In short, the coarse graining of matrices and graphs involves:
 * \olist
 *   \oli Retrieving a matrix or a graph matrix M from the
 *     problem.
 *   \oli Computing the eigenpairs of M to be preserved in the
 *     coarse-grained graph or matrix.
 *   \oli Setting some problem-specific constraints (e.g. dimension of
 *     the coarse-grained object).
 *   \oli Solving the constrained SCG problem, that is finding P'.
 *   \oli Computing from P' two semi-projectors L' and R'
 *     (e.g. following the method proposed in [1]).
 *   \oli Working out the product M'=L'MR'* and, if needed, defining
 *     from M' a coarse-grained graph.
 * \endolist
 * </para>
 * </section>
 *
 * <section id="functions-for-performing-scg"><title>Functions for performing SCG</title>
 * <para>
 * The main functions are \ref igraph_scg_adjacency(), \ref
 * igraph_scg_laplacian() and \ref igraph_scg_stochastic().
 * These functions handle all the steps involved in the
 * Spectral Coarse Graining (SCG) of some particular matrices and graphs
 * as described above and in reference [1]. In more details,
 * they compute some prescribed eigenpairs of a matrix or a
 * graph matrix, (for now adjacency, Laplacian and stochastic matrices are
 * available), work out an optimal partition to preserve the eigenpairs,
 * and finally output a coarse-grained matrix or graph along with other
 * useful information.
 * </para>
 *
 * <para>
 * These steps can also be carried out independently: (1) Use
 * \ref igraph_get_adjacency(), \ref igraph_get_sparsemat(),
 * \ref igraph_laplacian(), \ref igraph_get_stochastic() or \ref
 * igraph_get_stochastic_sparsemat() to compute a matrix M.
 * (2) Work out some prescribed eigenpairs of M e.g. by
 * means of \ref igraph_arpack_rssolve() or \ref
 * igraph_arpack_rnsolve(). (3) Invoke one the four
 * algorithms of the function \ref igraph_scg_grouping() to get a
 * partition that will preserve the eigenpairs in the coarse-grained
 * matrix. (4) Compute the semi-projectors L and R using
 * \ref igraph_scg_semiprojectors() and from there the coarse-grained
 * matrix M'=LMR*. If necessary, construct a coarse-grained graph from
 * M' (e.g. as in [1]).
 * </para>
 * </section>
 *
 * <section id="scg-references"><title>References</title>
 * <para>
 * [1] D. Morton de Lachapelle, D. Gfeller, and P. De Los Rios,
 * Shrinking Matrices while Preserving their Eigenpairs with Application
 * to the Spectral Coarse Graining of Graphs. Submitted to
 * <emphasis>SIAM Journal on Matrix Analysis and
 * Applications</emphasis>, 2008.
 * http://people.epfl.ch/david.morton
 * </para>
 * <para>
 * [2] D. Gfeller, and P. De Los Rios, Spectral Coarse Graining and
 * Synchronization in Oscillator Networks.
 * <emphasis>Physical Review Letters</emphasis>,
 * <emphasis role="strong">100</emphasis>(17), 2008.
 * http://arxiv.org/abs/0708.2055
 * </para>
 * <para>
 * [3] D. Gfeller, and P. De Los Rios, Spectral Coarse Graining of Complex
 * Networks, <emphasis>Physical Review Letters</emphasis>,
 * <emphasis role="strong">99</emphasis>(3), 2007.
 * http://arxiv.org/abs/0706.0812
 * </para>
 * </section>
 */

/**
 * \function igraph_scg_grouping
 * \brief SCG problem solver.
 *
 * This function solves the Spectral Coarse Graining (SCG) problem;
 * either exactly, or approximately but faster.
 *
 * </para><para>
 * The algorithm \c IGRAPH_SCG_OPTIMUM solves the SCG problem exactly
 * for each eigenvector in \p V. The running time of this algorithm is
 * O(max(nt) m^2) for the symmetric and Laplacian matrix problems.
 * It is O(m^3) for the stochastic problem. Here m is the number
 * of rows in \p V. In all three cases, the memory usage is O(m^2).
 *
 * </para><para>
 * The algorithms \c IGRAPH_SCG_INTERV and \c IGRAPH_SCG_INTERV_KM solve
 * the SCG problem approximately by performing a (for now) constant
 * binning of the components of the eigenvectors, that is <code>nt_vec[i]</code>
 * constant-size bins are used to partition the <code>i</code>th eigenvector in \c V.
 * When \p algo is \c IGRAPH_SCG_INTERV_KM, the (Lloyd) k-means algorithm is
 * run on each partition obtained by \c IGRAPH_SCG_INTERV to improve
 * accuracy.
 *
 * </para><para>
 * Once a minimizing partition (either exact or approximate) has been
 * found for each eigenvector, the final grouping is worked out as
 * follows: two vertices are grouped together in the final partition if
 * they are grouped together in each minimizing partition. In general, the
 * size of the final partition is not known in advance when the number
 * of columns in \p V is larger than one.
 *
 * </para><para>
 * Finally, the algorithm \c IGRAPH_SCG_EXACT groups the vertices with
 * equal components in each eigenvector. The last three algorithms
 * essentially have linear running time and memory load.
 *
 * \param V The matrix of eigenvectors to be preserved by coarse
 *    graining, each column is an eigenvector.
 * \param groups Pointer to an initialized vector; the result of the
 *    SCG is stored here.
 * \param nt Positive integer. When \p algo is \c IGRAPH_SCG_OPTIMUM,
 *    it gives the number of groups to partition each eigenvector
 *    separately. When \p algo is \c IGRAPH_SCG_INTERV or \c
 *    IGRAPH_SCG_INTERV_KM, it gives the number of intervals to
 *    partition each eigenvector. This is ignored when \p algo is \c
 *    IGRAPH_SCG_EXACT.
 * \param nt_vec May be (1) a numeric vector of length one, or
 *    (2) a vector of the same length as the number of eigenvectors given in \p V, or
 *    (3) a \c NULL pointer.
 *    If not \c NULL, then this argument gives the number of
 *    groups or intervals, and \p nt is ignored. Different number of
 *    groups or intervals can be specified for each eigenvector.
 * \param mtype The type of semi-projectors used in the SCG. Possible
 *    values are \c IGRAPH_SCG_SYMMETRIC, \c IGRAPH_SCG_STOCHASTIC and
 *    \c IGRAPH_SCG_LAPLACIAN.
 * \param algo The algorithm to solve the SCG problem. Possible
 *    values: \c IGRAPH_SCG_OPTIMUM, \c IGRAPH_SCG_INTERV_KM, \c
 *    IGRAPH_SCG_INTERV and \c IGRAPH_SCG_EXACT. Please see the
 *    details about them above.
 * \param p A probability vector, or \c NULL. This argument must be
 *    given if \p mtype is \c IGRAPH_SCG_STOCHASTIC, but it is ignored
 *    otherwise. For the stochastic case it gives the stationary
 *    probability distribution of a Markov chain, the one specified by
 *    the graph/matrix under study.
 * \param maxiter A positive integer giving the number of iterations
 *    of the k-means algorithm when \p algo is \c
 *    IGRAPH_SCG_INTERV_KM. It is ignored in other cases. A reasonable
 *    (initial) value for this argument is 100.
 * \return Error code.
 *
 * Time complexity: see description above.
 *
 * \sa \ref igraph_scg_adjacency(), \ref igraph_scg_laplacian(), \ref
 * igraph_scg_stochastic().
 *
 * \example examples/simple/igraph_scg_grouping.c
 * \example examples/simple/igraph_scg_grouping2.c
 * \example examples/simple/igraph_scg_grouping3.c
 * \example examples/simple/igraph_scg_grouping4.c
 */

int igraph_scg_grouping(const igraph_matrix_t *V,
                        igraph_vector_t *groups,
                        igraph_integer_t nt,
                        const igraph_vector_t *nt_vec,
                        igraph_scg_matrix_t mtype,
                        igraph_scg_algorithm_t algo,
                        const igraph_vector_t *p,
                        igraph_integer_t maxiter) {

    int no_of_nodes = (int) igraph_matrix_nrow(V);
    int nev = (int) igraph_matrix_ncol(V);
    igraph_matrix_int_t gr_mat;
    int i;

    if (nt_vec && igraph_vector_size(nt_vec) != 1 &&
        igraph_vector_size(nt_vec) != nev) {
        IGRAPH_ERROR("Invalid length for interval specification", IGRAPH_EINVAL);
    }
    if (nt_vec && igraph_vector_size(nt_vec) == 1) {
        nt = (igraph_integer_t) VECTOR(*nt_vec)[0];
        nt_vec = 0;
    }

    if (!nt_vec && algo != IGRAPH_SCG_EXACT) {
        if (nt <= 1 || nt >= no_of_nodes) {
            IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
        }
    } else if (algo != IGRAPH_SCG_EXACT) {
        igraph_real_t min, max;
        igraph_vector_minmax(nt_vec, &min, &max);
        if (min <= 1 || max >= no_of_nodes) {
            IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
        }
    }

    if (mtype == IGRAPH_SCG_STOCHASTIC && !p) {
        IGRAPH_ERROR("The p vector must be given for the stochastic matrix case",
                     IGRAPH_EINVAL);
    }

    if (p) {
        if (igraph_vector_size(p) != no_of_nodes) {
            IGRAPH_ERROR("Invalid p vector size", IGRAPH_EINVAL);
        }

        if (igraph_vector_min(p) < 0) {
            IGRAPH_ERROR("The elements of the p vector must be non-negative", IGRAPH_EINVAL);
        }
    }

    IGRAPH_CHECK(igraph_vector_resize(groups, no_of_nodes));

#define INVEC(i) (nt_vec ? VECTOR(*nt_vec)[i] : nt)

    IGRAPH_CHECK(igraph_matrix_int_init(&gr_mat, no_of_nodes, nev));
    IGRAPH_FINALLY(igraph_matrix_int_destroy, &gr_mat);

    switch (algo) {
    case IGRAPH_SCG_OPTIMUM:
        for (i = 0; i < nev; i++) {
            IGRAPH_CHECK(igraph_i_optimal_partition(&MATRIX(*V, 0, i),
                                                    &MATRIX(gr_mat, 0, i),
                                                    no_of_nodes, (int) INVEC(i),
                                                    mtype,
                                                    p ? VECTOR(*p) : 0, 0));
        }
        break;
    case IGRAPH_SCG_INTERV_KM:
        for (i = 0; i < nev; i++) {
            igraph_vector_t tmpv;
            igraph_vector_view(&tmpv, &MATRIX(*V, 0, i), no_of_nodes);
            IGRAPH_CHECK(igraph_i_intervals_plus_kmeans(&tmpv,
                         &MATRIX(gr_mat, 0, i),
                         no_of_nodes, (int) INVEC(i),
                         maxiter));
        }
        break;
    case IGRAPH_SCG_INTERV:
        for (i = 0; i < nev; i++) {
            igraph_vector_t tmpv;
            igraph_vector_view(&tmpv, &MATRIX(*V, 0, i), no_of_nodes);
            IGRAPH_CHECK(igraph_i_intervals_method(&tmpv,
                                                   &MATRIX(gr_mat, 0, i),
                                                   no_of_nodes, (int) INVEC(i)));
        }
        break;
    case IGRAPH_SCG_EXACT:
        for (i = 0; i < nev; i++) {
            IGRAPH_CHECK(igraph_i_exact_coarse_graining(&MATRIX(*V, 0, i),
                         &MATRIX(gr_mat, 0, i),
                         no_of_nodes));
        }
        break;
    }

#undef INVEC

    if (nev == 1) {
        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*groups)[i] = MATRIX(gr_mat, i, 0);
        }
    } else {
        igraph_i_scg_groups_t *g;
        int gr_nb = 0;
        
        g = IGRAPH_CALLOC(no_of_nodes, igraph_i_scg_groups_t);
        IGRAPH_FINALLY(igraph_free, g);

        IGRAPH_CHECK(igraph_matrix_int_transpose(&gr_mat));
        for (i = 0; i < no_of_nodes; i++) {
            g[i].ind = i;
            g[i].n = nev;
            g[i].gr = &MATRIX(gr_mat, 0, i);
        }

        igraph_qsort(g, (size_t) no_of_nodes, sizeof(igraph_i_scg_groups_t),
              igraph_i_compare_groups);
        VECTOR(*groups)[g[0].ind] = gr_nb;
        for (i = 1; i < no_of_nodes; i++) {
            if (igraph_i_compare_groups(&g[i], &g[i - 1]) != 0) {
                gr_nb++;
            }
            VECTOR(*groups)[g[i].ind] = gr_nb;
        }
        
        IGRAPH_FREE(g);
        IGRAPH_FINALLY_CLEAN(1);
    }

    igraph_matrix_int_destroy(&gr_mat);
    IGRAPH_FINALLY_CLEAN(1);

    IGRAPH_CHECK(igraph_reindex_membership(groups, 0, 0));

    return 0;
}

static int igraph_i_scg_semiprojectors_sym(const igraph_vector_t *groups,
                                           igraph_matrix_t *L,
                                           igraph_matrix_t *R,
                                           igraph_sparsemat_t *Lsparse,
                                           igraph_sparsemat_t *Rsparse,
                                           int no_of_groups,
                                           int no_of_nodes) {

    igraph_vector_t tab;
    int i;

    IGRAPH_VECTOR_INIT_FINALLY(&tab, no_of_groups);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(tab)[ (int) VECTOR(*groups)[i] ] += 1;
    }
    for (i = 0; i < no_of_groups; i++) {
        VECTOR(tab)[i] = sqrt(VECTOR(tab)[i]);
    }

    if (L) {
        IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
        igraph_matrix_null(L);
        for (i = 0; i < no_of_nodes; i++) {
            int g = (int) VECTOR(*groups)[i];
            MATRIX(*L, g, i) = 1 / VECTOR(tab)[g];
        }
    }

    if (R) {
        if (L) {
            IGRAPH_CHECK(igraph_matrix_update(R, L));
        } else {
            IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
            igraph_matrix_null(R);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*R, g, i) = 1 / VECTOR(tab)[g];
            }
        }
    }

    if (Lsparse) {
        IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
                                           /* nzmax= */ no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            int g = (int) VECTOR(*groups)[i];
            IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i, 1 / VECTOR(tab)[g]));
        }
    }

    if (Rsparse) {
        IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
                                           /* nzmax= */ no_of_nodes));
        for (i = 0; i < no_of_nodes; i++) {
            int g = (int) VECTOR(*groups)[i];
            IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i, 1 / VECTOR(tab)[g]));
        }
    }

    igraph_vector_destroy(&tab);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_scg_semiprojectors_lap(const igraph_vector_t *groups,
                                           igraph_matrix_t *L,
                                           igraph_matrix_t *R,
                                           igraph_sparsemat_t *Lsparse,
                                           igraph_sparsemat_t *Rsparse,
                                           int no_of_groups,
                                           int no_of_nodes,
                                           igraph_scg_norm_t norm) {

    igraph_vector_t tab;
    int i;

    IGRAPH_VECTOR_INIT_FINALLY(&tab, no_of_groups);
    for (i = 0; i < no_of_nodes; i++) {
        VECTOR(tab)[ (int) VECTOR(*groups)[i] ] += 1;
    }
    for (i = 0; i < no_of_groups; i++) {
        VECTOR(tab)[i] = VECTOR(tab)[i];
    }

    if (norm == IGRAPH_SCG_NORM_ROW) {
        if (L) {
            IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
            igraph_matrix_null(L);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*L, g, i) = 1.0 / VECTOR(tab)[g];
            }
        }
        if (R) {
            IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
            igraph_matrix_null(R);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*R, g, i) = 1.0;
            }
        }
        if (Lsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i,
                                                    1.0 / VECTOR(tab)[g]));
            }
        }
        if (Rsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i, 1.0));
            }
        }
    } else {
        if (L) {
            IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
            igraph_matrix_null(L);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*L, g, i) = 1.0;
            }
        }
        if (R) {
            IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
            igraph_matrix_null(R);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*R, g, i) = 1.0 / VECTOR(tab)[g];
            }
        }
        if (Lsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i, 1.0));
            }
        }
        if (Rsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i,
                                                    1.0 / VECTOR(tab)[g]));
            }
        }

    }

    igraph_vector_destroy(&tab);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_scg_semiprojectors_sto(const igraph_vector_t *groups,
                                           igraph_matrix_t *L,
                                           igraph_matrix_t *R,
                                           igraph_sparsemat_t *Lsparse,
                                           igraph_sparsemat_t *Rsparse,
                                           int no_of_groups,
                                           int no_of_nodes,
                                           const igraph_vector_t *p,
                                           igraph_scg_norm_t norm) {

    igraph_vector_t pgr, pnormed;
    int i;

    IGRAPH_VECTOR_INIT_FINALLY(&pgr, no_of_groups);
    IGRAPH_VECTOR_INIT_FINALLY(&pnormed, no_of_nodes);
    for (i = 0; i < no_of_nodes; i++) {
        int g = (int) VECTOR(*groups)[i];
        VECTOR(pgr)[g] += VECTOR(*p)[i];
    }
    for (i = 0; i < no_of_nodes; i++) {
        int g = (int) VECTOR(*groups)[i];
        VECTOR(pnormed)[i] = VECTOR(*p)[i] / VECTOR(pgr)[g];
    }

    if (norm == IGRAPH_SCG_NORM_ROW) {
        if (L) {
            IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
            igraph_matrix_null(L);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*L, g, i) = VECTOR(pnormed)[i];
            }
        }
        if (R) {
            IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
            igraph_matrix_null(R);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*R, g, i) = 1.0;
            }
        }
        if (Lsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i,
                                                    VECTOR(pnormed)[i]));
            }
        }
        if (Rsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i, 1.0));
            }
        }
    } else {
        if (L) {
            IGRAPH_CHECK(igraph_matrix_resize(L, no_of_groups, no_of_nodes));
            igraph_matrix_null(L);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int ) VECTOR(*groups)[i];
                MATRIX(*L, g, i) = 1.0;
            }
        }
        if (R) {
            IGRAPH_CHECK(igraph_matrix_resize(R, no_of_groups, no_of_nodes));
            igraph_matrix_null(R);
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                MATRIX(*R, g, i) = VECTOR(pnormed)[i];
            }
        }
        if (Lsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Lsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Lsparse, g, i, 1.0));
            }
        }
        if (Rsparse) {
            IGRAPH_CHECK(igraph_sparsemat_init(Rsparse, no_of_groups, no_of_nodes,
                                               /* nzmax= */ no_of_nodes));
            for (i = 0; i < no_of_nodes; i++) {
                int g = (int) VECTOR(*groups)[i];
                IGRAPH_CHECK(igraph_sparsemat_entry(Rsparse, g, i,
                                                    VECTOR(pnormed)[i]));
            }
        }
    }


    igraph_vector_destroy(&pnormed);
    igraph_vector_destroy(&pgr);
    IGRAPH_FINALLY_CLEAN(2);

    return 0;
}

/**
 * \function igraph_scg_semiprojectors
 * \brief Compute SCG semi-projectors for a given partition.
 *
 * The three types of semi-projectors are defined as follows.
 * Let gamma(j) label the group of vertex j in a partition of all the
 * vertices.
 *
 * </para><para>
 * The symmetric semi-projectors are defined as
 * <blockquote><para><phrase role="math">
 *   L[alpha,j] = R[alpha,j] = 1/sqrt(|alpha|) delta[alpha,gamma(j)],
 * </phrase></para></blockquote>
 * the (row) Laplacian semi-projectors as
 * <blockquote><para><phrase role="math">
 *   L[alpha,j] = 1/|alpha| delta[alpha,gamma(j)]
 * </phrase></para></blockquote>
 * and
 * <blockquote><para><phrase role="math">
 *   R[alpha,j] = delta[alpha,gamma(j)],
 * </phrase></para></blockquote>
 * and the (row) stochastic semi-projectors as
 * <blockquote><para><phrase role="math">
 *     L[alpha,j] = p[1][j] / sum(p[1][k]; k in gamma(j))
 *     delta[alpha,gamma(j)]
 * </phrase></para></blockquote>
 * and
 * <blockquote><para><phrase role="math">
 *     R[alpha,j] = delta[alpha,gamma(j)],
 * </phrase></para></blockquote>
 * where p[1] is the (left) eigenvector associated with the
 * one-eigenvalue of the stochastic matrix. L and R are
 * defined in a symmetric way when \p norm is \c
 * IGRAPH_SCG_NORM_COL. All these semi-projectors verify various
 * properties described in the reference.
 * \param groups A vector of integers, giving the group label of every
 *    vertex in the partition. Group labels should start at zero and
 *    should be sequential.
 * \param mtype The type of semi-projectors. For now \c
 *    IGRAPH_SCG_SYMMETRIC, \c IGRAPH_SCG_STOCHASTIC and \c
 *    IGRAP_SCG_LAPLACIAN are supported.
 * \param L If not a \c NULL pointer, then it must be a pointer to
 *    an initialized matrix. The left semi-projector is stored here.
 * \param R If not a \c NULL pointer, then it must be a pointer to
 *    an initialized matrix. The right semi-projector is stored here.
 * \param Lsparse If not a \c NULL pointer, then it must be a pointer
 *    to an uninitialized sparse matrix. The left semi-projector is
 *    stored here.
 * \param Rsparse If not a \c NULL pointer, then it must be a pointer
 *    to an uninitialized sparse matrix. The right semi-projector is
 *    stored here.
 * \param p \c NULL, or a probability vector of the same length as \p
 *    groups. \p p is the stationary probability distribution of a
 *    Markov chain when \p mtype is \c IGRAPH_SCG_STOCHASTIC. This
 *    argument is ignored in all other cases.
 * \param norm Either \c IGRAPH_SCG_NORM_ROW or \c IGRAPH_SCG_NORM_COL.
 *    Specifies whether the rows or the columns of the Laplacian
 *    matrix sum up to zero, or whether the rows or the columns of the
 *    stochastic matrix sum up to one.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \sa \ref igraph_scg_adjacency(), \ref igraph_scg_stochastic() and
 * \ref igraph_scg_laplacian(), \ref igraph_scg_grouping().
 *
 * \example examples/simple/igraph_scg_semiprojectors.c
 * \example examples/simple/igraph_scg_semiprojectors2.c
 * \example examples/simple/igraph_scg_semiprojectors3.c
 */

int igraph_scg_semiprojectors(const igraph_vector_t *groups,
                              igraph_scg_matrix_t mtype,
                              igraph_matrix_t *L,
                              igraph_matrix_t *R,
                              igraph_sparsemat_t *Lsparse,
                              igraph_sparsemat_t *Rsparse,
                              const igraph_vector_t *p,
                              igraph_scg_norm_t norm) {

    int no_of_nodes = (int) igraph_vector_size(groups);
    int no_of_groups;
    igraph_real_t min, max;

    igraph_vector_minmax(groups, &min, &max);
    no_of_groups = (int) max + 1;

    if (min < 0 || max >= no_of_nodes) {
        IGRAPH_ERROR("Invalid membership vector", IGRAPH_EINVAL);
    }

    if (mtype == IGRAPH_SCG_STOCHASTIC && !p) {
        IGRAPH_ERROR("`p' must be given for the stochastic matrix case",
                     IGRAPH_EINVAL);
    }

    if (p && igraph_vector_size(p) != no_of_nodes) {
        IGRAPH_ERROR("Invalid `p' vector length, should match number of vertices",
                     IGRAPH_EINVAL);
    }

    switch (mtype) {
    case IGRAPH_SCG_SYMMETRIC:
        IGRAPH_CHECK(igraph_i_scg_semiprojectors_sym(groups, L, R, Lsparse,
                     Rsparse, no_of_groups,
                     no_of_nodes));
        break;

    case IGRAPH_SCG_LAPLACIAN:
        IGRAPH_CHECK(igraph_i_scg_semiprojectors_lap(groups, L, R, Lsparse,
                     Rsparse, no_of_groups,
                     no_of_nodes, norm));
        break;

    case IGRAPH_SCG_STOCHASTIC:
        IGRAPH_CHECK(igraph_i_scg_semiprojectors_sto(groups, L, R, Lsparse,
                     Rsparse, no_of_groups,
                     no_of_nodes, p, norm));
        break;
    }

    return 0;
}

/**
 * \function igraph_scg_norm_eps
 * \brief Calculate SCG residuals.
 *
 * Computes |v[i]-Pv[i]|, where v[i] is the i-th eigenvector in \p V
 * and P is the projector corresponding to the \p mtype argument.
 *
 * \param V The matrix of eigenvectors to be preserved by coarse
 *    graining, each column is an eigenvector.
 * \param groups A vector of integers, giving the group label of every
 *    vertex in the partition. Group labels should start at zero and
 *    should be sequential.
 * \param eps Pointer to a real value, the result is stored here.
 * \param mtype The type of semi-projectors. For now \c
 *    IGRAPH_SCG_SYMMETRIC, \c IGRAPH_SCG_STOCHASTIC and \c
 *    IGRAP_SCG_LAPLACIAN are supported.
 * \param p \c NULL, or a probability vector of the same length as \p
 *    groups. \p p is the stationary probability distribution of a
 *    Markov chain when \p mtype is \c IGRAPH_SCG_STOCHASTIC. This
 *    argument is ignored in all other cases.
 * \param norm Either \c IGRAPH_SCG_NORM_ROW or \c IGRAPH_SCG_NORM_COL.
 *    Specifies whether the rows or the columns of the Laplacian
 *    matrix sum up to zero, or whether the rows or the columns of the
 *    stochastic matrix sum up to one.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \sa \ref igraph_scg_adjacency(), \ref igraph_scg_stochastic() and
 * \ref igraph_scg_laplacian(), \ref igraph_scg_grouping(), \ref
 * igraph_scg_semiprojectors().
 */

int igraph_scg_norm_eps(const igraph_matrix_t *V,
                        const igraph_vector_t *groups,
                        igraph_vector_t *eps,
                        igraph_scg_matrix_t mtype,
                        const igraph_vector_t *p,
                        igraph_scg_norm_t norm) {

    int no_of_nodes = (int) igraph_vector_size(groups);
    int no_of_vectors = (int) igraph_matrix_ncol(V);
    igraph_real_t min, max;
    igraph_sparsemat_t Lsparse, Rsparse, Lsparse2, Rsparse2, Rsparse3, proj;
    igraph_vector_t x, res;
    int k, i;

    if (igraph_matrix_nrow(V) != no_of_nodes) {
        IGRAPH_ERROR("Eigenvector length and group vector length do not match",
                     IGRAPH_EINVAL);
    }

    igraph_vector_minmax(groups, &min, &max);

    if (min < 0 || max >= no_of_nodes) {
        IGRAPH_ERROR("Invalid membership vector", IGRAPH_EINVAL);
    }

    if (mtype == IGRAPH_SCG_STOCHASTIC && !p) {
        IGRAPH_ERROR("`p' must be given for the stochastic matrix case",
                     IGRAPH_EINVAL);
    }

    if (p && igraph_vector_size(p) != no_of_nodes) {
        IGRAPH_ERROR("Invalid `p' vector length, should match number of vertices",
                     IGRAPH_EINVAL);
    }

    IGRAPH_CHECK(igraph_scg_semiprojectors(groups, mtype, /* L= */ 0,
                                           /* R= */ 0, &Lsparse, &Rsparse, p,
                                           norm));

    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Lsparse);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse);

    IGRAPH_CHECK(igraph_sparsemat_compress(&Lsparse, &Lsparse2));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Lsparse2);
    IGRAPH_CHECK(igraph_sparsemat_compress(&Rsparse, &Rsparse2));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse2);
    IGRAPH_CHECK(igraph_sparsemat_transpose(&Rsparse2, &Rsparse3,
                                            /*values=*/ 1));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse3);

    IGRAPH_CHECK(igraph_sparsemat_multiply(&Rsparse3, &Lsparse2, &proj));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &proj);

    IGRAPH_VECTOR_INIT_FINALLY(&res, no_of_nodes);
    IGRAPH_CHECK(igraph_vector_resize(eps, no_of_vectors));

    for (k = 0; k < no_of_vectors; k++) {
        igraph_vector_view(&x, &MATRIX(*V, 0, k), no_of_nodes);
        igraph_vector_null(&res);
        IGRAPH_CHECK(igraph_sparsemat_gaxpy(&proj, &x, &res));
        VECTOR(*eps)[k] = 0.0;
        for (i = 0; i < no_of_nodes; i++) {
            igraph_real_t di = MATRIX(*V, i, k) - VECTOR(res)[i];
            VECTOR(*eps)[k] += di * di;
        }
        VECTOR(*eps)[k] = sqrt(VECTOR(*eps)[k]);
    }

    igraph_vector_destroy(&res);
    igraph_sparsemat_destroy(&proj);
    igraph_sparsemat_destroy(&Rsparse3);
    igraph_sparsemat_destroy(&Rsparse2);
    igraph_sparsemat_destroy(&Lsparse2);
    igraph_sparsemat_destroy(&Rsparse);
    igraph_sparsemat_destroy(&Lsparse);
    IGRAPH_FINALLY_CLEAN(7);

    return 0;
}

static int igraph_i_matrix_laplacian(const igraph_matrix_t *matrix,
                                     igraph_matrix_t *mymatrix,
                                     igraph_scg_norm_t norm) {

    igraph_vector_t degree;
    int i, j, n = (int) igraph_matrix_nrow(matrix);
    IGRAPH_CHECK(igraph_matrix_resize(mymatrix, n, n));

    IGRAPH_VECTOR_INIT_FINALLY(&degree, n);

    if (norm == IGRAPH_SCG_NORM_ROW) {
        IGRAPH_CHECK(igraph_matrix_rowsum(matrix, &degree));
    } else {
        IGRAPH_CHECK(igraph_matrix_colsum(matrix, &degree));
    }
    for (i = 0; i < n; i++) {
        VECTOR(degree)[i] -= MATRIX(*matrix, i, i);
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            MATRIX(*mymatrix, i, j) = - MATRIX(*matrix, i, j);
        }
        MATRIX(*mymatrix, i, i) = VECTOR(degree)[i];
    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_sparsemat_laplacian(const igraph_sparsemat_t *sparse,
                                        igraph_sparsemat_t *mysparse,
                                        igraph_scg_norm_t norm) {

    igraph_vector_t degree;
    int i, n = (int) igraph_sparsemat_nrow(sparse);
    int nzmax = igraph_sparsemat_nzmax(sparse);
    igraph_sparsemat_iterator_t it;

    IGRAPH_CHECK(igraph_sparsemat_init(mysparse, n, n, nzmax + n));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparse);
    igraph_sparsemat_iterator_init(&it, (igraph_sparsemat_t *) sparse);

    IGRAPH_VECTOR_INIT_FINALLY(&degree, n);
    for (igraph_sparsemat_iterator_reset(&it);
         !igraph_sparsemat_iterator_end(&it);
         igraph_sparsemat_iterator_next(&it)) {
        int row = igraph_sparsemat_iterator_row(&it);
        int col = igraph_sparsemat_iterator_col(&it);
        if (row != col) {
            igraph_real_t val = igraph_sparsemat_iterator_get(&it);
            if (norm == IGRAPH_SCG_NORM_ROW) {
                VECTOR(degree)[row] += val;
            } else {
                VECTOR(degree)[col] += val;
            }
        }
    }

    /* Diagonal */
    for (i = 0; i < n; i++) {
        igraph_sparsemat_entry(mysparse, i, i, VECTOR(degree)[i]);
    }

    /* And the rest, filter out diagonal elements */
    for (igraph_sparsemat_iterator_reset(&it);
         !igraph_sparsemat_iterator_end(&it);
         igraph_sparsemat_iterator_next(&it)) {
        int row = igraph_sparsemat_iterator_row(&it);
        int col = igraph_sparsemat_iterator_col(&it);
        if (row != col) {
            igraph_real_t val = igraph_sparsemat_iterator_get(&it);
            igraph_sparsemat_entry(mysparse, row, col, -val);
        }
    }

    igraph_vector_destroy(&degree);
    IGRAPH_FINALLY_CLEAN(2);  /* + mysparse */

    return 0;
}

static int igraph_i_matrix_stochastic(const igraph_matrix_t *matrix,
                                      igraph_matrix_t *mymatrix,
                                      igraph_scg_norm_t norm) {

    int i, j, n = (int) igraph_matrix_nrow(matrix);
    IGRAPH_CHECK(igraph_matrix_copy(mymatrix, matrix));

    if (norm == IGRAPH_SCG_NORM_ROW) {
        for (i = 0; i < n; i++) {
            igraph_real_t sum = 0.0;
            for (j = 0; j < n; j++) {
                sum += MATRIX(*matrix, i, j);
            }
            if (sum == 0) {
                IGRAPH_WARNING("Zero degree vertices");
            }
            for (j = 0; j < n; j++) {
                MATRIX(*mymatrix, i, j) = MATRIX(*matrix, i, j) / sum;
            }
        }
    } else {
        for (i = 0; i < n; i++) {
            igraph_real_t sum = 0.0;
            for (j = 0; j < n; j++) {
                sum += MATRIX(*matrix, j, i);
            }
            if (sum == 0) {
                IGRAPH_WARNING("Zero degree vertices");
            }
            for (j = 0; j < n; j++) {
                MATRIX(*mymatrix, j, i) = MATRIX(*matrix, j, i) / sum;
            }
        }
    }

    return 0;
}

static int igraph_i_sparsemat_stochastic(const igraph_sparsemat_t *sparse,
                                         igraph_sparsemat_t *mysparse,
                                         igraph_scg_norm_t norm) {

    IGRAPH_CHECK(igraph_sparsemat_copy(mysparse, sparse));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparse);
    IGRAPH_CHECK(igraph_i_normalize_sparsemat(mysparse,
                 norm == IGRAPH_SCG_NORM_COL));
    IGRAPH_FINALLY_CLEAN(1);

    return 0;
}

static int igraph_i_scg_get_result(igraph_scg_matrix_t type,
                                   const igraph_matrix_t *matrix,
                                   const igraph_sparsemat_t *sparsemat,
                                   const igraph_sparsemat_t *Lsparse,
                                   const igraph_sparsemat_t *Rsparse_t,
                                   igraph_t *scg_graph,
                                   igraph_matrix_t *scg_matrix,
                                   igraph_sparsemat_t *scg_sparsemat,
                                   igraph_bool_t directed) {

    /* We need to calculate either scg_matrix (if input is dense), or
       scg_sparsemat (if input is sparse). For the latter we might need
       to temporarily use another matrix. */


    if (matrix) {
        igraph_matrix_t *my_scg_matrix = scg_matrix, v_scg_matrix;
        igraph_matrix_t tmp;
        igraph_sparsemat_t *myLsparse = (igraph_sparsemat_t *) Lsparse, v_Lsparse;

        if (!scg_matrix) {
            my_scg_matrix = &v_scg_matrix;
            IGRAPH_CHECK(igraph_matrix_init(my_scg_matrix, 0, 0));
            IGRAPH_FINALLY(igraph_matrix_destroy, my_scg_matrix);
        }

        if (!igraph_sparsemat_is_cc(Lsparse)) {
            myLsparse = &v_Lsparse;
            IGRAPH_CHECK(igraph_sparsemat_compress(Lsparse, myLsparse));
            IGRAPH_FINALLY(igraph_sparsemat_destroy, myLsparse);
        }

        IGRAPH_CHECK(igraph_matrix_init(&tmp, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_destroy, &tmp);
        IGRAPH_CHECK(igraph_sparsemat_dense_multiply(matrix, Rsparse_t, &tmp));
        IGRAPH_CHECK(igraph_sparsemat_multiply_by_dense(myLsparse, &tmp,
                     my_scg_matrix));
        igraph_matrix_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);

        if (scg_sparsemat) {
            IGRAPH_CHECK(igraph_matrix_as_sparsemat(scg_sparsemat, my_scg_matrix,
                                                    /* tol= */ 0));
            IGRAPH_FINALLY(igraph_sparsemat_destroy, scg_sparsemat);
        }

        if (scg_graph) {
            if (type != IGRAPH_SCG_LAPLACIAN) {
                IGRAPH_CHECK(igraph_weighted_adjacency(scg_graph, my_scg_matrix,
                                                       directed ?
                                                       IGRAPH_ADJ_DIRECTED :
                                                       IGRAPH_ADJ_UNDIRECTED,
                                                       "weight", /*loops=*/ 1));
            } else {
                int i, j, n = (int) igraph_matrix_nrow(my_scg_matrix);
                igraph_matrix_t tmp;
                IGRAPH_MATRIX_INIT_FINALLY(&tmp, n, n);
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        MATRIX(tmp, i, j) = -MATRIX(*my_scg_matrix, i, j);
                    }
                    MATRIX(tmp, i, i) = 0;
                }
                IGRAPH_CHECK(igraph_weighted_adjacency(scg_graph, &tmp, directed ?
                                                       IGRAPH_ADJ_DIRECTED :
                                                       IGRAPH_ADJ_UNDIRECTED,
                                                       "weight", /*loops=*/ 0));
                igraph_matrix_destroy(&tmp);
                IGRAPH_FINALLY_CLEAN(1);
            }
            IGRAPH_FINALLY(igraph_destroy, scg_graph);
        }

        if (scg_graph)     {
            IGRAPH_FINALLY_CLEAN(1);
        }
        if (scg_sparsemat) {
            IGRAPH_FINALLY_CLEAN(1);
        }

        if (!igraph_sparsemat_is_cc(Lsparse)) {
            igraph_sparsemat_destroy(myLsparse);
            IGRAPH_FINALLY_CLEAN(1);
        }

        if (!scg_matrix) {
            igraph_matrix_destroy(my_scg_matrix);
            IGRAPH_FINALLY_CLEAN(1);
        }

    } else { /* sparsemat */
        igraph_sparsemat_t *my_scg_sparsemat = scg_sparsemat, v_scg_sparsemat;
        igraph_sparsemat_t tmp, *mysparsemat = (igraph_sparsemat_t *) sparsemat,
                                 v_sparsemat, *myLsparse = (igraph_sparsemat_t *) Lsparse, v_Lsparse;
        if (!scg_sparsemat) {
            my_scg_sparsemat = &v_scg_sparsemat;
        }
        if (!igraph_sparsemat_is_cc(sparsemat)) {
            mysparsemat = &v_sparsemat;
            IGRAPH_CHECK(igraph_sparsemat_compress(sparsemat, mysparsemat));
            IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
        }
        if (!igraph_sparsemat_is_cc(Lsparse)) {
            myLsparse = &v_Lsparse;
            IGRAPH_CHECK(igraph_sparsemat_compress(Lsparse, myLsparse));
            IGRAPH_FINALLY(igraph_sparsemat_destroy, myLsparse);
        }
        IGRAPH_CHECK(igraph_sparsemat_multiply(mysparsemat, Rsparse_t,
                                               &tmp));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
        IGRAPH_CHECK(igraph_sparsemat_multiply(myLsparse, &tmp,
                                               my_scg_sparsemat));
        igraph_sparsemat_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
        IGRAPH_FINALLY(igraph_sparsemat_destroy, my_scg_sparsemat);

        if (scg_matrix) {
            IGRAPH_CHECK(igraph_sparsemat_as_matrix(scg_matrix, my_scg_sparsemat));
        }
        if (scg_graph) {
            if (type != IGRAPH_SCG_LAPLACIAN) {
                IGRAPH_CHECK(igraph_weighted_sparsemat(scg_graph, my_scg_sparsemat,
                                                       directed, "weight",
                                                       /*loops=*/ 1));
            } else {
                igraph_sparsemat_t tmp;
                IGRAPH_CHECK(igraph_sparsemat_copy(&tmp, my_scg_sparsemat));
                IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmp);
                IGRAPH_CHECK(igraph_sparsemat_neg(&tmp));
                IGRAPH_CHECK(igraph_weighted_sparsemat(scg_graph, &tmp, directed,
                                                       "weight", /*loops=*/ 0));
                igraph_sparsemat_destroy(&tmp);
                IGRAPH_FINALLY_CLEAN(1);
            }
            IGRAPH_FINALLY(igraph_destroy, scg_graph);
        }

        if (scg_graph) {
            IGRAPH_FINALLY_CLEAN(1);
        }
        if (!scg_sparsemat) {
            igraph_sparsemat_destroy(my_scg_sparsemat);
        }
        IGRAPH_FINALLY_CLEAN(1);    /* my_scg_sparsemat */
        if (!igraph_sparsemat_is_cc(Lsparse)) {
            igraph_sparsemat_destroy(myLsparse);
            IGRAPH_FINALLY_CLEAN(1);
        }
        if (!igraph_sparsemat_is_cc(sparsemat)) {
            igraph_sparsemat_destroy(mysparsemat);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    return 0;
}

static int igraph_i_scg_common_checks(const igraph_t *graph,
                                      const igraph_matrix_t *matrix,
                                      const igraph_sparsemat_t *sparsemat,
                                      const igraph_vector_t *ev,
                                      igraph_integer_t nt,
                                      const igraph_vector_t *nt_vec,
                                      const igraph_matrix_t *vectors,
                                      const igraph_matrix_complex_t *vectors_cmplx,
                                      const igraph_vector_t *groups,
                                      const igraph_t *scg_graph,
                                      const igraph_matrix_t *scg_matrix,
                                      const igraph_sparsemat_t *scg_sparsemat,
                                      const igraph_vector_t *p,
                                      igraph_real_t *evmin, igraph_real_t *evmax) {

    int no_of_nodes = -1;
    igraph_real_t min, max;
    int no_of_ev = (int) igraph_vector_size(ev);

    if ( (graph ? 1 : 0) + (matrix ? 1 : 0) + (sparsemat ? 1 : 0) != 1 ) {
        IGRAPH_ERROR("Give exactly one of `graph', `matrix' and `sparsemat'",
                     IGRAPH_EINVAL);
    }

    if (graph) {
        no_of_nodes = igraph_vcount(graph);
    } else if (matrix) {
        no_of_nodes = (int) igraph_matrix_nrow(matrix);
    } else if (sparsemat) {
        no_of_nodes = (int) igraph_sparsemat_nrow(sparsemat);
    }

    if ((matrix && igraph_matrix_ncol(matrix) != no_of_nodes) ||
        (sparsemat && igraph_sparsemat_ncol(sparsemat) != no_of_nodes)) {
        IGRAPH_ERROR("Matrix must be square", IGRAPH_NONSQUARE);
    }

    igraph_vector_minmax(ev, evmin, evmax);
    if (*evmin < 0 || *evmax >= no_of_nodes) {
        IGRAPH_ERROR("Invalid eigenvectors given", IGRAPH_EINVAL);
    }

    if (!nt_vec && (nt <= 1 || nt >= no_of_nodes)) {
        IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
    }

    if (nt_vec) {
        if (igraph_vector_size(nt_vec) != 1 &&
            igraph_vector_size(nt_vec) != no_of_ev) {
            IGRAPH_ERROR("Invalid length for interval specification",
                         IGRAPH_EINVAL);
        }
        igraph_vector_minmax(nt_vec, &min, &max);
        if (min <= 1 || max >= no_of_nodes) {
            IGRAPH_ERROR("Invalid interval specification", IGRAPH_EINVAL);
        }
    }

    if (vectors && igraph_matrix_size(vectors) != 0 &&
        (igraph_matrix_ncol(vectors) != no_of_ev ||
         igraph_matrix_nrow(vectors) != no_of_nodes)) {
        IGRAPH_ERROR("Invalid eigenvector matrix size", IGRAPH_EINVAL);
    }

    if (vectors_cmplx && igraph_matrix_complex_size(vectors_cmplx) != 0 &&
        (igraph_matrix_complex_ncol(vectors_cmplx) != no_of_ev ||
         igraph_matrix_complex_nrow(vectors_cmplx) != no_of_nodes)) {
        IGRAPH_ERROR("Invalid eigenvector matrix size", IGRAPH_EINVAL);
    }

    if (groups && igraph_vector_size(groups) != 0 &&
        igraph_vector_size(groups) != no_of_nodes) {
        IGRAPH_ERROR("Invalid `groups' vector size", IGRAPH_EINVAL);
    }

    if ( (scg_graph != 0) + (scg_matrix != 0) + (scg_sparsemat != 0) == 0 ) {
        IGRAPH_ERROR("No output is requested, please give at least one of "
                     "`scg_graph', `scg_matrix' and `scg_sparsemat'",
                     IGRAPH_EINVAL);
    }

    if (p && igraph_vector_size(p) != 0 &&
        igraph_vector_size(p) != no_of_nodes) {
        IGRAPH_ERROR("Invalid `p' vector size", IGRAPH_EINVAL);
    }

    return 0;
}

/**
 * \function igraph_scg_adjacency
 * Spectral coarse graining, symmetric case.
 *
 * This function handles all the steps involved in the Spectral Coarse
 * Graining (SCG) of some matrices and graphs as described in the
 * reference below.
 *
 * \param graph The input graph. Exactly one of \p graph, \p matrix
 *    and \p sparsemat must be given, the other two must be \c NULL
 *    pointers.
 * \param matrix The input matrix. Exactly one of \p graph, \p matrix
 *    and \p sparsemat must be given, the other two must be \c NULL
 *    pointers.
 * \param sparsemat The input sparse matrix. Exactly one of \p graph,
 *    \p matrix and \p sparsemat must be given, the other two must be
 *    \c NULL pointers.
 * \param ev A vector of positive integers giving the indexes of the
 *   eigenpairs to be preserved. 1 designates the eigenvalue with
 *    largest algebraic value, 2 the one with second largest algebraic
 *    value, etc.
 * \param nt Positive integer. When \p algo is \c IGRAPH_SCG_OPTIMUM,
 *    it gives the number of groups to partition each eigenvector
 *    separately. When \p algo is \c IGRAPH_SCG_INTERV or \c
 *    IGRAPH_SCG_INTERV_KM, it gives the number of intervals to
 *    partition each eigenvector. This is ignored when \p algo is \c
 *    IGRAPH_SCG_EXACT.
 * \param nt_vec A numeric vector of length one or the length must
 *    match the number of eigenvectors given in \p V, or a \c NULL
 *    pointer. If not \c NULL, then this argument gives the number of
 *    groups or intervals, and \p nt is ignored. Different number of
 *    groups or intervals can be specified for each eigenvector.
 * \param algo The algorithm to solve the SCG problem. Possible
 *    values: \c IGRAPH_SCG_OPTIMUM, \c IGRAPH_SCG_INTERV_KM, \c
 *    IGRAPH_SCG_INTERV and \c IGRAPH_SCG_EXACT. Please see the
 *    details about them above.
 * \param values If this is not \c NULL and the eigenvectors are
 *    re-calculated, then the eigenvalues are stored here.
 * \param vectors If this is not \c NULL, and not a zero-length
 *    matrix, then it is interpreted as the eigenvectors to use for
 *    the coarse-graining. Otherwise the eigenvectors are
 *    re-calculated, and they are stored here. (If this is not \c NULL.)
 * \param groups If this is not \c NULL, and not a zero-length vector,
 *    then it is interpreted as the vector of group labels. (Group
 *    labels are integers from zero and are sequential.) Otherwise
 *    group labels are re-calculated and stored here, if this argument
 *    is not a null pointer.
 * \param use_arpack Whether to use ARPACK for solving the
 *    eigenproblem. Currently ARPACK is not implemented.
 * \param maxiter A positive integer giving the number of iterations
 *    of the k-means algorithm when \p algo is \c
 *    IGRAPH_SCG_INTERV_KM. It is ignored in other cases. A reasonable
 *    (initial) value for this argument is 100.
 * \param scg_graph If not a \c NULL pointer, then the coarse-grained
 *    graph is returned here.
 * \param scg_matrix If not a \c NULL pointer, then it must be an
 *    initialied matrix, and the coarse-grained matrix is returned
 *    here.
 * \param scg_sparsemat If not a \c NULL pointer, then the coarse
 *    grained matrix is returned here, in sparse matrix form.
 * \param L If not a \c NULL pointer, then it must be an initialized
 *    matrix and the left semi-projector is returned here.
 * \param R If not a \c NULL pointer, then it must be an initialized
 *    matrix and the right semi-projector is returned here.
 * \param Lsparse If not a \c NULL pointer, then the left
 *    semi-projector is returned here.
 * \param Rsparse If not a \c NULL pointer, then the right
 *    semi-projector is returned here.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \sa \ref igraph_scg_grouping(), \ref igraph_scg_semiprojectors(),
 * \ref igraph_scg_stochastic() and \ref igraph_scg_laplacian().
 *
 * \example examples/simple/scg.c
 */

int igraph_scg_adjacency(const igraph_t *graph,
                         const igraph_matrix_t *matrix,
                         const igraph_sparsemat_t *sparsemat,
                         const igraph_vector_t *ev,
                         igraph_integer_t nt,
                         const igraph_vector_t *nt_vec,
                         igraph_scg_algorithm_t algo,
                         igraph_vector_t *values,
                         igraph_matrix_t *vectors,
                         igraph_vector_t *groups,
                         igraph_bool_t use_arpack,
                         igraph_integer_t maxiter,
                         igraph_t *scg_graph,
                         igraph_matrix_t *scg_matrix,
                         igraph_sparsemat_t *scg_sparsemat,
                         igraph_matrix_t *L,
                         igraph_matrix_t *R,
                         igraph_sparsemat_t *Lsparse,
                         igraph_sparsemat_t *Rsparse) {

    igraph_sparsemat_t *mysparsemat = (igraph_sparsemat_t*) sparsemat,
                        real_sparsemat;
    int no_of_ev = (int) igraph_vector_size(ev);
    /* eigenvectors are calculated and returned */
    igraph_bool_t do_vectors = vectors && igraph_matrix_size(vectors) == 0;
    /* groups are calculated */
    igraph_bool_t do_groups = !groups || igraph_vector_size(groups) == 0;
    /* eigenvectors are not returned but must be calculated for groups */
    igraph_bool_t tmp_vectors = !do_vectors && do_groups;
    /* need temporary vector for groups */
    igraph_bool_t tmp_groups = !groups;
    igraph_matrix_t myvectors;
    igraph_vector_t mygroups;
    igraph_bool_t tmp_lsparse = !Lsparse, tmp_rsparse = !Rsparse;
    igraph_sparsemat_t myLsparse, myRsparse, tmpsparse, Rsparse_t;
    int no_of_nodes;
    igraph_real_t evmin, evmax;
    igraph_bool_t directed;

    /* --------------------------------------------------------------------*/
    /* Argument checks */

    IGRAPH_CHECK(igraph_i_scg_common_checks(graph, matrix, sparsemat,
                                            ev, nt, nt_vec,
                                            vectors, 0, groups, scg_graph,
                                            scg_matrix, scg_sparsemat,
                                            /*p=*/ 0, &evmin, &evmax));

    if (graph) {
        no_of_nodes = igraph_vcount(graph);
        directed = igraph_is_directed(graph);
    } else if (matrix) {
        no_of_nodes = (int) igraph_matrix_nrow(matrix);
        directed = !igraph_matrix_is_symmetric(matrix);
    } else {
        no_of_nodes = (int) igraph_sparsemat_nrow(sparsemat);
        directed = !igraph_sparsemat_is_symmetric(sparsemat);
    }

    /* -------------------------------------------------------------------- */
    /* Convert graph, if needed */

    if (graph) {
        mysparsemat = &real_sparsemat;
        IGRAPH_CHECK(igraph_get_sparsemat(graph, mysparsemat));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
    }

    /* -------------------------------------------------------------------- */
    /* Compute eigenpairs, if needed */
    if (tmp_vectors) {
        vectors = &myvectors;
        IGRAPH_MATRIX_INIT_FINALLY(vectors, no_of_nodes, no_of_ev);
    }

    if (do_vectors || tmp_vectors) {
        igraph_arpack_options_t options;
        igraph_eigen_which_t which;
        igraph_matrix_t tmp;
        igraph_vector_t tmpev;
        igraph_vector_t tmpeval;
        int i;

        which.pos = IGRAPH_EIGEN_SELECT;
        which.il = (int) (no_of_nodes - evmax + 1);
        which.iu = (int) (no_of_nodes - evmin + 1);

        if (values) {
            IGRAPH_VECTOR_INIT_FINALLY(&tmpeval, 0);
        }
        IGRAPH_CHECK(igraph_matrix_init(&tmp, no_of_nodes,
                                        which.iu - which.il + 1));
        IGRAPH_FINALLY(igraph_matrix_destroy, &tmp);
        IGRAPH_CHECK(igraph_eigen_matrix_symmetric(matrix, mysparsemat,
                     /* fun= */ 0, no_of_nodes,
                     /* extra= */ 0,
                     /* algorithm= */
                     use_arpack ?
                     IGRAPH_EIGEN_ARPACK :
                     IGRAPH_EIGEN_LAPACK, &which,
                     &options, /*storage=*/ 0,
                     values ? &tmpeval : 0,
                     &tmp));
        IGRAPH_VECTOR_INIT_FINALLY(&tmpev, no_of_ev);
        for (i = 0; i < no_of_ev; i++) {
            VECTOR(tmpev)[i] = evmax - VECTOR(*ev)[i];
        }
        if (values) {
            IGRAPH_CHECK(igraph_vector_index(&tmpeval, values, &tmpev));
        }
        IGRAPH_CHECK(igraph_matrix_select_cols(&tmp, vectors, &tmpev));
        igraph_vector_destroy(&tmpev);
        igraph_matrix_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(2);
        if (values) {
            igraph_vector_destroy(&tmpeval);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    /* -------------------------------------------------------------------- */
    /* Work out groups, if needed */
    if (tmp_groups) {
        groups = &mygroups;
        IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)groups, no_of_nodes);
    }
    if (do_groups) {
        IGRAPH_CHECK(igraph_scg_grouping(vectors, (igraph_vector_t*)groups,
                                         nt, nt_vec,
                                         IGRAPH_SCG_SYMMETRIC, algo,
                                         /*p=*/ 0, maxiter));
    }

    /* -------------------------------------------------------------------- */
    /* Perform coarse graining */
    if (tmp_lsparse) {
        Lsparse = &myLsparse;
    }
    if (tmp_rsparse) {
        Rsparse = &myRsparse;
    }
    IGRAPH_CHECK(igraph_scg_semiprojectors(groups, IGRAPH_SCG_SYMMETRIC,
                                           L, R, Lsparse, Rsparse, /*p=*/ 0,
                                           IGRAPH_SCG_NORM_ROW));
    if (tmp_groups) {
        igraph_vector_destroy((igraph_vector_t*) groups);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (tmp_vectors) {
        igraph_matrix_destroy(vectors);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (Rsparse) {
        IGRAPH_FINALLY(igraph_sparsemat_destroy, Rsparse);
    }
    if (Lsparse) {
        IGRAPH_FINALLY(igraph_sparsemat_destroy, Lsparse);
    }

    /* -------------------------------------------------------------------- */
    /* Compute coarse grained matrix/graph/sparse matrix */
    IGRAPH_CHECK(igraph_sparsemat_compress(Rsparse, &tmpsparse));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmpsparse);
    IGRAPH_CHECK(igraph_sparsemat_transpose(&tmpsparse, &Rsparse_t,
                                            /*values=*/ 1));
    igraph_sparsemat_destroy(&tmpsparse);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse_t);

    IGRAPH_CHECK(igraph_i_scg_get_result(IGRAPH_SCG_SYMMETRIC,
                                         matrix, mysparsemat,
                                         Lsparse, &Rsparse_t,
                                         scg_graph, scg_matrix,
                                         scg_sparsemat, directed));

    /* -------------------------------------------------------------------- */
    /* Clean up */

    igraph_sparsemat_destroy(&Rsparse_t);
    IGRAPH_FINALLY_CLEAN(1);
    if (Lsparse) {
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (Rsparse) {
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (graph) {
        igraph_sparsemat_destroy(mysparsemat);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_scg_stochastic
 * Spectral coarse graining, stochastic case.
 *
 * This function handles all the steps involved in the Spectral Coarse
 * Graining (SCG) of some matrices and graphs as described in the
 * reference below.
 *
 * \param graph The input graph. Exactly one of \p graph, \p matrix
 *    and \p sparsemat must be given, the other two must be \c NULL
 *    pointers.
 * \param matrix The input matrix. Exactly one of \p graph, \p matrix
 *    and \p sparsemat must be given, the other two must be \c NULL
 *    pointers.
 * \param sparsemat The input sparse matrix. Exactly one of \p graph,
 *    \p matrix and \p sparsemat must be given, the other two must be
 *    \c NULL pointers.
 * \param ev A vector of positive integers giving the indexes of the
 *   eigenpairs to be preserved. 1 designates the eigenvalue with
 *    largest magnitude, 2 the one with second largest magnitude, etc.
 * \param nt Positive integer. When \p algo is \c IGRAPH_SCG_OPTIMUM,
 *    it gives the number of groups to partition each eigenvector
 *    separately. When \p algo is \c IGRAPH_SCG_INTERV or \c
 *    IGRAPH_SCG_INTERV_KM, it gives the number of intervals to
 *    partition each eigenvector. This is ignored when \p algo is \c
 *    IGRAPH_SCG_EXACT.
 * \param nt_vec A numeric vector of length one or the length must
 *    match the number of eigenvectors given in \p V, or a \c NULL
 *    pointer. If not \c NULL, then this argument gives the number of
 *    groups or intervals, and \p nt is ignored. Different number of
 *    groups or intervals can be specified for each eigenvector.
 * \param algo The algorithm to solve the SCG problem. Possible
 *    values: \c IGRAPH_SCG_OPTIMUM, \c IGRAPH_SCG_INTERV_KM, \c
 *    IGRAPH_SCG_INTERV and \c IGRAPH_SCG_EXACT. Please see the
 *    details about them above.
 * \param norm Either \c IGRAPH_SCG_NORM_ROW or \c IGRAPH_SCG_NORM_COL.
 *    Specifies whether the rows or the columns of the
 *    stochastic matrix sum up to one.
 * \param values If this is not \c NULL and the eigenvectors are
 *    re-calculated, then the eigenvalues are stored here.
 * \param vectors If this is not \c NULL, and not a zero-length
 *    matrix, then it is interpreted as the eigenvectors to use for
 *    the coarse-graining. Otherwise the eigenvectors are
 *    re-calculated, and they are stored here. (If this is not \c NULL.)
 * \param groups If this is not \c NULL, and not a zero-length vector,
 *    then it is interpreted as the vector of group labels. (Group
 *    labels are integers from zero and are sequential.) Otherwise
 *    group labels are re-calculated and stored here, if this argument
 *    is not a null pointer.
 * \param p If this is not \c NULL, and not zero length, then it is
 *    interpreted as the stationary probability distribution of the
 *    Markov chain corresponding to the input matrix/graph. Its length
 *    must match the number of  vertices in the input graph (or number
 *    of rows in the input matrix). If not given, then the stationary
 *    distribution is calculated and stored here. (Unless this
 *    argument is a \c NULL pointer, in which case it is not stored.)
 * \param use_arpack Whether to use ARPACK for solving the
 *    eigenproblem. Currently ARPACK is not implemented.
 * \param maxiter A positive integer giving the number of iterations
 *    of the k-means algorithm when \p algo is \c
 *    IGRAPH_SCG_INTERV_KM. It is ignored in other cases. A reasonable
 *    (initial) value for this argument is 100.
 * \param scg_graph If not a \c NULL pointer, then the coarse-grained
 *    graph is returned here.
 * \param scg_matrix If not a \c NULL pointer, then it must be an
 *    initialied matrix, and the coarse-grained matrix is returned
 *    here.
 * \param scg_sparsemat If not a \c NULL pointer, then the coarse
 *    grained matrix is returned here, in sparse matrix form.
 * \param L If not a \c NULL pointer, then it must be an initialized
 *    matrix and the left semi-projector is returned here.
 * \param R If not a \c NULL pointer, then it must be an initialized
 *    matrix and the right semi-projector is returned here.
 * \param Lsparse If not a \c NULL pointer, then the left
 *    semi-projector is returned here.
 * \param Rsparse If not a \c NULL pointer, then the right
 *    semi-projector is returned here.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \sa \ref igraph_scg_grouping(), \ref igraph_scg_semiprojectors(),
 * \ref igraph_scg_adjacency() and \ref igraph_scg_laplacian().
 */

int igraph_scg_stochastic(const igraph_t *graph,
                          const igraph_matrix_t *matrix,
                          const igraph_sparsemat_t *sparsemat,
                          const igraph_vector_t *ev,
                          igraph_integer_t nt,
                          const igraph_vector_t *nt_vec,
                          igraph_scg_algorithm_t algo,
                          igraph_scg_norm_t norm,
                          igraph_vector_complex_t *values,
                          igraph_matrix_complex_t *vectors,
                          igraph_vector_t *groups,
                          igraph_vector_t *p,
                          igraph_bool_t use_arpack,
                          igraph_integer_t maxiter,
                          igraph_t *scg_graph,
                          igraph_matrix_t *scg_matrix,
                          igraph_sparsemat_t *scg_sparsemat,
                          igraph_matrix_t *L,
                          igraph_matrix_t *R,
                          igraph_sparsemat_t *Lsparse,
                          igraph_sparsemat_t *Rsparse) {

    igraph_matrix_t *mymatrix = (igraph_matrix_t*) matrix, real_matrix;
    igraph_sparsemat_t *mysparsemat = (igraph_sparsemat_t*) sparsemat,
                        real_sparsemat;
    int no_of_nodes;
    igraph_real_t evmin, evmax;
    igraph_arpack_options_t options;
    igraph_eigen_which_t which;
    /* eigenvectors are calculated and returned */
    igraph_bool_t do_vectors = vectors && igraph_matrix_complex_size(vectors) == 0;
    /* groups are calculated */
    igraph_bool_t do_groups = !groups || igraph_vector_size(groups) == 0;
    igraph_bool_t tmp_groups = !groups;
    /* eigenvectors are not returned but must be calculated for groups */
    igraph_bool_t tmp_vectors = !do_vectors && do_groups;
    igraph_matrix_complex_t myvectors;
    igraph_vector_t mygroups;
    igraph_bool_t do_p = !p || igraph_vector_size(p) == 0;
    igraph_vector_t *myp = (igraph_vector_t *) p, real_p;
    int no_of_ev = (int) igraph_vector_size(ev);
    igraph_bool_t tmp_lsparse = !Lsparse, tmp_rsparse = !Rsparse;
    igraph_sparsemat_t myLsparse, myRsparse, tmpsparse, Rsparse_t;

    /* --------------------------------------------------------------------*/
    /* Argument checks */

    IGRAPH_CHECK(igraph_i_scg_common_checks(graph, matrix, sparsemat,
                                            ev, nt, nt_vec,
                                            0, vectors, groups, scg_graph,
                                            scg_matrix, scg_sparsemat, p,
                                            &evmin, &evmax));

    if (graph) {
        no_of_nodes = igraph_vcount(graph);
    } else if (matrix) {
        no_of_nodes = (int) igraph_matrix_nrow(matrix);
    } else {
        no_of_nodes = (int) igraph_sparsemat_nrow(sparsemat);
    }

    /* -------------------------------------------------------------------- */
    /* Convert graph, if needed */

    if (graph) {
        mysparsemat = &real_sparsemat;
        IGRAPH_CHECK(igraph_get_stochastic_sparsemat(graph, mysparsemat,
                     norm == IGRAPH_SCG_NORM_COL));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
    } else if (matrix) {
        mymatrix = &real_matrix;
        IGRAPH_CHECK(igraph_i_matrix_stochastic(matrix, mymatrix, norm));
        IGRAPH_FINALLY(igraph_matrix_destroy, mymatrix);
    } else { /* sparsemat */
        mysparsemat = &real_sparsemat;
        IGRAPH_CHECK(igraph_i_sparsemat_stochastic(sparsemat, mysparsemat, norm));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
    }

    /* -------------------------------------------------------------------- */
    /* Compute eigenpairs, if needed */

    if (tmp_vectors) {
        vectors = &myvectors;
        IGRAPH_CHECK(igraph_matrix_complex_init(vectors, no_of_nodes, no_of_ev));
        IGRAPH_FINALLY(igraph_matrix_complex_destroy, vectors);
    }

    if (do_vectors || tmp_vectors) {
        igraph_matrix_complex_t tmp;
        igraph_vector_t tmpev;
        igraph_vector_complex_t tmpeval;
        int i;

        which.pos = IGRAPH_EIGEN_SELECT;
        which.il = (int) (no_of_nodes - evmax + 1);
        which.iu = (int) (no_of_nodes - evmin + 1);

        if (values) {
            IGRAPH_CHECK(igraph_vector_complex_init(&tmpeval, 0));
            IGRAPH_FINALLY(igraph_vector_complex_destroy, &tmpeval);
        }
        IGRAPH_CHECK(igraph_matrix_complex_init(&tmp, no_of_nodes,
                                                which.iu - which.il + 1));
        IGRAPH_FINALLY(igraph_matrix_complex_destroy, &tmp);
        IGRAPH_CHECK(igraph_eigen_matrix(mymatrix, mysparsemat, /*fun=*/ 0,
                                         no_of_nodes, /*extra=*/ 0, use_arpack ?
                                         IGRAPH_EIGEN_ARPACK :
                                         IGRAPH_EIGEN_LAPACK, &which, &options,
                                         /*storage=*/ 0,
                                         values ? &tmpeval : 0, &tmp));

        IGRAPH_VECTOR_INIT_FINALLY(&tmpev, no_of_ev);
        for (i = 0; i < no_of_ev; i++) {
            VECTOR(tmpev)[i] = evmax - VECTOR(*ev)[i];
        }
        if (values) {
            IGRAPH_CHECK(igraph_vector_complex_index(&tmpeval, values, &tmpev));
        }
        IGRAPH_CHECK(igraph_matrix_complex_select_cols(&tmp, vectors, &tmpev));
        igraph_vector_destroy(&tmpev);
        igraph_matrix_complex_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(2);
        if (values) {
            igraph_vector_complex_destroy(&tmpeval);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    /* Compute p if not supplied */
    if (do_p) {
        igraph_eigen_which_t w;
        igraph_matrix_complex_t tmp;
        igraph_arpack_options_t o;
        igraph_matrix_t trans, *mytrans = &trans;
        igraph_sparsemat_t sparse_trans, *mysparse_trans = &sparse_trans;
        int i;
        igraph_arpack_options_init(&o);
        if (!p) {
            IGRAPH_VECTOR_INIT_FINALLY(&real_p, no_of_nodes);
            myp = &real_p;
        } else {
            IGRAPH_CHECK(igraph_vector_resize(p, no_of_nodes));
        }
        IGRAPH_CHECK(igraph_matrix_complex_init(&tmp, 0, 0));
        IGRAPH_FINALLY(igraph_matrix_complex_destroy, &tmp);
        w.pos = IGRAPH_EIGEN_LR;
        w.howmany = 1;

        if (mymatrix) {
            IGRAPH_CHECK(igraph_matrix_copy(&trans, mymatrix));
            IGRAPH_FINALLY(igraph_matrix_destroy, &trans);
            IGRAPH_CHECK(igraph_matrix_transpose(&trans));
            mysparse_trans = 0;
        } else {
            IGRAPH_CHECK(igraph_sparsemat_transpose(mysparsemat, &sparse_trans,
                                                    /*values=*/ 1));
            IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparse_trans);
            mytrans = 0;
        }

        IGRAPH_CHECK(igraph_eigen_matrix(mytrans, mysparse_trans, /*fun=*/ 0,
                                         no_of_nodes, /*extra=*/ 0, /*algorith=*/
                                         use_arpack ?
                                         IGRAPH_EIGEN_ARPACK :
                                         IGRAPH_EIGEN_LAPACK, &w, &o,
                                         /*storage=*/ 0, /*values=*/ 0, &tmp));

        if (mymatrix) {
            igraph_matrix_destroy(&trans);
            IGRAPH_FINALLY_CLEAN(1);
        } else {
            igraph_sparsemat_destroy(mysparse_trans);
            IGRAPH_FINALLY_CLEAN(1);
        }

        for (i = 0; i < no_of_nodes; i++) {
            VECTOR(*myp)[i] = fabs(IGRAPH_REAL(MATRIX(tmp, i, 0)));
        }
        igraph_matrix_complex_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* -------------------------------------------------------------------- */
    /* Work out groups, if needed */
    /* TODO: use complex part as well */
    if (tmp_groups) {
        groups = &mygroups;
        IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)groups, no_of_nodes);
    }
    if (do_groups) {
        igraph_matrix_t tmp;
        IGRAPH_MATRIX_INIT_FINALLY(&tmp, 0, 0);
        IGRAPH_CHECK(igraph_matrix_complex_real(vectors, &tmp));
        IGRAPH_CHECK(igraph_scg_grouping(&tmp, (igraph_vector_t*)groups,
                                         nt, nt_vec,
                                         IGRAPH_SCG_STOCHASTIC, algo,
                                         myp, maxiter));
        igraph_matrix_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* -------------------------------------------------------------------- */
    /* Perform coarse graining */
    if (tmp_lsparse) {
        Lsparse = &myLsparse;
    }
    if (tmp_rsparse) {
        Rsparse = &myRsparse;
    }
    IGRAPH_CHECK(igraph_scg_semiprojectors(groups, IGRAPH_SCG_STOCHASTIC,
                                           L, R, Lsparse, Rsparse, myp, norm));
    if (tmp_groups) {
        igraph_vector_destroy((igraph_vector_t*) groups);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (!p && do_p) {
        igraph_vector_destroy(myp);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (tmp_vectors) {
        igraph_matrix_complex_destroy(vectors);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (Rsparse) {
        IGRAPH_FINALLY(igraph_sparsemat_destroy, Rsparse);
    }
    if (Lsparse) {
        IGRAPH_FINALLY(igraph_sparsemat_destroy, Lsparse);
    }

    /* -------------------------------------------------------------------- */
    /* Compute coarse grained matrix/graph/sparse matrix */
    IGRAPH_CHECK(igraph_sparsemat_compress(Rsparse, &tmpsparse));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmpsparse);
    IGRAPH_CHECK(igraph_sparsemat_transpose(&tmpsparse, &Rsparse_t,
                                            /*values=*/ 1));
    igraph_sparsemat_destroy(&tmpsparse);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse_t);

    IGRAPH_CHECK(igraph_i_scg_get_result(IGRAPH_SCG_STOCHASTIC,
                                         mymatrix, mysparsemat,
                                         Lsparse, &Rsparse_t,
                                         scg_graph, scg_matrix,
                                         scg_sparsemat, /*directed=*/ 1));

    /* -------------------------------------------------------------------- */
    /* Clean up */

    igraph_sparsemat_destroy(&Rsparse_t);
    IGRAPH_FINALLY_CLEAN(1);
    if (Lsparse) {
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (Rsparse) {
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (graph) {
        igraph_sparsemat_destroy(mysparsemat);
        IGRAPH_FINALLY_CLEAN(1);
    } else if (matrix) {
        igraph_matrix_destroy(mymatrix);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_sparsemat_destroy(mysparsemat);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}

/**
 * \function igraph_scg_laplacian
 * \brief Spectral coarse graining, Laplacian case.
 *
 * This function handles all the steps involved in the Spectral Coarse
 * Graining (SCG) of some matrices and graphs as described in the
 * reference below.
 *
 * \param graph The input graph. Exactly one of \p graph, \p matrix
 *    and \p sparsemat must be given, the other two must be \c NULL
 *    pointers.
 * \param matrix The input matrix. Exactly one of \p graph, \p matrix
 *    and \p sparsemat must be given, the other two must be \c NULL
 *    pointers.
 * \param sparsemat The input sparse matrix. Exactly one of \p graph,
 *    \p matrix and \p sparsemat must be given, the other two must be
 *    \c NULL pointers.
 * \param ev A vector of positive integers giving the indexes of the
 *   eigenpairs to be preserved. 1 designates the eigenvalue with
 *    largest magnitude, 2 the one with second largest magnitude, etc.
 * \param nt Positive integer. When \p algo is \c IGRAPH_SCG_OPTIMUM,
 *    it gives the number of groups to partition each eigenvector
 *    separately. When \p algo is \c IGRAPH_SCG_INTERV or \c
 *    IGRAPH_SCG_INTERV_KM, it gives the number of intervals to
 *    partition each eigenvector. This is ignored when \p algo is \c
 *    IGRAPH_SCG_EXACT.
 * \param nt_vec A numeric vector of length one or the length must
 *    match the number of eigenvectors given in \p V, or a \c NULL
 *    pointer. If not \c NULL, then this argument gives the number of
 *    groups or intervals, and \p nt is ignored. Different number of
 *    groups or intervals can be specified for each eigenvector.
 * \param algo The algorithm to solve the SCG problem. Possible
 *    values: \c IGRAPH_SCG_OPTIMUM, \c IGRAPH_SCG_INTERV_KM, \c
 *    IGRAPH_SCG_INTERV and \c IGRAPH_SCG_EXACT. Please see the
 *    details about them above.
 * \param norm Either \c IGRAPH_SCG_NORM_ROW or \c IGRAPH_SCG_NORM_COL.
 *    Specifies whether the rows or the columns of the Laplacian
 *    matrix sum up to zero.
 * \param direction Whether to work with left or right eigenvectors.
 *    Possible values: \c IGRAPH_SCG_DIRECTION_DEFAULT, \c
 *    IGRAPH_SCG_DIRECTION_LEFT, \c IGRAPH_SCG_DIRECTION_RIGHT. This
 *    argument is currently ignored and right eigenvectors are always
 *    used.
 * \param values If this is not \c NULL and the eigenvectors are
 *    re-calculated, then the eigenvalues are stored here.
 * \param vectors If this is not \c NULL, and not a zero-length
 *    matrix, then it is interpreted as the eigenvectors to use for
 *    the coarse-graining. Otherwise the eigenvectors are
 *    re-calculated, and they are stored here. (If this is not \c NULL.)
 * \param groups If this is not \c NULL, and not a zero-length vector,
 *    then it is interpreted as the vector of group labels. (Group
 *    labels are integers from zero and are sequential.) Otherwise
 *    group labels are re-calculated and stored here, if this argument
 *    is not a null pointer.
 * \param use_arpack Whether to use ARPACK for solving the
 *    eigenproblem. Currently ARPACK is not implemented.
 * \param maxiter A positive integer giving the number of iterations
 *    of the k-means algorithm when \p algo is \c
 *    IGRAPH_SCG_INTERV_KM. It is ignored in other cases. A reasonable
 *    (initial) value for this argument is 100.
 * \param scg_graph If not a \c NULL pointer, then the coarse-grained
 *    graph is returned here.
 * \param scg_matrix If not a \c NULL pointer, then it must be an
 *    initialied matrix, and the coarse-grained matrix is returned
 *    here.
 * \param scg_sparsemat If not a \c NULL pointer, then the coarse
 *    grained matrix is returned here, in sparse matrix form.
 * \param L If not a \c NULL pointer, then it must be an initialized
 *    matrix and the left semi-projector is returned here.
 * \param R If not a \c NULL pointer, then it must be an initialized
 *    matrix and the right semi-projector is returned here.
 * \param Lsparse If not a \c NULL pointer, then the left
 *    semi-projector is returned here.
 * \param Rsparse If not a \c NULL pointer, then the right
 *    semi-projector is returned here.
 * \return Error code.
 *
 * Time complexity: TODO.
 *
 * \sa \ref igraph_scg_grouping(), \ref igraph_scg_semiprojectors(),
 * \ref igraph_scg_stochastic() and \ref igraph_scg_adjacency().
 */

int igraph_scg_laplacian(const igraph_t *graph,
                         const igraph_matrix_t *matrix,
                         const igraph_sparsemat_t *sparsemat,
                         const igraph_vector_t *ev,
                         igraph_integer_t nt,
                         const igraph_vector_t *nt_vec,
                         igraph_scg_algorithm_t algo,
                         igraph_scg_norm_t norm,
                         igraph_scg_direction_t direction,
                         igraph_vector_complex_t *values,
                         igraph_matrix_complex_t *vectors,
                         igraph_vector_t *groups,
                         igraph_bool_t use_arpack,
                         igraph_integer_t maxiter,
                         igraph_t *scg_graph,
                         igraph_matrix_t *scg_matrix,
                         igraph_sparsemat_t *scg_sparsemat,
                         igraph_matrix_t *L,
                         igraph_matrix_t *R,
                         igraph_sparsemat_t *Lsparse,
                         igraph_sparsemat_t *Rsparse) {

    igraph_matrix_t *mymatrix = (igraph_matrix_t*) matrix, real_matrix;
    igraph_sparsemat_t *mysparsemat = (igraph_sparsemat_t*) sparsemat,
                        real_sparsemat;
    int no_of_nodes;
    igraph_real_t evmin, evmax;
    igraph_arpack_options_t options;
    igraph_eigen_which_t which;
    /* eigenvectors are calculated and returned */
    igraph_bool_t do_vectors = vectors && igraph_matrix_complex_size(vectors) == 0;
    /* groups are calculated */
    igraph_bool_t do_groups = !groups || igraph_vector_size(groups) == 0;
    igraph_bool_t tmp_groups = !groups;
    /* eigenvectors are not returned but must be calculated for groups */
    igraph_bool_t tmp_vectors = !do_vectors && do_groups;
    igraph_matrix_complex_t myvectors;
    igraph_vector_t mygroups;
    int no_of_ev = (int) igraph_vector_size(ev);
    igraph_bool_t tmp_lsparse = !Lsparse, tmp_rsparse = !Rsparse;
    igraph_sparsemat_t myLsparse, myRsparse, tmpsparse, Rsparse_t;

    IGRAPH_UNUSED(direction);

    /* --------------------------------------------------------------------*/
    /* Argument checks */

    IGRAPH_CHECK(igraph_i_scg_common_checks(graph, matrix, sparsemat,
                                            ev, nt, nt_vec,
                                            0, vectors, groups, scg_graph,
                                            scg_matrix, scg_sparsemat, /*p=*/ 0,
                                            &evmin, &evmax));

    if (graph) {
        no_of_nodes = igraph_vcount(graph);
    } else if (matrix) {
        no_of_nodes = (int) igraph_matrix_nrow(matrix);
    } else {
        no_of_nodes = (int) igraph_sparsemat_nrow(sparsemat);
    }

    /* -------------------------------------------------------------------- */
    /* Convert graph, if needed, get Laplacian matrix */

    if (graph) {
        mysparsemat = &real_sparsemat;
        IGRAPH_CHECK(igraph_sparsemat_init(mysparsemat, 0, 0, 0));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
        IGRAPH_CHECK(igraph_laplacian(graph, 0, mysparsemat, /*normalized=*/ 0,
                                      /*weights=*/ 0));
    } else if (matrix) {
        mymatrix = &real_matrix;
        IGRAPH_MATRIX_INIT_FINALLY(mymatrix, no_of_nodes, no_of_nodes);
        IGRAPH_CHECK(igraph_i_matrix_laplacian(matrix, mymatrix, norm));
    } else { /* sparsemat */
        mysparsemat = &real_sparsemat;
        IGRAPH_CHECK(igraph_i_sparsemat_laplacian(sparsemat, mysparsemat,
                     norm == IGRAPH_SCG_NORM_COL));
        IGRAPH_FINALLY(igraph_sparsemat_destroy, mysparsemat);
    }

    /* -------------------------------------------------------------------- */
    /* Compute eigenpairs, if needed */

    if (tmp_vectors) {
        vectors = &myvectors;
        IGRAPH_CHECK(igraph_matrix_complex_init(vectors, no_of_nodes, no_of_ev));
        IGRAPH_FINALLY(igraph_matrix_complex_destroy, vectors);
    }

    if (do_vectors || tmp_vectors) {
        igraph_matrix_complex_t tmp;
        igraph_vector_t tmpev;
        igraph_vector_complex_t tmpeval;
        int i;

        which.pos = IGRAPH_EIGEN_SELECT;
        which.il = (int) (no_of_nodes - evmax + 1);
        which.iu = (int) (no_of_nodes - evmin + 1);

        if (values) {
            IGRAPH_CHECK(igraph_vector_complex_init(&tmpeval, 0));
            IGRAPH_FINALLY(igraph_vector_complex_destroy, &tmpeval);
        }
        IGRAPH_CHECK(igraph_matrix_complex_init(&tmp, no_of_nodes,
                                                which.iu - which.il + 1));
        IGRAPH_FINALLY(igraph_matrix_complex_destroy, &tmp);
        IGRAPH_CHECK(igraph_eigen_matrix(mymatrix, mysparsemat, /*fun=*/ 0,
                                         no_of_nodes, /*extra=*/ 0, use_arpack ?
                                         IGRAPH_EIGEN_ARPACK :
                                         IGRAPH_EIGEN_LAPACK, &which, &options,
                                         /*storage=*/ 0,
                                         values ? &tmpeval : 0, &tmp));

        IGRAPH_VECTOR_INIT_FINALLY(&tmpev, no_of_ev);
        for (i = 0; i < no_of_ev; i++) {
            VECTOR(tmpev)[i] = evmax - VECTOR(*ev)[i];
        }
        if (values) {
            IGRAPH_CHECK(igraph_vector_complex_index(&tmpeval, values, &tmpev));
        }
        IGRAPH_CHECK(igraph_matrix_complex_select_cols(&tmp, vectors, &tmpev));
        igraph_vector_destroy(&tmpev);
        igraph_matrix_complex_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(2);
        if (values) {
            igraph_vector_complex_destroy(&tmpeval);
            IGRAPH_FINALLY_CLEAN(1);
        }
    }

    /* -------------------------------------------------------------------- */
    /* Work out groups, if needed */
    /* TODO: use complex part as well */
    if (tmp_groups) {
        groups = &mygroups;
        IGRAPH_VECTOR_INIT_FINALLY((igraph_vector_t*)groups, no_of_nodes);
    }
    if (do_groups) {
        igraph_matrix_t tmp;
        IGRAPH_MATRIX_INIT_FINALLY(&tmp, 0, 0);
        IGRAPH_CHECK(igraph_matrix_complex_real(vectors, &tmp));
        IGRAPH_CHECK(igraph_scg_grouping(&tmp, (igraph_vector_t*)groups,
                                         nt, nt_vec,
                                         IGRAPH_SCG_LAPLACIAN, algo,
                                         /*p=*/ 0, maxiter));
        igraph_matrix_destroy(&tmp);
        IGRAPH_FINALLY_CLEAN(1);
    }

    /* -------------------------------------------------------------------- */
    /* Perform coarse graining */
    if (tmp_lsparse) {
        Lsparse = &myLsparse;
    }
    if (tmp_rsparse) {
        Rsparse = &myRsparse;
    }
    IGRAPH_CHECK(igraph_scg_semiprojectors(groups, IGRAPH_SCG_LAPLACIAN,
                                           L, R, Lsparse, Rsparse, /*p=*/ 0,
                                           norm));
    if (tmp_groups) {
        igraph_vector_destroy((igraph_vector_t*) groups);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (tmp_vectors) {
        igraph_matrix_complex_destroy(vectors);
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (Rsparse) {
        IGRAPH_FINALLY(igraph_sparsemat_destroy, Rsparse);
    }
    if (Lsparse) {
        IGRAPH_FINALLY(igraph_sparsemat_destroy, Lsparse);
    }

    /* -------------------------------------------------------------------- */
    /* Compute coarse grained matrix/graph/sparse matrix */
    IGRAPH_CHECK(igraph_sparsemat_compress(Rsparse, &tmpsparse));
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &tmpsparse);
    IGRAPH_CHECK(igraph_sparsemat_transpose(&tmpsparse, &Rsparse_t,
                                            /*values=*/ 1));
    igraph_sparsemat_destroy(&tmpsparse);
    IGRAPH_FINALLY_CLEAN(1);
    IGRAPH_FINALLY(igraph_sparsemat_destroy, &Rsparse_t);

    IGRAPH_CHECK(igraph_i_scg_get_result(IGRAPH_SCG_LAPLACIAN,
                                         mymatrix, mysparsemat,
                                         Lsparse, &Rsparse_t,
                                         scg_graph, scg_matrix,
                                         scg_sparsemat, /*directed=*/ 1));

    /* -------------------------------------------------------------------- */
    /* Clean up */

    igraph_sparsemat_destroy(&Rsparse_t);
    IGRAPH_FINALLY_CLEAN(1);
    if (Lsparse) {
        IGRAPH_FINALLY_CLEAN(1);
    }
    if (Rsparse) {
        IGRAPH_FINALLY_CLEAN(1);
    }

    if (graph) {
        igraph_sparsemat_destroy(mysparsemat);
        IGRAPH_FINALLY_CLEAN(1);
    } else if (matrix) {
        igraph_matrix_destroy(mymatrix);
        IGRAPH_FINALLY_CLEAN(1);
    } else {
        igraph_sparsemat_destroy(mysparsemat);
        IGRAPH_FINALLY_CLEAN(1);
    }

    return 0;
}
