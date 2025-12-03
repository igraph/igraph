#include <igraph.h>

/**
 * This function implements the maximum cardinality search plus algorithm.
 * It computes minimal chordal completions of graphs.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Anne Berry, Jean R. S. Blair, Pinar Heggernes and Barry W. Peyton:
 * Maximum Cardinality Search for Computing Minimal Triangulations of Graphs.
 * Algorithmica 39, 287–298 (2004)
 * https://doi.org/10.1007/s00453-004-1084-3
 * 
 * </para><para>
 * Anne Berry, Romain Pogorelcnik and Geneviève Simonet:
 * An Introduction to Clique Minimal Separator Decomposition.
 * Algorithms 2010, 3(2), 197-215
 * https://doi.org/10.3390/a3020197
 *
 * \param g The input graph. Edge directions will be ignored.
 * \param alpha A minimal elimination ordering on the vertex set.
 * \param minimal_chordal A minimal chordal completionof graph.
 * \param min_sep_gen The set of vertices which generate a minimal separator of minimal_chordal.
 *
 */

void MCS_M_plus(const igraph_t *g,
                    igraph_vector_int_t *alpha,
                    igraph_vector_int_t *min_sep_gen,
                    igraph_t *minimal_chordal) {

    // the number of vertices in the graph
    const igraph_integer_t n = igraph_vcount(g);

    // init
    igraph_copy(minimal_chordal, g);
    igraph_vector_int_init(alpha, n);
    igraph_vector_int_init(min_sep_gen, 0);

    igraph_vector_int_t F; // to store the fill_in edges
    igraph_vector_int_init(&F, 0);    

    igraph_integer_t s = -1;

    // initialise the labels of all vertices as 0
    // label = -1 if the vertex has been chosen
    igraph_vector_int_t label;
    igraph_vector_int_init(&label, n);

    // for loop
    igraph_integer_t x, y, z, min_label, size_Y, size_Z, index;
    igraph_vector_int_t Y, Z, reached;
    igraph_vector_int_list_t reach;
    igraph_bool_t x_y_adj;
    for (int i = n-1; i >= 0; i--) {
        // choose a vertex x of g of maximal label
        igraph_vector_int_which_minmax(&label, &min_label, &x);

        // Y <- N_{g}(x)
        igraph_vector_int_init(&Y, 0);
        igraph_neighbors(g, &Y, x, IGRAPH_ALL);

        // remove from Y the vertices with label -1
        size_Y = igraph_vector_int_size(&Y);
        index = 0;
        for (int j = 0; j < size_Y; j++) {
            if (VECTOR(label)[VECTOR(Y)[j]] > -1) {
                VECTOR(Y)[index] = VECTOR(Y)[j];
                index ++;
            }
        }
        igraph_vector_int_resize(&Y, index);
        
        // if loop
        if (VECTOR(label)[x] <= s) {
            igraph_vector_int_push_back(min_sep_gen, x);
        }

        // s <- label(x)
        s = VECTOR(label)[x];

        // mark x reached and mark all other vertices of g unreached
        igraph_vector_int_init(&reached, n);
        VECTOR(reached)[x] = 1;

        // for loop, initialise reach as a list of vectors
        igraph_vector_int_list_init(&reach, n);
        for (int j = 0; j < n; j++) {
            igraph_vector_int_init(igraph_vector_int_list_get_ptr(&reach, j), 0);
        }

        // for loop, neighbours of x
        size_Y = igraph_vector_int_size(&Y);
        for (int j = 0; j < size_Y; j++) {
            y = VECTOR(Y)[j];

            // make y reached
            VECTOR(reached)[y] = 1;

            // add y to reach(label(y))
            igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&reach, VECTOR(label)[y]), y);      
        }

        // for loop, other vertices
        for (int j = 0; j < n; j++) {

            // while reach(j) \neq empty
            while (igraph_vector_int_size(igraph_vector_int_list_get_ptr(&reach, j)) > 0) {

                // remove a vertex y from reach(j)
                y = igraph_vector_int_pop_back(igraph_vector_int_list_get_ptr(&reach, j));

                // for loop, neighbours of y, Z is the neighbours of y in g
                igraph_vector_int_init(&Z, 0);
                igraph_neighbors(g, &Z, y, IGRAPH_ALL);
                size_Z = igraph_vector_int_size(&Z);
                for (int k = 0; k < size_Z; k++) {
                    z = VECTOR(Z)[k];

                    // if z is unreached
                    if (VECTOR(label)[z] > -1 && VECTOR(reached)[z] == 0) {

                        // mark z reached
                        VECTOR(reached)[z] = 1;

                        // if label(z) > j
                        if (VECTOR(label)[z] > j) {

                            // Y <- Y + {z}
                            igraph_vector_int_push_back(&Y, z);

                            // add z to reach(label(z))
                            igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&reach, VECTOR(label)[z]), z);

                        } else {
                            // add z to reach(j)
                            igraph_vector_int_push_back(igraph_vector_int_list_get_ptr(&reach, j), z);
                        }
                    }
                }
            }
        }
        // for loop
        size_Y = igraph_vector_int_size(&Y);
        for (int j = 0; j < size_Y; j++) {
            y = VECTOR(Y)[j];

            // F <- F + {xy} if {x,y} is not an edge in g, here one can also add everything then use igraph_simplify() to remove multiple edges
            igraph_are_adjacent(g, x, y, &x_y_adj);
            if (!x_y_adj) {
                igraph_vector_int_push_back(&F, x);
                igraph_vector_int_push_back(&F, y);
            }

            // label(y) <- label(y) + 1
            VECTOR(label)[y] ++;
        }

        // alpha(i) <- x
        VECTOR(*alpha)[i] = x;

        // remove x from g, i.e. set the label -1
        VECTOR(label)[x] = -1;      
    }

    // add the edges of F to minimal_chordal
    igraph_add_edges(minimal_chordal, &F, NULL);
}

/**
 * This function implements Atoms algorithm.
 * It computes the atoms and minimal separators of graphs.
 *
 * </para><para>
 * References:
 *
 * </para><para>
 * Anne Berry, Romain Pogorelcnik and Geneviève Simonet:
 * An Introduction to Clique Minimal Separator Decomposition.
 * Algorithms 2010, 3(2), 197-215
 * https://doi.org/10.3390/a3020197
 *
 * \param g The input graph. Edge directions will be ignored.
 * \param A The set of atoms of g.
 * \param S_C The set of cliques minimal seaparators of g. 
 * 
 */

void Atoms(const igraph_t *g,
            igraph_vector_int_list_t *A,
            igraph_vector_int_list_t *S_C) {

    // the number of vertices in the graph
    const igraph_integer_t n = igraph_vcount(g);

    igraph_t minimal_chordal;
    igraph_vector_int_t alpha, min_sep_gen;
    igraph_vector_int_list_t S_H;

    igraph_vector_int_init(&alpha,n);
    igraph_vector_int_init(&min_sep_gen,0);
    igraph_vector_int_list_init(A, 0);
    igraph_vector_int_list_init(&S_H, 0);
    igraph_vector_int_list_init(S_C, 0);

    MCS_M_plus(g, &alpha, &min_sep_gen, &minimal_chordal);

    // label = 1 if the vertex has been chosen, 0 otherwise (for minimal_chordal)
    igraph_vector_int_t label;
    igraph_vector_int_init(&label, n);

    // g1 is G'
    igraph_t g1;
    igraph_copy(&g1, g);

    // g1_v is the vertices of g1
    igraph_vector_int_t g1_v;
    igraph_vector_int_init_range(&g1_v, 0, n);
    igraph_vector_int_t map, invmap;
    igraph_vector_int_init(&map, 0);
    igraph_vector_int_init(&invmap, 0);

    // for loop
    igraph_integer_t x, size_S, size_C, index;
    igraph_vector_int_t S,g1_minus_S_v, C, g1_minus_C_v;
    igraph_bool_t S_is_clique;
    igraph_vs_t S_vids, g1_minus_S_vids, g1_minus_C_vids;
    igraph_t g1_minus_S, g2;
    for (int i = 0; i < n; i++) {

        // x <- alpha(i)
        x = VECTOR(alpha)[i];

        // if x in X
        if (igraph_vector_int_contains(&min_sep_gen, x)) {

            // S <- N_{H'}(x)
            //// get the neighbours of x in H
            igraph_vector_int_init(&S, 0);
            igraph_neighbors(&minimal_chordal, &S, x, IGRAPH_ALL);
            //// remove those neighbours not in H', i.e., with label != 0
            size_S = igraph_vector_int_size(&S);
            index = 0;
            for (int j = 0; j < size_S; j++) {
                if (VECTOR(label)[VECTOR(S)[j]] == 0) {
                    VECTOR(S)[index] = VECTOR(S)[j];
                    index ++;
                }
            }
            igraph_vector_int_resize(&S, index);
            //// make a copy of S
            igraph_vector_int_t S_copy;
            igraph_vector_int_init_copy(&S_copy, &S);
            
            // S_H append S
            if (igraph_vector_int_size(&S_copy) > 0) {
                igraph_vector_int_list_push_back(&S_H, &S_copy);
            }

            // if S is a clique in G
            igraph_vs_vector(&S_vids, &S_copy);
            igraph_is_clique(g, S_vids, 0, &S_is_clique);
            if (S_is_clique == 1) {

                // S_C append S
                if (igraph_vector_int_size(&S_copy) > 0) {
                    igraph_vector_int_list_push_back(S_C, &S_copy);
                }

                // C <- the connected component of G' - S containing x
                //// g1_minus_S is the graph G'-S
                igraph_vector_int_init(&g1_minus_S_v, 0);
                igraph_vector_int_difference_sorted(&g1_v, &S, &g1_minus_S_v);
                igraph_vs_vector(&g1_minus_S_vids, &g1_minus_S_v);
                igraph_induced_subgraph_map(g, &g1_minus_S, g1_minus_S_vids, IGRAPH_SUBGRAPH_AUTO, &map, &invmap);
                //// C is the component containing x
                igraph_vector_int_init(&C, 0);
                igraph_subcomponent(&g1_minus_S, &C, VECTOR(map)[x] - 1, IGRAPH_ALL);
                size_C = igraph_vector_int_size(&C);
                for (int k = 0; k < size_C; k++) {
                    VECTOR(C)[k] = VECTOR(invmap)[VECTOR(C)[k]];
                }
                igraph_vector_int_sort(&C);
                
                // A append S union C
                igraph_vector_int_append(&S, &C);
                igraph_vector_int_list_push_back(A, &S);

                // G' <- G' - C
                igraph_vector_int_init (&g1_minus_C_v, 0);
                igraph_vector_int_difference_sorted(&g1_v, &C, &g1_minus_C_v);
                igraph_vs_vector(&g1_minus_C_vids, &g1_minus_C_v);
                igraph_induced_subgraph_map(g, &g2, g1_minus_C_vids, IGRAPH_SUBGRAPH_AUTO, &map, &invmap);
                g1_v = g1_minus_C_v;
                g1 = g2;
            }
        }
        //remove x from H', i.e., set the label to 1
        VECTOR(label)[x] = 1;
    }
    // A append g1_v
    igraph_vector_int_list_push_back(A, &g1_v);
}