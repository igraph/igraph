/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>

int main() {

    igraph_t g;
    igraph_real_t flow;
    igraph_vector_t capacity;
    igraph_integer_t source, target;
    FILE *infile;
    igraph_maxflow_stats_t stats;

    igraph_vector_init(&capacity, 0);

    /***************/
    infile = fopen("ak-4102.max", "r");
    igraph_read_graph_dimacs(&g, infile, 0, 0, &source, &target, &capacity,
                             IGRAPH_DIRECTED);
    fclose(infile);

    igraph_maxflow_value(&g, &flow, source, target, &capacity, &stats);

    if (flow != 8207) {
        return 1;
    }
    igraph_destroy(&g);
    /***************/

    /*   /\***************\/ */
    /*   infile=fopen("ak-8198.max", "r"); */
    /*   igraph_read_graph_dimacs(&g, infile, 0, 0, &source, &target, &capacity, */
    /*             IGRAPH_DIRECTED); */
    /*   fclose(infile); */

    /*   t=timer(); */
    /*   igraph_maxflow_value(&g, &flow, source, target, &capacity, &stats); */
    /*   t=timer()-t; */
    /*   printf("8198: %g (time %.10f)\n", flow, t); */
    /*   igraph_destroy(&g); */
    /*   /\***************\/ */

    /*   /\***************\/ */
    /*   infile=fopen("ak-16390.max", "r"); */
    /*   igraph_read_graph_dimacs(&g, infile, 0, 0, &source, &target, &capacity, */
    /*             IGRAPH_DIRECTED); */
    /*   fclose(infile); */

    /*   t=timer(); */
    /*   igraph_maxflow_value(&g, &flow, source, target, &capacity, &stats); */
    /*   t=timer()-t; */
    /*   printf("16390: %g (time %.10f)\n", flow, t); */
    /*   igraph_destroy(&g); */
    /*   /\***************\/ */

    /*   /\***************\/ */
    /*   infile=fopen("ak-32774.max", "r"); */
    /*   igraph_read_graph_dimacs(&g, infile, 0, 0, &source, &target, &capacity, */
    /*             IGRAPH_DIRECTED); */
    /*   fclose(infile); */

    /*   t=timer(); */
    /*   igraph_maxflow_value(&g, &flow, source, target, &capacity, &stats); */
    /*   t=timer()-t; */
    /*   printf("32774: %g (time %.10f)\n", flow, t); */
    /*   igraph_destroy(&g); */
    /*   /\***************\/ */

    /*   /\***************\/ */
    /*   infile=fopen("ak-65542.max", "r"); */
    /*   igraph_read_graph_dimacs(&g, infile, 0, 0, &source, &target, &capacity, */
    /*             IGRAPH_DIRECTED); */
    /*   fclose(infile); */

    /*   t=timer(); */
    /*   igraph_maxflow_value(&g, &flow, source, target, &capacity, &stats); */
    /*   t=timer()-t; */
    /*   printf("65542: %g (time %.10f)\n", flow, t); */
    /*   igraph_destroy(&g); */
    /*   /\***************\/ */

    igraph_vector_destroy(&capacity);

    return 0;
}
