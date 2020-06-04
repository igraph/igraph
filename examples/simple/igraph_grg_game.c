/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard st, Cambridge MA, 02139 USA

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
#include <math.h>

int main() {

    igraph_t g;
    /* long int i; */
    /* struct tms time; */
    /* clock_t current_time,start_time; */
    /* long int runs=100, n=10000; */
    /* igraph_real_t r=0.01; */

    /* Empty graph */
    igraph_grg_game(&g, 100, 0, 0, 0, 0);
    if (igraph_ecount(&g) != 0) {
        return 1;
    }
    igraph_destroy(&g);

    /* Full graph */
    igraph_grg_game(&g, 10, sqrt(2.0) / 2, 1, 0, 0);
    if (igraph_ecount(&g) != igraph_vcount(&g) * (igraph_vcount(&g) - 1) / 2) {
        return 2;
    }
    igraph_destroy(&g);

    /* Measure running time */
    /*   tps=sysconf(_SC_CLK_TCK); // clock ticks per second  */
    /*   times(&time); start_time=time.tms_utime; */
    /*   for (i=0; i<runs; i++) { */
    /*     igraph_grg_game2(&g, n, r, 1);  */
    /*     igraph_destroy(&g); */
    /*   } */
    /*   times(&time); current_time=time.tms_utime; */
    /*   fprintf(stdout,"    sorted: time=%.3fs\n",(current_time-start_time)/(double)tps); */

    /*   tps=sysconf(_SC_CLK_TCK); // clock ticks per second  */
    /*   times(&time); start_time=time.tms_utime; */
    /*   for (i=0; i<runs; i++) { */
    /*     igraph_grg_game(&g, n, r, 1); */
    /*     igraph_destroy(&g); */
    /*   } */
    /*   times(&time); current_time=time.tms_utime; */
    /*   fprintf(stdout,"non-sorted: time=%.3fs\n", */
    /*    (current_time-start_time)/(double)tps); */


    return 0;
}
