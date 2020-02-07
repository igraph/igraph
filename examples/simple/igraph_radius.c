/* -*- mode: C -*-  */
/* vim:set ts=4 sts=4 sw=4 et: */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

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

#include <igraph.h>

int main() {

    igraph_t g;
    igraph_real_t radius;

    igraph_star(&g, 10, IGRAPH_STAR_UNDIRECTED, 0);
    igraph_radius(&g, &radius, IGRAPH_OUT);
    if (radius != 1) {
        return 1;
    }
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_radius(&g, &radius, IGRAPH_ALL);
    if (radius != 1) {
        return 2;
    }
    igraph_destroy(&g);

    igraph_star(&g, 10, IGRAPH_STAR_OUT, 0);
    igraph_radius(&g, &radius, IGRAPH_OUT);
    if (radius != 0) {
        return 3;
    }
    igraph_destroy(&g);

    return 0;
}
