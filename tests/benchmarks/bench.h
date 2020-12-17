/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2013  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_BENCH_H
#define IGRAPH_BENCH_H

static inline void igraph_get_cpu_time(igraph_real_t *data) {

    struct rusage self, children;
    getrusage(RUSAGE_SELF, &self);
    getrusage(RUSAGE_CHILDREN, &children);
    data[0] = (double) self.ru_utime.tv_sec +
              1e-3 * (self.ru_utime.tv_usec / 1000);
    data[1] = (double) self.ru_stime.tv_sec +
              1e-3 * (self.ru_stime.tv_usec / 1000);
    data[2] = (double) children.ru_utime.tv_sec +
              1e-3 * (children.ru_utime.tv_usec / 1000);
    data[3] = (double) children.ru_stime.tv_sec +
              1e-3 * (children.ru_stime.tv_usec / 1000);
}

#define BENCH(NAME, ...)    do {                                                         \
        double start[4], stop[4];                                                                \
        igraph_get_cpu_time(start);                                                          \
        { __VA_ARGS__; };                                                                                \
        igraph_get_cpu_time(stop);                                                           \
        printf("%s %.3gs\n", NAME,                                                           \
               stop[0]+stop[1]+stop[2]+stop[3] -                               \
               start[0]-start[1]-start[2]-start[3]);                       \
    } while (0)

#endif

