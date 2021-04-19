/* -*- mode: C -*-  */
/* vim:set ts=4 sw=4 sts=4 et: */
/*
   IGraph R package.
   Copyright (C) 2014  Gabor Csardi <csardi.gabor@gmail.com>
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
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <math.h>

#include "layout/layout_internal.h"

#include "test_utilities.inc"

int intersect() {

    float negative[][8] = {
        { 1, 2, 2, 2, 1, 1, 2, 1 }, /* 1 */
        { 1, 2, 1, 1, 2, 2, 2, 1 }, /* 2 */
        { 1, 0, 0, 1, 2, 0, 3, 1 }, /* 3 */
        { 1, 0, 1, 1, 0, 2, 2, 2 }, /* 4 */
        { 1, 0, 1, 2, 3, 1, 3, 3 }, /* 5 */
        { 0, 0, 0, 2, 1, 1, 1, 2 }, /* 6 */
        { 0, 1, 1, 1, 2, 0, 2, 3 }, /* 7 */
        { 0, 0, 5, 0, 2, 1, 4, 3 }, /* 8 */
        { 0, 0, 5, 5, 3, 2, 3, 2 }  /* 9 */
    };

    float positive[][8] = {
        { 0, 1, 2, 1, 1, 0, 1, 2 }, /* 10 */
        { 0, 2, 5, 2, 1, 1, 4, 3 }, /* 11 */
        { 0, 0, 0, 3, 0, 1, 5, 1 }, /* 12 */
        { 0, 4, 2, 6, 0, 4, 2, 2 }  /* 13 */
    };
    /* { 1,1,1,1, 1,1,0,0 }, /\* 14 *\/ */
    /* { 0,0,1,1, 1,1,1,1 }, /\* 15 *\/ */
    /* { 0,0,2,2, 1,1,1,1 }}; /\* 16 *\/ */

    int no_neg = sizeof(negative) / sizeof(float) / 8;
    int no_pos = sizeof(positive) / sizeof(float) / 8;
    int i;

    for (i = 0; i < no_neg; i++) {
        float *co = negative[i];
        if (igraph_i_layout_segments_intersect(co[0], co[1], co[2], co[3],
                                        co[4], co[5], co[6], co[7])) {
            return i + 1;
        }
    }

    for (i = 0; i < no_pos; i++) {
        float *co = positive[i];
        if (!igraph_i_layout_segments_intersect(co[0], co[1], co[2], co[3],
                                         co[4], co[5], co[6], co[7])) {
            return no_neg + i + 1;
        }
    }

    return 0;
}

int distance() {

    float configs[][7] = {
        { 1, 1, 2, 0, 2, 3, 1.0 }, /* 1 */
        { 1, 1, 1, 0, 1, 3, 0.0 }, /* 2 */
        { 1, 1, 0, 1, 1, 0, 0.5 }, /* 3 */
        { 1, 2, 0, 0, 0, 1, 2.0 }, /* 4 */
        { 1, 0, 0, 1, 0, 2, 2.0 }, /* 5 */
        { 0, 0, 1, 1, 1, 2, 2.0 }, /* 6 */
        { 0, 3, 1, 1, 1, 2, 2.0 }  /* 7 */
    };

    int no = sizeof(configs) / sizeof(float) / 8;
    int i;

    for (i = 0; i < no; i++) {
        float *co = configs[i];
        float res = igraph_i_layout_point_segment_dist2(co[0], co[1],
                    co[2], co[3], co[4], co[5]);
        if (fabsf(res - co[6]) > 1e-12) {
            printf("%g\n", (double) res);
            return i + 1;
        }
    }

    return 0;
}

int main() {
    int res1, res2;

    res1 = intersect();
    if (res1 != 0) {
        printf("Unexpected result from intersect(), config %d.\n", res1);
        return res1;
    }
    res2 = distance() ;
    if (res2 != 0) {
        printf("Unexpected result from distance(), config %d.\n", res2);
        return res2;
    }

    VERIFY_FINALLY_STACK();

    return 0;
}
