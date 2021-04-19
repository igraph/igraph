/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2008-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_version.h"

#include <stdio.h>

static const char *igraph_version_string = IGRAPH_VERSION;

/**
 * \function igraph_version
 * Return the version of the igraph C library
 *
 * \param version_string Pointer to a string pointer. If not null, it
 *    is set to the igraph version string, e.g. "0.6" or "0.5.3". This
 *    string should not be modified or deallocated.
 * \param major If not a null pointer, then it is set to the major
 *    igraph version. E.g. for version "0.5.3" this is 0.
 * \param minor If not a null pointer, then it is set to the minor
 *    igraph version. E.g. for version "0.5.3" this is 5.
 * \param subminor If not a null pointer, then it is set to the
 *    subminor igraph version. E.g. for version "0.5.3" this is 3.
 * \return Error code.
 *
 * Time complexity: O(1).
 *
 * \example examples/simple/igraph_version.c
 */

int igraph_version(const char **version_string,
                   int *major,
                   int *minor,
                   int *subminor) {
    int i1, i2, i3;
    int *p1 = major ? major : &i1,
         *p2 = minor ? minor : &i2,
          *p3 = subminor ? subminor : &i3;

    if (version_string) {
        *version_string = igraph_version_string;
    }

    *p1 = *p2 = *p3 = 0;
    sscanf(IGRAPH_VERSION, "%i.%i.%i", p1, p2, p3);

    return 0;
}
