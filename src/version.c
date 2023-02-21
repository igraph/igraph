/*
   IGraph library.
   Copyright (C) 2008-2022  The igraph development team <igraph@igraph.org>

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

#include "igraph_version.h"

#include <stdio.h>

static const char *igraph_version_string = IGRAPH_VERSION;

/**
 * \function igraph_version
 * \brief The version of the igraph C library.
 *
 * \param version_string Pointer to a string pointer. If not null, it
 *    is set to the igraph version string, e.g. "0.9.11" or "0.10.0". This
 *    string must not be modified or deallocated.
 * \param major If not a null pointer, then it is set to the major
 *    igraph version. E.g. for version "0.9.11" this is 0.
 * \param minor If not a null pointer, then it is set to the minor
 *    igraph version. E.g. for version "0.9.11" this is 11.
 * \param subminor If not a null pointer, then it is set to the
 *    subminor igraph version. E.g. for version "0.9.11" this is 11.
 *
 * \example examples/simple/igraph_version.c
 */

void igraph_version(const char **version_string,
                    int *major,
                    int *minor,
                    int *subminor) {
    int i1, i2, i3;
    int *p1 = major ? major : &i1;
    int *p2 = minor ? minor : &i2;
    int *p3 = subminor ? subminor : &i3;

    if (version_string) {
        *version_string = igraph_version_string;
    }

    *p1 = *p2 = *p3 = 0;
    sscanf(IGRAPH_VERSION, "%i.%i.%i", p1, p2, p3);
}
