/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2010-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "internal/hacks.h"

#include <string.h>
#include <stdlib.h>

/* These are implementations of common C functions that may be missing from some compilers. */

/**
 * Drop-in replacement for strdup.
 * Used only in compilers that do not have strdup or _strdup
 */
char *igraph_i_strdup(const char *s) {
    size_t n = strlen(s) + 1;
    char *result = malloc(sizeof(char) * n);
    if (result) {
        memcpy(result, s, n);
    }
    return result;
}

/**
 * Drop-in replacement for strndup.
 * Used only in compilers that do not have strndup or _strndup
 */
char *igraph_i_strndup(const char *s1, size_t n) {
    size_t i;
    /* We need to check if the string is shorter than n characters.
     * We could use strlen, but that would do more work for long s1 and small n.
     * TODO: Maybe memchr would be nicer here.
     */
    for (i = 0; s1[i] != '\0' && i < n; i++) {}
    n = i;
    char *result = malloc(sizeof(char) * (n + 1));
    if (result) {
        memcpy(result, s1, n);
        result[n] = '\0';
    }
    return result;
}
