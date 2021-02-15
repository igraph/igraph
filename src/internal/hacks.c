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

/* These are implementations of common C functions that may be missing from some
 * compilers; for instance, icc does not provide stpcpy so we implement it
 * here. */

/**
 * Drop-in replacement for strdup.
 * Used only in compilers that do not have strdup or _strdup
 */
char* igraph_i_strdup(const char *s) {
    size_t n = strlen(s) + 1;
    char* result = (char*)malloc(sizeof(char) * n);
    if (result) {
        memcpy(result, s, n);
    }
    return result;
}

/**
 * Drop-in replacement for stpcpy.
 * Used only in compilers that do not have stpcpy
 */
char* igraph_i_stpcpy(char* s1, const char* s2) {
    char* result = strcpy(s1, s2);
    return result + strlen(s1);
}

