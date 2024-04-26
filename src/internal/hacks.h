/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#ifndef IGRAPH_HACKS_INTERNAL_H
#define IGRAPH_HACKS_INTERNAL_H

#include "igraph_decls.h"

#include "config.h"

/* The CMake feature test looks for strcasecmp/strncasecmp in strings.h */
#if defined(HAVE_STRCASECMP) || defined(HAVE_STRNCASECMP)
#include <strings.h>
#endif

#include <stdlib.h>

__BEGIN_DECLS

#ifndef HAVE_STRDUP
    #define strdup igraph_i_strdup
    char* igraph_i_strdup(const char *s);
#endif

#ifndef HAVE_STRNDUP
    #define strndup igraph_i_strndup
    char* igraph_i_strndup(const char *s, size_t n);
#endif

#ifndef HAVE_STRCASECMP
    #ifdef HAVE__STRICMP
        #define strcasecmp _stricmp
    #else
        #error "igraph needs strcasecmp() or _stricmp()"
    #endif
#endif

#ifndef HAVE_STRNCASECMP
    #ifdef HAVE__STRNICMP
        #define strncasecmp _strnicmp
    #else
        #error "igraph needs strncasecmp() or _strnicmp()"
    #endif
#endif

/* Magic macro to fail the build if certain condition does not hold. See:
 * https://stackoverflow.com/questions/4079243/how-can-i-use-sizeof-in-a-preprocessor-macro
 */
#define IGRAPH_STATIC_ASSERT(condition) ((void)sizeof(char[1 - 2*!(condition)]))

__END_DECLS

#endif
