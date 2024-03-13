/*
   IGraph library.
   Copyright (C) 2003-2024  The igraph development team <igraph@igraph.org>

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

/* Semi-public header; not included from igraph.h. This is intentional.
 * The macros defined in this header are usually not useful to the end user,
 * with some exceptions like the source of the R interface. You need to
 * include this header explicitly if needed. */

#ifndef IGRAPM_MISC_H
#define IGRAPM_MISC_H

#include "igraph_decls.h"

__BEGIN_DECLS

/* Magic macro to fail the build if certain condition does not hold. See:
 * https://stackoverflow.com/questions/4079243/how-can-i-use-sizeof-in-a-preprocessor-macro
 */
#define IGRAPH_STATIC_ASSERT(condition) ((void)sizeof(char[1 - 2*!(condition)]))

__END_DECLS

#endif
