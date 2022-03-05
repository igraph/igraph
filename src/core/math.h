/*
   IGraph library.
   Copyright (C) 2022  The igraph development team <igraph@igraph.org>

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

#ifndef IGRAPH_CORE_MATH_H
#define IGRAPH_CORE_MATH_H

/* Use math constants with MSVC */
#if defined(_MSC_VER) && !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES
#endif

#include <math.h>


/* Math constants are not part of standard C */
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
    #define M_PI_2 1.57079632679489661923
#endif
#ifndef M_LN2
    #define M_LN2 0.69314718055994530942
#endif
#ifndef M_SQRT2
    #define M_SQRT2 1.4142135623730950488016887
#endif
#ifndef M_LN_SQRT_2PI
    #define M_LN_SQRT_2PI 0.918938533204672741780329736406 /* log(sqrt(2*pi)) == log(2*pi)/2 */
#endif

#endif /* IGRAPH_CORE_MATH_H */
