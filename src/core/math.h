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

#include <math.h>

/* Math constants are not part of standard C */

/* The following definitions contain enough precision for
 * an IEEE-754 quadruple-precision floating point format. */

#ifndef M_E
#define M_E         2.71828182845904523536028747135266250
#endif

#ifndef M_LOG2E
#define M_LOG2E     1.44269504088896340735992468100189214
#endif

#ifndef M_LOG10E
#define M_LOG10E    0.434294481903251827651128918916605082
#endif

#ifndef M_LN2
#define M_LN2       0.693147180559945309417232121458176568
#endif

#ifndef M_LN10
#define M_LN10      2.30258509299404568401799145468436421
#endif

#ifndef M_PI
#define M_PI        3.14159265358979323846264338327950288
#endif

#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923132169163975144
#endif

#ifndef M_PI_4
#define M_PI_4      0.785398163397448309615660845819875721
#endif

#ifndef M_1_PI
#define M_1_PI      0.318309886183790671537767526745028724
#endif

#ifndef M_2_PI
#define M_2_PI      0.636619772367581343075535053490057448
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257389615890312154517
#endif

#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880168872420969808
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2   0.707106781186547524400844362104849039
#endif

#endif /* IGRAPH_CORE_MATH_H */
