/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2012  Gabor Csardi <csardi.gabor@gmail.com>
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

#include "igraph_types.h"

#include <float.h>

#ifdef _MSC_VER
    #define snprintf _snprintf
#endif

/* The number of digits chosen here will be used in all places where
 * igraph_real_fprintf_precise() is used, including all textual graph
 * formats such as GML, GraphML, Pajek, etc. DBL_DIG digits are sufficient
 * to preserve the decimal reprensentation during a
 * decimal (textual) -> binary -> decimal (textual) round-trip conversion.
 * This many digits are however not sufficient for a lossless
 * binary -> decimal -> binary conversion. Thus, writing numerical attributes
 * to a file and reading them back in may cause a tiny change in the last
 * binary digit of numbers. This change is minute, always smaller than 10^-15
 * times the original number, thus acceptable.
 *
 * We could output more digits, but that would come with its own problem:
 * It would sometimes cause a change in the decimal representation originally
 * input by users, which is surprising and confusing. For example,
 *
 * printf("%.17g\n", 100.1)
 *
 * outputs 100.09999999999999 instead of 100.1. We can prevent this by
 * using DBL_DIG == 15 digits instead of 17, which would be required
 * for a lossless binary -> decimal -> binary round-tripping.
 *
 * This justifies using DBL_DIG digits, and not more, in all places.
 */
#ifdef DBL_DIG
    /* Use DBL_DIG to determine the maximum precision used for %g */
    #define STRINGIFY_HELPER(x) #x
    #define STRINGIFY(x) STRINGIFY_HELPER(x)
    #define IGRAPH_REAL_PRINTF_PRECISE_FORMAT "%." STRINGIFY(DBL_DIG) "g"
#else
    /* Assume a precision of 15 digits for %g */
    #define IGRAPH_REAL_PRINTF_PRECISE_FORMAT "%.15g"
#endif

#ifndef USING_R
int igraph_real_printf(igraph_real_t val) {
    if (igraph_finite(val)) {
        return printf("%g", val);
    } else if (igraph_is_nan(val)) {
        return printf("NaN");
    } else if (igraph_is_inf(val)) {
        if (val < 0) {
            return printf("-Inf");
        } else {
            return printf("Inf");
        }
    } else {
        /* fallback */
        return printf("%g", val);
    }
}
#endif

int igraph_real_fprintf(FILE *file, igraph_real_t val) {
    if (igraph_finite(val)) {
        return fprintf(file, "%g", val);
    } else if (igraph_is_nan(val)) {
        return fprintf(file, "NaN");
    } else if (igraph_is_inf(val)) {
        if (val < 0) {
            return fprintf(file, "-Inf");
        } else {
            return fprintf(file, "Inf");
        }
    } else {
        /* fallback */
        return fprintf(file, "%g", val);
    }
}

int igraph_real_snprintf(char* str, size_t size, igraph_real_t val) {
    if (igraph_finite(val)) {
        return snprintf(str, size, "%g", val);
    } else if (igraph_is_nan(val)) {
        return snprintf(str, size, "NaN");
    } else if (igraph_is_inf(val)) {
        if (val < 0) {
            return snprintf(str, size, "-Inf");
        } else {
            return snprintf(str, size, "Inf");
        }
    } else {
        /* fallback */
        return snprintf(str, size, "%g", val);
    }
}

#ifndef USING_R
int igraph_real_printf_precise(igraph_real_t val) {
    if (igraph_finite(val)) {
        return printf(IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    } else if (igraph_is_nan(val)) {
        return printf("NaN");
    } else if (igraph_is_inf(val)) {
        if (val < 0) {
            return printf("-Inf");
        } else {
            return printf("Inf");
        }
    } else {
        /* fallback */
        return printf(IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    }
}
#endif

int igraph_real_fprintf_precise(FILE *file, igraph_real_t val) {
    if (igraph_finite(val)) {
        return fprintf(file, IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    } else if (igraph_is_nan(val)) {
        return fprintf(file, "NaN");
    } else if (igraph_is_inf(val)) {
        if (val < 0) {
            return fprintf(file, "-Inf");
        } else {
            return fprintf(file, "Inf");
        }
    } else {
        /* fallback */
        return fprintf(file, IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    }
}

int igraph_real_snprintf_precise(char* str, size_t size, igraph_real_t val) {
    if (igraph_finite(val)) {
        return snprintf(str, size, IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    } else if (igraph_is_nan(val)) {
        return snprintf(str, size, "NaN");
    } else if (igraph_is_inf(val)) {
        if (val < 0) {
            return snprintf(str, size, "-Inf");
        } else {
            return snprintf(str, size, "Inf");
        }
    } else {
        /* fallback */
        return snprintf(str, size, IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    }
}
