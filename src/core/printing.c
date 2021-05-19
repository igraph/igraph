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

#ifdef DBL_DIG
    /* Use DBL_DIG to determine the maximum precision used for %g */
    #define STRINGIFY_HELPER(x) #x
    #define STRINGIFY(x) STRINGIFY_HELPER(x)
    #define IGRAPH_REAL_PRINTF_PRECISE_FORMAT "%." STRINGIFY(DBL_DIG) "g"
#else
    /* Assume a precision of 10 digits for %g */
    #define IGRAPH_REAL_PRINTF_PRECISE_FORMAT "%.10g"
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
