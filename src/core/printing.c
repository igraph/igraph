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

#include "igraph_complex.h"
#include "igraph_error.h"
#include "igraph_types.h"

#include <float.h>

/* The number of digits chosen here will be used in all places where
 * igraph_real_fprintf_precise() is used, including all textual graph
 * formats such as GML, GraphML, Pajek, etc. DBL_DIG digits are sufficient
 * to preserve the decimal representation during a
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
    #define IGRAPH_REAL_PRINTF_PRECISE_FORMAT "%." IGRAPH_I_STRINGIFY(DBL_DIG) "g"
#else
    /* Assume a precision of 15 digits for %g, which is what IEEE-754 doubles require. */
    #define IGRAPH_REAL_PRINTF_PRECISE_FORMAT "%.15g"
#endif

int igraph_real_fprintf(FILE *file, igraph_real_t val) {
    if (isfinite(val)) {
        return fprintf(file, "%g", val);
    } else if (isnan(val)) {
        return fprintf(file, "NaN");
    } else if (isinf(val)) {
        if (val < 0) {
            return fprintf(file, "-Inf");
        } else {
            return fprintf(file, "Inf");
        }
    }
    IGRAPH_FATAL("Value is not finite, not infinite and not NaN either!");  /* LCOV_EXCL_LINE */
}

#ifndef USING_R
int igraph_real_printf(igraph_real_t val) {
    return igraph_real_fprintf(stdout, val);
}
#endif

int igraph_real_fprintf_aligned(FILE *file, int width, igraph_real_t val) {
    if (isfinite(val)) {
        return fprintf(file, "%*g", width, val);
    } else if (isnan(val)) {
        return fprintf(file, "%*s", width, "NaN");
    } else if (isinf(val)) {
        if (val < 0) {
            return fprintf(file, "%*s", width, "-Inf");
        } else {
            return fprintf(file, "%*s", width, "Inf");
        }
    }
    IGRAPH_FATAL("Value is not finite, not infinite and not NaN either!");  /* LCOV_EXCL_LINE */
}

#ifndef USING_R
int igraph_real_printf_aligned(int width, igraph_real_t val) {
    return igraph_real_fprintf_aligned(stdout, width, val);
}
#endif

int igraph_real_snprintf(char *str, size_t size, igraph_real_t val) {
    if (isfinite(val)) {
        return snprintf(str, size, "%g", val);
    } else if (isnan(val)) {
        return snprintf(str, size, "NaN");
    } else if (isinf(val)) {
        if (val < 0) {
            return snprintf(str, size, "-Inf");
        } else {
            return snprintf(str, size, "Inf");
        }
    }
    IGRAPH_FATAL("Value is not finite, not infinite and not NaN either!");  /* LCOV_EXCL_LINE */
}

int igraph_real_fprintf_precise(FILE *file, igraph_real_t val) {
    if (isfinite(val)) {
        return fprintf(file, IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    } else if (isnan(val)) {
        return fprintf(file, "NaN");
    } else if (isinf(val)) {
        if (val < 0) {
            return fprintf(file, "-Inf");
        } else {
            return fprintf(file, "Inf");
        }
    }
    IGRAPH_FATAL("Value is not finite, not infinite and not NaN either!");  /* LCOV_EXCL_LINE */
}

#ifndef USING_R
int igraph_real_printf_precise(igraph_real_t val) {
    return igraph_real_fprintf_precise(stdout, val);
}
#endif

int igraph_real_snprintf_precise(char *str, size_t size, igraph_real_t val) {
    if (isfinite(val)) {
        return snprintf(str, size, IGRAPH_REAL_PRINTF_PRECISE_FORMAT, val);
    } else if (isnan(val)) {
        return snprintf(str, size, "NaN");
    } else if (isinf(val)) {
        if (val < 0) {
            return snprintf(str, size, "-Inf");
        } else {
            return snprintf(str, size, "Inf");
        }
    }
    IGRAPH_FATAL("Value is not finite, not infinite and not NaN either!");  /* LCOV_EXCL_LINE */
}

#define PROPAGATE() \
    do { \
        if (res < 0) { \
            return -1; \
        } \
        cnt += res; \
    } while (0)

int igraph_complex_fprintf(FILE *file, igraph_complex_t val) {
    int res, cnt = 0;
    igraph_real_t re = IGRAPH_REAL(val), im = IGRAPH_IMAG(val);
    res = igraph_real_fprintf(file, re);
    PROPAGATE();
    if (! signbit(im)) {
        res = fprintf(file, "+");
        PROPAGATE();
    }
    res = igraph_real_fprintf(file, im);
    PROPAGATE();
    res = fprintf(file, "i");
    PROPAGATE();
    return cnt;
}

#undef PROPAGATE

#ifndef USING_R
int igraph_complex_printf(igraph_complex_t val) {
    return igraph_complex_fprintf(stdout, val);
}
#endif

#define PROPAGATE() \
    do { \
        if (res < 0) { \
            return -1; \
        } \
        cnt += res; \
        /* remember that 'size' is unsigned, can't check if size - res < 0! */ \
        if (size > res) size -= res; \
        else size = 0; \
        if (size == 0) str = NULL; else str += res; \
    } while (0)

int igraph_complex_snprintf(char *str, size_t size, igraph_complex_t val) {
    int res, cnt = 0;
    igraph_real_t re = IGRAPH_REAL(val), im = IGRAPH_IMAG(val);
    res = igraph_real_snprintf(str, size, re);
    PROPAGATE();
    if (! signbit(im)) {
        res = snprintf(str, size, "+");
        PROPAGATE();
    }
    res = igraph_real_snprintf(str, size, im);
    PROPAGATE();
    res = snprintf(str, size, "i");
    PROPAGATE();
    return cnt;
}

int igraph_complex_fprintf_aligned(FILE *file, int width, igraph_complex_t val) {
    /* Most characters produces by %g is 13, so including 'i' and null terminator we
     * need up to 13 + 13 + 1 + 1 = 28 characters in total. */
    char buf[28];

    if (igraph_complex_snprintf(buf, sizeof(buf) / sizeof(buf[0]), val) < 0) {
        return -1;
    }
    return fprintf(file, "%*s", width, buf);
}

#ifndef USING_R
int igraph_complex_printf_aligned(int width, igraph_complex_t val) {
    return igraph_complex_fprintf_aligned(stdout, width, val);
}
#endif

#undef PROPAGATE
