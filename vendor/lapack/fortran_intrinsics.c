/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-12  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

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

#include <float.h>

double digitsdbl_(double *x) {
    return (double) DBL_MANT_DIG;
}

double epsilondbl_(double *x) {
    return DBL_EPSILON;
}

double hugedbl_(double *x) {
    return DBL_MAX;
}

double tinydbl_(double *x) {
    return DBL_MIN;
}

int maxexponentdbl_(double *x) {
    return DBL_MAX_EXP;
}

int minexponentdbl_(double *x) {
    return DBL_MIN_EXP;
}

double radixdbl_(double *x) {
    return (double) FLT_RADIX;
}

