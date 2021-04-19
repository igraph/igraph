/* kolmogorov.c
 *
 * Copyright (C) 2010-2011 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <math.h>
#include "kolmogorov.h"

double plfit_kolmogorov(double z) {
    const double fj[4] = { -2, -8, -18, -32 };
    const double w = 2.50662827;
    const double c1 = -1.2337005501361697;   /* -pi^2 / 8 */
    const double c2 = -11.103304951225528;   /*  9*c1 */
    const double c3 = -30.842513753404244;   /* 25*c1 */

    double u = fabs(z);
    double v;

    if (u < 0.2)
        return 1;

    if (u < 0.755) {
        v = 1.0 / (u*u);
        return 1 - w * (exp(c1*v) + exp(c2*v) + exp(c3*v)) / u;
    }

    if (u < 6.8116) {
        double r[4] = { 0, 0, 0, 0 };
        long int maxj = (long int)(3.0 / u + 0.5);
        long int j;

        if (maxj < 1)
            maxj = 1;

        v = u*u;
        for (j = 0; j < maxj; j++) {
            r[j] = exp(fj[j] * v);
        }

        return 2*(r[0] - r[1] + r[2] - r[3]);
    }

    return 0;
}

double plfit_ks_test_one_sample_p(double d, size_t n) {
    return plfit_kolmogorov(d * sqrt(n));
}

double plfit_ks_test_two_sample_p(double d, size_t n1, size_t n2) {
    return plfit_kolmogorov(d * sqrt(n1*n2 / ((double)(n1+n2))));
}
