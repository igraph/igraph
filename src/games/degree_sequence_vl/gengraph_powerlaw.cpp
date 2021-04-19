/*
 *
 * gengraph - generation of random simple connected graphs with prescribed
 *            degree sequence
 *
 * Copyright (C) 2006  Fabien Viger
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
// Pascalou ...
#ifdef pascalou
    #define my_random() random()
    #define MY_RAND_MAX 0x7FFFFFFF
#else
    #include "gengraph_definitions.h"
#endif

#include "gengraph_powerlaw.h"
#include <cstdio>
#include <cmath>
#include <cassert>

#include "igraph_error.h"

namespace gengraph {

// Destructor
powerlaw::~powerlaw() {
    delete[] table;
    if (dt != NULL) {
        delete[] dt;
    }
}

// Constructor
powerlaw::powerlaw(double _alpha, int _mini, int _maxi) {
    alpha = _alpha;
    mini = _mini;
    maxi = _maxi;
    if (alpha <= 2.0 && maxi < 0)
        igraph_warningf("powerlaw exponent %f should be > 2 when no "
                        "Maximum is specified", IGRAPH_FILE_BASENAME, __LINE__, -1, alpha);
    if (alpha <= 1.0 && maxi >= 0)
        igraph_warningf("powerlaw exponent %f should be > 1", IGRAPH_FILE_BASENAME, __LINE__,
                        -1, alpha);
    if (maxi >= 0 && mini > maxi)
        igraph_warningf("powerlaw max %d should be greater than min %d",
                        IGRAPH_FILE_BASENAME, __LINE__, -1, maxi, mini);
    table = new int[POWERLAW_TABLE];
    tabulated = 0;
    dt = NULL;
}

// Sample
int powerlaw::sample() {
    if (proba_big != 0 && test_proba(proba_big)) {
        return int(floor(0.5 + big_sample(random_float())));
    }
    int r = my_random();
    // table[] contains integer from MY_RAND_MAX downto 0, in blocks. Search block...
    if (r > (MY_RAND_MAX >> max_dt)) {
        return mini;
    }
    int k = 0;
    while (k < max_dt) {
        r <<= 1;
        r += random_bit();
        k++;
    };
    int a = 0;
    int b;
    while ((b = dt[k++]) < 0 || r < table[b]) {
        if (b >= 0) {
            a = b + 1;
            if (a == tabulated - 1) {
                break;
            }
            r <<= 1;
            r += random_bit();
        }
    }

    // Now that we found the good block, run a dichotomy on this block [a,b]
    while (a < b) {
        int c = (a + b) / 2;
        if (r < table[c]) {
            a = c + 1;
        } else {
            b = c;
        }
    }
    return mini + a;
}

// Proba
double powerlaw::proba(int k) {
    if (k < mini || (maxi >= 0 && k > maxi)) {
        return 0.0;
    }
    if (k >= mini + tabulated) {
        return proba_big * (big_inv_sample(double(k) - 0.5) - big_inv_sample(double(k) + 0.5));
    } else {
        double div = table_mul;
        int prev_pos_in_table = k - mini - 1;
        if (prev_pos_in_table < 0) {
            return (double(MY_RAND_MAX) + 1.0 - double(table[0] >> max_dt)) * div;
        }
        // what block are we in ?
        int k1 = 0;
        while (k1 < max_dt) {
            div *= 0.5;
            k1++;
        };
        while (dt[k1] < 0 || dt[k1] < prev_pos_in_table) {
            k1++;
            div *= 0.5;
        };
        double prob2 = double(table[prev_pos_in_table + 1]);
        if (dt[k1] == prev_pos_in_table) do {
                prob2 *= 0.5;
            } while (dt[++k1] < 0);
        return (double(table[prev_pos_in_table]) - prob2) * div;
    }
}

// Relative Error
double powerlaw::error() {
    return 1.0 / (double(tabulated) * double(tabulated));
}

// Mean
double powerlaw::mean() {
    double sum = 0.0;
    for (int i = mini + tabulated; --i >= mini; ) {
        sum += double(i) * proba(i);
    }
    // add proba_big * integral(big_sample(t),t=0..1)
    if (proba_big != 0) {
        sum += proba_big * ((pow(_a + _b, _exp + 1.0) - pow(_b, _exp + 1.0)) / (_a * (_exp + 1.0)) + double(mini) - offset - sum);
    }
    return sum;
}

// Median. Returns integer Med such that P(X<=Med) >= 1/2
int powerlaw::median() {
    if (proba_big > 0.5) {
        return int(floor(0.5 + big_sample(1.0 - 0.5 / proba_big)));
    }
    double sum = 0.0;
    int i = mini;
    while (sum < 0.5) {
        sum += proba(i++);
    }
    return i - 1;
}

void powerlaw::init_to_offset(double _offset, int _tabulated) {
    offset = _offset;
    tabulated = _tabulated;
    if (maxi >= 0 && tabulated > maxi - mini) {
        tabulated = maxi - mini + 1;
    }
    double sum = 0.0;
    double item = double(tabulated) + offset;
    // Compute sum of tabulated probabilities
    for (int i = tabulated; i--; ) {
        sum += pow(item -= 1.0, -alpha);
    }
    // Compute others parameters : proba_big, table_mul, _a, _b, _exp
    if (maxi > 0 && maxi <= mini + tabulated - 1) {
        proba_big = 0;
        table_mul = inv_RANDMAX;
    } else {
        if (maxi < 0) {
            _b = 0.0;
        } else {
            _b = pow(double(maxi - mini) + 0.5 + offset, 1.0 - alpha);
        }
        _a = pow(double(tabulated) - 0.5 + offset, 1.0 - alpha) - _b;
        _exp = 1.0 / (1.0 - alpha);
        double sum_big = _a * (-_exp);
        proba_big = sum_big / (sum + sum_big);
        table_mul = inv_RANDMAX * sum / (sum + sum_big);
    }
    // How many delimiters will be necessary for the table ?
    max_dt = max(0, int(floor(alpha * log(double(tabulated)) / log(2.0))) - 6);
    if (dt != NULL) {
        delete[] dt;
    }
    dt = new int[max_dt + 1];
    // Create table as decreasing integers from MY_RAND_MAX+1 (in virtual position -1) down to 0
    // Every time the index crosses a delimiter, numbers get doubled.
    double ssum = 0;
    double mul = (double(MY_RAND_MAX) + 1.0) * pow(2.0, max_dt) / sum;
    item = double(tabulated) + offset;
    int k = max_dt;
    dt[k--] = tabulated - 1;
    for (int i = tabulated; --i > 0; ) {
        table[i] = int(floor(0.5 + ssum));
        ssum += mul * pow(item -= 1.0, -alpha);
        if (ssum > double(MY_RAND_MAX / 2) && k >= 0) {
            while ((ssum *= 0.5) > double(MY_RAND_MAX / 2)) {
                mul *= 0.5;
                dt[k--] = -1;
            };
            mul *= 0.5; dt[k--] = i - 1;
        }
    }
    table[0] = int(floor(0.5 + ssum));
    max_dt = k + 1;
}

void powerlaw::adjust_offset_mean(double _mean, double err, double factor) {
    // Set two bounds for offset
    double ol = offset;
    double oh = offset;
    if (mean() < _mean) {
        do {
            ol = oh;
            oh *= factor;
            init_to_offset(oh, tabulated);
        } while (mean() < _mean);
    } else {
        do {
            oh = ol;
            ol /= factor;
            init_to_offset(ol, tabulated);
        } while (mean() > _mean);
    }
    // Now, dichotomy
    while (fabs(oh - ol) > err * ol) {
        double oc = sqrt(oh * ol);
        init_to_offset(oc, tabulated);
        if (mean() < _mean) {
            ol = oc;
        } else {
            oh = oc;
        }
    }
    init_to_offset(sqrt(ol * oh), tabulated);
}

double powerlaw::init_to_mean(double _mean) {
    if (maxi >= 0 && _mean >= 0.5 * double((mini + maxi))) {
        /* Cannot use IGRAPH_ERRORF() as this function does not
         * return an igraph error code. */
        igraph_errorf("Fatal error in powerlaw::init_to_mean(%f): "
                      "Mean must be in ]min, (min+max)/2[ = ]%d, %d[",
                      IGRAPH_FILE_BASENAME, __LINE__, IGRAPH_EINVAL,
                      _mean, mini, (mini + maxi) / 2);
        return (-1.0);
    }
    init_to_offset(_mean - double(mini), 100);
    adjust_offset_mean(_mean, 0.01, 2);
    init_to_offset(offset, POWERLAW_TABLE);
    double eps = 1.0 / (double(POWERLAW_TABLE));
    adjust_offset_mean(_mean, eps * eps, 1.01);
    return offset;
}

} // namespace gengraph
