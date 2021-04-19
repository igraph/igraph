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
#ifndef _POWERLAW_H
#define _POWERLAW_H

// pascalou
#ifndef pascalou
    #include "gengraph_definitions.h"
#endif

// Discrete integer power-law : P(X=min+k) is proportionnal to (k+k0)^-alpha
// - possibility to determine a range [Min, Max] of possible samples
// - possibility to automatically compute k0 to obtain a given mean z

namespace gengraph {

#define POWERLAW_TABLE 10000

class powerlaw {
private:
    double alpha;  // Exponent
    int mini; // Minimum sample
    int maxi; // Maximum sample
    double offset; // Offset
    int tabulated; // Number of values to tabulate
    int *table;    // Table containing cumulative distribution for k=mini..mini+tabulated-1
    int *dt;        // Table delimiters
    int max_dt;     // number of delimiters - 1
    double proba_big;   // Probability to take a non-tabulated value
    double table_mul;   // equal to (1-proba_big)/(RAND_MAX+1)

    // Sample a non-tabulated value >= mini+tabulated
    inline double big_sample(double randomfloat) {
        return double(mini) + pow(_a * randomfloat + _b, _exp) - offset;
    }
    inline double big_inv_sample(double s) {
        return (pow(s - double(mini) + offset, 1.0 / _exp) - _b) / _a;
    }
    double _exp, _a, _b; // Cached values used by big_sample();

    // Dichotomic adjust of offset, so that to_adjust() returns value with
    // a precision of eps. Note that to_adjust() must be an increasing function of offset.
    void adjust_offset_mean(double value, double eps, double fac);

public:
    int sample();      // Return a random integer
    double proba(int); // Return probability to return integer
    double error();    // Returns relative numerical error done by this class
    double mean();     // Returns mean of the sampler
    int median();      // Returns median of the sampler

    // Initialize the power-law sampler.
    void init_to_offset(double, int);
    // Same, but also returns the offset found
    double init_to_mean(double);
    double init_to_median(double);

    inline void init() {
        init_to_offset(double(mini), POWERLAW_TABLE);
    };

    ~powerlaw();
    powerlaw(double exponent, int mini, int maxi = -1);
};

} // namespace gengraph

#endif //_POWERLAW_H
