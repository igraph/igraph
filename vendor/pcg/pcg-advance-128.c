/*
 * PCG Random Number Generation for C.
 *
 * Copyright 2014-2019 Melissa O'Neill <oneill@pcg-random.org>,
 *                     and the PCG Project contributors.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 * Licensed under the Apache License, Version 2.0 (provided in
 * LICENSE-APACHE.txt and at http://www.apache.org/licenses/LICENSE-2.0)
 * or under the MIT license (provided in LICENSE-MIT.txt and at
 * http://opensource.org/licenses/MIT), at your option. This file may not
 * be copied, modified, or distributed except according to those terms.
 *
 * Distributed on an "AS IS" BASIS, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied.  See your chosen license for details.
 *
 * For additional information about the PCG random number generation scheme,
 * visit http://www.pcg-random.org/.
 */

/*
 * This code is derived from the canonical C++ PCG implementation, which
 * has many additional features and is preferable if you can use C++ in
 * your project.
 *
 * Repetative C code is derived using C preprocessor metaprogramming
 * techniques.
 */

#include "pcg_variants.h"

/* Multi-step advance functions (jump-ahead, jump-back)
 *
 * The method used here is based on Brown, "Random Number Generation
 * with Arbitrary Stride,", Transactions of the American Nuclear
 * Society (Nov. 1994).  The algorithm is very similar to fast
 * exponentiation.
 *
 * Even though delta is an unsigned integer, we can pass a
 * signed integer to go backwards, it just goes "the long way round".
 */

#if PCG_HAS_128BIT_OPS
pcg128_t pcg_advance_lcg_128(pcg128_t state, pcg128_t delta, pcg128_t cur_mult,
                             pcg128_t cur_plus)
{
    pcg128_t acc_mult = 1u;
    pcg128_t acc_plus = 0u;
    while (delta > 0) {
        if (delta & 1) {
            acc_mult *= cur_mult;
            acc_plus = acc_plus * cur_mult + cur_plus;
        }
        cur_plus = (cur_mult + 1) * cur_plus;
        cur_mult *= cur_mult;
        delta /= 2;
    }
    return acc_mult * state + acc_plus;
}
#endif
