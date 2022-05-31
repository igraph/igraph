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
 * The contents of this file were mechanically derived from pcg_variants.h
 * (every inline function defined there gets a generated extern declaration).
 */

#include "pcg_variants.h"

/*
 * Rotate helper functions.
 */

extern inline uint32_t pcg_rotr_32(uint32_t value, unsigned int rot);

/*
 * Output functions.  These are the core of the PCG generation scheme.
 */

/* XSH RS */

extern inline uint32_t pcg_output_xsh_rs_64_32(uint64_t state);

/* XSH RR */

extern inline uint32_t pcg_output_xsh_rr_64_32(uint64_t state);

/* RXS M XS */

extern inline uint32_t pcg_output_rxs_m_xs_32_32(uint32_t state);

/* RXS M */

extern inline uint32_t pcg_output_rxs_m_64_32(uint64_t state);

/* XSL RR (only defined for >= 64 bits) */

extern inline uint32_t pcg_output_xsl_rr_64_32(uint64_t state);

/* XSL RR RR (only defined for >= 64 bits) */
