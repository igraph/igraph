/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMATH_H_
#define INFOMATH_H_

#include <cmath>
#include <cstdlib>

namespace infomap {
namespace infomath {

  using std::log2;

  inline double plogp(double p)
  {
    return p > 0.0 ? p * log2(p) : 0.0;
  }

  inline double isEqual(double a, double b, double tol = 1e-8)
  {
    return std::abs(a - b) <= tol;
  }

  /**
   * Tsallis entropy S_q of a uniform probability distribution of length n
   */
  inline double tsallisEntropyUniform(double n, double q = 1)
  {
    if (isEqual(q, 1)) {
      return std::log2(n);
    }
    return 1 / (q - 1) * (1 - pow(n, (1 - q))) / std::log(2);
  }

  /**
   * Interpolate from linear (q = 0) to log (q = 1)
   * linlog(k, 0) = k
   * linlog(k, 1) = log2(k)
   */
  inline double linlog(double k, double q = 1)
  {
    double baseCorrection = q <= 1 ? (1 - q) * std::log(2) + q : 1;
    double offsetCorrection = q <= 1 ? 1 - q : 0;
    return tsallisEntropyUniform(k, q) * baseCorrection + offsetCorrection;
  }

} // namespace infomath
} // namespace infomap

#endif // INFOMATH_H_
