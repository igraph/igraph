/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "convert.h"
#include "format.h"

namespace infomap {
namespace io {

  std::string toPrecision(double value, unsigned int precision, bool fixed)
  {
    // {:.{}f} == std::fixed << std::setprecision(precision); {:.{}g} ==
    // std::setprecision(precision) with the default floatfield. Verified
    // byte-for-byte (including rounding and subnormal edge cases) against the
    // previous std::ostringstream implementation. Kept in a .cpp so the fmt
    // header stays out of the widely-included convert.h.
    return fixed ? fmt::format(FMT_STRING("{:.{}f}"), value, precision)
                 : fmt::format(FMT_STRING("{:.{}g}"), value, precision);
  }

} // namespace io
} // namespace infomap
