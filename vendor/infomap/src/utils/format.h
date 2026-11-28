/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_FORMAT_H_
#define INFOMAP_FORMAT_H_

// Single entry point for the full {fmt} API (fmt::format and friends). Always
// include this wrapper instead of <fmt/format.h> directly, so FMT_HEADER_ONLY is
// defined identically on every build surface (CLI, Python, R, JS) and we never
// get an ODR mismatch. Headers that only need fmt's typed-argument machinery
// (not the renderer) should include the lighter utils/format_core.h instead.
//
// Only the iostream-free formatting core is vendored (vendor/fmt). Use
// fmt::format to build strings, then route them through infomap::Log. Never use
// fmt's print/println family (they write to stdout, bypassing the Log sink and
// breaking the R/Rprintf redirect) and never include fmt/ostream.h. This is
// enforced by scripts/check_cpp_stream_policy.py.

#include "format_core.h" // FMT_HEADER_ONLY + <stdlib.h> (global ::malloc) + <fmt/base.h>

#include <fmt/format.h> // IWYU pragma: export

#endif // INFOMAP_FORMAT_H_
