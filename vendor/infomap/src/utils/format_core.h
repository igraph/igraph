/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_FORMAT_CORE_H_
#define INFOMAP_FORMAT_CORE_H_

// Lightweight {fmt} entry point: the core API only (format_string, string_view,
// make_format_args, formatter). Include this in widely-included headers (e.g.
// utils/Log.h) so they get fmt's typed-argument machinery WITHOUT pulling the
// heavy fmt/format.h into every translation unit; pair it with utils/format.h
// in the .cpp that actually renders. FMT_HEADER_ONLY is pinned identically to
// format.h to avoid ODR mismatches. Never include <fmt/...> directly elsewhere;
// see utils/format.h.
//
// The lightweight core is <fmt/base.h>. As of {fmt} 11 the old <fmt/core.h> is a
// compatibility shim that pulls in the full <fmt/format.h>, so including it here
// would defeat the purpose; include <fmt/base.h> directly.

// Header-only is mandatory: Infomap vendors only fmt's headers (vendor/fmt) and
// compiles no fmt source files, so a build surface that forces FMT_HEADER_ONLY
// off (or to some other value) would link-fail with undefined fmt symbols.
// Surface that as a clear compile error here rather than a cryptic link error,
// so every build surface (CLI, Python, R, JS) stays consistent by construction.
#if defined(FMT_HEADER_ONLY) && (FMT_HEADER_ONLY != 1)
#error "Infomap requires header-only {fmt}: FMT_HEADER_ONLY must be 1 (no fmt sources are compiled)."
#endif

#ifndef FMT_HEADER_ONLY
#define FMT_HEADER_ONLY 1
#endif

// Under FMT_HEADER_ONLY, <fmt/base.h> pulls in the full <fmt/format.h> at its
// end, whose detail::allocator calls *unqualified* malloc/free without including
// a header for them. That resolves on libstdc++ via transitive declarations but
// not on Emscripten's libc++ ("use of undeclared identifier 'malloc'"). Declare
// the global ::malloc/::free here, before fmt is pulled in, so every fmt-using
// TU is covered (this is the single header both format.h and Log.h include
// first). <cstdlib> is insufficient: it only guarantees the names in std::.
#include <stdlib.h> // NOLINT(modernize-deprecated-headers) -- need global ::malloc

#include <fmt/base.h> // IWYU pragma: export

#endif // INFOMAP_FORMAT_CORE_H_
