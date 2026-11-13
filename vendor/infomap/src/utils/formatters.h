/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_FORMATTERS_H_
#define INFOMAP_FORMATTERS_H_

// {fmt} formatter specializations for Infomap value types that already define
// operator<<(std::ostream&). Include this where such a type is passed to
// Log::print / fmt::format with a "{}" placeholder.
//
// fmt's own ostream bridge (fmt/ostream.h) is banned by the stream policy
// because it pulls <ostream>/FILE* into shared runtime code. Instead we reuse
// io::stringify (already iostream-based, confined to convert.*) once, behind a
// single OstreamFormatter adapter, so each type below is a one-line opt-in.

#include "format.h"
#include "convert.h"

#include "Date.h"
#include "Stopwatch.h"
#include "../io/Config.h" // FlowModel

namespace infomap {

template <typename T>
struct OstreamFormatter : fmt::formatter<std::string> {
  template <typename FormatContext>
  auto format(const T& value, FormatContext& ctx) const -> decltype(ctx.out())
  {
    return fmt::formatter<std::string>::format(io::stringify(value), ctx);
  }
};

} // namespace infomap

template <>
struct fmt::formatter<infomap::Date> : infomap::OstreamFormatter<infomap::Date> {};
template <>
struct fmt::formatter<infomap::ElapsedTime> : infomap::OstreamFormatter<infomap::ElapsedTime> {};
template <>
struct fmt::formatter<infomap::Stopwatch> : infomap::OstreamFormatter<infomap::Stopwatch> {};
template <>
struct fmt::formatter<infomap::FlowModel> : infomap::OstreamFormatter<infomap::FlowModel> {};

#endif // INFOMAP_FORMATTERS_H_
