/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "MemoryUsage.h"

#if defined(__APPLE__) || defined(__unix__) || defined(__unix)
#include <sys/resource.h>
#endif

namespace infomap {

namespace {

  unsigned long long peakRssBytes()
  {
#if defined(__APPLE__)
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
      return 0;
    }
    return static_cast<unsigned long long>(usage.ru_maxrss);
#elif defined(__unix__) || defined(__unix)
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
      return 0;
    }
    return static_cast<unsigned long long>(usage.ru_maxrss) * 1024ULL;
#else
    return 0;
#endif
  }

} // namespace

MemoryReport currentMemoryReport(unsigned long long numNodes, unsigned long long numLinks)
{
  MemoryReport report;
  const auto rssBytes = peakRssBytes();
  report.rssPeakMb = static_cast<double>(rssBytes) / (1024.0 * 1024.0);
  if (numNodes > 0) {
    report.haveBytesPerNode = true;
    report.bytesPerNode = static_cast<double>(rssBytes) / static_cast<double>(numNodes);
  }
  if (numLinks > 0) {
    report.haveBytesPerLink = true;
    report.bytesPerLink = static_cast<double>(rssBytes) / static_cast<double>(numLinks);
  }
  return report;
}

} // namespace infomap
