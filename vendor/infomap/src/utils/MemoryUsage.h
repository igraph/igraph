/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef MEMORYUSAGE_H_
#define MEMORYUSAGE_H_

namespace infomap {

struct MemoryReport {
  double rssPeakMb = 0.0;
  bool haveBytesPerNode = false;
  bool haveBytesPerLink = false;
  double bytesPerNode = 0.0;
  double bytesPerLink = 0.0;
};

MemoryReport currentMemoryReport(unsigned long long numNodes, unsigned long long numLinks);

} // namespace infomap

#endif // MEMORYUSAGE_H_
