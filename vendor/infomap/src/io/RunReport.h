/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef RUNREPORT_H_
#define RUNREPORT_H_

#include "../utils/MemoryUsage.h"
#include "../utils/TimingRegistry.h"

#include <string>
#include <utility>
#include <vector>

namespace infomap {

struct RunReportNetwork {
  unsigned long long nodes = 0;
  unsigned long long links = 0;
  bool directed = false;
};

struct RunSummaryReport {
  std::string version;
  double codelength = 0.0;
  unsigned int topModules = 0;
  unsigned int levels = 0;
  unsigned int trials = 0;
  unsigned int bestTrial = 0;
  bool autoStopped = false; // --converge: trials stopped early on a codelength plateau
  // Power-iteration outcome, emitted only when one ran. A failed iteration means the
  // flow -- and therefore every codelength below -- is not the network's stationary
  // distribution, which no automated consumer could see before (#899).
  bool haveFlowConvergence = false;
  bool flowConverged = true;
  unsigned int flowIterations = 0;
  double flowError = 0.0;
  std::vector<double> trialCodelengths;
  std::vector<unsigned int> trialTopModules;
};

struct RunTimingReport {
  std::string version;
  bool openmp = false;
  unsigned int threadsRequested = 1;
  unsigned int threadsUsed = 1;
  std::string threadSource = "serial";
  RunReportNetwork network;
  std::vector<std::pair<std::string, double>> timing;
  std::vector<TrialTimingRecord> trials;
  bool includeMemory = false;
  MemoryReport memory;
};

std::string runSummaryReportJson(const RunSummaryReport& report);
std::string runTimingReportJson(const RunTimingReport& report);
void writeJsonReport(const std::string& path, const std::string& json, bool overwrite = true);

} // namespace infomap

#endif // RUNREPORT_H_
