/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef TRIAL_RESULTS_H_
#define TRIAL_RESULTS_H_

#include <string>
#include <vector>

namespace infomap {

struct TrialResultEntry {
  unsigned int trial = 0; // global trial index (0-based)
  unsigned long seed = 0;
  double codelength = 0.0;
  unsigned int numTopModules = 0;
  unsigned int numLevels = 0;
  int thread = 0;
  double timeSec = 0.0;
};

struct TrialResultsFile {
  std::string networkFingerprint;
  std::string configFingerprint;
  std::string infomapVersion;
  unsigned long baseSeed = 0;
  unsigned int trialOffset = 0;
  unsigned int numTrials = 0;
  std::string bestTreeFile; // path relative to this results file's directory
  // The state-level companion for higher-order input, same relative form. Empty for
  // first-order networks. Without it a distributed higher-order run published only the
  // physical projection, where one node id appears in several modules, and the merged
  // .clu inherited that (#906).
  std::string bestStatesTreeFile;
  std::vector<TrialResultEntry> trials;
};

std::string serializeTrialResults(const TrialResultsFile& results);

} // namespace infomap

#endif // TRIAL_RESULTS_H_
