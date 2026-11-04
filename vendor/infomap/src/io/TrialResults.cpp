/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include <nlohmann/json.hpp>

#include "TrialResults.h"

#include <utility>

namespace infomap {

namespace {

  using Json = nlohmann::ordered_json;

} // namespace

// ---------------------------------------------------------------------------
// serializeTrialResults
// ---------------------------------------------------------------------------

std::string serializeTrialResults(const TrialResultsFile& r)
{
  Json json;
  json["network_fingerprint"] = r.networkFingerprint;
  json["config_fingerprint"] = r.configFingerprint;
  json["infomap_version"] = r.infomapVersion;
  json["base_seed"] = r.baseSeed;
  json["trial_offset"] = r.trialOffset;
  json["num_trials"] = r.numTrials;
  json["best_tree_file"] = r.bestTreeFile;
  if (!r.bestStatesTreeFile.empty())
    json["best_states_tree_file"] = r.bestStatesTreeFile;

  Json trials = Json::array();
  for (const auto& t : r.trials) {
    Json trial;
    trial["trial"] = t.trial;
    trial["seed"] = t.seed;
    trial["codelength"] = t.codelength;
    trial["num_top_modules"] = t.numTopModules;
    trial["num_levels"] = t.numLevels;
    trial["thread"] = t.thread;
    trial["time_s"] = t.timeSec;
    trials.push_back(std::move(trial));
  }
  json["trials"] = std::move(trials);

  return json.dump() + '\n';
}

} // namespace infomap
