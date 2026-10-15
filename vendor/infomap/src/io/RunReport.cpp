/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include <nlohmann/json.hpp>

#include "RunReport.h"
#include "SafeFile.h"
#include "../utils/Log.h"

#include <utility>

namespace {

using Json = nlohmann::ordered_json;

} // namespace

namespace infomap {

namespace {

  std::string dumpJsonLine(const Json& json) { return json.dump() + '\n'; }

} // namespace

std::string runSummaryReportJson(const RunSummaryReport& report)
{
  Json json;
  json["version"] = report.version;
  json["codelength"] = report.codelength;
  json["top_modules"] = report.topModules;
  json["levels"] = report.levels;
  json["trials"] = report.trials;
  json["best_trial"] = report.bestTrial;
  json["auto_stopped"] = report.autoStopped;
  if (report.haveFlowConvergence) {
    Json flow;
    flow["converged"] = report.flowConverged;
    flow["iterations"] = report.flowIterations;
    flow["error"] = report.flowError;
    json["flow"] = std::move(flow);
  }
  json["trial_codelengths"] = report.trialCodelengths;
  json["trial_top_modules"] = report.trialTopModules;
  return dumpJsonLine(json);
}

std::string runTimingReportJson(const RunTimingReport& report)
{
  Json json;
  json["version"] = report.version;
  json["openmp"] = report.openmp;
  json["threads_requested"] = report.threadsRequested;
  json["threads_used"] = report.threadsUsed;
  json["thread_source"] = report.threadSource;

  Json network;
  network["nodes"] = report.network.nodes;
  network["links"] = report.network.links;
  network["directed"] = report.network.directed;
  json["network"] = std::move(network);

  Json timing;
  for (const auto& item : report.timing) {
    timing[item.first] = item.second;
  }
  json["timing"] = std::move(timing);

  Json trials = Json::array();
  for (const auto& trial : report.trials) {
    if (!trial.valid) {
      continue;
    }
    Json item;
    item["trial"] = trial.trial;
    item["thread"] = trial.thread;
    item["seed"] = trial.seed;
    item["time_s"] = trial.timeSec;
    item["codelength"] = trial.codelength;
    item["top_modules"] = trial.topModules;
    item["num_levels"] = trial.numLevels;
    trials.push_back(std::move(item));
  }
  json["trials"] = std::move(trials);

  if (report.includeMemory) {
    Json memory;
    memory["rss_peak_mb"] = report.memory.rssPeakMb;
    if (report.memory.haveBytesPerNode) {
      memory["bytes_per_node"] = report.memory.bytesPerNode;
    }
    if (report.memory.haveBytesPerLink) {
      memory["bytes_per_link"] = report.memory.bytesPerLink;
    }
    json["memory"] = std::move(memory);
  }

  return dumpJsonLine(json);
}

void writeJsonReport(const std::string& path, const std::string& json, bool overwrite)
{
  if (path == "-") {
    Log::getOutputStream() << json;
    return;
  }

  SafeOutFile output(path, std::ios_base::out, overwrite);
  output << json;
  output.commit();
}

} // namespace infomap
