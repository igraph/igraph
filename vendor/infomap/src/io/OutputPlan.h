/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef OUTPUT_PLAN_H_
#define OUTPUT_PLAN_H_

#include "Config.h"
#include "OutputFormats.h"

#include <string>
#include <utility>
#include <vector>

namespace infomap {

class InfomapBase;
class Network;

enum class OutputPhase : std::uint8_t {
  BeforeFlow,
  AfterFlow,
  AfterPartition
};

struct OutputArtifact {
  std::string resultKey;
  std::string label;
  std::string filename;
  OutputPhase phase = OutputPhase::AfterPartition;
  OutputKind kind = OutputKind::Tree;
  bool states = false;
  int cluLevel = 1;
  bool writeLinks = false;
  bool printFlow = false;
};

std::string outputPlanBasename(const Config& config, int trial = -1);

std::vector<OutputArtifact> planOutputArtifacts(const Config& config, const std::string& basename, OutputPhase phase);

std::vector<OutputArtifact> planOutputArtifacts(const Config& config, OutputPhase phase, int trial = -1);

// JSON report sidecars (summary, timing, run manifest) as {resultKey, path} pairs.
// Stdout targets ("-") and unset paths are excluded. Shared by the run manifest
// and the no-overwrite pre-flight so the two never drift apart.
std::vector<std::pair<std::string, std::string>> planReportArtifacts(const Config& config);

// Whether the network the run is about to partition is higher-order, as
// configureNetworkMode() will classify it once the network is read. The pre-flight
// runs before that classification but has to plan the paths it produces, so the
// answer is passed in rather than read off Config::stateOutput, which is still
// false at that point (#1018).
enum class HigherOrderInput : std::uint8_t {
  No,
  Yes
};

// Every file path a CLI run would write: modular artifacts across all phases
// (expanded per trial when --print-all-trials writes separate files) plus the
// report sidecars. Used by the no-overwrite pre-flight.
//
// `higherOrder` stands in for Config::stateOutput, which the caller cannot have
// set yet. It applies to every phase: all three are written after
// configureNetworkMode() installs the flag.
std::vector<std::string> planAllOutputPaths(const Config& config, HigherOrderInput higherOrder);

// Throws InfomapError(OutputError) when a planned output path is one of the run's
// own input files, and additionally when --no-overwrite is set and any planned
// output path already exists. No-op when overwriting is allowed and nothing
// collides with an input.
void preflightOutputTargets(const Config& config, HigherOrderInput higherOrder);

void writeOutputArtifact(InfomapBase& infomap, Network& network, const OutputArtifact& artifact);

void writeOutputArtifacts(InfomapBase& infomap, Network& network, OutputPhase phase, int trial = -1);

} // namespace infomap

#endif // OUTPUT_PLAN_H_
