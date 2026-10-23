/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include <nlohmann/json.hpp>

#include "RunMetadata.h"
#include "Config.h"
#include "OutputPlan.h"
#include "../utils/format.h"
#include "../version.h"

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <utility>
#include <vector>
#include <algorithm>

namespace infomap {

namespace {

  constexpr unsigned long long FNV_OFFSET = 14695981039346656037ull;
  constexpr unsigned long long FNV_PRIME = 1099511628211ull;
  using Json = nlohmann::ordered_json;

  template <typename T>
  void addCanonicalNumber(Json& json, const std::string& key, const T& value)
  {
    json[key] = value;
  }

  void hashBytes(unsigned long long& hash, const char* data, std::streamsize size)
  {
    for (std::streamsize i = 0; i < size; ++i) {
      hash ^= static_cast<unsigned char>(data[i]);
      hash *= FNV_PRIME;
    }
  }

  std::string fnvHex(const std::string& value)
  {
    auto hash = FNV_OFFSET;
    hashBytes(hash, value.data(), static_cast<std::streamsize>(value.size()));
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
  }

  // FNV-1a 64 over the whole file, one streaming pass. It used to hash the size
  // plus the first and last 64 KiB, which was blind to a same-size edit anywhere
  // in the interior of a file larger than 128 KiB -- precisely the edit-and-rerun
  // workflow the fingerprint exists to catch (#1026). A full pass costs one read
  // of data the run is about to parse anyway.
  std::string fileContentFingerprint(const std::string& path)
  {
    std::ifstream input(path.c_str(), std::ios_base::binary);
    if (!input) {
      throw std::runtime_error(fmt::format(FMT_STRING("Cannot open input file '{}'. Check that the path points to a file and that you have read permissions."), path));
    }

    auto hash = FNV_OFFSET;
    constexpr std::streamsize CHUNK_SIZE = 65536;
    std::vector<char> buffer(static_cast<std::size_t>(CHUNK_SIZE));
    while (true) {
      input.read(buffer.data(), CHUNK_SIZE);
      const auto got = input.gcount();
      if (got > 0)
        hashBytes(hash, buffer.data(), got);
      if (got < CHUNK_SIZE)
        break;
    }
    // A short read is EOF *or* an I/O error, and the two must not be confused:
    // hashing only the readable prefix would publish a valid-looking identity
    // for content nobody read.
    if (input.bad()) {
      throw std::runtime_error(fmt::format(FMT_STRING("Cannot read input file '{}': the read failed partway through."), path));
    }

    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
  }

  std::string outputArtifactsJson(const Config& config)
  {
    Json json = Json::array();
    const auto appendArtifact = [&](const std::string& key, const std::string& path) {
      Json artifact;
      artifact["key"] = key;
      artifact["path"] = path;
      json.push_back(std::move(artifact));
    };
    for (auto phase : { OutputPhase::BeforeFlow, OutputPhase::AfterFlow, OutputPhase::AfterPartition }) {
      for (const auto& output : planOutputArtifacts(config, phase))
        appendArtifact(output.resultKey, output.filename);
    }
    for (const auto& report : planReportArtifacts(config))
      appendArtifact(report.first, report.second);
    return json.dump();
  }

} // namespace

std::string canonicalConfigJson(const Config& config)
{
  Json json;

  // Input identity (path, size, content) is captured separately by the input
  // fingerprints -- network, cluster data and metadata each get one in the
  // manifest and the artifact headers (#1026); the config fingerprint covers
  // only algorithm-affecting settings. That makes it stable across paths for
  // the *network* alone: the cluster-data and meta-data paths are still fields
  // below, so moving either of those unchanged files does change this
  // fingerprint. Their identities are published beside it, so a consumer that
  // wants path-independence can compare hashes instead.
  //
  // How the trial budget is *divided* is deliberately excluded -- numTrials and
  // trialOffset name which slice of one budget a run executed, not what it
  // executed. infomap.merge requires all shards of a run to share this
  // fingerprint, and shards differ in exactly those two fields, so including
  // them would reject every legitimate distributed run. The base seed is
  // included, and seeding is per-trial and offset-derived (see
  // RunSession::trialSeed), so shards that agree here really are slices of the
  // same sequence of trials.
  //
  // Output-only settings are excluded because they cannot change the result.
  // That held for trialResultsPath only after #905: it used to switch the
  // per-trial reseeding on, which made two runs with the same fingerprint
  // publish different partitions.
  addCanonicalNumber(json, "additional_input_count", config.additionalInput.size());
  json["flow_model"] = flowModelToString(config.flowModel);
  json["directed"] = config.directed;
  addCanonicalNumber(json, "seed", config.seedToRandomNumberGenerator);
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  json["lossy"] = config.lossy;
  addCanonicalNumber(json, "lossy_lambda", config.lossyLambda);
#endif
  json["two_level"] = config.twoLevel;
  json["no_infomap"] = config.noInfomap;
  json["regularized"] = config.regularized;
  addCanonicalNumber(json, "regularization_strength", config.regularizationStrength);
  json["recorded_teleportation"] = config.recordedTeleportation;
  addCanonicalNumber(json, "teleportation_probability", config.teleportationProbability);
  addCanonicalNumber(json, "markov_time", config.markovTime);
  json["variable_markov_time"] = config.variableMarkovTime;
  addCanonicalNumber(json, "variable_markov_damping", config.variableMarkovTimeDamping);
  addCanonicalNumber(json, "variable_markov_min_scale", config.variableMarkovTimeMinLocalScale);
  json["entropy_corrected"] = config.entropyBiasCorrection;
  addCanonicalNumber(json, "entropy_correction_strength", config.entropyBiasCorrectionMultiplier);
  json["use_node_weights_as_flow"] = config.useNodeWeightsAsFlow;
  json["teleport_to_nodes"] = config.teleportToNodes;
  addCanonicalNumber(json, "weight_threshold", config.weightThreshold);
  json["no_self_links"] = config.noSelfLinks;
  addCanonicalNumber(json, "node_limit", config.nodeLimit);
  addCanonicalNumber(json, "matchable_multilayer_ids", config.matchableMultilayerIds);
  json["bipartite"] = config.bipartite;
  json["bipartite_teleportation"] = config.bipartiteTeleportation;
  json["cluster_data"] = config.clusterDataFile;
  // Whether that initial partition is kept or optimized away, and where nodes the file
  // does not mention end up: both change the partition, and both were missing, so two
  // shards that disagreed on them fingerprinted identically and merged (#906).
  json["cluster_data_is_hard"] = config.clusterDataIsHard;
  json["assign_to_neighbouring_module"] = config.assignToNeighbouringModule;
  json["meta_data"] = config.metaDataFile;
  addCanonicalNumber(json, "meta_data_rate", config.metaDataRate);
  json["meta_data_unweighted"] = config.unweightedMetaData;
  addCanonicalNumber(json, "preferred_number_of_modules", config.preferredNumberOfModules);
  addCanonicalNumber(json, "preferred_number_of_levels", config.preferredNumberOfLevels);
  addCanonicalNumber(json, "preferred_number_of_levels_strength", config.preferredNumberOfLevelsStrength);
  json["parallel_trials"] = config.parallelTrials;
  json["converge_trials"] = config.convergeTrials;
  json["inner_parallelization"] = config.innerParallelization;
  addCanonicalNumber(json, "core_loop_limit", config.coreLoopLimit);
  addCanonicalNumber(json, "core_level_limit", config.levelAggregationLimit);
  addCanonicalNumber(json, "tune_iteration_limit", config.tuneIterationLimit);
  addCanonicalNumber(json, "core_loop_codelength_threshold", config.minimumCodelengthImprovement);
  addCanonicalNumber(json, "tune_iteration_relative_threshold", config.minimumRelativeTuneIterationImprovement);
  addCanonicalNumber(json, "max_flow_iterations", config.maxFlowIterations);
  addCanonicalNumber(json, "min_flow_iterations", config.minFlowIterations);
  addCanonicalNumber(json, "flow_tolerance", config.flowTolerance);
  json["prefer_modular_solution"] = config.preferModularSolution;
  addCanonicalNumber(json, "num_random_moves", config.numRandomMoves);
  addCanonicalNumber(json, "max_degree_for_random_moves", config.maxDegreeForRandomMoves);
  json["skip_adjust_bipartite_flow"] = config.skipAdjustBipartiteFlow;
  json["markov_time_no_self_links"] = config.markovTimeNoSelfLinks;
  addCanonicalNumber(json, "multilayer_relax_rate", config.multilayerRelaxRate);
  addCanonicalNumber(json, "multilayer_relax_limit", config.multilayerRelaxLimit);
  addCanonicalNumber(json, "multilayer_relax_limit_up", config.multilayerRelaxLimitUp);
  addCanonicalNumber(json, "multilayer_relax_limit_down", config.multilayerRelaxLimitDown);
  addCanonicalNumber(json, "multilayer_js_relax_rate", config.multilayerJSRelaxRate);
  json["multilayer_relax_to_self"] = config.multilayerRelaxToSelf;
  json["multilayer_relax_by_jsd"] = config.multilayerRelaxByJensenShannonDivergence;
  addCanonicalNumber(json, "multilayer_js_relax_limit", config.multilayerJSRelaxLimit);
  json["no_coarse_tune"] = config.noCoarseTune;
  json["only_super_modules"] = config.onlySuperModules;
  addCanonicalNumber(json, "fast_hierarchical_solution", config.fastHierarchicalSolution);
  json["randomize_core_loop_limit"] = config.randomizeCoreLoopLimit;
  addCanonicalNumber(json, "minimum_single_node_codelength_improvement", config.minimumSingleNodeCodelengthImprovement);

  return json.dump();
}

std::string configFingerprint(const Config& config)
{
  return fnvHex(canonicalConfigJson(config));
}

InputIdentity inputIdentity(const std::string& path)
{
  InputIdentity identity;
  if (path.empty())
    return identity;

  struct stat info;
  if (stat(path.c_str(), &info) != 0) {
    // Not an error here: a missing --cluster-data or network file is reported
    // by the reader that opens it, as a parse error the bindings classify
    // (NetworkParseError). Capturing identity first must not change that.
    return identity;
  }

  identity.path = path;
  identity.size = static_cast<unsigned long long>(info.st_size);
  identity.mtime = static_cast<long long>(info.st_mtime);
  identity.hash = fileContentFingerprint(path);
  return identity;
}

std::string inputIdentityJson(const InputIdentity& identity)
{
  if (!identity.known())
    return "null";
  Json json;
  json["path"] = identity.path;
  json["size"] = identity.size;
  json["mtime"] = identity.mtime;
  json["hash"] = identity.hash;
  return json.dump();
}

std::string inputFingerprintJson(const std::string& path)
{
  return inputIdentityJson(inputIdentity(path));
}

std::string networkFingerprint(const std::string& path)
{
  return inputIdentity(path).hash;
}

std::string runManifestJson(const Config& config, const RunIdentities& identities)
{
  const auto canonicalConfig = canonicalConfigJson(config);
  Json json;
  json["version"] = INFOMAP_VERSION;
  json["command"] = config.parsedString;
  json["num_trials"] = config.numTrials;
  json["config"] = Json::parse(canonicalConfig);
  json["config_fingerprint"] = identities.configFingerprint.empty() ? fnvHex(canonicalConfig) : identities.configFingerprint;
  // The identities the run captured when it started, not a re-hash of whatever
  // the paths point at now: re-reading could disagree if a file changed since,
  // and a network built through the bindings has no path in the config to
  // re-read at all, so the manifest used to report a null input for a run whose
  // headers named a file (#1026).
  json["input"] = Json::parse(inputIdentityJson(identities.input));
  // Every file that changes the published partition, not only the network:
  // a --cluster-data seed and a --meta-data file used to be recorded by path
  // alone, so two runs on different files could not be told apart.
  json["cluster_data"] = Json::parse(inputIdentityJson(identities.clusterData));
  json["meta_data"] = Json::parse(inputIdentityJson(identities.metaData));
  json["outputs"] = Json::parse(outputArtifactsJson(config));
  return json.dump() + '\n';
}

std::string runManifestJson(const Config& config)
{
  RunIdentities identities;
  identities.input = inputIdentity(config.networkFile);
  identities.clusterData = inputIdentity(config.clusterDataFile);
  identities.metaData = inputIdentity(config.metaDataFile);
  return runManifestJson(config, identities);
}

} // namespace infomap
