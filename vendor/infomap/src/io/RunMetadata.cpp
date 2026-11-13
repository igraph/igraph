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

#include <cstdio>
#include <memory>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <utility>
#include <vector>
#include <algorithm>
#include <cstdint>

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

  // XXH64 (seed 0), the 64-bit xxHash: https://github.com/Cyan4973/xxHash
  // It consumes the input eight bytes at a time in four independent lanes, where
  // byte-at-a-time FNV-1a cost about 7.5 instructions per input byte on every run
  // (#1167). The value is the one `xxhsum -H64` prints for the same file.
  constexpr std::uint64_t XXH_PRIME1 = 11400714785074694791ull;
  constexpr std::uint64_t XXH_PRIME2 = 14029467366897019727ull;
  constexpr std::uint64_t XXH_PRIME3 = 1609587929392839161ull;
  constexpr std::uint64_t XXH_PRIME4 = 9650029242287828579ull;
  constexpr std::uint64_t XXH_PRIME5 = 2870177450012600261ull;
  constexpr std::size_t XXH_STRIPE = 32;

  inline std::uint64_t rotl64(std::uint64_t x, int r) { return (x << r) | (x >> (64 - r)); }

  // Little-endian on every host; compilers fold this into a single load.
  inline std::uint64_t readLE64(const unsigned char* p)
  {
    std::uint64_t v = 0;
    for (int i = 7; i >= 0; --i)
      v = (v << 8) | p[i];
    return v;
  }

  inline std::uint32_t readLE32(const unsigned char* p)
  {
    return static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8) | (static_cast<std::uint32_t>(p[2]) << 16) | (static_cast<std::uint32_t>(p[3]) << 24);
  }

  inline std::uint64_t xxhRound(std::uint64_t acc, std::uint64_t input)
  {
    acc += input * XXH_PRIME2;
    acc = rotl64(acc, 31);
    return acc * XXH_PRIME1;
  }

  inline std::uint64_t xxhMergeRound(std::uint64_t acc, std::uint64_t val)
  {
    acc ^= xxhRound(0, val);
    return acc * XXH_PRIME1 + XXH_PRIME4;
  }

  class Xxh64 {
  public:
    // Whole stripes only; the caller passes any remainder to finish().
    void stripes(const unsigned char* p, std::size_t size)
    {
      for (const auto* end = p + size; p != end; p += XXH_STRIPE) {
        m_v1 = xxhRound(m_v1, readLE64(p));
        m_v2 = xxhRound(m_v2, readLE64(p + 8));
        m_v3 = xxhRound(m_v3, readLE64(p + 16));
        m_v4 = xxhRound(m_v4, readLE64(p + 24));
      }
      m_length += size;
    }

    std::uint64_t finish(const unsigned char* p, std::size_t size) const
    {
      std::uint64_t h;
      if (m_length >= XXH_STRIPE) {
        h = rotl64(m_v1, 1) + rotl64(m_v2, 7) + rotl64(m_v3, 12) + rotl64(m_v4, 18);
        h = xxhMergeRound(h, m_v1);
        h = xxhMergeRound(h, m_v2);
        h = xxhMergeRound(h, m_v3);
        h = xxhMergeRound(h, m_v4);
      } else {
        h = XXH_PRIME5;
      }
      h += m_length + size;

      for (; size >= 8; p += 8, size -= 8)
        h = rotl64(h ^ xxhRound(0, readLE64(p)), 27) * XXH_PRIME1 + XXH_PRIME4;
      if (size >= 4) {
        h = rotl64(h ^ (readLE32(p) * XXH_PRIME1), 23) * XXH_PRIME2 + XXH_PRIME3;
        p += 4;
        size -= 4;
      }
      for (; size > 0; ++p, --size)
        h = rotl64(h ^ (*p * XXH_PRIME5), 11) * XXH_PRIME1;

      h ^= h >> 33;
      h *= XXH_PRIME2;
      h ^= h >> 29;
      h *= XXH_PRIME3;
      h ^= h >> 32;
      return h;
    }

  private:
    std::uint64_t m_v1 = XXH_PRIME1 + XXH_PRIME2;
    std::uint64_t m_v2 = XXH_PRIME2;
    std::uint64_t m_v3 = 0;
    std::uint64_t m_v4 = 0 - XXH_PRIME1;
    std::uint64_t m_length = 0;
  };

  // A hash of the whole file, one streaming pass. It used to hash the size
  // plus the first and last 64 KiB, which was blind to a same-size edit anywhere
  // in the interior of a file larger than 128 KiB -- precisely the edit-and-rerun
  // workflow the fingerprint exists to catch (#1026).
  std::string fileContentFingerprint(const std::string& path)
  {
    // stdio, unbuffered: std::ifstream::read copied every byte through the
    // stream's own buffer, which cost more than the hash itself (#1167).
    const std::unique_ptr<std::FILE, int (*)(std::FILE*)> input(std::fopen(path.c_str(), "rb"), &std::fclose);
    if (!input) {
      throw std::runtime_error(fmt::format(FMT_STRING("Cannot open input file '{}'. Check that the path points to a file and that you have read permissions."), path));
    }
    std::setvbuf(input.get(), nullptr, _IONBF, 0);

    // A multiple of the stripe size, so only the last chunk leaves a remainder.
    constexpr std::size_t CHUNK_SIZE = 65536;
    static_assert(CHUNK_SIZE % XXH_STRIPE == 0, "chunks must hold whole stripes");
    std::vector<unsigned char> buffer(CHUNK_SIZE);
    Xxh64 hasher;
    std::size_t got = 0;
    std::size_t whole = 0;
    do {
      got = std::fread(buffer.data(), 1, CHUNK_SIZE, input.get());
      whole = got - got % XXH_STRIPE;
      hasher.stripes(buffer.data(), whole);
    } while (got == CHUNK_SIZE);
    // A short read is EOF *or* an I/O error, and the two must not be confused:
    // hashing only the readable prefix would publish a valid-looking identity
    // for content nobody read.
    if (std::ferror(input.get())) {
      throw std::runtime_error(fmt::format(FMT_STRING("Cannot read input file '{}': the read failed partway through."), path));
    }

    const auto hash = hasher.finish(buffer.data() + whole, got - whole);
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
  addCanonicalNumber(json, "intra_regularization_strength", config.intraRegularizationStrength);
  addCanonicalNumber(json, "inter_regularization_strength", config.interRegularizationStrength);
  json["multilayer_skip_absent_nodes"] = config.multilayerSkipAbsentNodes;
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
