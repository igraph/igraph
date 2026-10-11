/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "Config.h"
#include "InfomapError.h"
#include "OutputFormats.h"
#include "ParameterCatalog.h"
#include "ProgramInterface.h"
#include "SafeFile.h"
#include "../utils/FileURI.h"
#include "../utils/Log.h"
#include "../utils/convert.h"
#include "../utils/format.h"
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <sstream>
#include <vector>
#include <stdexcept>
#include <utility>

namespace infomap {

constexpr int FlowModel::undirected;
constexpr int FlowModel::directed;
constexpr int FlowModel::undirdir;
constexpr int FlowModel::outdirdir;
constexpr int FlowModel::rawdir;
constexpr int FlowModel::precomputed;

namespace {

  const std::vector<std::pair<std::string, FlowModel>>& flowModelMappings()
  {
    static const std::vector<std::pair<std::string, FlowModel>> mappings = {
      { "undirected", FlowModel::undirected },
      { "directed", FlowModel::directed },
      { "undirdir", FlowModel::undirdir },
      { "outdirdir", FlowModel::outdirdir },
      { "rawdir", FlowModel::rawdir },
      { "precomputed", FlowModel::precomputed },
    };
    return mappings;
  }

  void enableOutputFormat(Config& config, const OutputFormat& format)
  {
    switch (format.kind) {
    case OutputKind::Clu:
      config.printClu = true;
      break;
    case OutputKind::Tree:
      config.printTree = true;
      break;
    case OutputKind::FlowTree:
      config.printFlowTree = true;
      break;
    case OutputKind::Newick:
      config.printNewick = true;
      break;
    case OutputKind::Json:
      config.printJson = true;
      break;
    case OutputKind::Csv:
      config.printCsv = true;
      break;
    case OutputKind::PajekNetwork:
      config.printPajekNetwork = true;
      break;
    case OutputKind::StateNetwork:
      config.printStateNetwork = true;
      break;
    case OutputKind::FlowNetwork:
      config.printFlowNetwork = true;
      break;
    }
  }

  // Config invariants applied by adaptDefaults. These describe the coherent
  // state of a Config after fields have been set, regardless of whether the
  // fields were set by flag parsing or by direct mutation through bindings.

  void applyLibraryOutputDefaults(Config& config)
  {
    if (!config.isCLI && config.outDirectory.empty())
      config.noFileOutput = true;
  }

  void validateRequiredCliOutput(const Config& config)
  {
    if (!config.noFileOutput && config.outDirectory.empty() && config.isCLI) {
      throw std::runtime_error("Missing out_directory");
    }
  }

  // Whether the caller actually passed --num-trials/-N, as opposed to landing on
  // its default. parsedOptions records the long name whichever form was typed.
  bool numTrialsGivenExplicitly(const Config& config)
  {
    for (const auto& option : config.parsedOptions) {
      if (option.longName == "num-trials")
        return true;
    }
    return false;
  }

  void applyOptionInteractions(Config& config)
  {
    if (config.regularized) {
      config.recordedTeleportation = true;
    }

    if (config.noInfomap) {
      config.numTrials = 1;
    }

    // --converge reinterprets numTrials as a cap, and falls back to the default cap
    // when the caller did not give one. Both signals are needed to decide that:
    //   * the value, because a library caller mutates the struct and calls
    //     adaptDefaults(), leaving parsedOptions empty -- provenance alone would
    //     read a deliberate numTrials = 4 as "not given" and overwrite it;
    //   * the provenance, because numTrials == 1 was the sole sentinel before, and
    //     as the option's minimum it is also what an explicit `-N 1 --converge`
    //     produces, which was silently rewritten to the default cap. A cap of one
    //     makes --converge a no-op, but that is the caller's call to make.
    //
    // So the fix reaches the CLI only, and two paths keep the old behaviour:
    //   * a library caller who sets numTrials = 1 deliberately still gets the
    //     default cap, since nothing distinguishes that from leaving it untouched;
    //   * the Python/R bindings do parse an argument string, but their option
    //     renderer omits any value equal to its default, so `num_trials=1` with
    //     `converge=True` renders no --num-trials at all (measured: 6 trials on
    //     unbalanced_hierarchy either way).
    // That elision was assumed behaviour-preserving because the rendered defaults
    // come from this catalog -- this option is the counter-example, since
    // --converge makes "given as 1" mean something different from "not given".
    // Fixing it belongs in the generator (the residue recorded in #914).
    const bool trialCapWasChosen = numTrialsGivenExplicitly(config) || config.numTrials != 1;
    if (config.convergeTrials && !config.noInfomap && !trialCapWasChosen) {
      config.numTrials = Config::convergeDefaultMaxTrials;
    }
  }

  // The matchable state id is `physId << (ceil(log2(N)) + 1) | layerId`, so an N
  // this large makes the shift count 32 or more: undefined behaviour, which in
  // practice changed the flows even for single-digit node ids at exit 0. Bounding N
  // keeps the shift inside 32 bits; the companion check on the physical id lives in
  // StateNetwork::addStateNodeWithDeterministicId.
  void validateMatchableMultilayerIds(const Config& config)
  {
    // Local copy avoids odr-use of the constexpr member (no out-of-line definition).
    const unsigned int maxLargestLayerId = Config::maxMatchableMultilayerIds;
    if (config.matchableMultilayerIds > maxLargestLayerId) {
      throw std::runtime_error(fmt::format(FMT_STRING("--matchable-multilayer-ids {} is too large: the largest layer id must be at most {}, so the state id shift stays within 32 bits"), config.matchableMultilayerIds, maxLargestLayerId));
    }
  }

  void validateConvergeTrials(const Config& config)
  {
    if (!config.convergeTrials) {
      return;
    }
    if (config.parallelTrials) {
      throw std::runtime_error("--converge cannot be combined with --parallel-trials (auto-convergence needs each trial's result before deciding to continue)");
    }
    if (config.trialOffset > 0 || !config.trialResultsPath.empty()) {
      throw std::runtime_error("--converge cannot be combined with distributed sharding (--trial-offset / --trial-results); independent shards cannot coordinate a stop");
    }
  }

#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  void applyAndValidateLossyInteraction(Config& config)
  {
    if (!config.lossy)
      return;

    if (config.directed || (config.flowModelIsSet && config.flowModel != FlowModel::undirected))
      throw std::runtime_error("--lossy requires undirected flow");
    if (config.stateInput || config.multilayerInput || !config.additionalInput.empty())
      throw std::runtime_error("--lossy does not support memory or multilayer networks");
    if (config.haveMetaData())
      throw std::runtime_error("--lossy does not support meta data");
    if (config.recordedTeleportation || config.regularized)
      throw std::runtime_error("--lossy does not support recorded teleportation or regularization");
    if (config.variableMarkovTime || config.markovTime != 1.0)
      throw std::runtime_error("--lossy does not support Markov time scaling");

    config.twoLevel = true;
  }
#endif

  void applyFingerprintOnlyOutputInteraction(Config& config)
  {
    if (config.printConfigFingerprint) {
      config.noFileOutput = true;
    }
  }

  void applyRuntimeOutputInteractions(Config& config)
  {
    // Pretty is the only console rendering and is always on. --pretty/--no-pretty
    // are deprecated no-ops kept for backward compatibility; neutralize them here.
    config.prettyOutput = true;

    if (config.printAllTrials && config.numTrials < 2) {
      config.printAllTrials = false;
    }
  }

  void validateRunReportOutput(const Config& config)
  {
    if (config.memoryReport && config.timingJsonPath.empty()) {
      throw std::runtime_error("--memory-report requires --timing-json");
    }
    if (config.timingJsonPath == "-" && config.summaryJsonPath == "-") {
      throw std::runtime_error("--timing-json - and --summary-json - cannot both write to stdout");
    }
    if (config.timingJsonPath == "-" && config.runManifestPath == "-") {
      throw std::runtime_error("--timing-json - and --manifest-json - cannot both write to stdout");
    }
    if (config.summaryJsonPath == "-" && config.runManifestPath == "-") {
      throw std::runtime_error("--summary-json - and --manifest-json - cannot both write to stdout");
    }
    if (config.timingJsonPath == "-" && !config.silent) {
      throw std::runtime_error("--timing-json - requires --silent");
    }
    if (config.summaryJsonPath == "-" && !config.silent) {
      throw std::runtime_error("--summary-json - requires --silent");
    }
    if (config.runManifestPath == "-" && !config.silent) {
      throw std::runtime_error("--manifest-json - requires --silent");
    }
  }

  void normalizeOutputDirectory(Config& config)
  {
    if (!config.haveOutput() || config.outDirectory.empty())
      return;

    if (config.outDirectory.back() != '/')
      config.outDirectory.push_back('/');
  }

  void applyOutputNameDefault(Config& config)
  {
    if (config.outName.empty()) {
      config.outName = !config.networkFile.empty() ? FileURI(config.networkFile).getName() : "no-name";
    }
  }

  void applyThreadBudgetInteraction(Config& config)
  {
    if (config.numThreadsArg.empty() || config.numThreadsArg == "auto") {
      config.numThreads = 0;
      return;
    }
    const auto invalid = []() {
      throw std::runtime_error("--num-threads must be 'auto' or a positive integer");
    };
    long value = 0;
    try {
      std::size_t consumed = 0;
      value = std::stol(config.numThreadsArg, &consumed);
      if (consumed != config.numThreadsArg.size()) {
        invalid(); // trailing garbage like "4x"
      }
    } catch (const std::invalid_argument&) {
      invalid();
    } catch (const std::out_of_range&) {
      invalid();
    }
    if (value < 1 || static_cast<unsigned long>(value) > std::numeric_limits<unsigned int>::max()) {
      invalid();
    }
    config.numThreads = static_cast<unsigned int>(value);
  }

  // Lifecycle-only steps. These read staged parse state, touch the filesystem,
  // or mutate global state — they are not Config invariants and must not fire
  // when a library user calls adaptDefaults() on a mutated Config.

  void rejectDeprecatedAliases(const ParsedParameterSet& parsed)
  {
    if (parsed.deprecatedIncludeSelfLinks) {
      throw std::runtime_error("The --include-self-links flag is deprecated; self-links are included by default. Use --no-self-links to exclude them.");
    }
  }

  void applyOutputDirectory(Config& config, const ParsedParameterSet& parsed)
  {
    if (!parsed.optionalOutputDir.empty())
      config.outDirectory = parsed.optionalOutputDir[0];
  }

  void applyFlowModelSelection(Config& config, const ParsedParameterSet& parsed)
  {
    if (config.directed) {
      config.setFlowModel(FlowModel::directed);
      return;
    }

    if (parsed.flowModelArg.empty()) {
      return;
    }

    FlowModel flowModel = FlowModel::undirected;
    if (!parseFlowModel(parsed.flowModelArg, flowModel)) {
      throw std::runtime_error(fmt::format(FMT_STRING("Unrecognized flow model: '{}'"), parsed.flowModelArg));
    }
    config.setFlowModel(flowModel);
  }

  void validateOutputDirectory(const Config& config)
  {
    if (config.haveOutput()) {
      ensureDirectoryExists(config.outDirectory);
    }
    if (config.haveOutput() && !isDirectoryWritable(config.outDirectory))
      throw InfomapError(ExitCode::OutputError, fmt::format(FMT_STRING("Can't write to directory '{}'. Check that the directory exists and that you have write permissions."), config.outDirectory));
  }

  void initializeLogging(const Config& config)
  {
    Log::init(config.verbosity, config.silent, config.verboseNumberPrecision);
  }

  // What the flags say about silence, read before anything is parsed.
  //
  // Log visibility comes from the parsed Config, which does not exist yet while
  // the flags are being read -- so a rejected argument emitted the version banner
  // and the usage line through Log() at default verbosity
  // (ProgramInterface::exitWithError) even for a caller that asked for silence.
  // The Python API passes --silent on every call, so any rejected keyword put two
  // lines of English on the process's stdout, flushed at exit after everything the
  // Python program had written, against its documented quiet-by-default contract.
  //
  // Split on whitespace, the way tokenizeArgs does, and compared whole, so
  // "--silentish" does not match; there is no short form. The parser accepts the
  // --no-<option> negation and lets the last occurrence win (verified:
  // "--silent --no-silent" parses to silent = false and the reverse to true), so
  // this reads the last of either rather than the first --silent and agrees with
  // the Config the parse is about to produce. Defaults to Config::silent's own
  // default when neither appears.
  //
  // Not option-aware, and that is a deliberate limit rather than an oversight.
  // The parser decides option-versus-value by position -- parseLongOption consumes
  // the following token as an argument even when it starts with "--" -- so
  // `--out-name --silent` sets out_name to the literal "--silent" with silent
  // false, while this scan counts the token and silences the window. Telling the
  // two apart needs the per-option requireArgument knowledge that lives in
  // ProgramInterface, and mirroring it here would put the parser's rules in a
  // second place that can drift. The whole cost of being wrong is that a rejected
  // argument prints no banner in a run whose out-name is literally "--silent".
  bool flagsRequestSilence(const std::string& flags)
  {
    bool silent = false;
    std::istringstream stream(flags);
    std::string token;
    while (stream >> token) {
      if (token == "--silent")
        silent = true;
      else if (token == "--no-silent")
        silent = false;
    }
    return silent;
  }

  // Sets Log's silence for the window between the first token and
  // initializeLogging, from the current flags alone, and restores the previous
  // value on the way out.
  //
  // Unconditional on purpose. Engaging only for a silent request would leave a
  // non-silent parse reading whatever the last run in this process happened to
  // set -- a successful silent run leaves Log silent, so a later non-silent
  // rejection printed nothing. Parse-time output has to follow the flags in front
  // of it, not process history. Restoring also means a library caller that catches
  // the error is not left with a muted engine.
  class ScopedParseSilence {
  public:
    explicit ScopedParseSilence(bool silent)
        : m_previous(Log::isSilent())
    {
      Log::setSilent(silent);
    }

    ~ScopedParseSilence()
    {
      Log::setSilent(m_previous);
    }

  private:
    bool m_previous;
  };

  void buildConfigFromFlags(Config& config, const std::string& flags, bool isCLI)
  {
    config.parsedString = flags;

    ProgramInterface api("Infomap",
                         "Implementation of the Infomap clustering algorithm based on the Map Equation (see www.mapequation.org)",
                         INFOMAP_VERSION);

    api.setGroups({ "Input", "Algorithm", "Accuracy", "Output" });
    api.setJsonParametersProvider(parameterCatalogJson);

    ParsedParameterSet staging;
    registerCatalogWithProgramInterface(api, { config, staging }, isCLI);

    // Scoped to end before initializeLogging: the guard restores the previous
    // silence on destruction, so leaving it open would undo the real setting.
    {
      const ScopedParseSilence parseSilence(flagsRequestSilence(flags));

      api.parseArgs(flags);
      config.parsedOptions = api.getUsedOptionArguments();

      rejectDeprecatedAliases(staging);
      applyOutputDirectory(config, staging);
      applyFlowModelSelection(config, staging);

      config.adaptDefaults();

      validateOutputDirectory(config);
    }

    initializeLogging(config);
  }

} // namespace

const std::vector<std::string>& flowModelNames()
{
  static const std::vector<std::string> names = [] {
    const auto& mappings = flowModelMappings();
    std::vector<std::string> values;
    values.reserve(mappings.size());
    std::transform(mappings.begin(), mappings.end(), std::back_inserter(values), [](const auto& mapping) { return mapping.first; });
    return values;
  }();
  return names;
}

bool parseFlowModel(const std::string& name, FlowModel& flowModel)
{
  for (const auto& mapping : flowModelMappings()) {
    if (mapping.first == name) {
      flowModel = mapping.second;
      return true;
    }
  }
  return false;
}

const char* flowModelToString(const FlowModel& flowModel)
{
  for (const auto& mapping : flowModelMappings()) {
    if (mapping.second == flowModel) {
      return mapping.first.c_str();
    }
  }
  return "undirected";
}

Config::Config(const std::string& flags, bool isCLI) : isCLI(isCLI)
{
  buildConfigFromFlags(*this, flags, isCLI);
}

void Config::adaptDefaults()
{
  auto outputs = io::split(outputFormats, ',');
  for (std::string& o : outputs) {
    const auto* format = findOutputFormat(o);
    if (format == nullptr) {
      throw std::runtime_error(fmt::format(FMT_STRING("Unrecognized output format: '{}'."), o));
    }
    enableOutputFormat(*this, *format);
  }

  // Of no output format specified, use tree as default (if not used as a library).
  if (isCLI && !haveModularResultOutput()) {
    printTree = true;
  }

  // Cross-field invariants. These run whether construction was via flag parsing
  // or library mutation followed by adaptDefaults().
  applyFingerprintOnlyOutputInteraction(*this);
  applyLibraryOutputDefaults(*this);
  validateRequiredCliOutput(*this);
  applyOptionInteractions(*this);
  validateConvergeTrials(*this);
  validateMatchableMultilayerIds(*this);
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  applyAndValidateLossyInteraction(*this);
#endif
  applyThreadBudgetInteraction(*this);
  validateRunReportOutput(*this);
  normalizeOutputDirectory(*this);
  applyOutputNameDefault(*this);
  applyRuntimeOutputInteractions(*this);
}

std::ostream& operator<<(std::ostream& out, FlowModel f)
{
  return out << flowModelToString(f);
}

} // namespace infomap
