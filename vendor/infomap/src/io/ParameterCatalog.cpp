/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include <nlohmann/json.hpp>

#include "ParameterCatalog.h"
#include "Config.h"
#include "OutputFormats.h"
#include "ProgramInterface.h"
#include "../utils/convert.h"
#include "../utils/format.h"

#include <stdexcept>
#include <utility>

namespace infomap {

namespace {

  template <typename T>
  void registerOptionForTarget(ProgramInterface& api, T& target, const ParameterSpec& parameter);

  void registerOptionForTarget(ProgramInterface& api, bool& target, const ParameterSpec& parameter);

  void registerIncrementalOptionForTarget(ProgramInterface& api, unsigned int& target, const ParameterSpec& parameter);

  class SpecBuilder {
  public:
    SpecBuilder() = default;

    SpecBuilder& shortName(char value)
    {
      parameter_.shortName = value;
      return *this;
    }

    SpecBuilder& longName(std::string value)
    {
      parameter_.longName = std::move(value);
      return *this;
    }

    SpecBuilder& description(std::string value)
    {
      parameter_.description = std::move(value);
      return *this;
    }

    SpecBuilder& group(std::string value)
    {
      parameter_.group = std::move(value);
      return *this;
    }

    SpecBuilder& argument(std::string value)
    {
      parameter_.argumentName = std::move(value);
      parameter_.requireArgument = true;
      return *this;
    }

    SpecBuilder& defaultValue(std::string value)
    {
      parameter_.defaultValue = std::move(value);
      return *this;
    }

    SpecBuilder& advanced()
    {
      parameter_.isAdvanced = true;
      return *this;
    }

    SpecBuilder& hidden()
    {
      parameter_.hidden = true;
      return *this;
    }

    SpecBuilder& incremental()
    {
      parameter_.incrementalArgument = true;
      return *this;
    }

    SpecBuilder& choices(std::vector<std::string> values)
    {
      parameter_.choices = std::move(values);
      return *this;
    }

    SpecBuilder& min(std::string value)
    {
      parameter_.minValue = std::move(value);
      return *this;
    }

    SpecBuilder& range(std::string minValue, std::string maxValue)
    {
      parameter_.minValue = std::move(minValue);
      parameter_.maxValue = std::move(maxValue);
      return *this;
    }

    SpecBuilder& renderPolicy(std::string value)
    {
      parameter_.renderPolicy = std::move(value);
      return *this;
    }

    SpecBuilder& cliOnly()
    {
      parameter_.cliOnly = true;
      return *this;
    }

    SpecBuilder& libraryOnly()
    {
      parameter_.libraryOnly = true;
      return *this;
    }

    SpecBuilder& hideFromJson()
    {
      parameter_.includeInJson = false;
      return *this;
    }

    template <typename T>
    SpecBuilder& configTarget(T Config::* member)
    {
      parameter_.registrar = [member](ProgramInterface& api, ConfigParameterTargets& targets, const ParameterSpec& parameter) {
        registerOptionForTarget(api, targets.config.*member, parameter);
      };
      return *this;
    }

    SpecBuilder& incrementalTarget(unsigned int Config::* member)
    {
      parameter_.registrar = [member](ProgramInterface& api, ConfigParameterTargets& targets, const ParameterSpec& parameter) {
        registerIncrementalOptionForTarget(api, targets.config.*member, parameter);
      };
      return *this;
    }

    SpecBuilder& flowModelTarget()
    {
      parameter_.registrar = [](ProgramInterface& api, ConfigParameterTargets& targets, const ParameterSpec& parameter) {
        registerOptionForTarget(api, targets.parsed.flowModelArg, parameter);
      };
      return *this;
    }

    SpecBuilder& deprecatedIncludeSelfLinksTarget()
    {
      parameter_.registrar = [](ProgramInterface& api, ConfigParameterTargets& targets, const ParameterSpec& parameter) {
        registerOptionForTarget(api, targets.parsed.deprecatedIncludeSelfLinks, parameter);
      };
      return *this;
    }

    SpecBuilder& networkFileArgument()
    {
      parameter_.registrar = [](ProgramInterface& api, ConfigParameterTargets& targets, const ParameterSpec& parameter) {
        api.addNonOptionArgument(targets.config.networkFile, "network_file", parameter.description, parameter.group);
      };
      return *this;
    }

    SpecBuilder& outputDirectoryArgument()
    {
      parameter_.registrar = [](ProgramInterface& api, ConfigParameterTargets& targets, const ParameterSpec& parameter) {
        api.addOptionalNonOptionArguments(targets.parsed.optionalOutputDir, "out_directory", parameter.description, parameter.group);
      };
      return *this;
    }

    operator ParameterSpec() const { return parameter_; }

  private:
    ParameterSpec parameter_;
  };

  SpecBuilder param()
  {
    return SpecBuilder();
  }

  using Json = nlohmann::ordered_json;

  std::string stringDefault(const ParameterSpec& parameter)
  {
    return parameter.requireArgument ? parameter.defaultValue : "";
  }

  Option& addConfiguredOption(Option& option, const ParameterSpec& parameter)
  {
    if (parameter.hidden)
      option.setHidden(true);
    if (!parameter.choices.empty())
      option.setChoices(parameter.choices);
    return option;
  }

  Option& addFlagOption(ProgramInterface& api, bool& target, const ParameterSpec& parameter)
  {
    if (parameter.shortName == '\0') {
      return api.addOptionArgument(target, parameter.longName, parameter.description, parameter.group, parameter.isAdvanced);
    }
    return api.addOptionArgument(target, parameter.shortName, parameter.longName, parameter.description, parameter.group, parameter.isAdvanced);
  }

  template <typename T>
  Option& addValueOption(ProgramInterface& api, T& target, const ParameterSpec& parameter)
  {
    if (parameter.shortName == '\0') {
      return api.addOptionArgument(target, parameter.longName, parameter.description, parameter.argumentName, parameter.group, parameter.isAdvanced);
    }
    return api.addOptionArgument(target, parameter.shortName, parameter.longName, parameter.description, parameter.argumentName, parameter.group, parameter.isAdvanced);
  }

  template <typename T>
  T numericConstraint(const ParameterSpec& parameter, const std::string& value, const std::string& name)
  {
    T constraint;
    if (value.empty() || !io::stringToValue(value, constraint)) {
      throw std::runtime_error(fmt::format(FMT_STRING("Invalid {} for parameter --{}."), name, parameter.longName));
    }
    return constraint;
  }

  template <typename T>
  Option& addLowerBoundOption(ProgramInterface& api, T& target, const ParameterSpec& parameter)
  {
    const auto minValue = numericConstraint<T>(parameter, parameter.minValue, "minValue");
    if (parameter.shortName == '\0') {
      return api.addOptionArgument(target, parameter.longName, parameter.description, parameter.argumentName, parameter.group, minValue, parameter.isAdvanced);
    }
    return api.addOptionArgument(target, parameter.shortName, parameter.longName, parameter.description, parameter.argumentName, parameter.group, minValue, parameter.isAdvanced);
  }

  template <typename T>
  Option& addRangeOption(ProgramInterface& api, T& target, const ParameterSpec& parameter)
  {
    const auto minValue = numericConstraint<T>(parameter, parameter.minValue, "minValue");
    const auto maxValue = numericConstraint<T>(parameter, parameter.maxValue, "maxValue");
    if (parameter.shortName == '\0') {
      return api.addOptionArgument(target, parameter.longName, parameter.description, parameter.argumentName, parameter.group, minValue, maxValue, parameter.isAdvanced);
    }
    return api.addOptionArgument(target, parameter.shortName, parameter.longName, parameter.description, parameter.argumentName, parameter.group, minValue, maxValue, parameter.isAdvanced);
  }

  template <typename T>
  Option& addValueOptionForSpec(ProgramInterface& api, T& target, const ParameterSpec& parameter)
  {
    if (!parameter.minValue.empty() && !parameter.maxValue.empty()) {
      return addRangeOption(api, target, parameter);
    }
    if (!parameter.minValue.empty()) {
      return addLowerBoundOption(api, target, parameter);
    }
    return addValueOption(api, target, parameter);
  }

  void registerOptionForTarget(ProgramInterface& api, bool& target, const ParameterSpec& parameter)
  {
    addConfiguredOption(addFlagOption(api, target, parameter), parameter);
  }

  template <typename T>
  void registerOptionForTarget(ProgramInterface& api, T& target, const ParameterSpec& parameter)
  {
    addConfiguredOption(addValueOptionForSpec(api, target, parameter), parameter);
  }

  void registerIncrementalOptionForTarget(ProgramInterface& api, unsigned int& target, const ParameterSpec& parameter)
  {
    addConfiguredOption(api.addIncrementalOptionArgument(target, parameter.shortName, parameter.longName, parameter.description, parameter.group, parameter.isAdvanced), parameter);
  }

} // namespace

const std::vector<ParameterSpec>& parameterCatalog()
{
  static const std::vector<ParameterSpec> parameters = {
    param()
        .shortName('h')
        .longName("help")
        .description("Prints this help message. Use -hh to show advanced options.")
        .group("About")
        .incremental()
        .defaultValue("false")
        .renderPolicy("repeated_short"),
    param()
        .shortName('V')
        .longName("version")
        .description("Display program version information.")
        .group("About")
        .defaultValue("false"),
    param()
        .longName("completion")
        .description("Print shell completion script for bash or zsh.")
        .argument(ArgType::option)
        .group("About")
        .choices({ "bash", "zsh" })
        .hideFromJson(),
    param()
        .longName("network_file")
        .description("Network file to read. Infomap assumes link-list format unless the file starts with a Pajek heading.")
        .group("Input")
        .cliOnly()
        .hideFromJson()
        .networkFileArgument(),
    param()
        .longName("input")
        .description("Network file to read. Infomap assumes link-list format unless the file starts with a Pajek heading.")
        .argument(ArgType::path)
        .group("Input")
        .libraryOnly()
        .hideFromJson()
        .configTarget(&Config::networkFile),
    param()
        .longName("skip-adjust-bipartite-flow")
        .description("Keep flow on bipartite nodes instead of distributing it to primary nodes.")
        .group("Input")
        .advanced()
        .configTarget(&Config::skipAdjustBipartiteFlow),
    param()
        .longName("bipartite-teleportation")
        .description("Use bipartite teleportation instead of the default two-step unipartite teleportation.")
        .group("Input")
        .advanced()
        .configTarget(&Config::bipartiteTeleportation),
    param()
        .longName("weight-threshold")
        .description("Ignore input links with weight below this threshold.")
        .argument(ArgType::number)
        .group("Input")
        .advanced()
        .defaultValue("0")
        .min("0")
        .configTarget(&Config::weightThreshold),
    param()
        .shortName('k')
        .longName("include-self-links")
        .description("DEPRECATED. Self-links are included by default; use --no-self-links to exclude them.")
        .group("Input")
        .advanced()
        .hidden()
        .deprecatedIncludeSelfLinksTarget(),
    param()
        .longName("no-self-links")
        .description("Exclude self-links from the input network.")
        .group("Input")
        .advanced()
        .defaultValue("false")
        .renderPolicy("deprecated_alias_target")
        .configTarget(&Config::noSelfLinks),
    param()
        .longName("node-limit")
        .description("Read only nodes up to this node id and ignore links connected to higher node ids.")
        .argument(ArgType::integer)
        .group("Input")
        .advanced()
        .defaultValue("0")
        .min("1")
        .configTarget(&Config::nodeLimit),
    param()
        .longName("matchable-multilayer-ids")
        .description("Construct state ids from node ids and layer ids that stay comparable across networks. Set at least to the largest layer id among networks to match.")
        .argument(ArgType::integer)
        .group("Input")
        .advanced()
        .defaultValue("0")
        .min("1")
        .configTarget(&Config::matchableMultilayerIds),
    param()
        .shortName('c')
        .longName("cluster-data")
        .description("Read an initial partition from a clu file or a hierarchy from a tree/ftree file. Tree input may use physical or state nodes for higher-order networks.")
        .argument(ArgType::path)
        .group("Input")
        .configTarget(&Config::clusterDataFile),
    param()
        .longName("assign-to-neighbouring-module")
        .description("With --cluster-data, assign nodes missing module ids to a neighboring node's module when possible.")
        .group("Input")
        .advanced()
        .configTarget(&Config::assignToNeighbouringModule),
    param()
        .longName("meta-data")
        .description("Read metadata to encode from a clu-format file.")
        .argument(ArgType::path)
        .group("Input")
        .advanced()
        .configTarget(&Config::metaDataFile),
    param()
        .longName("meta-data-rate")
        .description("With --meta-data, set the metadata encoding rate. The default encodes metadata at each step.")
        .argument(ArgType::number)
        .group("Input")
        .advanced()
        .defaultValue("1")
        .min("0")
        .configTarget(&Config::metaDataRate),
    param()
        .longName("meta-data-unweighted")
        .description("With --meta-data, encode metadata without weighting by node flow.")
        .group("Input")
        .advanced()
        .configTarget(&Config::unweightedMetaData),
    param()
        .longName("no-infomap")
        .description("Skip optimization. Use this to calculate codelength for --cluster-data or to print non-modular statistics.")
        .group("Input")
        .configTarget(&Config::noInfomap),
    param()
        .longName("out-name")
        .description("Base name for output files, for example [out_directory]/[out-name].tree.")
        .argument(ArgType::string)
        .group("Output")
        .advanced()
        .configTarget(&Config::outName),
    param()
        .shortName('0')
        .longName("no-file-output")
        .description("Do not write output files.")
        .group("Output")
        .advanced()
        .configTarget(&Config::noFileOutput),
    param()
        .longName("tree")
        .description("Write the modular hierarchy to a tree file. Enabled by default when no other output format is selected.")
        .group("Output")
        .configTarget(&Config::printTree),
    param()
        .longName("ftree")
        .description("Write the modular hierarchy and aggregated links between nested modules to an ftree file. Used by Network Navigator.")
        .group("Output")
        .configTarget(&Config::printFlowTree),
    param()
        .longName("clu")
        .description("Write top-level module ids for each node to a clu file.")
        .group("Output")
        .configTarget(&Config::printClu),
    param()
        .longName("clu-level")
        .description("With --clu or --output clu, write module ids at this depth from the root. Use -1 for bottom-level modules.")
        .argument(ArgType::integer)
        .group("Output")
        .advanced()
        .defaultValue("1")
        .min("-1")
        .configTarget(&Config::cluLevel),
    param()
        .shortName('o')
        .longName("output")
        .description("Write selected output formats as a comma-separated list without spaces, e.g. -o clu,tree,ftree. Options: clu, tree, ftree, newick, json, csv, network, states, flow.")
        .argument(ArgType::list)
        .group("Output")
        .advanced()
        .choices(outputFormatNames())
        .renderPolicy("comma_list")
        .configTarget(&Config::outputFormats),
    param()
        .longName("hide-bipartite-nodes")
        .description("Hide bipartite nodes in output by projecting the solution to primary nodes.")
        .group("Output")
        .advanced()
        .configTarget(&Config::hideBipartiteNodes),
    param()
        .longName("print-all-trials")
        .description("Write each trial to separate output files. Has effect only when --num-trials is greater than 1.")
        .group("Output")
        .advanced()
        .configTarget(&Config::printAllTrials),
    param()
        .longName("no-overwrite")
        .description("Fail with an output error if any target output file already exists. By default existing files are replaced.")
        .group("Output")
        .advanced()
        .configTarget(&Config::noOverwriteOutput),
    param()
        .longName("print-config-fingerprint")
        .description("Print the canonical configuration fingerprint and exit.")
        .group("Output")
        .advanced()
        .configTarget(&Config::printConfigFingerprint),
    param()
        .longName("timing-json")
        .description("Write machine-readable run timing JSON to this path. Use - for stdout.")
        .argument(ArgType::path)
        .group("Output")
        .advanced()
        .configTarget(&Config::timingJsonPath),
    param()
        .longName("summary-json")
        .description("Write machine-readable final run summary JSON to this path. Use - for stdout.")
        .argument(ArgType::path)
        .group("Output")
        .advanced()
        .configTarget(&Config::summaryJsonPath),
    param()
        .longName("manifest-json")
        .description("Write a machine-readable run manifest JSON to this path. Use - for stdout.")
        .argument(ArgType::path)
        .group("Output")
        .advanced()
        .configTarget(&Config::runManifestPath),
    param()
        .longName("memory-report")
        .description("Include peak RSS and best-effort bytes per node/link estimates in timing JSON. Requires --timing-json.")
        .group("Output")
        .advanced()
        .configTarget(&Config::memoryReport),
    param()
        .shortName('2')
        .longName("two-level")
        .description("Optimize a two-level partition instead of the default multi-level hierarchy.")
        .group("Algorithm")
        .configTarget(&Config::twoLevel),
    param()
        .shortName('f')
        .longName("flow-model")
        .description("Choose how Infomap derives flow from the input links. Options: undirected, directed, undirdir, outdirdir, rawdir, precomputed.")
        .argument(ArgType::option)
        .group("Algorithm")
        .choices(flowModelNames())
        .flowModelTarget(),
    param()
        .shortName('d')
        .longName("directed")
        .description("Treat input links as directed. Shorthand for --flow-model directed.")
        .group("Algorithm")
        .defaultValue("false")
        .renderPolicy("directed_alias")
        .configTarget(&Config::directed),
    param()
        .shortName('e')
        .longName("recorded-teleportation")
        .description("When teleportation is used to calculate flow, also record teleportation steps in the codelength.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::recordedTeleportation),
    param()
        .longName("use-node-weights-as-flow")
        .description("Use node weights from the API or Pajek node records as normalized node flow.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::useNodeWeightsAsFlow),
    param()
        .longName("to-nodes")
        .description("Teleport to nodes instead of links. Uses uniform node weights unless node weights are provided.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::teleportToNodes),
    param()
        .shortName('p')
        .longName("teleportation-probability")
        .description("Set the probability of teleporting to a random node or link when calculating flow.")
        .argument(ArgType::probability)
        .group("Algorithm")
        .advanced()
        .defaultValue("0.15")
        .range("0", "1")
        .configTarget(&Config::teleportationProbability),
    param()
        .longName("max-flow-iterations")
        .description("Limit the power iteration used to calculate flow (directed and regularized flow models) to this many iterations.")
        .argument(ArgType::integer)
        .group("Algorithm")
        .advanced()
        .defaultValue("400")
        .min("1")
        .configTarget(&Config::maxFlowIterations),
    param()
        .longName("min-flow-iterations")
        .description("Require at least this many power iterations before the flow calculation can converge, even if --flow-tolerance is already met.")
        .argument(ArgType::integer)
        .group("Algorithm")
        .advanced()
        .defaultValue("50")
        .min("0")
        .configTarget(&Config::minFlowIterations),
    param()
        .longName("flow-tolerance")
        .description("Convergence tolerance for the power iteration used to calculate flow. Iteration stops once the per-iteration change in flow drops to or below this value, after --min-flow-iterations have run.")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1e-15")
        .min("0")
        .configTarget(&Config::flowTolerance),
    param()
        .longName("regularized")
        .description("Add a fully connected Bayesian prior network to reduce overfitting to missing links. Activates --recorded-teleportation.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::regularized),
    param()
        .longName("regularization-strength")
        .description("Scale the relative strength of the Bayesian prior network used by --regularized.")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1")
        .min("0")
        .configTarget(&Config::regularizationStrength),
    param()
        .longName("entropy-corrected")
        .description("Correct for negative entropy bias in small samples, especially solutions with many modules.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::entropyBiasCorrection),
    param()
        .longName("entropy-correction-strength")
        .description("Scale the default correction used by --entropy-corrected.")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1")
        .min("0")
        .configTarget(&Config::entropyBiasCorrectionMultiplier),
    param()
        .longName("markov-time")
        .description("Scale link flow to change the cost of moving between modules. Higher values result in fewer modules.")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1")
        .min("0")
        .configTarget(&Config::markovTime),
    param()
        .longName("variable-markov-time")
        .description("Vary Markov time locally to reduce overpartitioning in sparse areas while keeping higher resolution in dense areas.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::variableMarkovTime),
    param()
        .longName("variable-markov-damping")
        .description("With --variable-markov-time, set damping between local effective degree (0) and local entropy (1).")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1")
        .range("0", "1")
        .configTarget(&Config::variableMarkovTimeDamping),
    param()
        .longName("variable-markov-min-scale")
        .description("With --variable-markov-time, set the minimum local scale for zero-entropy nodes. Local Markov time is max scale divided by local scale.")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1")
        .configTarget(&Config::variableMarkovTimeMinLocalScale),
    param()
        .longName("preferred-number-of-modules")
        .description("Penalize solutions by how far their number of modules differs from this value.")
        .argument(ArgType::integer)
        .group("Algorithm")
        .advanced()
        .defaultValue("0")
        .min("1")
        .configTarget(&Config::preferredNumberOfModules),
    param()
        .longName("preferred-number-of-levels")
        .description("Soft preference for the depth of the hierarchy. Steering to a shallower depth is reliable at a small codelength cost; deeper is best-effort, bounded by what the optimizer proposes. No-op with --two-level or strength 0.")
        .argument(ArgType::integer)
        .group("Algorithm")
        .advanced()
        .defaultValue("0")
        .min("1")
        .configTarget(&Config::preferredNumberOfLevels),
    param()
        .longName("preferred-number-of-levels-strength")
        .description("Scale the strength of --preferred-number-of-levels. 0 disables the preference; larger values increase the cost of deviating from the preferred depth.")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1")
        .min("0")
        .configTarget(&Config::preferredNumberOfLevelsStrength),
    param()
        .longName("multilayer-relax-rate")
        .description("Set the probability of relaxing from a state node to neighboring layers instead of staying in the current layer.")
        .argument(ArgType::probability)
        .group("Algorithm")
        .advanced()
        .defaultValue("0.15")
        .range("0", "1")
        .configTarget(&Config::multilayerRelaxRate),
    param()
        .longName("multilayer-relax-limit")
        .description("Limit relaxation to this many neighboring layer ids in each direction. Use a negative value to allow relaxation to any layer.")
        .argument(ArgType::integer)
        .group("Algorithm")
        .advanced()
        .defaultValue("-1")
        .configTarget(&Config::multilayerRelaxLimit),
    param()
        .longName("multilayer-relax-limit-up")
        .description("Limit relaxation upward to this many higher neighboring layer ids. Use a negative value to allow relaxation to any higher layer.")
        .argument(ArgType::integer)
        .group("Algorithm")
        .advanced()
        .defaultValue("-1")
        .configTarget(&Config::multilayerRelaxLimitUp),
    param()
        .longName("multilayer-relax-limit-down")
        .description("Limit relaxation downward to this many lower neighboring layer ids. Use a negative value to allow relaxation to any lower layer.")
        .argument(ArgType::integer)
        .group("Algorithm")
        .advanced()
        .defaultValue("-1")
        .configTarget(&Config::multilayerRelaxLimitDown),
    param()
        .longName("multilayer-relax-by-jsd")
        .description("Weight multilayer relaxation by out-link similarity measured with Jensen-Shannon divergence.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::multilayerRelaxByJensenShannonDivergence),
    param()
        .longName("multilayer-relax-to-self")
        .description("On relaxation, link a state node to its own physical node in the target layer instead of spreading to its out-neighbors. Builds a smaller state network with the same flow as the default.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::multilayerRelaxToSelf),
    param()
        .shortName('s')
        .longName("seed")
        .description("Set the random number generator seed for reproducible results.")
        .argument(ArgType::integer)
        .group("Accuracy")
        .defaultValue("123")
        .min("1")
        .configTarget(&Config::seedToRandomNumberGenerator),
    param()
        .shortName('N')
        .longName("num-trials")
        .description("Run this many independent trials and keep the best solution.")
        .argument(ArgType::integer)
        .group("Accuracy")
        .defaultValue("1")
        .min("1")
        .configTarget(&Config::numTrials),
    param()
        .shortName('M')
        .longName("core-loop-limit")
        .description("Limit how many core loops try to move each node to the best module.")
        .argument(ArgType::integer)
        .group("Accuracy")
        .advanced()
        .defaultValue("10")
        .min("1")
        .configTarget(&Config::coreLoopLimit),
    param()
        .shortName('L')
        .longName("core-level-limit")
        .description("Limit how many times core loops are reapplied to the aggregated modular network to find larger structures. 0 means no limit.")
        .argument(ArgType::integer)
        .group("Accuracy")
        .advanced()
        .defaultValue("0")
        .min("0")
        .configTarget(&Config::levelAggregationLimit),
    param()
        .shortName('T')
        .longName("tune-iteration-limit")
        .description("Limit the main iterations in the two-level partition algorithm. 0 means no limit.")
        .argument(ArgType::integer)
        .group("Accuracy")
        .advanced()
        .defaultValue("0")
        .min("0")
        .configTarget(&Config::tuneIterationLimit),
    param()
        .longName("core-loop-codelength-threshold")
        .description("Require at least this codelength improvement to accept a new solution in a core loop.")
        .argument(ArgType::number)
        .group("Accuracy")
        .advanced()
        .defaultValue("1e-10")
        .min("0")
        .configTarget(&Config::minimumCodelengthImprovement),
    param()
        .longName("tune-iteration-relative-threshold")
        .description("Require each tune iteration to improve codelength by this fraction of the initial two-level codelength.")
        .argument(ArgType::number)
        .group("Accuracy")
        .advanced()
        .defaultValue("1e-05")
        .min("0")
        .configTarget(&Config::minimumRelativeTuneIterationImprovement),
    param()
        .shortName('F')
        .longName("fast-hierarchical-solution")
        .description("Find top modules quickly. Use -FF to keep all fast levels. Use -FFF to skip recursive refinement.")
        .group("Accuracy")
        .advanced()
        .incremental()
        .defaultValue("false")
        .renderPolicy("repeated_short")
        .incrementalTarget(&Config::fastHierarchicalSolution),
    param()
        .longName("inner-parallelization")
        .description("Experimental: use batched parallel node moves for coarse optimization. Performance gains are workload-dependent, often require a relaxed core-loop-codelength-threshold and low tune-iteration-limit, and may produce a different partition than serial optimization.")
        .group("Accuracy")
        .advanced()
        .configTarget(&Config::innerParallelization),
    param()
        .longName("parallel-trials")
        .description("Run independent trials in parallel with OpenMP. --num-trials remains the total number of trials; the number of parallel workers follows the OpenMP thread count (e.g. OMP_NUM_THREADS), clamped to --num-trials. Peak memory scales with the worker count. Nested OpenMP and --inner-parallelization are disabled inside workers.")
        .group("Accuracy")
        .advanced()
        .configTarget(&Config::parallelTrials),
    param()
        .longName("converge")
        .description("Treat the trial count as a cap and stop early once the best codelength has plateaued (no meaningful improvement over several consecutive trials). Runs trials serially; cannot be combined with parallel trials or distributed sharding. With no explicit trial count, a default cap is used.")
        .group("Accuracy")
        .configTarget(&Config::convergeTrials),
    param()
        .longName("num-threads")
        .description("Effective thread budget: 'auto' (resolve from --num-threads > INFOMAP_NUM_THREADS > SLURM_CPUS_PER_TASK > OMP_NUM_THREADS > cpuset > hardware), or a positive integer. 1 forces fully serial. Governs the recursive partition, parallel trials, and inner parallelization.")
        .argument(ArgType::string)
        .group("Accuracy")
        .advanced()
        .defaultValue("auto")
        .configTarget(&Config::numThreadsArg),
    param()
        .longName("threads")
        .description("Alias for --num-threads.")
        .argument(ArgType::string)
        .group("Accuracy")
        .advanced()
        .defaultValue("auto")
        .configTarget(&Config::numThreadsArg),
    param()
        .longName("trial-offset")
        .description("Global index of the first trial this process runs; trial i uses seed = base_seed + (trial_offset + i). Default 0 (single-process behavior).")
        .argument(ArgType::integer)
        .group("Output")
        .advanced()
        .defaultValue("0")
        .min("0")
        .configTarget(&Config::trialOffset),
    param()
        .longName("trial-results")
        .description("Write this shard's per-trial results (codelengths, seeds, best-tree reference, fingerprints) as JSON to this path, for deterministic merging of distributed shard runs into a final solution.")
        .argument(ArgType::path)
        .group("Output")
        .advanced()
        .configTarget(&Config::trialResultsPath),
    param()
        .longName("no-final-output")
        .description("Skip writing this process's aggregate best result. Per-trial outputs and --trial-results are still written.")
        .group("Output")
        .advanced()
        .configTarget(&Config::noFinalOutput),
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
    param()
        .longName("lossy")
        .description("Use the lossy map equation: modules whose naming overhead exceeds lambda times their entropy rate become noise modules with one shared visit codeword. Implies --two-level. Requires undirected flow.")
        .group("Algorithm")
        .advanced()
        .configTarget(&Config::lossy),
    param()
        .longName("lambda")
        .description("With --lossy, the distortion price: one bit of forgone dynamics costs lambda transmitted bits. Low values lump more structure into noise modules; high values approach the standard map equation.")
        .argument(ArgType::number)
        .group("Algorithm")
        .advanced()
        .defaultValue("1")
        .min("0")
        .configTarget(&Config::lossyLambda),
#endif
#if INFOMAP_FEATURE_TEST_FEATURE
    param()
        .longName("test-feature")
        .description("Enable the internal feature flag canary.")
        .group("Accuracy")
        .advanced()
        .configTarget(&Config::testFeature),
#endif
    param()
        .longName("prefer-modular-solution")
        .description("Prefer a modular solution even when one module gives a lower codelength.")
        .group("Accuracy")
        .advanced()
        .configTarget(&Config::preferModularSolution),
    param()
        .longName("num-random-moves")
        .description("Try this many random moves in each core loop to merge weakly connected nodes.")
        .argument(ArgType::integer)
        .group("Accuracy")
        .advanced()
        .defaultValue("5")
        .min("0")
        .configTarget(&Config::numRandomMoves),
    param()
        .longName("max-degree-for-random-moves")
        .description("Try random moves only for nodes with degree at most this value.")
        .argument(ArgType::integer)
        .group("Accuracy")
        .advanced()
        .defaultValue("2")
        .min("0")
        .configTarget(&Config::maxDegreeForRandomMoves),
    param()
        .longName("out_directory")
        .description("Directory where output files are written.")
        .group("Output")
        .hideFromJson()
        .outputDirectoryArgument(),
    param()
        .shortName('v')
        .longName("verbose")
        .description("Increase console verbosity. Add more v flags to increase verbosity up to -vvv.")
        .group("Output")
        .incremental()
        .defaultValue("false")
        .renderPolicy("repeated_short")
        .incrementalTarget(&Config::verbosity),
    param()
        .longName("silent")
        .description("Suppress console output.")
        .group("Output")
        .configTarget(&Config::silent),
    param()
        .longName("pretty")
        .description("Deprecated no-op: pretty console output is always on. Kept for backward compatibility.")
        .group("Output")
        .hidden()
        .configTarget(&Config::prettyOutput),
  };
  return parameters;
}

void registerCatalogWithProgramInterface(ProgramInterface& api, ConfigParameterTargets targets, bool isCli)
{
  for (const auto& parameter : parameterCatalog()) {
    if (!parameter.registrar) {
      continue;
    }
    if (isCli && parameter.libraryOnly) {
      continue;
    }
    if (!isCli && parameter.cliOnly) {
      continue;
    }
    parameter.registrar(api, targets, parameter);
  }
}

void registerConfigParameters(ProgramInterface& api, ConfigParameterTargets targets, bool isCli)
{
  registerCatalogWithProgramInterface(api, targets, isCli);
}

std::string parameterCatalogJson()
{
  Json json;
  Json parameters = Json::array();
  for (const auto& parameter : parameterCatalog()) {
    if (parameter.hidden || !parameter.includeInJson) {
      continue;
    }
    Json item;
    item["long"] = "--" + parameter.longName;
    item["short"] = parameter.shortName == '\0' ? "" : std::string("-") + parameter.shortName;
    item["description"] = parameter.description;
    item["group"] = parameter.group;
    item["required"] = parameter.requireArgument;
    item["advanced"] = parameter.isAdvanced;
    item["incremental"] = parameter.incrementalArgument;
    if (parameter.requireArgument) {
      item["longType"] = parameter.argumentName;
      item["shortType"] = std::string(1, ArgType::toShort.at(parameter.argumentName));
      item["default"] = stringDefault(parameter);
    } else {
      item["default"] = false;
    }
    if (!parameter.choices.empty()) {
      item["choices"] = parameter.choices;
    }
    if (!parameter.minValue.empty()) {
      item["min"] = parameter.minValue;
    }
    if (!parameter.maxValue.empty()) {
      item["max"] = parameter.maxValue;
    }
    if (!parameter.renderPolicy.empty()) {
      item["renderPolicy"] = parameter.renderPolicy;
    }
    parameters.push_back(std::move(item));
  }
  json["parameters"] = std::move(parameters);
  return json.dump(2);
}

} // namespace infomap
