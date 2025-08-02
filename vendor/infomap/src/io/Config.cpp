/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "Config.h"
#include "ProgramInterface.h"
#include "SafeFile.h"
#include "../utils/FileURI.h"
#include "../utils/Log.h"
#include <vector>
#include <stdexcept>

namespace infomap {

constexpr int FlowModel::undirected;
constexpr int FlowModel::directed;
constexpr int FlowModel::undirdir;
constexpr int FlowModel::outdirdir;
constexpr int FlowModel::rawdir;
constexpr int FlowModel::precomputed;

Config::Config(const std::string& flags, bool isCLI) : isCLI(isCLI)
{
  ProgramInterface api("Infomap",
                       "Implementation of the Infomap clustering algorithm based on the Map Equation (see www.mapequation.org)",
                       INFOMAP_VERSION);

  api.setGroups({ "Input", "Algorithm", "Accuracy", "Output" });

  std::vector<std::string> optionalOutputDir; // Used if !isCLI
  // --------------------- Input options ---------------------
  if (isCLI) {
    api.addNonOptionArgument(networkFile, "network_file", "File containing the network data. Assumes a link list format if no Pajek formatted heading.", "Input");
  } else {
    api.addOptionArgument(networkFile, "input", "File containing the network data. Assumes a link list format if no Pajek formatted heading.", ArgType::path, "Input");
  }

  api.addOptionArgument(skipAdjustBipartiteFlow, "skip-adjust-bipartite-flow", "Skip distributing all flow from the bipartite nodes to the primary nodes.", "Input", true);

  api.addOptionArgument(bipartiteTeleportation, "bipartite-teleportation", "Teleport like the bipartite flow instead of two-step (unipartite) teleportation.", "Input", true);

  api.addOptionArgument(weightThreshold, "weight-threshold", "Limit the number of links to read from the network. Ignore links with less weight than the threshold.", ArgType::number, "Input", 0.0, true);

  bool deprecated_includeSelfLinks = false;
  api.addOptionArgument(deprecated_includeSelfLinks, 'k', "include-self-links", "DEPRECATED. Include self links by default now, exclude with --no-self-links.", "Input", true).setHidden(true);

  api.addOptionArgument(noSelfLinks, "no-self-links", "Exclude self links in the input network.", "Input", true);

  api.addOptionArgument(nodeLimit, "node-limit", "Limit the number of nodes to read from the network. Ignore links connected to ignored nodes.", ArgType::integer, "Input", 1u, true);

  api.addOptionArgument(matchableMultilayerIds, "matchable-multilayer-ids", "Construct state ids from node and layer ids that are consistent across networks for the same max number of layers. Set to at least the largest layer id among networks to match.", ArgType::integer, "Input", 1u, true);

  api.addOptionArgument(clusterDataFile, 'c', "cluster-data", "Provide an initial two-level (clu format) or multi-layer (tree format) solution.", ArgType::path, "Input");

  api.addOptionArgument(assignToNeighbouringModule, "assign-to-neighbouring-module", "Assign nodes without module assignments (from --cluster-data) to the module assignment of a neighbouring node if possible.", "Input", true);

  api.addOptionArgument(metaDataFile, "meta-data", "Provide meta data (clu format) that should be encoded.", ArgType::path, "Input", true);

  api.addOptionArgument(metaDataRate, "meta-data-rate", "Metadata encoding rate. Default is to encode each step.", ArgType::number, "Input", 0.0, true);

  api.addOptionArgument(unweightedMetaData, "meta-data-unweighted", "Don't weight meta data by node flow.", "Input", true);

  api.addOptionArgument(noInfomap, "no-infomap", "Don't run the optimizer. Useful to calculate codelength of provided cluster data or to print non-modular statistics.", "Input");

  // --------------------- Output options ---------------------

  api.addOptionArgument(outName, "out-name", "Name for the output files, e.g. [output_directory]/[out-name].tree", ArgType::string, "Output", true);

  api.addOptionArgument(noFileOutput, '0', "no-file-output", "Don't write output to file.", "Output", true);

  api.addOptionArgument(printTree, "tree", "Write a tree file with the modular hierarchy. Automatically enabled if no other output is specified.", "Output");

  api.addOptionArgument(printFlowTree, "ftree", "Write a ftree file with the modular hierarchy including aggregated links between (nested) modules. (Used by Network Navigator)", "Output");

  api.addOptionArgument(printClu, "clu", "Write a clu file with the top cluster ids for each node.", "Output");

  api.addOptionArgument(cluLevel, "clu-level", "For clu output, print modules at specified depth from root. Use -1 for bottom level modules.", ArgType::integer, "Output", -1, true);

  api.addOptionArgument(outputFormats, 'o', "output", "Comma-separated output formats without spaces, e.g. -o clu,tree,ftree. Options: clu, tree, ftree, newick, json, csv, network, states, flow.", ArgType::list, "Output", true);

  api.addOptionArgument(hideBipartiteNodes, "hide-bipartite-nodes", "Project bipartite solution to unipartite.", "Output", true);

  api.addOptionArgument(printAllTrials, "print-all-trials", "Print all trials to separate files.", "Output", true);

  // --------------------- Core algorithm options ---------------------
  api.addOptionArgument(twoLevel, '2', "two-level", "Optimize a two-level partition of the network. Default is multi-level.", "Algorithm");

  std::string flowModelArg;

  api.addOptionArgument(flowModelArg, 'f', "flow-model", "Specify flow model. Options: undirected, directed, undirdir, outdirdir, rawdir, precomputed.", ArgType::option, "Algorithm");

  api.addOptionArgument(directed, 'd', "directed", "Assume directed links. Shorthand for '--flow-model directed'.", "Algorithm");

  api.addOptionArgument(recordedTeleportation, 'e', "recorded-teleportation", "If teleportation is used to calculate the flow, also record it when minimizing codelength.", "Algorithm", true);

  api.addOptionArgument(useNodeWeightsAsFlow, "use-node-weights-as-flow", "Use node weights (from api or after names in Pajek format) as flow, normalized to sum to 1", "Algorithm", true);

  api.addOptionArgument(teleportToNodes, "to-nodes", "Teleport to nodes instead of to links, assuming uniform node weights if no such input data.", "Algorithm", true);

  api.addOptionArgument(teleportationProbability, 'p', "teleportation-probability", "Probability of teleporting to a random node or link.", ArgType::probability, "Algorithm", 0.0, 1.0, true);

  api.addOptionArgument(regularized, "regularized", "Effectively add a fully connected Bayesian prior network to not overfit due to missing links. Implies recorded teleportation", "Algorithm", true);

  api.addOptionArgument(regularizationStrength, "regularization-strength", "Adjust relative strength of Bayesian prior network with this multiplier.", ArgType::number, "Algorithm", 0.0, true);

  api.addOptionArgument(entropyBiasCorrection, "entropy-corrected", "Correct for negative entropy bias in small samples (many modules).", "Algorithm", true);

  api.addOptionArgument(entropyBiasCorrectionMultiplier, "entropy-correction-strength", "Increase or decrease the default entropy correction with this factor.", ArgType::number, "Algorithm", true);

  api.addOptionArgument(markovTime, "markov-time", "Scale link flow to change the cost of moving between modules. Higher values results in fewer modules.", ArgType::number, "Algorithm", 0.0, true);

  api.addOptionArgument(variableMarkovTime, "variable-markov-time", "Increase Markov time locally to level out link flow. Reduces risk of overpartitioning sparse areas while keeping high resolution in dense areas.", "Algorithm", true);

  api.addOptionArgument(variableMarkovTimeDamping, "variable-markov-damping", "Damping parameter for variable Markov time, to scale with local effective degree (0) or local entropy (1).", ArgType::number, "Algorithm", true);

  api.addOptionArgument(variableMarkovTimeMinLocalScale, "variable-markov-min-scale", "Minimum local scale for nodes with zero entropy to avoid division by zero. Local Markov time is max scale divided by local scale.", ArgType::number, "Algorithm", true);

  // api.addOptionArgument(markovTimeNoSelfLinks, "markov-time-no-self-links", "For testing.", "Algorithm", true);

  api.addOptionArgument(preferredNumberOfModules, "preferred-number-of-modules", "Penalize solutions the more they differ from this number.", ArgType::integer, "Algorithm", 1u, true);

  api.addOptionArgument(multilayerRelaxRate, "multilayer-relax-rate", "Probability to relax the constraint to move only in the current layer.", ArgType::probability, "Algorithm", 0.0, 1.0, true);

  api.addOptionArgument(multilayerRelaxLimit, "multilayer-relax-limit", "Number of neighboring layers in each direction to relax to. If negative, relax to any layer.", ArgType::integer, "Algorithm", -1, true);

  api.addOptionArgument(multilayerRelaxLimitUp, "multilayer-relax-limit-up", "Number of neighboring layers with higher id to relax to. If negative, relax to any layer.", ArgType::integer, "Algorithm", -1, true);

  api.addOptionArgument(multilayerRelaxLimitDown, "multilayer-relax-limit-down", "Number of neighboring layers with lower id to relax to. If negative, relax to any layer.", ArgType::integer, "Algorithm", -1, true);

  api.addOptionArgument(multilayerRelaxByJensenShannonDivergence, "multilayer-relax-by-jsd", "Relax proportional to the out-link similarity measured by the Jensen-Shannon divergence.", "Algorithm", true);

  // --------------------- Performance and accuracy options ---------------------
  // api.addOptionArgument(seedToRandomNumberGenerator, 's', "seed", "A seed (integer) to the random number generator for reproducible results.", ArgType::integer, "Accuracy", 1ul);

  api.addOptionArgument(numTrials, 'N', "num-trials", "Number of outer-most loops to run before picking the best solution.", ArgType::integer, "Accuracy", 1u);

  api.addOptionArgument(coreLoopLimit, 'M', "core-loop-limit", "Limit the number of loops that tries to move each node into the best possible module.", ArgType::integer, "Accuracy", 1u, true);

  api.addOptionArgument(levelAggregationLimit, 'L', "core-level-limit", "Limit the number of times the core loops are reapplied on existing modular network to search bigger structures.", ArgType::integer, "Accuracy", 1u, true);

  api.addOptionArgument(tuneIterationLimit, 'T', "tune-iteration-limit", "Limit the number of main iterations in the two-level partition algorithm. 0 means no limit.", ArgType::integer, "Accuracy", 1u, true);

  api.addOptionArgument(minimumCodelengthImprovement, "core-loop-codelength-threshold", "Minimum codelength threshold for accepting a new solution in core loop.", ArgType::number, "Accuracy", 0.0, true);

  api.addOptionArgument(minimumRelativeTuneIterationImprovement, "tune-iteration-relative-threshold", "Set codelength improvement threshold of each new tune iteration to 'f' times the initial two-level codelength.", ArgType::number, "Accuracy", 0.0, true);

  api.addIncrementalOptionArgument(fastHierarchicalSolution, 'F', "fast-hierarchical-solution", "Find top modules fast. Use -FF to keep all fast levels. Use -FFF to skip recursive part.", "Accuracy", true);

  api.addOptionArgument(preferModularSolution, "prefer-modular-solution", "Prefer modular solutions even if they are worse than putting all nodes in one module.", "Accuracy", true);

  api.addOptionArgument(innerParallelization, "inner-parallelization", "Parallelize the inner-most loop for greater speed. This may give some accuracy tradeoff.", "Accuracy", true);

  api.addOptionalNonOptionArguments(optionalOutputDir, "out_directory", "Directory to write the results to.", "Output");

  api.addIncrementalOptionArgument(verbosity, 'v', "verbose", "Verbose output on the console. Add additional 'v' flags to increase verbosity up to -vvv.", "Output");

  api.addOptionArgument(silent, "silent", "No output on the console.", "Output");

  api.parseArgs(flags);

  if (deprecated_includeSelfLinks) {
    throw std::runtime_error("The --include-self-links flag is deprecated to include self links by default. Use --no-self-links to exclude.");
  }

  if (!optionalOutputDir.empty())
    outDirectory = optionalOutputDir[0];

  if (!isCLI && outDirectory.empty())
    noFileOutput = true;

  if (!noFileOutput && outDirectory.empty() && isCLI) {
    throw std::runtime_error("Missing out_directory");
  }

  if (flowModelArg == "directed" || directed) {
    setFlowModel(FlowModel::directed);
  } else if (flowModelArg == "undirected") {
    setFlowModel(FlowModel::undirected);
  } else if (flowModelArg == "undirdir") {
    setFlowModel(FlowModel::undirdir);
  } else if (flowModelArg == "outdirdir") {
    setFlowModel(FlowModel::outdirdir);
  } else if (flowModelArg == "rawdir") {
    setFlowModel(FlowModel::rawdir);
  } else if (flowModelArg == "precomputed") {
    setFlowModel(FlowModel::precomputed);
  } else if (!flowModelArg.empty()) {
    throw std::runtime_error(io::Str() << "Unrecognized flow model: '" << flowModelArg << "'");
  }

  if (regularized) {
    recordedTeleportation = true;
  }

  if (*--outDirectory.end() != '/')
    outDirectory.append("/");

  if (haveOutput() && !isDirectoryWritable(outDirectory))
    throw std::runtime_error(io::Str() << "Can't write to directory '" << outDirectory << "'. Check that the directory exists and that you have write permissions.");

  if (outName.empty()) {
    outName = !networkFile.empty() ? FileURI(networkFile).getName() : "no-name";
  }

  if (noInfomap) {
    numTrials = 1;
  }

  parsedString = flags;
  parsedOptions = api.getUsedOptionArguments();

  if (printAllTrials && numTrials < 2) {
    printAllTrials = false;
  }

  adaptDefaults();

  Log::init(verbosity, silent, verboseNumberPrecision);
}

void Config::adaptDefaults()
{
  auto outputs = io::split(outputFormats, ',');
  for (std::string& o : outputs) {
    if (o == "clu") {
      printClu = true;
    } else if (o == "tree") {
      printTree = true;
    } else if (o == "ftree") {
      printFlowTree = true;
    } else if (o == "newick") {
      printNewick = true;
    } else if (o == "json") {
      printJson = true;
    } else if (o == "csv") {
      printCsv = true;
    } else if (o == "network") {
      printPajekNetwork = true;
    } else if (o == "flow") {
      printFlowNetwork = true;
    } else if (o == "states") {
      printStateNetwork = true;
    } else {
      throw std::runtime_error(io::Str() << "Unrecognized output format: '" << o << "'.");
    }
  }

  // Of no output format specified, use tree as default (if not used as a library).
  if (isCLI && !haveModularResultOutput()) {
    printTree = true;
  }
}

std::ostream& operator<<(std::ostream& out, FlowModel f)
{
  return out << flowModelToString(f);
}

} // namespace infomap
