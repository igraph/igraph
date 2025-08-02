/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef CONFIG_H_
#define CONFIG_H_

#include "../utils/Date.h"
#include "../version.h"
#include "ProgramInterface.h"

#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <limits>

namespace infomap {

struct FlowModel {
  static constexpr int undirected = 0;
  static constexpr int directed = 1;
  static constexpr int undirdir = 2;
  static constexpr int outdirdir = 3;
  static constexpr int rawdir = 4;
  static constexpr int precomputed = 5;

  int value = 0;

  FlowModel(int val) : value(val) { }
  FlowModel& operator=(int val)
  {
    value = val;
    return *this;
  }

  operator int&() { return value; }
  operator int() const { return value; }
};

std::ostream& operator<<(std::ostream& out, FlowModel f);

inline const char* flowModelToString(const FlowModel& flowModel)
{
  switch (flowModel) {
  case FlowModel::directed:
    return "directed";
  case FlowModel::undirdir:
    return "undirdir";
  case FlowModel::outdirdir:
    return "outdirdir";
  case FlowModel::rawdir:
    return "rawdir";
  case FlowModel::precomputed:
    return "precomputed";
  case FlowModel::undirected:
  default:
    return "undirected";
  }
}

struct Config {
  // Input
  bool isCLI = false;
  std::string networkFile;
  std::vector<std::string> additionalInput;
  bool stateInput = false;
  bool stateOutput = false;
  bool multilayerInput = false;
  double weightThreshold = 0.0;
  bool bipartite = false;
  bool skipAdjustBipartiteFlow = false;
  bool bipartiteTeleportation = false;
  bool noSelfLinks = false; // Replaces includeSelfLinks
  unsigned int nodeLimit = 0;
  unsigned int matchableMultilayerIds = 0;
  std::string clusterDataFile;
  std::string metaDataFile;
  double metaDataRate = 1.0;
  bool unweightedMetaData = false;
  unsigned int numMetaDataDimensions = 0;
  bool clusterDataIsHard = false; // FIXME Not used
  bool assignToNeighbouringModule = false;
  bool noInfomap = false;

  FlowModel flowModel = FlowModel::undirected;
  bool flowModelIsSet = false;
  bool directed = false;
  bool useNodeWeightsAsFlow = false;
  bool teleportToNodes = false;
  double markovTime = 1.0;
  bool variableMarkovTime = false;
  double variableMarkovTimeDamping = 1.0; // 0 for linear scaling, 1 for log scaled.
  double variableMarkovTimeMinLocalScale = 1; // Correspond to two links in undirected unweighted networks. Avoids division by zero.
  bool markovTimeNoSelfLinks = false;
  double multilayerRelaxRate = 0.15;
  int multilayerRelaxLimit = -1; // Amount of layers allowed to jump up or down
  int multilayerRelaxLimitUp = -1; // One-sided limit to higher layers
  int multilayerRelaxLimitDown = -1; // One-sided limit to lower layers
  double multilayerJSRelaxRate = 0.15;
  bool multilayerRelaxByJensenShannonDivergence = false;
  int multilayerJSRelaxLimit = -1;

  // Clustering
  bool twoLevel = false;
  bool noCoarseTune = false;
  bool recordedTeleportation = false;
  bool regularized = false; // Add a Bayesian prior network with recorded teleportation (sets recordedTeleportation and teleportToNodes to true)
  double regularizationStrength = 1.0; // Scale Bayesian prior constant ln(N)/N with this factor
  double teleportationProbability = 0.15;
  unsigned int preferredNumberOfModules = 0;
  bool entropyBiasCorrection = false;
  double entropyBiasCorrectionMultiplier = 1;
  /* unsigned long seedToRandomNumberGenerator = 123; */

  // Performance and accuracy
  unsigned int numTrials = 1;
  double minimumCodelengthImprovement = 1e-10;
  double minimumSingleNodeCodelengthImprovement = 1e-16;
  bool randomizeCoreLoopLimit = false;
  unsigned int coreLoopLimit = 10;
  unsigned int levelAggregationLimit = 0;
  unsigned int tuneIterationLimit = 0; // Iterations of fine-tune/coarse-tune in two-level partition
  double minimumRelativeTuneIterationImprovement = 1e-5;
  bool onlySuperModules = false;
  unsigned int fastHierarchicalSolution = 0;
  bool preferModularSolution = false;
  bool innerParallelization = false;

  // Output
  std::string outDirectory;
  std::string outName;
  std::string outputFormats;
  bool printTree = false;
  bool printFlowTree = false;
  bool printNewick = false;
  bool printJson = false;
  bool printCsv = false;
  bool printClu = false;
  bool printAllTrials = false;
  int cluLevel = 1; // Write modules at specified depth from root. 1, 2, ... or -1 for bottom level
  bool printFlowNetwork = false;
  bool printPajekNetwork = false;
  bool printStateNetwork = false;
  bool noFileOutput = false;
  unsigned int verbosity = 0;
  unsigned int verboseNumberPrecision = 9;
  bool silent = false;
  bool hideBipartiteNodes = false;

  // Other
  Date startDate;
  std::string version = INFOMAP_VERSION;
  std::string parsedString;
  std::vector<ParsedOption> parsedOptions;
  infomap::interruptionHandlerFn *interruptionHandler = NULL;

  Config() = default;

  explicit Config(const std::string& flags, bool isCLI = false);

  Config& cloneAsNonMain(const Config& other)
  {
    isCLI = other.isCLI;
    networkFile = other.networkFile;
    additionalInput = other.additionalInput;
    stateInput = other.stateInput;
    stateOutput = other.stateOutput;
    multilayerInput = other.multilayerInput;
    weightThreshold = other.weightThreshold;
    bipartite = other.bipartite;
    skipAdjustBipartiteFlow = other.skipAdjustBipartiteFlow;
    bipartiteTeleportation = other.bipartiteTeleportation;
    noSelfLinks = other.noSelfLinks;
    nodeLimit = other.nodeLimit;
    matchableMultilayerIds = other.matchableMultilayerIds;
    metaDataRate = other.metaDataRate;
    unweightedMetaData = other.unweightedMetaData;
    numMetaDataDimensions = other.numMetaDataDimensions;
    assignToNeighbouringModule = other.assignToNeighbouringModule;
    noInfomap = other.noInfomap;
    flowModel = other.flowModel;
    flowModelIsSet = other.flowModelIsSet;
    directed = other.directed;
    useNodeWeightsAsFlow = other.useNodeWeightsAsFlow;
    teleportToNodes = other.teleportToNodes;
    markovTime = other.markovTime;
    variableMarkovTime = other.variableMarkovTime;
    variableMarkovTimeDamping = other.variableMarkovTimeDamping;
    markovTimeNoSelfLinks = other.markovTimeNoSelfLinks;
    multilayerRelaxRate = other.multilayerRelaxRate;
    multilayerRelaxLimit = other.multilayerRelaxLimit;
    multilayerRelaxLimitUp = other.multilayerRelaxLimitUp;
    multilayerRelaxLimitDown = other.multilayerRelaxLimitDown;
    multilayerJSRelaxRate = other.multilayerJSRelaxRate;
    multilayerRelaxByJensenShannonDivergence = other.multilayerRelaxByJensenShannonDivergence;
    multilayerJSRelaxLimit = other.multilayerJSRelaxLimit;
    twoLevel = other.twoLevel;
    noCoarseTune = other.noCoarseTune;
    recordedTeleportation = other.recordedTeleportation;
    regularized = other.regularized;
    regularizationStrength = other.regularizationStrength;
    teleportationProbability = other.teleportationProbability;
    entropyBiasCorrection = other.entropyBiasCorrection;
    entropyBiasCorrectionMultiplier = other.entropyBiasCorrectionMultiplier;
    // seedToRandomNumberGenerator = other.seedToRandomNumberGenerator;
    minimumCodelengthImprovement = other.minimumCodelengthImprovement;
    minimumSingleNodeCodelengthImprovement = other.minimumSingleNodeCodelengthImprovement;
    randomizeCoreLoopLimit = other.randomizeCoreLoopLimit;
    minimumRelativeTuneIterationImprovement = other.minimumRelativeTuneIterationImprovement;
    preferModularSolution = other.preferModularSolution;
    innerParallelization = other.innerParallelization;
    outDirectory = other.outDirectory;
    outName = other.outName;
    outputFormats = other.outputFormats;
    verbosity = other.verbosity;
    verboseNumberPrecision = other.verboseNumberPrecision;
    startDate = other.startDate;
    version = other.version;
    return *this;
  }

  void adaptDefaults();

  void setStateInput() { stateInput = true; }

  void setStateOutput() { stateOutput = true; }

  void setMultilayerInput() { multilayerInput = true; }

  void setFlowModel(FlowModel value)
  {
    flowModel = value;
    flowModelIsSet = true;
  }

  bool isUndirectedClustering() const { return flowModel == FlowModel::undirected; }

  bool isUndirectedFlow() const { return flowModel == FlowModel::undirected || flowModel == FlowModel::undirdir; }

  bool printAsUndirected() const { return isUndirectedClustering(); }

  bool isMultilayerNetwork() const { return multilayerInput || !additionalInput.empty(); }
  bool isBipartite() const { return bipartite; }

  bool haveMemory() const { return stateInput; }
  bool printStates() const { return stateOutput; }

  bool haveMetaData() const { return !metaDataFile.empty() || numMetaDataDimensions != 0; }

  bool haveOutput() const { return !noFileOutput; }

  bool haveModularResultOutput() const
  {
    return printTree || printFlowTree || printNewick || printJson || printCsv || printClu;
  }
};

} // namespace infomap

#endif // CONFIG_H_
