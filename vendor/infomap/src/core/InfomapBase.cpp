/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "InfomapBase.h"
#include "InfomapConfig.h"
#include "InfoNode.h"
#include "FlowData.h"
#include "BiasedMapEquation.h"
#include "MemMapEquation.h"
#include "MetaMapEquation.h"
#include "InfomapOptimizer.h"
#include "../io/SafeFile.h"
#include "../utils/FileURI.h"
#include "../utils/FlowCalculator.h"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace infomap {

std::map<unsigned int, std::vector<unsigned int>> InfomapBase::getMultilevelModules(bool states)
{
  if (haveMemory() && !states) {
    throw std::runtime_error("Cannot get multilevel modules on higher-order network without states.");
  }
  unsigned int maxDepth = maxTreeDepth();
  unsigned int numModuleLevels = maxDepth - 1;
  std::map<unsigned int, std::vector<unsigned int>> modules;
  if (maxDepth < 2) return modules;
  for (unsigned int level = 1; level <= numModuleLevels; ++level) {
    for (auto it(iterLeafNodes(static_cast<int>(level))); !it.isEnd(); ++it) {
      auto& node = *it;
      auto nodeId = states ? node.stateId : node.physicalId;
      modules[nodeId].push_back(it.moduleId());
    }
  }
  return modules;
}

// ===================================================
// IO
// ===================================================

std::ostream& operator<<(std::ostream& out, const InfomapBase& infomap)
{
  return infomap.toString(out);
}

// ===================================================
// Getters
// ===================================================

unsigned int InfomapBase::numLevels() const
{
  // TODO: Make sure this is not called unless tree is guaranteed to have even depth!
  unsigned int depth = 0;
  InfoNode* n = m_root.firstChild;
  while (n != nullptr) {
    ++depth;
    n = n->firstChild;
  }
  return depth;
}

unsigned int InfomapBase::maxTreeDepth() const
{
  unsigned int maxDepth = 0;
  for (auto it = root().begin_tree(); !it.isEnd(); ++it) {
    const InfoNode& node = *it;
    if (node.isLeaf()) {
      if (it.depth() > maxDepth) {
        maxDepth = it.depth();
      }
    }
  }
  return maxDepth;
}

// ===================================================
// Run
// ===================================================

void InfomapBase::run(const std::string& parameters)
{
  if (!isMainInfomap()) {
    runPartition();
    return;
  }

  m_elapsedTime = Stopwatch(true);
  m_startDate = Date();

  std::string currentParameters = io::Str() << m_initialParameters << (parameters.empty() ? "" : " ") << parameters;
  if (currentParameters != m_currentParameters) {
    m_currentParameters = currentParameters;
    setConfig(Config(m_currentParameters, isCLI));
    m_network.setConfig(*this);
  }

  Log::init(verbosity, silent, verboseNumberPrecision);

  Log() << "=======================================================\n";
  Log() << "  Infomap v" << INFOMAP_VERSION << " starts at " << m_startDate << "\n";
  Log() << "  -> Input network: " << networkFile << "\n";
  if (noFileOutput)
    Log() << "  -> No file output!\n";
  else
    Log() << "  -> Output path:   " << outDirectory << "\n";
  if (!parsedOptions.empty()) {
    for (unsigned int i = 0; i < parsedOptions.size(); ++i)
      Log() << (i == 0 ? "  -> Configuration: " : "                    ") << parsedOptions[i] << "\n";
  }
  if (!m_initialPartition.empty()) {
    Log() << "  -> " << m_initialPartition.size() << " initial module ids provided\n";
  }
  Log() << "=======================================================\n";
#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
  {
    Log() << "  OpenMP " << _OPENMP << " detected with " << omp_get_num_threads() << " threads...\n";
  }
#endif

  m_network.postProcessInputData();
  if (m_network.numNodes() == 0) {
    m_network.readInputData(networkFile);
  }

  if (!metaDataFile.empty()) {
    initMetaData(metaDataFile);
  }

  run(m_network);

  Log() << "===================================================\n";
  Log() << "  Infomap ends at " << m_endDate << "\n";
  Log() << "  (Elapsed time: " << m_elapsedTime << ")\n";
  Log() << "===================================================\n";
}

void InfomapBase::run(Network& network)
{
  if (!isMainInfomap())
    throw std::logic_error("Can't run a non-main Infomap with an input network");

  if (network.numNodes() == 0) {
    network.postProcessInputData();
    if (network.numNodes() == 0) {
      throw std::domain_error("Network is empty");
    }
  }

  if (printStateNetwork) {
    std::string filename = outDirectory + outName + "_states.net";
    Log() << "Writing state network to '" << filename << "'... ";
    network.writeStateNetwork(filename);
    Log() << "done!\n";
  }

  if (printPajekNetwork) {
    std::string filename;
    if (printStates()) {
      filename = outDirectory + outName + "_states_as_physical.net";
      Log() << "Writing state network as first order Pajek network to '" << filename << "'... ";
    } else {
      // Non-memory input
      filename = outDirectory + outName + ".net";
      Log() << "Writing Pajek network to '" << filename << "'... ";
    }
    network.writePajekNetwork(filename);
    Log() << "done!\n";
  }

  if (network.haveMemoryInput()) {
    Log() << "  -> Found higher order network input, using the Map Equation for higher order network flows\n";
    setStateInput();
    setStateOutput();

    if (network.isMultilayerNetwork() && !isMultilayerNetwork()) {
      setMultilayerInput();
    }
  } else {
    if (haveMemory() || network.higherOrderInputMethodCalled()) {
      Log() << "  -> Warning: Higher order network specified but no higher order input found.\n";
      // Use state output anyway for consistency even in the special case when input is first order
      setStateOutput();
    }
    Log() << "  -> Ordinary network input, using the Map Equation for first order network flows\n";
  }

  if (network.haveDirectedInput() && isUndirectedFlow()) {
    Log() << "  -> Notice: Directed input found, changing flow model from '" << flowModel << "' to '" << FlowModel(FlowModel::directed) << "'\n";
    flowModel = FlowModel::directed;
  }
  network.setConfig(*this);

  calculateFlow(network, *this);

  if (network.isBipartite()) {
    bipartite = true;
  }

  initNetwork(network);

  if (numLeafNodes() == 0)
    throw std::domain_error("No nodes to partition");

  if (printFlowNetwork) {
    std::string filename;
    if (printStates()) {
      filename = outDirectory + outName + "_states_as_physical_flow.net";
      Log() << "Writing flow state network as first order Pajek network to '" << filename << "'... ";
    } else {
      // Non-memory input
      filename = outDirectory + outName + "_flow.net";
      Log() << "Writing flow network to '" << filename << "'... ";
    }
    network.writePajekNetwork(filename, true);
    Log() << "done!\n";
  }

  // If used as a library, we may want to reuse the network instance, else clear to use less memory
  // TODO: May have to use some meta data for output?
  if (isCLI) {
    network.clearLinks();
  }

  if (haveMemory())
    Log(2) << "Run Infomap with memory...\n";
  else
    Log(2) << "Run Infomap...\n";

  std::ostringstream bestSolutionStatistics;
  unsigned int bestNumLevels = 0;
  double bestHierarchicalCodelength = std::numeric_limits<double>::max();
  m_codelengths.clear();
  NodePaths bestTree(numLeafNodes());
  unsigned int bestTrialIndex = 0;

  unsigned int numTrials = this->numTrials;

  for (unsigned int i = 0; i < numTrials; ++i) {
    removeModules();
    auto startDate = Date();
    Stopwatch timer(true);

    if (isMainInfomap()) {
      Log() << "\n";
      Log() << "================================================\n";
      Log() << "Trial " << (i + 1) << "/" << numTrials << " starting at " << startDate << "\n";
      Log() << "================================================\n";

      if (!clusterDataFile.empty())
        initPartition(clusterDataFile, clusterDataIsHard, &network);
      else if (!m_initialPartition.empty())
        initPartition(m_initialPartition, clusterDataIsHard);
    }

    if (!noInfomap)
      runPartition();
    else
      m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);

    if (haveHardPartition())
      restoreHardPartition();

    if (isMainInfomap()) {
      Log() << "\n=> Trial " << (i + 1) << "/" << numTrials << " finished in " << timer.getElapsedTimeInSec() << "s with codelength " << m_hierarchicalCodelength << "\n";
      m_codelengths.push_back(m_hierarchicalCodelength);

      if (printAllTrials && numTrials > 1) {
        writeResult(static_cast<int>(i + 1));
      }

      if (m_hierarchicalCodelength < bestHierarchicalCodelength - 1e-10) {
        bestSolutionStatistics.clear();
        bestSolutionStatistics.str("");
        bestNumLevels = printPerLevelCodelength(root(), bestSolutionStatistics);
        bestHierarchicalCodelength = m_hierarchicalCodelength;
        bestTrialIndex = i;
        root().sortChildrenOnFlow();
        writeResult();
        if (numTrials > 1) {
          bestTree.clear();
          for (auto it(iterLeafNodes()); !it.isEnd(); ++it) {
            bestTree.emplace_back(it->stateId, it.path());
          }
        }
      }
    }
  }
  if (isMainInfomap()) {
    m_elapsedTime.stop();
    m_endDate = Date();
    Log() << "\n\n";
    Log() << "================================================\n";
    Log() << "Summary after " << numTrials << (numTrials > 1 ? " trials\n" : " trial\n");
    Log() << "================================================\n";
    if (numTrials > 1) {
      if (bestTrialIndex < numTrials - 1) {
        // Restore Infomap tree to best solution
        initTree(bestTree);
        writeResult(); // Overwrite result to get total elapsed time in output file header
      }

      Log() << std::fixed << std::setprecision(9);
      double averageCodelength = 0.0;
      double minCodelength = m_codelengths[0];
      double maxCodelength = m_codelengths[0];
      Log() << "Codelengths: [";
      for (auto codelength : m_codelengths) {
        Log() << codelength << ", ";
        averageCodelength += codelength;
        minCodelength = std::min(minCodelength, codelength);
        maxCodelength = std::max(maxCodelength, codelength);
      }
      averageCodelength /= numTrials;
      Log() << "\b\b]\n";
      Log() << "[min, average, max] codelength: [" << minCodelength << ", " << averageCodelength << ", " << maxCodelength << "]\n\n";
    }

    Log() << "Best end modular solution in " << bestNumLevels << " levels";
    if (bestHierarchicalCodelength > m_oneLevelCodelength)
      Log() << " (warning: worse than one-level solution)";
    Log() << ":\n";
    Log() << bestSolutionStatistics.str() << '\n';
  }
}

// ===================================================
// Run: Init: *
// ===================================================

InfomapBase& InfomapBase::initMetaData(const std::string& metaDataFile)
{
  m_network.readMetaData(metaDataFile);
  return *this;
}

InfomapBase& InfomapBase::initNetwork(Network& network)
{
  if (network.numNodes() == 0)
    throw std::domain_error("No nodes in network");
  if (m_root.childDegree() > 0) {
    m_root.deleteChildren();
    m_leafNodes.clear();
  }
  generateSubNetwork(network);

  initOptimizer();

  init();
  return *this;
}

InfomapBase& InfomapBase::initNetwork(InfoNode& parent, bool asSuperNetwork)
{
  generateSubNetwork(parent);

  if (asSuperNetwork)
    transformNodeFlowToEnterFlow(m_root);

  init();
  return *this;
}

InfomapBase& InfomapBase::initPartition(const std::string& clusterDataFile, bool hard, const Network* network)
{
  FileURI file(clusterDataFile);
  ClusterMap clusterMap;
  if (this->isMultilayerNetwork() && network != nullptr) {
    auto map = network->layerNodeToStateId();
    clusterMap.readClusterData(clusterDataFile, false, &map);
  } else {
    clusterMap.readClusterData(clusterDataFile);
  }

  Log() << "Init partition from file '" << clusterDataFile << "'... ";

  const auto& ext = clusterMap.extension();

  if (ext == "tree" || ext == "ftree") {
    initTree(clusterMap.nodePaths());
  } else if (ext == "clu") {
    initPartition(clusterMap.clusterIds(), hard);
  }

  Log() << "done! Generated " << numLevels() << " levels with codelength " << getIndexCodelength() << " + " << (m_hierarchicalCodelength - getIndexCodelength()) << " = " << io::toPrecision(m_hierarchicalCodelength) << '\n';

  return *this;
}

InfomapBase& InfomapBase::initTree(const NodePaths& tree)
{
  Log(4) << "Init tree... " << std::setprecision(9);
  auto maxDepth = 2;
  std::map<unsigned int, unsigned int> nodeIdToIndex;
  auto leafIndex = 0;
  for (auto& leafNode : m_leafNodes) {
    // Also detach leaf nodes to delete all modules, safe to call multiple times
    leafNode->parent->releaseChildren();
    // Set parent to nullptr to be able to collect any orphaned nodes in the end
    leafNode->parent = nullptr;
    nodeIdToIndex[leafNode->stateId] = leafIndex;
    ++leafIndex;
  }
  m_root.deleteChildren();

  auto numNodesFound = 0;
  auto numNodesNotInNetwork = 0;
  for (const auto& nodePath : tree) {
    ++numNodesFound;
    InfoNode* node = &root();
    auto depth = 0;
    const auto& path = nodePath.second;
    const auto nodeId = nodePath.first;
    InfoNode* leafNode = nullptr;

    try {
      auto nodeIndex = nodeIdToIndex.at(nodeId);
      leafNode = m_leafNodes[nodeIndex];
    } catch (std::exception& e) {
      ++numNodesNotInNetwork;
      continue;
    }

    for (size_t i = 0; i < path.size(); ++i) {
      auto childNumber = path[i]; // 1-based indexing

      // Create new node if path doesn't exist
      // TODO: Check correct tree indexing?
      if (node->childDegree() < childNumber) {
        InfoNode* child = leafNode;
        if (i + 1 < path.size()) {
          child = new InfoNode();
        }
        node->addChild(child);
      }

      node = node->lastChild;
      ++depth;
    }
    maxDepth = std::max(maxDepth, depth);
  }

  if (numNodesNotInNetwork > 0) {
    Log(1) << "\n -> " << numNodesNotInNetwork << "/" << numNodesFound << " nodes in tree not found in network.";
  }

  auto numNodesAddedToNeighbouringModules = 0;
  auto numNodesWithoutClusterInfo = 0;

  // Set orphaned nodes to their own or neighbouring module
  for (auto& leafNode : m_leafNodes) {
    if (leafNode->parent != nullptr) {
      continue;
    }

    ++numNodesWithoutClusterInfo;
    if (assignToNeighbouringModule) {
      // Take first neighbour that has a module assigned
      for (auto& edge : leafNode->outEdges()) {
        if (edge->target->parent != nullptr) {
          edge->target->parent->addChild(leafNode);
          ++numNodesAddedToNeighbouringModules;
          break;
        }
      }
      // Check incoming links if still orphan
      if (leafNode->parent == nullptr) {
        for (auto& edge : leafNode->inEdges()) {
          if (edge->source->parent != nullptr) {
            edge->source->parent->addChild(leafNode);
            ++numNodesAddedToNeighbouringModules;
            break;
          }
        }
      }
      // Set to own module if no neighbour module available
      if (leafNode->parent == nullptr) {
        auto module = new InfoNode();
        root().addChild(module);
        module->addChild(leafNode);
      }
    } else {
      // Set to own module if no neighbour module available
      auto module = new InfoNode();
      root().addChild(module);
      module->addChild(leafNode);
    }
  }
  if (numNodesWithoutClusterInfo > 0) {
    if (assignToNeighbouringModule) {
      Log() << "\n -> " << numNodesWithoutClusterInfo << " nodes not found in tree, " << numNodesAddedToNeighbouringModules << " are put into neighbouring modules and " << (numNodesWithoutClusterInfo - numNodesAddedToNeighbouringModules) << " in separate.";
    } else {
      Log() << "\n -> " << numNodesWithoutClusterInfo << " nodes not found in tree are put into separate modules.";
    }
  }

  aggregateFlowValuesFromLeafToRoot();
  initTree();
  initNetwork();
  m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();
  // Init partition to calculate indexCodelength and moduleCodelength
  if (haveModules())
    setActiveNetworkFromChildrenOfRoot();
  else
    setActiveNetworkFromLeafs();
  initPartition();

  m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
  Log() << "\n -> Initiated to codelength " << m_hierarchicalCodelength << " in " << maxDepth << " levels with " << numTopModules() << " top modules.\n";
  Log() << std::setprecision(6);
  return *this;
}

/**
 * clusterIds is a map from nodeId -> clusterId
 */
InfomapBase& InfomapBase::initPartition(const std::map<unsigned int, unsigned int>& clusterIds, bool hard)
{
  // Generate map from node id to index in leaf node vector
  std::map<unsigned int, unsigned int> nodeIdToIndex;
  unsigned int leafIndex = 0;
  for (auto& nodePtr : m_leafNodes) {
    auto& leafNode = *nodePtr;
    nodeIdToIndex[leafNode.stateId] = leafIndex;
    ++leafIndex;
  }

  // Remap clusterIds to vector from leaf node index to module index instead of nodeId -> clusterId
  std::map<unsigned int, std::set<unsigned int>> clusterIdToNodeIds;
  for (auto it : clusterIds) {
    clusterIdToNodeIds[it.second].insert(it.first);
  }
  unsigned int numNodes = numLeafNodes();
  std::vector<unsigned int> modules(numNodes);
  std::vector<unsigned int> selectedNodes(numNodes, 0);
  unsigned int moduleIndex = 0;
  for (const auto& it : clusterIdToNodeIds) {
    auto& nodes = it.second;
    for (auto& nodeId : nodes) {
      auto nodeIndex = nodeIdToIndex[nodeId];
      modules[nodeIndex] = moduleIndex;
      ++selectedNodes[nodeIndex];
    }
    ++moduleIndex;
  }

  unsigned int numNodesWithoutClusterInfo = 0;
  for (auto& count : selectedNodes) {
    if (count == 0) {
      ++numNodesWithoutClusterInfo;
    }
  }

  if (numNodesWithoutClusterInfo == 0) {
    return initPartition(modules, hard);
  }

  if (assignToNeighbouringModule) {
    Log() << "\n -> " << numNodesWithoutClusterInfo << " nodes not found in cluster file are put into neighbouring modules if possible. ";
    for (unsigned int i = 0; i < numNodes; ++i) {
      if (selectedNodes[i] == 0) {
        // Check out edges greedily for connected modules
        auto& node = *m_leafNodes[i];
        for (auto& e : node.outEdges()) {
          auto& edge = *e;
          auto targetNodeIndex = nodeIdToIndex[edge.target->stateId];
          if (selectedNodes[targetNodeIndex] != 0) {
            selectedNodes[i] = 1;
            modules[i] = modules[targetNodeIndex];
            break;
          }
        }
        if (selectedNodes[i] == 0) {
          // Check in edges greedily for connected modules
          for (auto& e : node.inEdges()) {
            auto& edge = *e;
            auto sourceNodeIndex = nodeIdToIndex[edge.source->stateId];
            if (selectedNodes[sourceNodeIndex] != 0) {
              selectedNodes[i] = 1;
              modules[i] = modules[sourceNodeIndex];
              break;
            }
          }
        }
        if (selectedNodes[i] == 0) {
          // No connected node with a module index, fall back to create new module
          selectedNodes[i] = 1;
          modules[i] = moduleIndex;
          ++moduleIndex;
        }
      }
    }
  } else {
    Log() << "\n -> " << numNodesWithoutClusterInfo << " nodes not found in cluster file are put into separate modules. ";
    // Put non-selected nodes in its own module
    for (unsigned int i = 0; i < numNodes; ++i) {
      if (selectedNodes[i] == 0) {
        modules[i] = moduleIndex;
        ++moduleIndex;
      }
    }
  }

  return initPartition(modules);
}

InfomapBase& InfomapBase::initPartition(std::vector<std::vector<unsigned int>>& clusters, bool hard)
{
  Log(3) << "\n -> Init " << (hard ? "hard" : "soft") << " partition... " << std::flush;
  unsigned int numNodes = numLeafNodes();
  // Get module indices from the merge structure
  std::vector<unsigned int> modules(numNodes);
  std::vector<unsigned int> selectedNodes(numNodes, 0);
  unsigned int moduleIndex = 0;

  for (auto& cluster : clusters) {
    if (cluster.empty())
      continue;
    for (auto id : cluster) {
      ++selectedNodes[id];
      modules[id] = moduleIndex;
    }
    ++moduleIndex;
  }

  // Put non-selected nodes in its own module
  unsigned int numNodesWithoutClusterInfo = 0;
  for (unsigned int i = 0; i < numNodes; ++i) {
    if (selectedNodes[i] == 0) {
      modules[i] = moduleIndex;
      ++moduleIndex;
      ++numNodesWithoutClusterInfo;
    }
  }
  if (numNodesWithoutClusterInfo > 0)
    Log() << "\n -> " << numNodesWithoutClusterInfo << " nodes not found in cluster file are put into separate modules. ";

  return initPartition(modules, hard);
}

InfomapBase& InfomapBase::initPartition(std::vector<unsigned int>& modules, bool hard)
{
  removeModules();
  setActiveNetworkFromLeafs();
  initPartition(); // TODO: confusing same name, should be able to init default without arguments here too
  moveActiveNodesToPredefinedModules(modules);
  consolidateModules(false);

  if (hard) {
    // Save the original network
    m_originalLeafNodes.swap(m_leafNodes);

    // Use the pre-partitioned network as the new leaf network
    unsigned int numNodesInNewNetwork = numTopModules();
    m_leafNodes.resize(numNodesInNewNetwork);
    unsigned int nodeIndex = 0;
    unsigned int numLinksInNewNetwork = 0;
    for (InfoNode& node : m_root) {
      m_leafNodes[nodeIndex] = &node;
      // Collapse children
      node.collapseChildren();
      numLinksInNewNetwork += node.outDegree();
      ++nodeIndex;
    }

    Log(1) << "\n -> Hard-partitioned the network to " << numNodesInNewNetwork << " nodes and " << numLinksInNewNetwork << " links with codelength " << *this << '\n';
  } else {
    Log(1) << "\n -> Initiated to codelength " << *this << " in " << numTopModules() << " top modules.\n";
  }
  m_hierarchicalCodelength = getCodelength();

  return *this;
}

void InfomapBase::generateSubNetwork(Network& network)
{
  unsigned int numNodes = network.numNodes();
  auto& metaData = network.metaData();
  numMetaDataDimensions = network.numMetaDataColumns();

  Log() << "Build internal network with " << numNodes << " nodes and " << network.numLinks() << " links...\n";
  if (!metaData.empty()) {
    Log() << "and " << metaData.size() << " meta-data records in " << numMetaDataDimensions << " dimensions\n";
  }

  m_leafNodes.reserve(numNodes);
  double sumNodeFlow = 0.0;
  double sumTeleFlow = 0.0;
  std::map<unsigned int, unsigned int> nodeIndexMap;
  for (auto& nodeIt : network.nodes()) {
    auto& networkNode = nodeIt.second;
    auto* node = new InfoNode(networkNode.flow, networkNode.id, networkNode.physicalId, networkNode.layerId);
    node->data.teleportWeight = networkNode.weight;
    node->data.teleportFlow = networkNode.teleFlow;
    node->data.exitFlow = networkNode.exitFlow;
    node->data.enterFlow = networkNode.enterFlow;
    if (haveMetaData()) {
      auto meta = metaData.find(networkNode.id);
      if (meta != metaData.end()) {
        node->metaData = meta->second;
      } else {
        node->metaData = std::vector<int>(numMetaDataDimensions, -1);
      }
    }
    sumNodeFlow += networkNode.flow;
    sumTeleFlow += networkNode.teleFlow;
    m_root.addChild(node);
    nodeIndexMap[networkNode.id] = m_leafNodes.size();
    m_leafNodes.push_back(node);
  }
  m_root.data.flow = sumNodeFlow;
  m_root.data.teleportFlow = sumTeleFlow;
  // m_calculateEnterExitFlow = true; //TODO: Implement always in flow calculation
  if (!this->regularized) {
    m_calculateEnterExitFlow = true;
  }
  BiasedMapEquation::setNetworkProperties(network);

  if (std::abs(sumNodeFlow - 1.0) > 1e-10)
    Log() << "Warning, total flow on nodes differs from 1.0 by " << sumNodeFlow - 1.0 << ".\n";

  bool changeMarkovTime = std::abs(markovTime - 1) > 1e-3;
  if (changeMarkovTime) {
    Log() << "  -> Rescale link flow with global Markov time " << markovTime << "\n";
  }

  for (auto& linkIt : network.nodeLinkMap()) {
    unsigned int linkSourceId = linkIt.first.id;
    unsigned int sourceIndex = nodeIndexMap[linkSourceId];
    const auto& subLinks = linkIt.second;
    for (auto& subIt : subLinks) {
      unsigned int linkTargetId = subIt.first.id;
      unsigned int targetIndex = nodeIndexMap[linkTargetId];
      // Ignore self-links in optimization as it doesn't change enter/exit flow on modular level
      if (sourceIndex != targetIndex) {
        auto& linkData = subIt.second;
        m_leafNodes[sourceIndex]->addOutEdge(*m_leafNodes[targetIndex], linkData.weight, linkData.flow * markovTime);
      }
    }
  }

  if (variableMarkovTime) {
    Log() << "  -> Rescale link flow with variable Markov time\n";
    if (std::abs(variableMarkovTimeDamping - 1) > 1e-9) {
      Log() << "  -> Use variable Markov time strength " << variableMarkovTimeDamping << "\n";
    }
  }

  double maxEntropy = 0.0;
  double maxFlow = 0.0;
  double entropyRate = 0.0;
  unsigned int maxDegree = 0;
  unsigned int maxOutDegree = 0;
  unsigned int maxInDegree = 0;
  double totDegree = network.sumDegree();
  std::vector<double> entropies(numNodes, 0);

  for (unsigned i = 0; i < numNodes; ++i) {
    InfoNode& node = *m_leafNodes[i];
    maxDegree = std::max(maxDegree, node.degree());
    maxOutDegree = std::max(maxOutDegree, node.outDegree());
    maxInDegree = std::max(maxInDegree, node.inDegree());
    double entropy = 0;
    double sumOut = 0;
    for (InfoEdge* e : node.outEdges()) {
      entropy -= infomath::plogp(e->data.weight);
      sumOut += e->data.weight;
    }
    if (isUndirectedFlow()) {
      for (InfoEdge* e : node.inEdges()) {
        entropy -= infomath::plogp(e->data.weight);
        sumOut += e->data.weight;
      }
    }
    entropy = sumOut > 1e-9 ? (entropy + infomath::plogp(sumOut)) / sumOut : 0;
    maxEntropy = std::max(maxEntropy, entropy);
    entropyRate += m_leafNodes[i]->data.flow * entropy;
    maxFlow = std::max(maxFlow, node.data.flow);
    entropies[i] = entropy; // Store for undirected networks
  }

  m_entropyRate = entropyRate;
  m_maxEntropy = maxEntropy;
  m_maxFlow = maxFlow;

  double minLocalScale = variableMarkovTimeMinLocalScale;
  double damping = variableMarkovTimeDamping;

  double maxScale = infomath::linlog(maxFlow * totDegree, damping);

  if (variableMarkovTime) {
    if (damping < 0) {
      maxScale = infomath::linlog(pow(2.0, maxEntropy), -damping);
    }
    for (unsigned i = 0; i < numNodes; ++i) {
      InfoNode& node = *m_leafNodes[i];
      double localScale = damping < 0 ? infomath::linlog(pow(2.0, entropies[i]), -damping) : infomath::linlog(std::max(minLocalScale, node.data.flow * totDegree), damping);
      for (InfoEdge* e : node.outEdges()) {
        if (isUndirectedFlow()) {
          double oppositeLocalScale = damping < 0 ? infomath::linlog(pow(2.0, entropies[nodeIndexMap[e->target->stateId]]), -damping) : infomath::linlog(std::max(minLocalScale, e->target->data.flow * totDegree), damping);
          localScale = std::max(localScale, oppositeLocalScale);
        }
        double localMarkovTimeScale = maxScale / std::max(minLocalScale, localScale);
        e->data.flow *= localMarkovTimeScale;
        network.nodeLinkMap()[e->source->stateId][e->target->stateId].flow = e->data.flow;
      }
    }
  }
}

void InfomapBase::generateSubNetwork(InfoNode& parent)
{
  // Store parent module
  m_root.owner = &parent;
  m_root.data = parent.data;

  unsigned int numNodes = parent.childDegree();
  m_leafNodes.resize(numNodes);

  Log(1) << "Generate sub network with " << numNodes << " nodes...\n";

  unsigned int childIndex = 0;
  for (InfoNode& node : parent) {
    auto* clonedNode = new InfoNode(node);
    clonedNode->initClean();
    m_root.addChild(clonedNode);
    node.index = childIndex; // Set index to its place in this subnetwork to be able to find edge target below
    m_leafNodes[childIndex] = clonedNode;
    ++childIndex;
  }

  InfoNode* parentPtr = &parent;
  // Clone edges
  for (InfoNode& node : parent) {
    for (InfoEdge* e : node.outEdges()) {
      InfoEdge& edge = *e;
      // If neighbour node is within the same module, add the link to this subnetwork.
      if (edge.target->parent == parentPtr) {
        m_leafNodes[edge.source->index]->addOutEdge(*m_leafNodes[edge.target->index], edge.data.weight, edge.data.flow);
      }
    }
  }
}

void InfomapBase::init()
{
  Log(3) << "InfomapBase::init()...\n";

  if (m_calculateEnterExitFlow && isMainInfomap())
    initEnterExitFlow();

  initNetwork();

  m_oneLevelCodelength = calcCodelength(m_root);
  Log() << "  -> One-level codelength: " << m_oneLevelCodelength << '\n';
}

// ===================================================
// Run: *
// ===================================================

void InfomapBase::hierarchicalPartition()
{
  Log(1) << "Hierarchical partition...\n";

  const auto depth = maxTreeDepth();
  if (depth > 2) {
    Log(1) << "Continuing from a tree with " << depth << " levels...\n";

    if (fastHierarchicalSolution == 0) {
      Log(1) << "Removing sub modules...\n";
      removeSubModules(true);
      m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
    }

    else if (fastHierarchicalSolution == 1) {
      Log(1) << "Fine-tune bottom modules... ";
      bool isSilent = isMainInfomap() && Log::isSilent();

      double codelengthBefore = 0.0;
      double codelengthAfter = 0.0;

      for (InfoNode& module : m_root.infomapTree()) {
        if (!module.isLeaf() && module.firstChild->isLeafModule()) {
          codelengthBefore += module.codelength;
          auto numLeafs = 0;
          for (auto& leafModule : module) {
            numLeafs += leafModule.childDegree();
          }

          std::vector<unsigned int> modules(numLeafs);
          std::vector<InfoNode*> leafs(numLeafs);

          auto i = 0;
          for (auto it = module.begin_tree(); !it.isEnd(); ++it) {
            if (it->isLeaf()) {
              modules[i] = it.moduleIndex();
              leafs[i] = it.current();
              ++i;
            }
          }

          module.replaceChildrenWithGrandChildren();

          auto& subInfomap = getSubInfomap(module).initNetwork(module);

          subInfomap.initPartition(modules);

          // Run two-level partition + find hierarchically super modules (skip recursion)
          subInfomap.setOnlySuperModules(true).run();

          // Collect sub Infomap modules
          i = 0;
          for (auto& subLeafPtr : subInfomap.leafNodes()) {
            modules[i] = subLeafPtr->index;
            ++i;
          }

          // Create new sub modules
          std::vector<InfoNode*> subModules(numLeafs, nullptr);
          module.releaseChildren();

          for (auto j = 0; j < numLeafs; ++j) {
            InfoNode* leaf = leafs[j];
            unsigned int moduleIndex = modules[j];
            if (subModules[moduleIndex] == nullptr) {
              subModules[moduleIndex] = new InfoNode(subInfomap.leafNodes()[j]->parent->data);
              subModules[moduleIndex]->index = moduleIndex;
              module.addChild(subModules[moduleIndex]);
            }
            subModules[moduleIndex]->addChild(leaf);
          }
          module.disposeInfomap();
          module.codelength = calcCodelength(module);
          codelengthAfter += module.codelength;
        }
      }

      if (isMainInfomap())
        Log::setSilent(isSilent);

      const double diffCodelength = codelengthBefore - codelengthAfter;
      Log() << "done! Codelength improvement " << (diffCodelength / codelengthBefore) * 100 << "% to codelength " << codelengthAfter << "\n";
    }

    recursivePartition();
    return;
  }

  partition();

  if (numTopModules() == 1 || numTopModules() == numLeafNodes()) {
    Log(1) << "Trivial partition, skip search for hierarchical solution.\n";
    return;
  }

  if (numTopModules() > preferredNumberOfModules) {
    findHierarchicalSuperModules();
  }

  if (onlySuperModules) {
    removeSubModules(true);
    m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
    return;
  }

  if (fastHierarchicalSolution >= 2) {
    // FIXME Calculate individual module codelengths and store on the modules?
    return;
  }

  if (fastHierarchicalSolution == 0) {
    removeSubModules(true);
    m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
  }

  recursivePartition();
}

void InfomapBase::partition()
{
  auto oldPrecision = Log::precision();
  Log(0, 0) << "Two-level compression: " << std::setprecision(2) << std::flush;
  Log(1) << "Trying to find modular structure... \n";
  double initialCodelength = m_oneLevelCodelength;
  double oldCodelength = initialCodelength;

  m_tuneIterationIndex = 0;
  findTopModulesRepeatedly(levelAggregationLimit);

  double newCodelength = getCodelength();
  double compression = oldCodelength < 1e-16 ? 0.0 : (oldCodelength - newCodelength) / oldCodelength;
  Log(0, 0) << (compression * 100) << "% " << std::flush;
  oldCodelength = newCodelength;

  bool doFineTune = true;
  bool coarseTuned = false;
  while (numTopModules() > 1 && (m_tuneIterationIndex + 1) != tuneIterationLimit) {
    ++m_tuneIterationIndex;
    if (doFineTune) {
      Log(3) << "\n";
      Log(1, 3) << "Fine tune... " << std::flush;
      unsigned int numEffectiveLoops = fineTune();
      if (numEffectiveLoops > 0) {
        Log(1, 3) << " -> done in " << numEffectiveLoops << " effective loops to codelength " << *this << " in " << numTopModules() << " modules.\n";
        // Continue to merge fine-tuned modules
        findTopModulesRepeatedly(levelAggregationLimit);
      } else {
        Log(1, 3) << "no improvement.\n";
      }
    } else {
      coarseTuned = true;
      if (!noCoarseTune) {
        Log(3) << "\n";
        Log(1, 3) << "Coarse tune... " << std::flush;
        unsigned int numEffectiveLoops = coarseTune();
        if (numEffectiveLoops > 0) {
          Log(1, 3) << "done in " << numEffectiveLoops << " effective loops to codelength " << *this << " in " << numTopModules() << " modules.\n";
          // Continue to merge fine-tuned modules
          findTopModulesRepeatedly(levelAggregationLimit);
        } else {
          Log(1, 3) << "no improvement.\n";
        }
      }
    }
    newCodelength = getCodelength();
    compression = oldCodelength < 1e-16 ? 0.0 : (oldCodelength - newCodelength) / oldCodelength;
    bool isImprovement = newCodelength <= oldCodelength - minimumCodelengthImprovement && newCodelength < oldCodelength - initialCodelength * minimumRelativeTuneIterationImprovement;
    if (!isImprovement) {
      // Make sure coarse-tuning have been tried
      if (coarseTuned)
        break;
    } else {
      oldCodelength = newCodelength;
    }
    Log(0, 0) << (compression * 100) << "% " << std::flush;
    doFineTune = !doFineTune;
  }

  Log(0, 0) << std::setprecision(oldPrecision) << '\n';
  Log() << "Partitioned to codelength " << *this << " in " << numTopModules();
  if (m_numNonTrivialTopModules != numTopModules())
    Log() << " (" << m_numNonTrivialTopModules << " non-trivial)";
  Log() << " modules.\n";

  if (!preferModularSolution && preferredNumberOfModules == 0 && haveNonTrivialModules() && getCodelength() > getOneLevelCodelength()) {
    Log() << "Worse codelength than one-level codelength, putting all nodes in one module... ";

    // Create new single module between modules and root
    auto& module = root().replaceChildrenWithOneNode();
    module.data = m_root.data;
    module.index = 0;
    for (auto& node : module) {
      node.index = 0;
    }
    module.codelength = getOneLevelCodelength();
    m_hierarchicalCodelength = getOneLevelCodelength();

  } else {
    // Set consolidated cluster index on nodes and modules
    for (auto clusterIt = m_root.begin_tree(); !clusterIt.isEnd(); ++clusterIt) {
      clusterIt->index = clusterIt.moduleIndex();
    }

    // Calculate individual module codelengths and store on the modules
    calcCodelengthOnTree(root(), false);
    m_root.codelength = getIndexCodelength();
    m_hierarchicalCodelength = getCodelength();
  }
  m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);
}

void InfomapBase::restoreHardPartition()
{
  // Collect all leaf nodes in a separate sequence to be able to iterate independent of tree structure
  std::vector<InfoNode*> leafNodes(numLeafNodes());
  unsigned int leafIndex = 0;
  for (InfoNode& node : m_root.infomapTree()) {
    if (node.isLeaf()) {
      leafNodes[leafIndex] = &node;
      ++leafIndex;
    }
  }
  // Expand all children
  unsigned int numExpandedChildren = 0;
  unsigned int numExpandedNodes = 0;
  for (InfoNode* node : leafNodes) {
    ++numExpandedNodes;
    numExpandedChildren += node->expandChildren();
    node->replaceWithChildren();
  }

  // Swap back original leaf nodes
  m_originalLeafNodes.swap(m_leafNodes);

  Log(1) << "Expanded " << numExpandedNodes << " hard modules to " << numExpandedChildren << " original nodes.\n";
}

// ===================================================
// Run: Init: *
// ===================================================

void InfomapBase::initEnterExitFlow()
{
  // Not done in Bayesian
  // TODO: Skip this, always add enter/exit/tele flow from flow calculator
  // Calculate enter/exit
  double alpha = teleportationProbability;
  double beta = 1.0 - alpha;

  for (auto* n : m_leafNodes) {
    n->data.enterFlow = n->data.exitFlow = 0.0;
  }
  if (!isUndirectedClustering()) {
    for (auto* n : m_leafNodes) {
      auto& node = *n;
      node.data.teleportFlow = alpha * node.data.flow;
      node.data.teleportSourceFlow = node.data.flow;
      if (node.isDangling()) {
        node.data.teleportFlow = node.data.flow;
        node.data.danglingFlow = node.data.flow;
        m_sumDanglingFlow += node.data.flow;
      }
    }
    for (auto* n : m_leafNodes) {
      auto& node = *n;
      for (InfoEdge* e : node.outEdges()) {
        InfoEdge& edge = *e;
        // Self-links not included here, should not add to enter and exit flow in its enclosing module
        edge.source->data.exitFlow += edge.data.flow;
        edge.target->data.enterFlow += edge.data.flow;
      }
      if (recordedTeleportation) {
        // Don't let self-teleportation add to the enter/exit flow (i.e. multiply with (1.0 - node.data.teleportWeight))
        node.data.exitFlow += (alpha * node.data.flow + beta * node.data.danglingFlow) * (1.0 - node.data.teleportWeight);
        node.data.enterFlow += (alpha * (1.0 - node.data.flow) + beta * (m_sumDanglingFlow - node.data.danglingFlow)) * node.data.teleportWeight;
      }
    }
  } else {
    for (auto* n : m_leafNodes) {
      auto& node = *n;
      node.data.teleportFlow = alpha * node.data.flow;
      node.data.teleportSourceFlow = node.data.flow;
      if (node.isDangling()) {
        node.data.teleportFlow = node.data.flow;
        node.data.danglingFlow = node.data.flow;
      }
      for (InfoEdge* e : node.outEdges()) {
        InfoEdge& edge = *e;
        // Self-links not included here, should not add to enter and exit flow in its enclosing module
        double halfFlow = edge.data.flow / 2;
        edge.source->data.exitFlow += halfFlow;
        edge.target->data.exitFlow += halfFlow;
        edge.source->data.enterFlow += halfFlow;
        edge.target->data.enterFlow += halfFlow;
      }
    }
  }
}

// Aggregate node and enter/exit flow to all tree nodes
void InfomapBase::aggregateFlowValuesFromLeafToRoot()
{
  // Aggregate flow from leaf nodes to root node
  unsigned int numLevels = 0;
  root().data.flow = 0.0;
  for (auto it = root().begin_post_depth_first(); !it.isEnd(); ++it) {
    auto& node = *it;
    if (!node.isRoot())
      node.parent->data += node.data;
    // Don't aggregate enter and exit flow
    if (!node.isLeaf()) {
      node.index = it.depth(); // Use index to store the depth on modules
      node.data.enterFlow = 0.0;
      node.data.exitFlow = 0.0;
    } else
      numLevels = std::max(numLevels, it.depth());
  }

  if (std::abs(root().data.flow - 1.0) > 1e-10) {
    Log() << "Warning, aggregated flow is not exactly 1.0, but " << std::setprecision(10) << root().data.flow << std::setprecision(9) << ".\n";
  }

  // Aggregate enter and exit flow between modules
  for (auto& leafNode : m_leafNodes) {
    auto& leafNodeSource = *leafNode;
    for (InfoEdge* e : leafNodeSource.outEdges()) {
      auto& edge = *e;
      auto& leafNodeTarget = edge.target;
      double linkFlow = edge.data.flow;
      double halfFlow = linkFlow / 2;

      InfoNode* node1 = leafNodeSource.parent;
      InfoNode* node2 = leafNodeTarget->parent;

      if (node1 == node2)
        continue;

      // First aggregate link flow until equal depth
      while (node1->index > node2->index) {
        if (isUndirectedClustering()) {
          node1->data.exitFlow += halfFlow;
          node1->data.enterFlow += halfFlow;
        } else {
          node1->data.exitFlow += linkFlow;
        }
        node1 = node1->parent;
      }
      while (node2->index > node1->index) {
        if (isUndirectedClustering()) {
          node2->data.enterFlow += halfFlow;
          node2->data.exitFlow += halfFlow;
        } else {
          node2->data.enterFlow += linkFlow;
        }
        node2 = node2->parent;
      }

      // Then aggregate link flow until equal parent
      while (node1 != node2) {
        if (isUndirectedClustering()) {
          node1->data.exitFlow += halfFlow;
          node1->data.enterFlow += halfFlow;
          node2->data.enterFlow += halfFlow;
          node2->data.exitFlow += halfFlow;
        } else {
          node1->data.exitFlow += linkFlow;
          node2->data.enterFlow += linkFlow;
        }
        node1 = node1->parent;
        node2 = node2->parent;
      }
    }
  }
  if (recordedTeleportation) {
    Log() << "\n\nAggregating enter/exit flow for recorded teleportation, sum teleFlow: " << m_root.data.teleportFlow << "\n";

    for (auto& node : m_root.infomapTree()) {
      if (!node.isLeaf()) {
        // Don't code self-teleportation

        node.data.enterFlow += (m_root.data.teleportFlow - node.data.teleportFlow) * node.data.teleportWeight;
        node.data.exitFlow += node.data.teleportFlow * (1.0 - node.data.teleportWeight);
      }
    }
  }
}

double InfomapBase::calcCodelengthOnTree(InfoNode& root, bool includeRoot) const
{
  double totalCodelength = 0.0;
  // Calculate enter/exit flow (assuming 0 by default)
  for (auto& node : root.infomapTree()) {
    if (node.isLeaf() || (!includeRoot && node.isRoot()))
      continue;
    node.codelength = calcCodelength(node);
    totalCodelength += node.codelength;
  }
  return totalCodelength;
}

// ===================================================
// Run: Partition: *
// ===================================================

void InfomapBase::setActiveNetworkFromChildrenOfRoot()
{
  m_moduleNodes.resize(m_root.childDegree());
  unsigned int i = 0;
  for (auto& node : m_root)
    m_moduleNodes[i++] = &node;
  m_activeNetwork = &m_moduleNodes;
}

void InfomapBase::findTopModulesRepeatedly(unsigned int maxLevels)
{
  Log(1, 2) << "Iteration " << (m_tuneIterationIndex + 1) << ", moving ";
  Log(3) << "\nIteration " << (m_tuneIterationIndex + 1) << ":\n";
  m_aggregationLevel = 0;
  unsigned int numLevelsConsolidated = numLevels() - 1;
  if (maxLevels == 0)
    maxLevels = std::numeric_limits<unsigned int>::max();

  std::string initialCodelength;

  // Reapply core algorithm on modular network, replacing modules with super modules
  while (numTopModules() > 1 && numLevelsConsolidated != maxLevels) {
    if (haveModules())
      setActiveNetworkFromChildrenOfRoot();
    else
      setActiveNetworkFromLeafs();
    initPartition();
    if (m_aggregationLevel == 0) {
      initialCodelength = io::Str() << "" << *this;
    }

    Log(1, 2) << activeNetwork().size() << "*" << std::flush;
    Log(3) << "Level " << (numLevelsConsolidated + 1) << " (codelength: " << *this << "): Moving " << activeNetwork().size() << " nodes... " << std::flush;
    // Core loop, merging modules
    unsigned int numOptimizationLoops = optimizeActiveNetwork();

    Log(1, 2) << numOptimizationLoops << ", " << std::flush;
    Log(3) << "done! -> codelength " << *this << " in " << numActiveModules() << " modules.\n";

    // If no improvement, revert codelength terms to the last consolidated state
    if (haveModules() && restoreConsolidatedOptimizationPointIfNoImprovement()) {
      Log(3) << "-> Restored to codelength " << *this << " in " << numTopModules() << " modules.\n";
      break;
    }

    // Consolidate modules
    bool replaceExistingModules = haveModules();
    consolidateModules(replaceExistingModules);
    ++numLevelsConsolidated;
    ++m_aggregationLevel;
  }
  m_aggregationLevel = 0;

  m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();

  Log(1, 2) << (m_isCoarseTune ? "modules" : "nodes") << "*loops from codelength " << initialCodelength << " to codelength " << *this << " in " << numTopModules() << " modules. (" << m_numNonTrivialTopModules << " non-trivial modules)\n";
}

unsigned int InfomapBase::fineTune()
{
  if (numLevels() != 2)
    throw std::logic_error("InfomapBase::fineTune() called but numLevels != 2");

  setActiveNetworkFromLeafs();
  initPartition();

  // Make sure module nodes have correct index //TODO: Assume always correct here?
  unsigned int moduleIndex = 0;
  for (auto& module : m_root) {
    module.index = moduleIndex++;
  }

  // Put leaf modules in existing modules
  std::vector<unsigned int> modules(numLeafNodes());
  for (unsigned int i = 0; i < numLeafNodes(); ++i) {
    modules[i] = m_leafNodes[i]->parent->index;
  }
  moveActiveNodesToPredefinedModules(modules);

  Log(3) << " -> moved to codelength " << *this << " in " << numActiveModules() << " existing modules. Try tuning...\n";

  // Continue to optimize from there to tune leaf nodes
  unsigned int numEffectiveLoops = optimizeActiveNetwork();
  if (numEffectiveLoops == 0) {
    restoreConsolidatedOptimizationPointIfNoImprovement();
    Log(4) << "Fine-tune didn't improve solution, restoring last.\n";
  } else {
    // Delete existing modules and consolidate fine-tuned modules
    root().replaceChildrenWithGrandChildren();
    consolidateModules(false);
    Log(4) << "Fine-tune done in " << numEffectiveLoops << " effective loops to codelength " << *this << " in " << numTopModules() << " modules.\n";
  }

  return numEffectiveLoops;
}

unsigned int InfomapBase::coarseTune()
{
  if (numLevels() != 2)
    throw std::logic_error("InfomapBase::coarseTune() called but numLevels != 2");

  Log(4) << "Coarse-tune...\nPartition each module in sub-modules for coarse tune...\n";

  bool isSilent = false;
  if (isMainInfomap()) {
    isSilent = Log::isSilent();
    Log::setSilent(true);
  }

  unsigned int moduleIndexOffset = 0;
  for (auto& node : m_root) {
    // Don't search for sub-modules in too small modules
    if (node.childDegree() < 2) {
      for (auto& child : node) {
        child.index = moduleIndexOffset;
      }
      moduleIndexOffset += 1;
      continue;
    } else {
      InfomapBase& subInfomap = getSubInfomap(node)
                                    .setTwoLevel(true)
                                    .setTuneIterationLimit(1);
      subInfomap.initNetwork(node).run();

      auto originalLeafIt = node.begin_child();
      for (auto& subLeafPtr : subInfomap.leafNodes()) {
        InfoNode& subLeaf = *subLeafPtr;
        originalLeafIt->index = subLeaf.index + moduleIndexOffset;
        ++originalLeafIt;
      }
      moduleIndexOffset += subInfomap.numTopModules();

      node.disposeInfomap();
    }
  }

  if (isMainInfomap())
    Log::setSilent(isSilent);

  Log(4) << "Move leaf nodes to " << moduleIndexOffset << " sub-modules... \n";
  // Put leaf modules in the calculated sub-modules
  std::vector<unsigned int> subModules(numLeafNodes());
  for (unsigned int i = 0; i < numLeafNodes(); ++i) {
    subModules[i] = m_leafNodes[i]->index;
  }

  setActiveNetworkFromLeafs();
  initPartition();
  moveActiveNodesToPredefinedModules(subModules);

  Log(4) << "Moved to " << numActiveModules() << " modules...\n";

  // Replace existing modules but store that structure on the sub-modules
  consolidateModules(true);

  Log(4) << "Consolidated " << numTopModules() << " sub-modules to codelength " << *this << '\n';

  Log(4) << "Move sub-modules to former modular structure... \n";
  // Put sub-modules in the former modular structure
  std::vector<unsigned int> modules(numTopModules());
  unsigned int iModule = 0;
  for (auto& subModule : m_root) {
    modules[iModule++] = subModule.index;
  }

  setActiveNetworkFromChildrenOfRoot();
  initPartition();
  // Move sub-modules to former modular structure
  moveActiveNodesToPredefinedModules(modules);

  Log(4) << "Tune sub-modules from codelength " << *this << " in " << numActiveModules() << " modules... \n";
  // Continue to optimize from there to tune sub-modules
  unsigned int numEffectiveLoops = optimizeActiveNetwork();
  Log(4) << "Tuned sub-modules in " << numEffectiveLoops << " effective loops to codelength " << *this << " in " << numActiveModules() << " modules.\n";
  consolidateModules(true);
  Log(4) << "Consolidated tuned sub-modules to codelength " << *this << " in " << numTopModules() << " modules.\n";
  return numEffectiveLoops;
}

unsigned int InfomapBase::calculateNumNonTrivialTopModules() const
{
  if (root().childDegree() == 1)
    return 0;
  unsigned int numNonTrivialTopModules = 0;
  for (auto& module : m_root) {
    if (module.childDegree() > 1)
      ++numNonTrivialTopModules;
  }
  return numNonTrivialTopModules;
}

unsigned int InfomapBase::calculateMaxDepth() const
{
  unsigned int maxDepth = 0;
  for (auto it(m_root.begin_tree()); !it.isEnd(); ++it) {
    if (it->isLeaf()) {
      maxDepth = std::max(maxDepth, it.depth());
    }
  }
  return maxDepth;
}

unsigned int InfomapBase::findHierarchicalSuperModulesFast(unsigned int superLevelLimit)
{
  if (superLevelLimit == 0)
    return 0;

  unsigned int numLevelsCreated = 0;
  double oldIndexLength = getIndexCodelength();
  double hierarchicalCodelength = getCodelength();
  double workingHierarchicalCodelength = hierarchicalCodelength;

  Log(1) << "\nFind hierarchical super modules iteratively...\n";
  Log(0, 0) << "Fast super-level compression: " << std::setprecision(2) << std::flush;

  // Add index codebooks as long as the code gets shorter
  do {
    Log(1) << "Iteration " << numLevelsCreated + 1 << ", finding super modules fast to " << numTopModules() << (haveModules() ? " modules" : " nodes") << "... \n";

    if (haveModules()) {
      setActiveNetworkFromChildrenOfRoot();
      transformNodeFlowToEnterFlow(m_root);
      initSuperNetwork();
    } else {
      setActiveNetworkFromLeafs();
    }

    initPartition();

    unsigned int numEffectiveLoops = optimizeActiveNetwork();

    double codelength = getCodelength();
    double indexCodelength = getIndexCodelength();
    unsigned int numSuperModules = numActiveModules();
    bool trivialSolution = numSuperModules == 1 || numSuperModules == activeNetwork().size();

    bool acceptSolution = !trivialSolution && codelength < oldIndexLength - minimumCodelengthImprovement;
    // Force at least one modular level!
    bool acceptByForce = !acceptSolution && !haveModules();
    if (acceptByForce)
      acceptSolution = true;

    workingHierarchicalCodelength += codelength - oldIndexLength;

    Log(0, 0) << hideIf(!acceptSolution) << ((hierarchicalCodelength - workingHierarchicalCodelength) / hierarchicalCodelength * 100) << "% " << std::flush;
    Log(1) << "  -> Found " << numActiveModules() << " super modules in " << numEffectiveLoops << " effective loops with hierarchical codelength " << indexCodelength << " + " << (workingHierarchicalCodelength - indexCodelength) << " = " << workingHierarchicalCodelength << (acceptSolution ? "\n" : ", discarding the solution.\n") << std::flush;
    Log(1) << oldIndexLength << " -> " << *this << "\n";

    if (!acceptSolution) {
      restoreConsolidatedOptimizationPointIfNoImprovement(true);
      break;
    }

    // Consolidate the dynamic modules without replacing any existing ones.
    consolidateModules(false);

    hierarchicalCodelength = workingHierarchicalCodelength;
    oldIndexLength = indexCodelength;

    ++numLevelsCreated;

    m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();

  } while (m_numNonTrivialTopModules > 1 && numLevelsCreated != superLevelLimit);

  // Restore the temporary transformation of flow to enter flow on super modules
  resetFlowOnModules();

  Log(0, 0) << "to codelength " << io::toPrecision(hierarchicalCodelength) << " in " << numTopModules() << " top modules.\n";
  Log(1) << "Finding super modules done! Added " << numLevelsCreated << " levels with hierarchical codelength " << io::toPrecision(hierarchicalCodelength) << " in " << numTopModules() << " top modules.\n";

  // Calculate and store hierarchical codelengths
  m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);

  return numLevelsCreated;
}

unsigned int InfomapBase::findHierarchicalSuperModules(unsigned int superLevelLimit)
{
  if (superLevelLimit == 0)
    return 0;

  unsigned int numLevelsCreated = 0;
  double oldIndexLength = getIndexCodelength();
  double hierarchicalCodelength = getCodelength();
  double workingHierarchicalCodelength = hierarchicalCodelength;

  if (!haveModules())
    throw std::logic_error("Trying to find hierarchical super modules without any modules");

  Log(1) << "\nFind hierarchical super modules iteratively...\n";
  Log(0, 0) << "Super-level compression: " << std::setprecision(2) << std::flush;

  // Add index codebooks as long as the code gets shorter
  do {
    Log(1) << "Iteration " << numLevelsCreated + 1 << ", finding super modules to " << numTopModules() << " modules... \n";
    bool isSilent = false;
    if (isMainInfomap()) {
      isSilent = Log::isSilent();
      Log::setSilent(true);
    }
    InfoNode tmp(m_root.data); // Temporary owner for the super Infomap instance
    auto& superInfomap = getSuperInfomap(tmp)
                             .setTwoLevel(true);
    superInfomap.initNetwork(m_root, true); //.initSuperNetwork();
    superInfomap.run();
    if (isMainInfomap()) {
      Log::setSilent(isSilent);
    }
    double superCodelength = superInfomap.getCodelength();
    double superIndexCodelength = superInfomap.getIndexCodelength();

    unsigned int numSuperModules = superInfomap.numTopModules();
    bool trivialSolution = numSuperModules == 1 || numSuperModules == numTopModules();

    if (trivialSolution) {
      Log(1) << "failed to find non-trivial super modules.\n";
      break;
    }

    bool improvedCodelength = superCodelength < oldIndexLength - minimumCodelengthImprovement;
    if (!improvedCodelength) {
      Log(1) << "two-level index codebook not improved over one-level.\n";
      break;
    }

    workingHierarchicalCodelength += superCodelength - oldIndexLength;

    Log(0, 0) << ((hierarchicalCodelength - workingHierarchicalCodelength)
                  / hierarchicalCodelength * 100)
              << "% " << std::flush;
    Log(1) << "  -> Found " << numSuperModules << " super modules in with hierarchical codelength " << superIndexCodelength << " + " << (superCodelength - superIndexCodelength) << " + " << (workingHierarchicalCodelength - superCodelength) << " = " << workingHierarchicalCodelength << '\n';

    // Merge current top modules to two-level solution found in the super Infomap instance.
    std::vector<unsigned int> modules(numTopModules());
    unsigned int iModule = 0;
    for (InfoNode* n : superInfomap.leafNodes()) {
      modules[iModule++] = n->index;
    }

    setActiveNetworkFromChildrenOfRoot();
    initPartition();
    // Move top modules to super modules
    moveActiveNodesToPredefinedModules(modules);

    // Consolidate the dynamic modules without replacing any existing ones.
    consolidateModules(false);

    Log(1) << "Consolidated " << numTopModules() << " super modules. Index codelength: " << oldIndexLength << " -> " << *this << "\n";

    hierarchicalCodelength = workingHierarchicalCodelength;
    oldIndexLength = superIndexCodelength;

    ++numLevelsCreated;

    m_numNonTrivialTopModules = calculateNumNonTrivialTopModules();

  } while (m_numNonTrivialTopModules > 1 && numLevelsCreated != superLevelLimit);

  // Restore the temporary transformation of flow to enter flow on super modules
  resetFlowOnModules();

  Log(0, 0) << "to codelength " << io::toPrecision(hierarchicalCodelength) << " in " << numTopModules() << " top modules.\n";
  Log(1) << "Finding super modules done! Added " << numLevelsCreated << " levels with hierarchical codelength " << io::toPrecision(hierarchicalCodelength) << " in " << numTopModules() << " top modules.\n";

  // Calculate and store hierarchical codelengths
  m_hierarchicalCodelength = calcCodelengthOnTree(root(), true);

  return numLevelsCreated;
}

void InfomapBase::transformNodeFlowToEnterFlow(InfoNode& parent)
{
  double sumFlow = 0.0;
  for (auto& module : m_root) {
    module.data.flow = module.data.enterFlow;
    sumFlow += module.data.flow;
  }
  parent.data.flow = sumFlow;
}

void InfomapBase::resetFlowOnModules()
{
  // Reset flow on all modules
  for (auto& module : m_root.infomapTree()) {
    if (!module.isLeaf())
      module.data.flow = 0.0;
  }
  // Aggregate flow from leaf nodes up in the tree
  for (InfoNode* n : m_leafNodes) {
    double leafNodeFlow = n->data.flow;
    InfoNode* module = n->parent;
    do {
      module->data.flow += leafNodeFlow;
      module = module->parent;
    } while (module != nullptr);
  }
}

unsigned int InfomapBase::removeModules()
{
  unsigned int numLevelsDeleted = 0;
  while (numLevels() > 1) {
    m_root.replaceChildrenWithGrandChildren();
    ++numLevelsDeleted;
  }
  if (numLevelsDeleted > 0) {
    Log(2) << "Removed " << numLevelsDeleted << " levels of modules.\n";
  }
  return numLevelsDeleted;
}

unsigned int InfomapBase::removeSubModules(bool recalculateCodelengthOnTree)
{
  unsigned int numLevelsDeleted = 0;
  while (numLevels() > 2) {
    for (InfoNode& module : m_root)
      module.replaceChildrenWithGrandChildren();
    ++numLevelsDeleted;
  }
  if (numLevelsDeleted > 0) {
    Log(2) << "Removed " << numLevelsDeleted << " levels of sub-modules.\n";
    if (recalculateCodelengthOnTree)
      calcCodelengthOnTree(root(), false);
  }
  return numLevelsDeleted;
}

unsigned int InfomapBase::recursivePartition()
{
  double indexCodelength = getIndexCodelength();
  double hierarchicalCodelength = m_hierarchicalCodelength;

  Log(0, 0) << "\nRecursive sub-structure compression: " << std::flush;

  PartitionQueue partitionQueue;
  if (fastHierarchicalSolution > 0) {
    queueLeafModules(partitionQueue);
    Log(1) << "\nFind sub modules recursively from " << partitionQueue.size() << " sub modules on level " << numLevels() - 2;
  } else {
    queueTopModules(partitionQueue);
    Log(1) << "\nFind sub modules recursively from " << partitionQueue.size() << " top modules";
    double moduleCodelength = 0.0;
    for (auto& module : m_root)
      moduleCodelength += module.codelength;
    hierarchicalCodelength = indexCodelength + moduleCodelength;
  }
  Log(1) << " with codelength: " << indexCodelength << " + " << (hierarchicalCodelength - indexCodelength) << " = " << io::toPrecision(hierarchicalCodelength) << '\n';

  double sumConsolidatedCodelength = hierarchicalCodelength - partitionQueue.moduleCodelength;

  bool isSilent = false;
  if (isMainInfomap()) {
    isSilent = Log::isSilent();
  }

  while (partitionQueue.size() > 0) {
    Log(1) << "Level " << partitionQueue.level << ": " << (partitionQueue.flow * 100) << "% of the flow in " << partitionQueue.size() << " modules. Partitioning... " << std::setprecision(6) << std::flush;

    if (isMainInfomap())
      Log::setSilent(true);

    // Partition all modules in the queue and fill up the next level queue
    PartitionQueue nextLevelQueue;
    processPartitionQueue(partitionQueue, nextLevelQueue);

    if (isMainInfomap())
      Log::setSilent(isSilent);

    double leftToImprove = partitionQueue.moduleCodelength;
    sumConsolidatedCodelength += partitionQueue.indexCodelength + partitionQueue.leafCodelength;
    double limitCodelength = sumConsolidatedCodelength + leftToImprove;

    Log(0, 0) << ((hierarchicalCodelength - limitCodelength) / hierarchicalCodelength) * 100 << "% " << std::flush;
    Log(1) << "done! Codelength: " << partitionQueue.indexCodelength << " + " << partitionQueue.leafCodelength << " (+ " << leftToImprove << " left to improve)"
           << " -> limit: " << io::toPrecision(limitCodelength) << " bits.\n";

    hierarchicalCodelength = limitCodelength;

    partitionQueue.swap(nextLevelQueue);
  }

  // Store resulting hierarchical codelength
  m_hierarchicalCodelength = hierarchicalCodelength;

  Log(0, 0) << ". Found " << partitionQueue.level << " levels with codelength " << io::toPrecision(hierarchicalCodelength) << "\n";
  Log(1) << "  -> Found " << partitionQueue.level << " levels with codelength " << io::toPrecision(hierarchicalCodelength) << "\n";

  return partitionQueue.level;
}

void InfomapBase::queueTopModules(PartitionQueue& partitionQueue)
{
  // Add modules to partition queue
  unsigned int numNonTrivialModules = 0;
  partitionQueue.resize(numTopModules());
  double sumFlow = 0.0;
  double sumNonTrivialFlow = 0.0;
  double sumModuleCodelength = 0.0;
  unsigned int moduleIndex = 0;
  for (auto& module : m_root) {
    partitionQueue[moduleIndex] = &module;
    sumFlow += module.data.flow;
    sumModuleCodelength += module.codelength;
    if (module.childDegree() > 1) {
      ++numNonTrivialModules;
      sumNonTrivialFlow += module.data.flow;
    }
    ++moduleIndex;
  }
  partitionQueue.flow = sumFlow;
  partitionQueue.numNonTrivialModules = numNonTrivialModules;
  partitionQueue.nonTrivialFlow = sumNonTrivialFlow;
  partitionQueue.indexCodelength = getIndexCodelength();
  partitionQueue.moduleCodelength = sumModuleCodelength;
  partitionQueue.level = 1;
}

void InfomapBase::queueLeafModules(PartitionQueue& partitionQueue)
{
  unsigned int numLeafModules = 0;
  for (auto& node : m_root.infomapTree()) {
    if (node.isLeafModule())
      ++numLeafModules;
  }

  // Add modules to partition queue
  partitionQueue.resize(numLeafModules);
  unsigned int numNonTrivialModules = 0;
  double sumFlow = 0.0;
  double sumNonTrivialFlow = 0.0;
  double sumModuleCodelength = 0.0;
  unsigned int moduleIndex = 0;
  unsigned int maxDepth = 0;
  for (auto it = m_root.begin_tree(); !it.isEnd(); ++it) {
    if (it->isLeafModule()) {
      auto& module = *it;
      partitionQueue[moduleIndex] = it.current();
      sumFlow += module.data.flow;
      sumModuleCodelength += module.codelength;
      if (module.childDegree() > 1) {
        ++numNonTrivialModules;
        sumNonTrivialFlow += module.data.flow;
      }
      maxDepth = std::max(maxDepth, it.depth());
      ++moduleIndex;
    }
  }
  partitionQueue.flow = sumFlow;
  partitionQueue.numNonTrivialModules = numNonTrivialModules;
  partitionQueue.nonTrivialFlow = sumNonTrivialFlow;
  partitionQueue.indexCodelength = getIndexCodelength();
  partitionQueue.moduleCodelength = sumModuleCodelength;
  partitionQueue.level = maxDepth;
}

bool InfomapBase::processPartitionQueue(PartitionQueue& queue, PartitionQueue& nextLevelQueue) const
{
  PartitionQueue::size_t numModules = queue.size();
  std::vector<double> indexCodelengths(numModules, 0.0);
  std::vector<double> moduleCodelengths(numModules, 0.0);
  std::vector<double> leafCodelengths(numModules, 0.0);
  std::vector<PartitionQueue> subQueues(numModules);

#pragma omp parallel for schedule(dynamic)
  for (PartitionQueue::size_t moduleIndex = 0; moduleIndex < numModules; ++moduleIndex) {
    InfoNode& module = *queue[moduleIndex];

    module.codelength = calcCodelength(module);
    // Delete former sub-structure if exists
    if (module.disposeInfomap())
      module.codelength = calcCodelength(module);

    // If only trivial substructure is to be found, no need to create infomap instance to find sub-module structures.
    if (module.childDegree() <= 2) {
      module.codelength = calcCodelength(module);
      leafCodelengths[moduleIndex] = module.codelength;
      continue;
    }

    double oldModuleCodelength = module.codelength;
    PartitionQueue& subQueue = subQueues[moduleIndex];
    subQueue.level = queue.level + 1;

    auto& subInfomap = getSubInfomap(module)
                           .initNetwork(module);
    // Run two-level partition + find hierarchically super modules (skip recursion)
    subInfomap.setOnlySuperModules(true).run();

    double subCodelength = subInfomap.getHierarchicalCodelength();
    double subIndexCodelength = subInfomap.root().codelength;
    double subModuleCodelength = subCodelength - subIndexCodelength;
    InfoNode& subRoot = *module.getInfomapRoot();
    unsigned int numSubModules = subRoot.childDegree();
    bool trivialSubPartition = numSubModules == 1 || numSubModules == module.childDegree();
    bool improvedCodelength = subCodelength < oldModuleCodelength - minimumCodelengthImprovement;

    if (trivialSubPartition || !improvedCodelength) {
      Log(1) << "Disposing unaccepted sub Infomap instance.\n";
      module.disposeInfomap();
      module.codelength = oldModuleCodelength;
      subQueue.skip = true;
      leafCodelengths[moduleIndex] = module.codelength;
    } else {
      // Improvement
      subInfomap.queueTopModules(subQueue);
      indexCodelengths[moduleIndex] = subIndexCodelength;
      moduleCodelengths[moduleIndex] = subModuleCodelength;
    }
  }

  double sumLeafCodelength = 0.0;
  double sumIndexCodelength = 0.0;
  double sumModuleCodelengths = 0.0;
  PartitionQueue::size_t nextLevelSize = 0;
  for (PartitionQueue::size_t moduleIndex = 0; moduleIndex < numModules; ++moduleIndex) {
    nextLevelSize += subQueues[moduleIndex].skip ? 0 : subQueues[moduleIndex].size();
    sumLeafCodelength += leafCodelengths[moduleIndex];
    sumIndexCodelength += indexCodelengths[moduleIndex];
    sumModuleCodelengths += moduleCodelengths[moduleIndex];
  }

  queue.indexCodelength = sumIndexCodelength;
  queue.leafCodelength = sumLeafCodelength;
  queue.moduleCodelength = sumModuleCodelengths;

  // Collect the sub-queues and build the next-level queue
  nextLevelQueue.level = queue.level + 1;
  nextLevelQueue.resize(nextLevelSize);
  PartitionQueue::size_t nextLevelIndex = 0;
  for (PartitionQueue::size_t moduleIndex = 0; moduleIndex < numModules; ++moduleIndex) {
    PartitionQueue& subQueue = subQueues[moduleIndex];
    if (!subQueue.skip) {
      for (PartitionQueue::size_t subIndex = 0; subIndex < subQueue.size(); ++subIndex) {
        nextLevelQueue[nextLevelIndex++] = subQueue[subIndex];
      }
      nextLevelQueue.flow += subQueue.flow;
      nextLevelQueue.nonTrivialFlow += subQueue.nonTrivialFlow;
      nextLevelQueue.numNonTrivialModules += subQueue.numNonTrivialModules;
    }
  }

  return nextLevelSize > 0;
}

// ===================================================
// Write output
// ===================================================

void InfomapBase::writeResult(int trial)
{
  if (noFileOutput)
    return;

  io::Str s;
  s << outDirectory + outName;

  if (printAllTrials && trial != -1 && numTrials > 1) {
    s << "_trial_" << trial;
  }

  std::string basename = s;

  if (printTree) {
    std::string filename = basename + ".tree";

    if (!printStates()) {
      Log() << "Write tree to " << filename << "... ";
      writeTree(filename);
      Log() << "done!\n";
    } else {
      // Write both physical and state level
      Log() << "Write physical tree to " << filename << "... ";
      writeTree(filename);
      Log() << "done!\n";
      std::string filenameStates = basename + "_states.tree";
      Log() << "Write state tree to " << filenameStates << "... ";
      writeTree(filenameStates, true);
      Log() << "done!\n";
    }
  }

  if (printFlowTree) {
    std::string filename = basename + ".ftree";

    if (!printStates()) {
      Log() << "Write flow tree to " << filename << "... ";
      writeFlowTree(filename);
      Log() << "done!\n";
    } else {
      // Write both physical and state level
      Log() << "Write physical flow tree to " << filename << "... ";
      writeFlowTree(filename, false);
      Log() << "done!\n";
      std::string filenameStates = basename + "_states.ftree";
      Log() << "Write state flow tree to " << filenameStates << "... ";
      writeFlowTree(filenameStates, true);
      Log() << "done!\n";
    }
  }

  if (printNewick) {
    std::string filename = basename + ".nwk";

    if (!printStates()) {
      Log() << "Write Newick tree to " << filename << "... ";
      writeNewickTree(filename);
      Log() << "done!\n";
    } else {
      // Write both physical and state level
      Log() << "Write physical Newick tree to " << filename << "... ";
      writeNewickTree(filename, false);
      Log() << "done!\n";
      std::string filenameStates = basename + "_states.nwk";
      Log() << "Write state Newick tree to " << filenameStates << "... ";
      writeNewickTree(filenameStates, true);
      Log() << "done!\n";
    }
  }

  if (printJson) {
    std::string filename = basename + ".json";
    const bool writeLinks = false;

    if (!printStates()) {
      Log() << "Write JSON tree to " << filename << "... ";
      writeJsonTree(filename, false, writeLinks);
      Log() << "done!\n";
    } else {
      // Write both physical and state level
      Log() << "Write physical JSON tree to " << filename << "... ";
      writeJsonTree(filename, false, writeLinks);
      Log() << "done!\n";
      std::string filenameStates = basename + "_states.json";
      Log() << "Write state JSON tree to " << filenameStates << "... ";
      writeJsonTree(filenameStates, true, writeLinks);
      Log() << "done!\n";
    }
  }

  if (printCsv) {
    std::string filename = basename + ".csv";

    if (!printStates()) {
      Log() << "Write CSV tree to " << filename << "... ";
      writeCsvTree(filename);
      Log() << "done!\n";
    } else {
      // Write both physical and state level
      Log() << "Write physical CSV tree to " << filename << "... ";
      writeCsvTree(filename, false);
      Log() << "done!\n";
      std::string filenameStates = basename + "_states.csv";
      Log() << "Write state CSV tree to " << filenameStates << "... ";
      writeCsvTree(filenameStates, true);
      Log() << "done!\n";
    }
  }

  if (printClu) {
    std::string filename = basename + ".clu";
    if (!printStates()) {
      Log() << "Write node modules to " << filename << "... ";
      writeClu(filename, false, cluLevel);
      Log() << "done!\n";
    } else {
      // Write both physical and state level
      Log() << "Write physical node modules to " << filename << "... ";
      writeClu(filename, false, cluLevel);
      Log() << "done!\n";
      std::string filenameStates = basename + "_states.clu";
      Log() << "Write state node modules to " << filenameStates << "... ";
      writeClu(filenameStates, true, cluLevel);
      Log() << "done!\n";
    }
  }
}

unsigned int printPerLevelCodelength(const InfoNode& parent, std::ostream& out)
{
  std::vector<detail::PerLevelStat> perLevelStats;
  aggregatePerLevelCodelength(parent, perLevelStats);

  unsigned int numLevels = perLevelStats.size();

  out << "Per level number of modules:         [";
  for (unsigned int i = 0; i < numLevels - 1; ++i) {
    out << io::padValue(perLevelStats[i].numModules, 11) << ", ";
  }
  out << io::padValue(perLevelStats[numLevels - 1].numModules, 11) << "]";
  unsigned int sumNumModules = 0;
  for (unsigned int i = 0; i < numLevels; ++i)
    sumNumModules += perLevelStats[i].numModules;
  out << " (sum: " << sumNumModules << ")\n";

  out << "Per level number of leaf nodes:      [";
  for (unsigned int i = 0; i < numLevels - 1; ++i) {
    out << io::padValue(perLevelStats[i].numLeafNodes, 11) << ", ";
  }
  out << io::padValue(perLevelStats[numLevels - 1].numLeafNodes, 11) << "]";
  unsigned int sumNumLeafNodes = 0;
  for (unsigned int i = 0; i < numLevels; ++i)
    sumNumLeafNodes += perLevelStats[i].numLeafNodes;
  out << " (sum: " << sumNumLeafNodes << ")\n";

  out << "Per level average child degree:      [";
  double childDegree = perLevelStats[0].numNodes();
  double sumAverageChildDegree = childDegree * childDegree;
  if (numLevels > 1) {
    out << io::padValue(perLevelStats[0].numModules, 11) << ", ";
  }
  for (unsigned int i = 1; i < numLevels - 1; ++i) {
    childDegree = perLevelStats[i].numNodes() * 1.0 / perLevelStats[i - 1].numModules;
    sumAverageChildDegree += childDegree * perLevelStats[i].numNodes();
    out << io::padValue(childDegree, 11) << ", ";
  }
  if (numLevels > 1) {
    childDegree = perLevelStats[numLevels - 1].numNodes() * 1.0 / perLevelStats[numLevels - 2].numModules;
    sumAverageChildDegree += childDegree * perLevelStats[numLevels - 1].numNodes();
  }
  out << io::padValue(childDegree, 11) << "]";
  out << " (average: " << sumAverageChildDegree / (sumNumModules + sumNumLeafNodes) << ")\n";

  out << std::fixed << std::setprecision(9);
  out << "Per level codelength for modules:    [";
  for (unsigned int i = 0; i < numLevels - 1; ++i) {
    out << perLevelStats[i].indexLength << ", ";
  }
  out << perLevelStats[numLevels - 1].indexLength << "]";
  double sumIndexLengths = 0.0;
  for (unsigned int i = 0; i < numLevels; ++i)
    sumIndexLengths += perLevelStats[i].indexLength;
  out << " (sum: " << sumIndexLengths << ")\n";

  out << "Per level codelength for leaf nodes: [";
  for (unsigned int i = 0; i < numLevels - 1; ++i) {
    out << perLevelStats[i].leafLength << ", ";
  }
  out << perLevelStats[numLevels - 1].leafLength << "]";

  double sumLeafLengths = 0.0;
  for (unsigned int i = 0; i < numLevels; ++i)
    sumLeafLengths += perLevelStats[i].leafLength;
  out << " (sum: " << sumLeafLengths << ")\n";

  out << "Per level codelength total:          [";
  for (unsigned int i = 0; i < numLevels - 1; ++i) {
    out << perLevelStats[i].codelength() << ", ";
  }
  out << perLevelStats[numLevels - 1].codelength() << "]";

  double sumCodelengths = 0.0;
  for (unsigned int i = 0; i < numLevels; ++i)
    sumCodelengths += perLevelStats[i].codelength();
  out << " (sum: " << sumCodelengths << ")\n";

  return numLevels;
}

void aggregatePerLevelCodelength(const InfoNode& parent, std::vector<detail::PerLevelStat>& perLevelStat, unsigned int level)
{
  if (perLevelStat.size() < level + 1)
    perLevelStat.resize(level + 1);

  if (parent.firstChild->isLeaf()) {
    perLevelStat[level].numLeafNodes += parent.childDegree();
    perLevelStat[level].leafLength += parent.codelength;
    return;
  }

  perLevelStat[level].numModules += parent.childDegree();
  perLevelStat[level].indexLength += parent.codelength;

  for (auto& module : parent) {
    if (module.getInfomapRoot() != nullptr)
      aggregatePerLevelCodelength(*module.getInfomapRoot(), perLevelStat, level + 1);
    else
      aggregatePerLevelCodelength(module, perLevelStat, level + 1);
  }
}

void InfomapBase::initOptimizer(bool forceNoMemory)
{
  if (haveMetaData()) {
    m_optimizer = std::make_unique<InfomapOptimizer<MetaMapEquation>>();
  } else if (haveMemory() && !forceNoMemory) {
    m_optimizer = std::make_unique<InfomapOptimizer<MemMapEquation>>();
  } else {
    m_optimizer = std::make_unique<InfomapOptimizer<BiasedMapEquation>>();
  }
  m_optimizer->init(this);
}

} // namespace infomap
