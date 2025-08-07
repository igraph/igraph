/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFOMAP_BASE_H_
#define INFOMAP_BASE_H_

#include "InfomapConfig.h"
#include "InfoEdge.h"
#include "InfoNode.h"
#include "InfomapOptimizerBase.h"
#include "iterators/InfomapIterator.h"
#include "../io/ClusterMap.h"
#include "../io/Network.h"
#include "../io/Output.h"
#include "../utils/Log.h"
#include "../utils/Date.h"
#include "../utils/Stopwatch.h"

#include <vector>
#include <deque>
#include <map>
#include <limits>
#include <string>
#include <iostream>
#include <sstream>

namespace infomap {

namespace detail {
  class PartitionQueue;
  struct PerLevelStat;
} // namespace detail

class InfomapBase : public InfomapConfig<InfomapBase> {
  template <typename Objective>
  friend class InfomapOptimizer;

  void initOptimizer(bool forceNoMemory = false);

public:
  using PartitionQueue = detail::PartitionQueue;

  InfomapBase() : InfomapConfig<InfomapBase>() { initOptimizer(); }

  explicit InfomapBase(const Config& conf) : InfomapConfig<InfomapBase>(conf), m_network(conf) { initOptimizer(); }

  explicit InfomapBase(const std::string& flags, bool isCli = false) : InfomapConfig<InfomapBase>(flags, isCli)
  {
    initOptimizer();
    m_network.setConfig(*this);
    m_initialParameters = m_currentParameters = flags;
  }

  virtual ~InfomapBase() = default;

  // ===================================================
  // Iterators
  // ===================================================

  InfomapIterator iterTree(int maxClusterLevel = 1) { return { &root(), maxClusterLevel }; }

  InfomapIteratorPhysical iterTreePhysical(int maxClusterLevel = 1) { return { &root(), maxClusterLevel }; }

  InfomapModuleIterator iterModules(int maxClusterLevel = 1) { return { &root(), maxClusterLevel }; }

  InfomapLeafModuleIterator iterLeafModules(int maxClusterLevel = 1) { return { &root(), maxClusterLevel }; }

  InfomapLeafIterator iterLeafNodes(int maxClusterLevel = 1) { return { &root(), maxClusterLevel }; }

  InfomapLeafIteratorPhysical iterLeafNodesPhysical(int maxClusterLevel = 1) { return { &root(), maxClusterLevel }; }

  InfomapIterator begin(int maxClusterLevel = 1) { return { &root(), maxClusterLevel }; }

  InfomapIterator end() const { return InfomapIterator(nullptr); }

  // ===================================================
  // Getters
  // ===================================================

  Network& network() { return m_network; }
  const Network& network() const { return m_network; }

  InfoNode& root() { return m_root; }
  const InfoNode& root() const { return m_root; }

  unsigned int numLeafNodes() const { return m_leafNodes.size(); }

  const std::vector<InfoNode*>& leafNodes() const { return m_leafNodes; }

  unsigned int numTopModules() const { return m_root.childDegree(); }

  unsigned int numActiveModules() const { return m_optimizer->numActiveModules(); }

  unsigned int numNonTrivialTopModules() const { return m_numNonTrivialTopModules; }

  bool haveModules() const { return !m_root.isLeaf() && !m_root.firstChild->isLeaf(); }

  bool haveNonTrivialModules() const { return numNonTrivialTopModules() > 0; }

  /**
   * Number of node levels below the root in current Infomap instance, 1 if no modules
   */
  unsigned int numLevels() const;

  /**
   * Get maximum depth of any child in the tree, following possible sub Infomap instances
   */
  unsigned int maxTreeDepth() const;

  double getCodelength() const { return m_optimizer->getCodelength(); }

  double getMetaCodelength(bool unweighted = false) const { return m_optimizer->getMetaCodelength(unweighted); }

  double codelength() const { return m_hierarchicalCodelength; }

  const std::vector<double>& codelengths() const { return m_codelengths; }

  double getIndexCodelength() const { return m_optimizer->getIndexCodelength(); }

  double getModuleCodelength() const { return m_hierarchicalCodelength - m_optimizer->getIndexCodelength(); }

  double getHierarchicalCodelength() const { return m_hierarchicalCodelength; }

  double getOneLevelCodelength() const { return m_oneLevelCodelength; }

  double getRelativeCodelengthSavings() const
  {
    auto oneLevelCodelength = getOneLevelCodelength();
    return oneLevelCodelength < 1e-16 ? 0 : 1.0 - codelength() / oneLevelCodelength;
  }

  double getEntropyRate() { return m_entropyRate; }
  double getMaxEntropy() { return m_maxEntropy; }
  double getMaxFlow() { return m_maxFlow; }

  const Date& getStartDate() const { return m_startDate; }
  const Stopwatch& getElapsedTime() const { return m_elapsedTime; }

  std::vector<InfoNode*>& activeNetwork() const { return *m_activeNetwork; }

  std::map<unsigned int, std::vector<unsigned int>> getMultilevelModules(bool states = false);

  // ===================================================
  // IO
  // ===================================================

  std::ostream& toString(std::ostream& out) const { return m_optimizer->toString(out); }

  // ===================================================
  // Run
  // ===================================================

  using InitialPartition = std::map<unsigned int, unsigned int>;

  const InitialPartition& getInitialPartition() const { return m_initialPartition; }

  InfomapBase& setInitialPartition(const InitialPartition& moduleIds)
  {
    m_initialPartition = moduleIds;
    return *this;
  }

  void run(const std::string& parameters = "");

  void run(Network& network);

private:
  bool isFullNetwork() const { return m_isMain && m_aggregationLevel == 0; }
  bool isFirstLoop() const { return m_tuneIterationIndex == 0 && isFullNetwork(); }

  InfomapBase* getNewInfomapInstance() const { return new InfomapBase(getConfig()); }
  InfomapBase* getNewInfomapInstanceWithoutMemory() const
  {
    auto im = new InfomapBase();
    im->initOptimizer(true);
    return im;
  }

  InfomapBase& getSubInfomap(InfoNode& node) const
  {
    return node.setInfomap(getNewInfomapInstance())
        .setIsMain(false)
        .setSubLevel(m_subLevel + 1)
        .setNonMainConfig(*this);
  }

  InfomapBase& getSuperInfomap(InfoNode& node) const
  {
    return node.setInfomap(getNewInfomapInstanceWithoutMemory())
        .setIsMain(false)
        .setSubLevel(m_subLevel + SUPER_LEVEL_ADDITION)
        .setNonMainConfig(*this);
  }

  /**
   * Only the main infomap reads an external cluster file if exist
   */
  InfomapBase& setIsMain(bool isMain)
  {
    m_isMain = isMain;
    return *this;
  }

  InfomapBase& setSubLevel(unsigned int level)
  {
    m_subLevel = level;
    return *this;
  }

  bool isTopLevel() const { return (m_subLevel & (SUPER_LEVEL_ADDITION - 1)) == 0; }

  bool isSuperLevelOnTopLevel() const { return m_subLevel == SUPER_LEVEL_ADDITION; }

  bool isMainInfomap() const { return m_isMain; }

  bool haveHardPartition() const { return !m_originalLeafNodes.empty(); }

  // ===================================================
  // Run: *
  // ===================================================

  InfomapBase& initNetwork(Network& network);
  InfomapBase& initNetwork(InfoNode& parent, bool asSuperNetwork = false);

  void generateSubNetwork(Network& network);
  void generateSubNetwork(InfoNode& parent);

  /**
   * Init categorical meta data on all nodes from a file with the following format:
   * # nodeId metaData
   * 1 1
   * 2 1
   * 3 2
   * 4 2
   * 5 3
   *
   */
  InfomapBase& initMetaData(const std::string& metaDataFile);

  /**
   * Provide an initial partition of the network.
   *
   * @param clusterDataFile A .clu file containing cluster data.
   * @param hard If true, the provided clusters will not be splitted. This reduces the
   * effective network size during the optimization phase but the hard partitions are
   * after that replaced by the original nodes.
   */
  InfomapBase& initPartition(const std::string& clusterDataFile, bool hard = false, const Network* network = nullptr);

  /**
   * Provide an initial partition of the network.
   *
   * @param clusterIds map from nodeId to clusterId, doesn't have to be complete
   * @param hard If true, the provided clusters will not be splitted. This reduces the
   * effective network size during the optimization phase but the hard partitions are
   * after that replaced by the original nodes.
   */
  InfomapBase& initPartition(const std::map<unsigned int, unsigned int>& clusterIds, bool hard = false);

  /**
   * Provide an initial partition of the network.
   *
   * @param clusters Each sub-vector contain node IDs for all nodes that should be merged.
   * @param hard If true, the provided clusters will not be splitted. This reduces the
   * effective network size during the optimization phase but the hard partitions are
   * after that replaced by the original nodes.
   */
  InfomapBase& initPartition(std::vector<std::vector<unsigned int>>& clusters, bool hard = false);

  /**
   * Provide an initial partition of the network.
   *
   * @param modules Module indices for each node
   */
  InfomapBase& initPartition(std::vector<unsigned int>& modules, bool hard = false);

  /**
   * Provide an initial hierarchical partition of the network
   *
   * @param tree A tree path for each node
   */
  InfomapBase& initTree(const NodePaths& tree);

  void init();

  void runPartition()
  {
    if (twoLevel)
      partition();
    else
      hierarchicalPartition();
  }

  void restoreHardPartition();

  void writeResult(int trial = -1);

  // ===================================================
  // runPartition: *
  // ===================================================

  void hierarchicalPartition();

  void partition();

  // ===================================================
  // runPartition: init: *
  // ===================================================

  /**
   * Done in network?
   */
  void initEnterExitFlow();

  void aggregateFlowValuesFromLeafToRoot();

  // Init terms that is constant for the whole network
  void initTree() { return m_optimizer->initTree(); }

  void initNetwork() { return m_optimizer->initNetwork(); }

  void initSuperNetwork() { return m_optimizer->initSuperNetwork(); }

  double calcCodelength(const InfoNode& parent) const { return m_optimizer->calcCodelength(parent); }

  /**
   * Calculate and store codelength on all modules in the tree
   * @param includeRoot Also calculate the codelength on the root node
   * @return the hierarchical codelength
   */
  double calcCodelengthOnTree(InfoNode& root, bool includeRoot = true) const;

  // ===================================================
  // Run: Partition: *
  // ===================================================

  void setActiveNetworkFromLeafs() { m_activeNetwork = &m_leafNodes; }

  void setActiveNetworkFromChildrenOfRoot();

  void initPartition() { return m_optimizer->initPartition(); }

  void findTopModulesRepeatedly(unsigned int maxLevels);

  unsigned int fineTune();

  unsigned int coarseTune();

  /**
   * Return the number of effective core loops, i.e. not the last if not at coreLoopLimit
   */
  unsigned int optimizeActiveNetwork() { return m_optimizer->optimizeActiveNetwork(); }

  void moveActiveNodesToPredefinedModules(std::vector<unsigned int>& modules)
  {
    return m_optimizer->moveActiveNodesToPredefinedModules(modules);
  }

  void consolidateModules(bool replaceExistingModules = true)
  {
    return m_optimizer->consolidateModules(replaceExistingModules);
  }

  unsigned int calculateNumNonTrivialTopModules() const;

  unsigned int calculateMaxDepth() const;

  // ===================================================
  // Partition: findTopModulesRepeatedly: *
  // ===================================================

  /**
   * Return true if restored to consolidated optimization state
   */
  bool restoreConsolidatedOptimizationPointIfNoImprovement(bool forceRestore = false)
  {
    return m_optimizer->restoreConsolidatedOptimizationPointIfNoImprovement(forceRestore);
  }

  // ===================================================
  // Run: Hierarchical Partition: *
  // ===================================================

  /**
   * Find super modules applying the whole two-level algorithm on the
   * top modules iteratively
   * @param levelLimit The maximum number of super module levels allowed
   * @return number of levels created
   */
  unsigned int findHierarchicalSuperModules(unsigned int superLevelLimit = std::numeric_limits<unsigned int>::max());

  /**
   * Find super modules fast by merge and consolidate top modules iteratively
   * @param levelLimit The maximum number of super module levels allowed
   * @return number of levels created
   */
  unsigned int findHierarchicalSuperModulesFast(unsigned int superLevelLimit = std::numeric_limits<unsigned int>::max());

  void transformNodeFlowToEnterFlow(InfoNode& parent);

  void resetFlowOnModules();

  unsigned int removeModules();

  unsigned int removeSubModules(bool recalculateCodelengthOnTree);

  unsigned int recursivePartition();

  void queueTopModules(PartitionQueue& partitionQueue);

  void queueLeafModules(PartitionQueue& partitionQueue);

  bool processPartitionQueue(PartitionQueue& queue, PartitionQueue& nextLevel) const;

public:
  // ===================================================
  // Output: *
  // ===================================================

  /**
   * Write tree to a .tree file.
   * @param filename the filename for the output file. If empty, use default
   * based on output directory and input file name
   * @param states if memory network, print the state-level network without merging physical nodes within modules
   * @return the filename written to
   */
  std::string writeTree(const std::string& filename = "", bool states = false) { return infomap::writeTree(*this, m_network, filename, states); }

  /**
   * Write flow tree to a .ftree file.
   * This is the same as a .tree file but appended with links aggregated
   * within modules on all levels in the tree
   * @param filename the filename for the output file. If empty, use default
   * based on output directory and input file name
   * @param states if memory network, print the state-level network without merging physical nodes within modules
   * @return the filename written to
   */
  std::string writeFlowTree(const std::string& filename = "", bool states = false) { return infomap::writeFlowTree(*this, m_network, filename, states); }

  /**
   * Write Newick tree to a .tre file.
   * @param filename the filename for the output file. If empty, use default
   * based on output directory and input file name
   * @param states if memory network, print the state-level network without merging physical nodes within modules
   * @return the filename written to
   */
  std::string writeNewickTree(const std::string& filename = "", bool states = false) { return infomap::writeNewickTree(*this, filename, states); }

  std::string writeJsonTree(const std::string& filename = "", bool states = false, bool writeLinks = false) { return infomap::writeJsonTree(*this, m_network, filename, states, writeLinks); }

  std::string writeCsvTree(const std::string& filename = "", bool states = false) { return infomap::writeCsvTree(*this, m_network, filename, states); }

  /**
   * Write tree to a .clu file.
   * @param filename the filename for the output file. If empty, use default
   * based on output directory and input file name
   * @param states if memory network, print the state-level network without merging physical nodes within modules
   * @param moduleIndexLevel the depth from the root on which to advance module index.
   * Value 1 (default) will give the module index on the coarsest level, 2 the level below and so on.
   * Value -1 will give the module index for the lowest level, i.e. the finest modular structure.
   * @return the filename written to
   */
  std::string writeClu(const std::string& filename = "", bool states = false, int moduleIndexLevel = 1) { return infomap::writeClu(*this, m_network, filename, states, moduleIndexLevel); }

private:
  // ===================================================
  // Debug: *
  // ===================================================

  void printDebug() const { return m_optimizer->printDebug(); }

  // ===================================================
  // Members
  // ===================================================

protected:
  InfoNode m_root;
  std::vector<InfoNode*> m_leafNodes;
  std::vector<InfoNode*> m_moduleNodes;
  std::vector<InfoNode*>* m_activeNetwork = nullptr;

  std::vector<InfoNode*> m_originalLeafNodes;

  Network m_network;
  InitialPartition m_initialPartition = {}; // nodeId -> moduleId

  const unsigned int SUPER_LEVEL_ADDITION = 1 << 20;
  bool m_isMain = true;
  unsigned int m_subLevel = 0;

  bool m_calculateEnterExitFlow = false;

  double m_oneLevelCodelength = 0.0;
  unsigned int m_numNonTrivialTopModules = 0;
  unsigned int m_tuneIterationIndex = 0;
  bool m_isCoarseTune = false;
  unsigned int m_aggregationLevel = 0;

  double m_hierarchicalCodelength = 0.0;
  std::vector<double> m_codelengths;
  double m_entropyRate = 0.0;
  double m_maxEntropy = 0.0;
  double m_maxFlow = 0.0;

  double m_sumDanglingFlow = 0.0;

  Date m_startDate;
  Date m_endDate;
  Stopwatch m_elapsedTime = Stopwatch(false);
  std::string m_initialParameters;
  std::string m_currentParameters;

  std::unique_ptr<InfomapOptimizerBase> m_optimizer;
};

/**
 * Print per level statistics
 */
unsigned int printPerLevelCodelength(const InfoNode& parent, std::ostream& out);

void aggregatePerLevelCodelength(const InfoNode& parent, std::vector<detail::PerLevelStat>& perLevelStat, unsigned int level = 0);

namespace detail {

  struct PerLevelStat {
    double codelength() const { return indexLength + leafLength; }

    unsigned int numNodes() const { return numModules + numLeafNodes; }

    unsigned int numModules = 0;
    unsigned int numLeafNodes = 0;
    double indexLength = 0.0;
    double leafLength = 0.0;
  };

  class PartitionQueue {
    using PendingModule = InfoNode*;

    std::deque<PendingModule> m_queue;

  public:
    unsigned int level = 1;
    unsigned int numNonTrivialModules = 0;
    double flow = 0.0;
    double nonTrivialFlow = 0.0;
    bool skip = false;
    double indexCodelength = 0.0; // Consolidated
    double leafCodelength = 0.0; // Consolidated
    double moduleCodelength = 0.0; // Left to improve on next level

    using size_t = std::deque<PendingModule>::size_type;

    void swap(PartitionQueue& other) noexcept
    {
      std::swap(level, other.level);
      std::swap(numNonTrivialModules, other.numNonTrivialModules);
      std::swap(flow, other.flow);
      std::swap(nonTrivialFlow, other.nonTrivialFlow);
      std::swap(skip, other.skip);
      std::swap(indexCodelength, other.indexCodelength);
      std::swap(leafCodelength, other.leafCodelength);
      std::swap(moduleCodelength, other.moduleCodelength);
      m_queue.swap(other.m_queue);
    }

    size_t size() const { return m_queue.size(); }

    void resize(size_t size) { m_queue.resize(size); }

    PendingModule& operator[](size_t i) { return m_queue[i]; }
  };

} // namespace detail

} // namespace infomap

#endif // INFOMAP_BASE_H_
