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
#include "ObjectPool.h"
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
#include <array>
#include <atomic>
#include <stdexcept>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <ostream>

namespace infomap {

namespace detail {
  class PartitionQueue;
  struct PartitionTaskRecord;
  struct PerLevelStat;
} // namespace detail

// Cooperative cancellation hook (issue #412). Return true to stop the run at the
// next checkpoint. Invoked only on the owner thread, so it is safe to enter
// host-language APIs (R_CheckUserInterrupt, PyErr_CheckSignals) from it.
using InterruptCallback = bool (*)(void* userData);

// Thrown out of run() on cancellation; std::runtime_error so existing
// catch(std::exception) unwinding still applies.
class InterruptionError : public std::runtime_error {
public:
  InterruptionError() : std::runtime_error("Infomap run interrupted.") {}
};

class InfomapBase : public InfomapConfig<InfomapBase> {
  template <typename Objective>
  friend class InfomapOptimizer;

  void initOptimizer(bool forceNoMemory = false);

public:
  using PartitionQueue = detail::PartitionQueue;

  InfomapBase() { initOptimizer(); }

  explicit InfomapBase(const Config& conf) : InfomapConfig<InfomapBase>(conf), m_network(conf) { initOptimizer(); }

  explicit InfomapBase(const std::string& flags, bool isCli = false) : InfomapConfig<InfomapBase>(flags, isCli)
  {
    initOptimizer();
    m_network.setConfig(*this);
    m_initialParameters = m_currentParameters = flags;
  }

  ~InfomapBase() override = default;

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

#ifndef SWIG
  // One-level reference reported to the user. In lossy mode this is the lossless
  // L_1 = H(p_alpha) (a lambda-independent upper bound), not the lambda-dependent
  // corrected one-module objective that m_oneLevelCodelength holds for the optimizer.
  // Guarded from SWIG so it is not exposed as a new binding; callers below use it.
  double getReferenceOneLevelCodelength() const
  {
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
    if (lossy)
      return m_optimizer->getLossyOneLevelLossless();
#endif
    return m_oneLevelCodelength;
  }
#endif

  double getRelativeCodelengthSavings() const
  {
    auto oneLevelCodelength = getReferenceOneLevelCodelength();
    return oneLevelCodelength < 1e-16 ? 0 : 1.0 - codelength() / oneLevelCodelength;
  }

#ifndef SWIG
  double getRelativeCodelengthSavings(double codelength) const
  {
    auto oneLevelCodelength = getReferenceOneLevelCodelength();
    return oneLevelCodelength < 1e-16 ? 0 : 1.0 - codelength / oneLevelCodelength;
  }
#endif

  double getEntropyRate() { return m_entropyRate; }
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  double getLossyRate() const { return m_optimizer->getLossyRate(); }
  double getLossyDistortion() const { return m_optimizer->getLossyDistortion(); }
  double getLossyOneLevelLossless() const { return m_optimizer->getLossyOneLevelLossless(); }
  // 1-based indices of the top modules that are noise modules (l_i > lambda * H_i).
  std::vector<unsigned int> noiseTopModules() const;
#endif
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
    m_multilayerInitialPartition.clear();
    return *this;
  }

  // Initial partition for a multilayer network keyed by physical identity, given
  // as parallel (layer_id, node_id, module_id) arrays. Resolution to state ids
  // is deferred until the network is built (the state ids don't exist before
  // then); see initTrialPartition. (issue #616)
  InfomapBase& setMultilayerInitialPartition(const std::vector<unsigned int>& layerIds,
                                             const std::vector<unsigned int>& nodeIds,
                                             const std::vector<unsigned int>& moduleIds)
  {
    if (layerIds.size() != nodeIds.size() || layerIds.size() != moduleIds.size())
      throw std::invalid_argument("setMultilayerInitialPartition: layer/node/module arrays must have equal length");
    m_multilayerInitialPartition.assign(layerIds.size(), {});
    for (std::size_t i = 0; i < layerIds.size(); ++i)
      m_multilayerInitialPartition[i] = { layerIds[i], nodeIds[i], moduleIds[i] };
    m_initialPartition.clear();
    return *this;
  }

  void run(const std::string& parameters = "");

  void run(Network& network);

  // ===================================================
  // Cooperative interruption (issue #412)
  // ===================================================

  // Polled on the owner thread during run(); return true to abort with an
  // InterruptionError. No handler = unchanged behavior.
  InfomapBase& setInterruptHandler(InterruptCallback callback, void* userData = nullptr)
  {
    m_interruptCallback = callback;
    m_interruptUserData = userData;
    return *this;
  }

  InfomapBase& clearInterruptHandler()
  {
    m_interruptCallback = nullptr;
    m_interruptUserData = nullptr;
    return *this;
  }

  bool hasInterruptHandler() const { return m_interruptCallback != nullptr; }

  // Thread-safe push to cancel an in-progress run, observed at the next
  // checkpoint. A request made before run() is cleared at run entry (so it does
  // not pre-cancel a fresh run). Relaxed (no companion data) — ordering comes
  // from the OpenMP joins.
  void requestInterrupt() noexcept { m_cancel->store(true, std::memory_order_relaxed); }

  bool interruptRequested() const noexcept { return m_cancel->load(std::memory_order_relaxed); }

private:
  class RunSession;

  // Allocate a tree node from this instance's pool and stamp its owning pool
  // back-pointer so node-local teardown (deleteChildren / remove /
  // replaceWithChildren) can reclaim it without an InfomapBase handle.
  template <typename... Args>
  InfoNode& allocNode(Args&&... args)
  {
    InfoNode& node = *m_nodePool.alloc(std::forward<Args>(args)...);
    node.m_pool = &m_nodePool;
    node.m_edgePool = &m_edgePool;
    return node;
  }

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
    auto& subInfomap = node.setInfomap(getNewInfomapInstance())
                           .setIsMain(false)
                           .setSubLevel(m_subLevel + 1)
                           .setNonMainConfig(*this)
                           .inheritRuntimeContext(*this);
    // Carry the full-network properties down so entropy bias correction stays consistent
    // across the hierarchy (this used to be a shared static).
    subInfomap.m_optimizer->inheritNetworkPropertiesFrom(*m_optimizer);
    return subInfomap;
  }

  InfomapBase& getSuperInfomap(InfoNode& node) const
  {
    auto& superInfomap = node.setInfomap(getNewInfomapInstanceWithoutMemory())
                             .setIsMain(false)
                             .setSubLevel(m_subLevel + SUPER_LEVEL_ADDITION)
                             .setNonMainConfig(*this)
                             .inheritRuntimeContext(*this);
    superInfomap.m_optimizer->inheritNetworkPropertiesFrom(*m_optimizer);
    // getNewInfomapInstanceWithoutMemory builds the objective from a default Config, so unlike a
    // sub-Infomap this one has to be told what to minimise. findHierarchicalSuperModules compares
    // its codelength directly against this instance's index codelength, and the two are only
    // comparable when both objectives are configured alike (#904).
    superInfomap.m_optimizer->inheritObjectiveParametersFrom(*m_optimizer);
    return superInfomap;
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

  // Share the cancel flag (only) with a sub/worker Infomap. The callback is NOT
  // inherited — only the main instance, on the owner thread, ever invokes it.
  InfomapBase& inheritRuntimeContext(const InfomapBase& other)
  {
    m_cancel = other.m_cancel;
    return *this;
  }

  // Flag-only, any thread. Call only where a throw can unwind legally — NOT
  // inside an OpenMP structured block.
  void checkCancelled() const
  {
    if (m_cancel->load(std::memory_order_relaxed))
      throw InterruptionError();
  }

  // The thread-id guard keeps the callback off worker threads even when this
  // object's methods run inside the task graph (off-thread it is a flag check).
  void pollInterrupt()
  {
    if (m_interruptCallback != nullptr
        && std::this_thread::get_id() == m_ownerThreadId
        && m_interruptCallback(m_interruptUserData))
      m_cancel->store(true, std::memory_order_relaxed);
    checkCancelled();
  }

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
   * @param clusterDataFile A .clu, .tree or .ftree file containing cluster data.
   * Tree input may be physical or state-level on higher-order networks.
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

  /**
   * Normalize tree leaf ids so the leaf id always identifies a state node in
   * m_leafNodes. For higher-order networks, physical-id rows are expanded by
   * mapping each physical id to all state ids that share it.
   *
   * @param tree Tree paths as parsed from a .tree/.ftree cluster file.
   * @param numNodesNotInNetwork Output count of physical ids that don't match any state node.
   */
  NodePaths normalizeTreePaths(const TreePaths& tree, unsigned int& numNodesNotInNetwork) const;

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

  void partitionModuleRecursively(InfoNode& module, unsigned int level, detail::PartitionTaskRecord& record) const;

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
  // Declared before m_root so the tree (whose nodes live in m_nodePool) is
  // reclaimed by ~m_root before the pool storage is dropped. Destruction runs
  // in reverse declaration order.
  NodePool m_nodePool;
  EdgePool m_edgePool;
  InfoNode m_root;
  std::vector<InfoNode*> m_leafNodes;
  std::vector<InfoNode*> m_moduleNodes;
  std::vector<InfoNode*>* m_activeNetwork = nullptr;

  std::vector<InfoNode*> m_originalLeafNodes;

  Network m_network;
  InitialPartition m_initialPartition; // nodeId -> moduleId
  // Pending physical multilayer partition: {layer_id, node_id, module_id} rows,
  // resolved to state ids in initTrialPartition once the network is built.
  std::vector<std::array<unsigned int, 3>> m_multilayerInitialPartition;

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
  std::vector<unsigned int> m_numTopModules;
  double m_entropyRate = 0.0;
  double m_maxEntropy = 0.0;
  double m_maxFlow = 0.0;

  double m_sumDanglingFlow = 0.0;

  Date m_startDate;
  Date m_endDate;
  Stopwatch m_elapsedTime = Stopwatch(false);
  std::string m_initialParameters;
  std::string m_currentParameters;

  // Cooperative interruption (issue #412). m_cancel points at m_cancelRequested
  // (self) or, for sub/worker instances, at the main's flag (inheritRuntimeContext),
  // so the run tree shares one flag. m_ownerThreadId is stamped in run(Network&).
  std::atomic<bool> m_cancelRequested { false };
  std::atomic<bool>* m_cancel = &m_cancelRequested;
  InterruptCallback m_interruptCallback = nullptr;
  void* m_interruptUserData = nullptr;
  std::thread::id m_ownerThreadId;

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
