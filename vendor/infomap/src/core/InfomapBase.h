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
#include "../io/RunMetadata.h"
#include "../utils/Log.h"
#include "../utils/Date.h"
#include "../utils/Stopwatch.h"

#include <algorithm>
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
   * Number of node levels below the root, counted over every branch, 1 if no modules.
   * Alias of maxTreeDepth. This is the depth the summary, the tree writer and the JSON
   * output report, so it is what any other report of a depth should use (#1036).
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
  // The seed the run was configured with, which is not the same as reading
  // seedToRandomNumberGenerator: the serial loop moves that field to the current
  // trial's seed (base + trialOffset + index) and only gives it back when the loop
  // ends, while updateBestResult writes the aggregate artifact from inside the loop.
  // A header reading the live field therefore recorded the last trial's seed
  // whenever the best trial was the last one, since restoreBestResult then performs
  // no rewrite. Falls back to the live value until a run has captured it.
  //
  // Guarded from SWIG so it is not exposed as a new binding, the same way
  // getReferenceOneLevelCodelength above is: the output writer is its only caller.
  unsigned long baseSeed() const { return m_haveBaseSeed ? m_baseSeed : seedToRandomNumberGenerator; }

  // Identity of every result-affecting input file and of the effective config,
  // captured once when a run starts so every artifact the run writes -- the
  // .tree/.ftree/.clu headers, the JSON tree, the manifest -- carries the same
  // values (#1026). Unknown (empty hash) when the network was built in memory
  // or the file was not used. Guarded from SWIG like baseSeed(); the bindings
  // read all of it through provenanceJson().
  const InputIdentity& inputIdentity() const { return m_inputIdentity; }
  const InputIdentity& clusterDataIdentity() const { return m_clusterDataIdentity; }
  const InputIdentity& metaDataIdentity() const { return m_metaDataIdentity; }
  const std::string& runConfigFingerprint() const { return m_configFingerprint; }
#endif

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

  /**
   * Depth of the leftmost branch: follows firstChild to the first leaf. O(depth) where
   * numLevels is O(nodes), which is why the search uses it -- consolidateModules alone
   * calls it ~270k times per web-NotreDame trial. It equals numLevels only on a
   * uniform-depth tree, which every one of its call sites has by construction: they run
   * inside the per-module optimization, below the point where sub-Infomap results are
   * stitched in and make the tree ragged. Measured over the test suite and the benchmark
   * networks, the two never disagreed at any of those sites. Never report this value --
   * that is the bug in #1036.
   */
  unsigned int depthOfFirstLeaf() const;

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

  // A sub-Infomap partitions one module, and its modules are never the tree's top modules, so it
  // runs without --preferred-number-of-modules. The tree charges that term once, on its top
  // modules (calcCodelengthOnTree); a sub-search charging it on its own module count kept or
  // rejected splits on a cost the tree does not have (#1292). The objective reads the
  // preference when the instance is built, so it is cleared here, not in cloneAsNonMain.
  InfomapBase* getNewInfomapInstance() const
  {
    Config subConfig = getConfig();
    subConfig.preferredNumberOfModules = 0;
    return new InfomapBase(subConfig);
  }
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
                           .inheritRuntimeContext(*this)
                           .inheritInterruptHandler(*this);
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
                             .inheritRuntimeContext(*this)
                             .inheritInterruptHandler(*this);
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

  // Share the cancel flag and the largest-team count with a sub/worker Infomap.
  // The callback is NOT inherited here — a parallel-trial worker only reads the flag.
  InfomapBase& inheritRuntimeContext(const InfomapBase& other)
  {
    m_cancel = other.m_cancel;
    m_largestTeam = other.m_largestTeam;
    return *this;
  }

  // A nested sub/super Infomap also takes the callback and the owner thread, so a
  // signal is still polled while the main instance waits inside its run(). The
  // owner-thread guard in pollInterrupt() keeps the callback on that thread (#1134).
  InfomapBase& inheritInterruptHandler(const InfomapBase& other)
  {
    m_interruptCallback = other.m_interruptCallback;
    m_interruptUserData = other.m_interruptUserData;
    m_ownerThreadId = other.m_ownerThreadId;
    return *this;
  }

  // Any thread, from inside a parallel region: record the size of its team. The
  // largest one is the timing report's threads_used, so that field counts the
  // threads the run actually started rather than predicting them (#1145).
  void recordTeamSize(unsigned int threads) noexcept
  {
    auto largest = m_largestTeam->load(std::memory_order_relaxed);
    while (threads > largest && !m_largestTeam->compare_exchange_weak(largest, threads, std::memory_order_relaxed)) {
    }
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

#ifndef SWIG
  /**
   * Apply cluster data that has already been parsed. A run parses --cluster-data
   * once and applies it from here in every trial, so the file is read exactly
   * once however many trials there are (#1072). Not binding API: the bindings
   * set clusterDataFile and let the run do this.
   */
  InfomapBase& initPartition(const ClusterMap& clusterData, bool hard = false);
#endif

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

  //! Reject cluster data that gives a module both a leaf and a sub-module as children.
  //! Call this where a tree is parsed from input; initTree does not, because it also
  //! materializes engine-generated trees, which cannot have that shape (#990).
  static void validateClusterDataTreeShape(const NodePaths& tree);

#ifndef SWIG
  //! Apply a tree parsed from input: a .tree/.ftree file or JSON node paths. A hard
  //! partition from a tree must be one level deep, and is applied as from a .clu file.
  InfomapBase& initTreeInput(const NodePaths& tree, bool hard);
#endif

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

  // The scale of enter and exit flow. They are built from link flow scaled by the Markov time, so
  // their rounding grows with it, and so does the tolerance of every check on them (#1192).
  double rateScale() const { return std::max(1.0, markovTime); }

  // Flow invariants of the tree under `root`: total flow, no negative flow, enter and exit
  // flow that balance at every module, and module flow that equals the sum of its children.
  // Throws InfomapError naming `phase` on the first violation. Compiled in only when
  // INFOMAP_CHECK_INVARIANTS is set (debug, sanitizer and fuzz builds); a no-op otherwise.
#if INFOMAP_CHECK_INVARIANTS
  void checkInvariants(const InfoNode& root, const char* phase) const;
  // Throws the InfomapOptimizer::checkMoveDelta violation, formatted here rather than in the header.
  [[noreturn]] void failMoveDelta(const char* where, double predictedDelta, double actualDelta) const;
  // The codelength the search tracked must be the one the objective prices on the partition it
  // ended on, and the codelength that recursivePartition compared for a module must be the one
  // of the subtree it kept (#1279). Throws InfomapError naming `phase` on a mismatch.
  void checkCodelength(double tracked, double priced, const char* phase) const;
  void checkPartitionCodelength() const;
#else
  void checkInvariants(const InfoNode& /*root*/, const char* /*phase*/) const {}
#endif

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

  // The objective's network-level terms describe whichever network was active when they
  // were last initialised -- root's children, which are the leaves only sometimes. Record
  // which, so a later leaf-level partition can restore them instead of scoring the leaves
  // against a module network's terms (#831; see initPartition(std::vector<unsigned int>&)).
  void initNetwork()
  {
    m_optimizer->initNetwork();
    m_objectiveTermsAreForLeafNetwork = m_root.firstChild != nullptr && m_root.firstChild->isLeaf();
  }

  // Re-init the network-level terms for the leaf network if they describe a module network.
  // Call it with the leaves as root's children, after removeModules().
  void restoreLeafNetworkTerms()
  {
    if (!m_objectiveTermsAreForLeafNetwork)
      initNetwork();
  }

  void initSuperNetwork()
  {
    m_optimizer->initSuperNetwork();
    m_objectiveTermsAreForLeafNetwork = false;
  }

  double calcCodelength(const InfoNode& parent) const { return m_optimizer->calcCodelength(parent); }

  double calcTreeCodelengthCost(unsigned int numTopModules) const { return m_optimizer->calcTreeCodelengthCost(numTopModules); }

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

  /**
   * Dissolve every intermediate module whose index codebook no longer pays for
   * itself: its children move into its parent, one level up. Best-first over the
   * whole tree; each gain is priced through the objective on the spliced tree.
   * @return the number of modules dissolved
   */
  unsigned int dissolveUnprofitableModules();

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
  unsigned long m_baseSeed = 0;
  bool m_haveBaseSeed = false;
  // See inputIdentity(). The wrapper's readInputData records the files it read
  // in m_readInputPaths, since the bindings never set networkFile.
  InputIdentity m_inputIdentity;
  InputIdentity m_clusterDataIdentity;
  InputIdentity m_metaDataIdentity;
  std::string m_configFingerprint;
  std::vector<std::string> m_readInputPaths;
  // Hashed when the file was read, not when the run starts: a file edited in
  // between must not be hashed in the name of a partition built from the old
  // content. Dropped by any in-memory build call (InfomapWrapper).
  InputIdentity m_inputIdentityFromRead;
  bool m_networkHasInMemoryContent = false;
  unsigned int m_numNonTrivialTopModules = 0;
  unsigned int m_tuneIterationIndex = 0;
  bool m_isCoarseTune = false;
  unsigned int m_aggregationLevel = 0;
  //! Whether the objective's network-level terms currently describe the leaf network (#831).
  bool m_objectiveTermsAreForLeafNetwork = false;

  double m_hierarchicalCodelength = 0.0;
  // Whether the leaf flow of this run sums to 1, for checkInvariants(). Set by the run once the
  // flow is calculated; false until then, which skips that one check.
  bool m_flowSumsToOne = false;
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
  // The largest OpenMP team the run started, shared through the run tree the same
  // way. Reset in run(Network&); see recordTeamSize.
  std::atomic<unsigned int> m_largestTeamSize { 1 };
  std::atomic<unsigned int>* m_largestTeam = &m_largestTeamSize;
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
