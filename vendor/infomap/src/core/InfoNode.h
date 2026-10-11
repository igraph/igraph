/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef INFONODE_H_
#define INFONODE_H_

#include "FlowData.h"
#include "InfoEdge.h"
#ifndef SWIG
#include "ObjectPool.h"
#endif
#include "iterators/infomapIterators.h"
#include "iterators/IterWrapper.h"
#include "../utils/MetaCollection.h"

#include <stdexcept>
#include <memory>
#include <ostream>
#include <vector>
#include <limits>

namespace infomap {

class InfomapBase;

#ifndef SWIG
// Internal allocation plumbing; kept out of the SWIG bindings so the pool type
// does not leak into the generated wrappers.
class InfoNode;

template <typename T>
class ObjectPool;
using NodePool = ObjectPool<InfoNode>;

// Feature-only InfoNode fields, allocated out-of-line only when meta-data or
// regularized-multilayer flow is active. Ordinary and higher-order (state)
// networks never allocate this, so a node pays only the owning pointer.
struct InfoNodeExtras {
  std::vector<int> metaData; // Categorical value for each meta data dimension
  std::vector<LayerTeleFlowData> layerTeleFlowData; // For regularized multilayer flow
  std::unique_ptr<MetaCollection> metaCollection; // For modules; meta data only

  InfoNodeExtras() = default;
  InfoNodeExtras(const InfoNodeExtras& other)
      : metaData(other.metaData),
        layerTeleFlowData(other.layerTeleFlowData),
        metaCollection(other.metaCollection ? std::make_unique<MetaCollection>(*other.metaCollection) : nullptr) {}
};
#endif

/**
 * Tree node with raw-pointer ownership.
 *
 * An InfoNode owns the active child chain reachable through firstChild/lastChild
 * and the collapsed child chain reachable through collapsedFirstChild/
 * collapsedLastChild. Destruction deletes both chains. releaseChildren() only
 * detaches the active chain from this parent; the caller must reattach or delete
 * those nodes. Reparenting helpers detach children before deleting the removed
 * intermediate node. An InfoNode also owns its sub-Infomap pointer and outgoing
 * InfoEdge objects; incoming edge pointers are non-owning back-references.
 */
class InfoNode {
public:
  using child_iterator = ChildIterator<InfoNode*>;
  using const_child_iterator = ChildIterator<InfoNode const*>;
  using infomap_child_iterator = InfomapChildIterator<InfoNode*>;
  using const_infomap_child_iterator = InfomapChildIterator<InfoNode const*>;

  using tree_iterator = TreeIterator<InfoNode*>;
  using const_tree_iterator = TreeIterator<InfoNode const*>;

  using leaf_node_iterator = LeafNodeIterator<InfoNode*>;
  using const_leaf_node_iterator = LeafNodeIterator<InfoNode const*>;
  using leaf_module_iterator = LeafModuleIterator<InfoNode*>;
  using const_leaf_module_iterator = LeafModuleIterator<InfoNode const*>;

  using post_depth_first_iterator = DepthFirstIterator<InfoNode*, false>;
  using const_post_depth_first_iterator = DepthFirstIterator<InfoNode const*, false>;

  using edge_iterator = std::vector<InfoEdge*>::iterator;
  using const_edge_iterator = std::vector<InfoEdge*>::const_iterator;

  using edge_iterator_wrapper = IterWrapper<edge_iterator>;
  using const_edge_iterator_wrapper = IterWrapper<const_edge_iterator>;

  using infomap_iterator_wrapper = IterWrapper<tree_iterator>;
  using const_infomap_iterator_wrapper = IterWrapper<const_tree_iterator>;

  using child_iterator_wrapper = IterWrapper<child_iterator>;
  using const_child_iterator_wrapper = IterWrapper<const_child_iterator>;

  using infomap_child_iterator_wrapper = IterWrapper<infomap_child_iterator>;
  using const_infomap_child_iterator_wrapper = IterWrapper<const_infomap_child_iterator>;

public:
  FlowData data;
  unsigned int index = 0; // Temporary index used in finding best module
  unsigned int stateId = 0; // Unique state node id for the leaf nodes
  unsigned int physicalId = 0; // Physical id equals stateId for first order networks, otherwise can be non-unique
  unsigned int layerId = 0; // Layer id for multilayer networks

  InfoNode* owner = nullptr; // Infomap owner (if this is an Infomap root)
  InfoNode* parent = nullptr;
  InfoNode* previous = nullptr; // sibling
  InfoNode* next = nullptr; // sibling
  InfoNode* firstChild = nullptr;
  InfoNode* lastChild = nullptr;
  InfoNode* collapsedFirstChild = nullptr;
  InfoNode* collapsedLastChild = nullptr;
  double codelength = 0.0; // TODO: Better design for hierarchical stuff!?
  bool dirty = false;

  std::vector<PhysData> physicalNodes;
  std::vector<unsigned int> stateNodes; // For physically aggregated nodes

#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
  // Lossy map equation leaf aggregates, additive over contained leaf nodes:
  // sum of p_alpha * h_alpha (share of the Markov entropy rate) and sum of plogp(p_alpha).
  double lossyEntropy = 0.0;
  double lossyFlowLogFlow = 0.0;
#endif

private:
  unsigned int m_childDegree = 0;
  bool m_childrenChanged = false;
  unsigned int m_numLeafMembers = 0;

  std::vector<InfoEdge*> m_outEdges;
  std::vector<InfoEdge*> m_inEdges;

  InfomapBase* m_infomap = nullptr;

#ifndef SWIG
  // Owning pool. Set by InfomapBase::allocNode for every pool-allocated node;
  // stays nullptr only for the value-member root, which is never pool-freed.
  NodePool* m_pool = nullptr;
  // Pool that out-edges of this node are allocated from (same instance's edge
  // pool). nullptr for standalone nodes, whose edges fall back to new/delete.
  EdgePool* m_edgePool = nullptr;
  // Out-of-line feature-only fields (meta data / regularized multilayer);
  // nullptr for the common case (ordinary + higher-order networks).
  std::unique_ptr<InfoNodeExtras> m_extras;
  friend class InfomapBase;
#endif

public:
  InfoNode(const FlowData& flowData)
      : data(flowData) {};

  // For first order nodes, physicalId equals stateId
  InfoNode(const FlowData& flowData, unsigned int stateId)
      : data(flowData), stateId(stateId), physicalId(stateId) {};

  InfoNode(const FlowData& flowData, unsigned int stateId, unsigned int physicalId)
      : data(flowData), stateId(stateId), physicalId(physicalId) {};

  InfoNode(const FlowData& flowData, unsigned int stateId, unsigned int physicalId, unsigned int layerId)
      : data(flowData), stateId(stateId), physicalId(physicalId), layerId(layerId) {};

  InfoNode() = default;

  InfoNode(const InfoNode& other)
      : data(other.data),
        index(other.index),
        stateId(other.stateId),
        physicalId(other.physicalId),
        layerId(other.layerId),
        parent(other.parent),
        previous(other.previous),
        next(other.next),
        firstChild(other.firstChild),
        lastChild(other.lastChild),
        collapsedFirstChild(other.collapsedFirstChild),
        collapsedLastChild(other.collapsedLastChild),
        codelength(other.codelength),
        dirty(other.dirty),
        physicalNodes(other.physicalNodes),
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
        lossyEntropy(other.lossyEntropy),
        lossyFlowLogFlow(other.lossyFlowLogFlow),
#endif
        m_childDegree(other.m_childDegree),
        m_childrenChanged(other.m_childrenChanged),
        m_numLeafMembers(other.m_numLeafMembers),
        m_extras(other.m_extras ? std::make_unique<InfoNodeExtras>(*other.m_extras) : nullptr)
  {
  }

  ~InfoNode() noexcept;

  // Copy-assignment is required: InfomapIteratorPhysical holds a
  // std::map<unsigned int, InfoNode> and its defaulted operator= copy-assigns
  // the mapped InfoNode values. Mirror the copy constructor member-for-member
  // so the two copy paths stay consistent (the previous hand-written version
  // had silently dropped physicalNodes). stateNodes is intentionally omitted in
  // both, matching the constructor's long-standing behavior.
  InfoNode& operator=(const InfoNode& other)
  {
    if (this == &other)
      return *this;
    data = other.data;
    index = other.index;
    stateId = other.stateId;
    physicalId = other.physicalId;
    layerId = other.layerId;
    parent = other.parent;
    previous = other.previous;
    next = other.next;
    firstChild = other.firstChild;
    lastChild = other.lastChild;
    collapsedFirstChild = other.collapsedFirstChild;
    collapsedLastChild = other.collapsedLastChild;
    codelength = other.codelength;
    dirty = other.dirty;
    physicalNodes = other.physicalNodes;
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
    lossyEntropy = other.lossyEntropy;
    lossyFlowLogFlow = other.lossyFlowLogFlow;
#endif
    m_childDegree = other.m_childDegree;
    m_childrenChanged = other.m_childrenChanged;
    m_numLeafMembers = other.m_numLeafMembers;
    m_extras = other.m_extras ? std::make_unique<InfoNodeExtras>(*other.m_extras) : nullptr;
    return *this;
  }

  // ---------------------------- Getters ----------------------------

  // Body lives in InfoNode.cpp so the SWIG-parsed header carries only the
  // declaration (metaData now lives in the SWIG-hidden extras struct).
  unsigned int getMetaData(unsigned int dimension = 0) noexcept;

#ifndef SWIG
  // ---------------------------- Out-of-line extras ----------------------------
  static const std::vector<int>& emptyMetaData() noexcept
  {
    static const std::vector<int> empty;
    return empty;
  }
  static const std::vector<LayerTeleFlowData>& emptyLayerTeleFlowData() noexcept
  {
    static const std::vector<LayerTeleFlowData> empty;
    return empty;
  }

  const std::vector<int>& metaData() const noexcept { return m_extras ? m_extras->metaData : emptyMetaData(); }
  const std::vector<LayerTeleFlowData>& layerTeleFlowData() const noexcept { return m_extras ? m_extras->layerTeleFlowData : emptyLayerTeleFlowData(); }

  bool hasExtras() const noexcept { return m_extras != nullptr; }
  InfoNodeExtras& ensureExtras()
  {
    if (!m_extras)
      m_extras = std::make_unique<InfoNodeExtras>();
    return *m_extras;
  }

  // ---------------------------- Meta collection ----------------------------
  // Allocated only for nodes that carry meta data, via the extras struct.
  bool hasMetaCollection() const noexcept { return m_extras && m_extras->metaCollection != nullptr; }
  MetaCollection& metaCollection() noexcept { return *m_extras->metaCollection; } // precondition: hasMetaCollection()
  const MetaCollection& metaCollection() const noexcept { return *m_extras->metaCollection; }
  const MetaCollection* metaCollectionPtr() const noexcept { return m_extras ? m_extras->metaCollection.get() : nullptr; }

  MetaCollection& ensureMetaCollection()
  {
    auto& extras = ensureExtras();
    if (!extras.metaCollection)
      extras.metaCollection = std::make_unique<MetaCollection>();
    return *extras.metaCollection;
  }
#endif

  // ---------------------------- Infomap ----------------------------
  InfomapBase& getInfomap();

  const InfomapBase& getInfomap() const;

  InfomapBase& setInfomap(InfomapBase*);

  InfoNode* getInfomapRoot() noexcept;

  InfoNode const* getInfomapRoot() const noexcept;

  /**
   * Dispose the Infomap instance if it exists
   * @return true if an existing Infomap instance was deleted
   */
  bool disposeInfomap() noexcept;

  /**
   * Number of physical nodes in memory nodes
   */
  unsigned int numPhysicalNodes() const noexcept { return physicalNodes.size(); }

  // ---------------------------- Tree iterators ----------------------------

  // Default iteration on children
  child_iterator begin() noexcept { return { this }; }

  child_iterator end() noexcept { return { nullptr }; }

  const_child_iterator begin() const noexcept { return { this }; }

  const_child_iterator end() const noexcept { return { nullptr }; }

  child_iterator begin_child() noexcept { return { this }; }

  child_iterator end_child() noexcept { return { nullptr }; }

  const_child_iterator begin_child() const noexcept { return { this }; }

  const_child_iterator end_child() const noexcept { return { nullptr }; }

  child_iterator_wrapper children() noexcept { return { { this }, { nullptr } }; }

  const_child_iterator_wrapper children() const noexcept { return { { this }, { nullptr } }; }

  infomap_child_iterator_wrapper infomap_children() noexcept { return { { this }, { nullptr } }; }

  const_infomap_child_iterator_wrapper infomap_children() const noexcept { return { { this }, { nullptr } }; }

  post_depth_first_iterator begin_post_depth_first() noexcept { return { this }; }

  leaf_node_iterator begin_leaf_nodes() noexcept { return { this }; }

  leaf_module_iterator begin_leaf_modules() noexcept { return { this }; }

  tree_iterator begin_tree(unsigned int maxClusterLevel = std::numeric_limits<unsigned int>::max()) noexcept { return { this, static_cast<int>(maxClusterLevel) }; }

  tree_iterator end_tree() noexcept { return { nullptr }; }

  const_tree_iterator begin_tree(unsigned int maxClusterLevel = std::numeric_limits<unsigned int>::max()) const noexcept { return { this, static_cast<int>(maxClusterLevel) }; }

  const_tree_iterator end_tree() const noexcept { return { nullptr }; }

  infomap_iterator_wrapper infomapTree(unsigned int maxClusterLevel = std::numeric_limits<unsigned int>::max()) noexcept { return { { this, static_cast<int>(maxClusterLevel) }, { nullptr } }; }

  const_infomap_iterator_wrapper infomapTree(unsigned int maxClusterLevel = std::numeric_limits<unsigned int>::max()) const noexcept { return { { this, static_cast<int>(maxClusterLevel) }, { nullptr } }; }

  // ---------------------------- Graph iterators ----------------------------

  edge_iterator begin_outEdge() noexcept { return m_outEdges.begin(); }

  edge_iterator end_outEdge() noexcept { return m_outEdges.end(); }

  edge_iterator begin_inEdge() noexcept { return m_inEdges.begin(); }

  edge_iterator end_inEdge() noexcept { return m_inEdges.end(); }

  edge_iterator_wrapper outEdges() noexcept { return { m_outEdges }; }

  edge_iterator_wrapper inEdges() noexcept { return { m_inEdges }; }

  // ---------------------------- Capacity ----------------------------

  unsigned int childDegree() const noexcept { return m_childDegree; }

  bool isLeaf() const noexcept { return firstChild == nullptr; }

  // TODO: Safe to assume all children are leaves if first child is leaf?
  bool isLeafModule() const noexcept { return m_infomap == nullptr && firstChild != nullptr && firstChild->firstChild == nullptr; }

  bool isRoot() const noexcept { return parent == nullptr; }

  unsigned int depth() const noexcept;

  unsigned int firstDepthBelow() const noexcept;

  unsigned int numLeafMembers() const noexcept { return m_numLeafMembers; }

  bool isDangling() const noexcept { return m_outEdges.empty(); }

  unsigned int outDegree() const noexcept { return m_outEdges.size(); }

  unsigned int inDegree() const noexcept { return m_inEdges.size(); }

  unsigned int degree() const noexcept { return outDegree() + inDegree(); }

  // ---------------------------- Order ----------------------------
  bool isFirst() const noexcept { return !parent || parent->firstChild == this; }

  bool isLast() const noexcept { return !parent || parent->lastChild == this; }

  unsigned int childIndex() const noexcept;

  // Generate 1-based tree path
  std::vector<unsigned int> calculatePath() const noexcept;

  unsigned int infomapChildDegree() const noexcept;

  unsigned int id() const noexcept { return stateId; }

  // ---------------------------- Operators ----------------------------

  bool operator==(const InfoNode& rhs) const noexcept { return this == &rhs; }

  bool operator!=(const InfoNode& rhs) const noexcept { return this != &rhs; }

  friend std::ostream& operator<<(std::ostream& out, const InfoNode& node) noexcept
  {
    if (node.isLeaf())
      out << "[" << node.physicalId << "]";
    else
      out << "[module]";
    return out;
  }

  // ---------------------------- Mutators ----------------------------

  /**
   * Clear a cloned node to initial state
   */
  void initClean() noexcept;

  void sortChildrenOnFlow(bool recurse = true) noexcept;

  /**
   * Move the active child chain to the collapsed child chain without deleting
   * it. The node owns the collapsed chain until expandChildren() restores it or
   * deleteChildren()/destruction deletes it.
   * @return the number of children collapsed
   */
  unsigned int collapseChildren() noexcept;

  /**
   * Move the collapsed child chain back to the active child chain.
   * @return the number of collapsed children expanded
   */
  unsigned int expandChildren();

  // ------ OLD -----

  // After change, set the child degree if known instead of lazily computing it by traversing the linked list
  void setChildDegree(unsigned int value) noexcept;

  void setNumLeafNodes(unsigned int value) noexcept { m_numLeafMembers = value; }

  /**
   * Append child to this node's active child chain and transfer tree ownership
   * to this parent. The child must not still be owned by another active chain.
   */
  void addChild(InfoNode* child) noexcept;

  /**
   * Detach this node from its active child chain without deleting any child.
   * Child parent/sibling links are not rewritten; callers must reattach or
   * delete the detached chain before those stale links can be observed.
   */
  void releaseChildren() noexcept;

  /**
   * If not already having a single child, replace children
   * with a single new node, assuming grandchildren.
   * @return the single child
   */
  InfoNode& replaceChildrenWithOneNode();

  /**
   * @return 1 if the node is removed, otherwise 0
   */
  unsigned int replaceWithChildren() noexcept;

  void replaceWithChildrenDebug() noexcept;

  /**
   * @return The number of children removed
   */
  unsigned int replaceChildrenWithGrandChildren() noexcept;

  void replaceChildrenWithGrandChildrenDebug() noexcept;

  /**
   * Delete this node.
   *
   * removeChildren=false leaves firstChild intact, so the destructor deletes the
   * active child chain. removeChildren=true clears firstChild before deletion,
   * leaving the previous child chain for the caller to reattach or delete.
   */
  void remove(bool removeChildren) noexcept;

  /**
   * Delete active and collapsed child chains owned by this node.
   */
  void deleteChildren() noexcept;

  void addOutEdge(InfoNode& target, double weight, double flow = 0.0) noexcept
  {
    auto* edge = m_edgePool != nullptr
        ? m_edgePool->alloc(*this, target, weight, flow)
        : new InfoEdge(*this, target, weight, flow);
    m_outEdges.push_back(edge);
    target.m_inEdges.push_back(edge);
  }

private:
  void calcChildDegree() noexcept;

  // Destroy a node through its owning pool, or via plain delete if it was not
  // pool-allocated (e.g. a standalone node built in a test). Lets InfoNode stay
  // usable without an InfomapBase pool while still recycling pooled nodes.
  static void destroyNode(InfoNode* node) noexcept;
};

// Guards the InfoNode shrink (feature-only fields moved out-of-line into
// InfoNodeExtras). The common case (ordinary + higher-order) pays only the
// extras pointer. The lossy feature adds two doubles inline.
// Achieved on this toolchain: 328 -> 296 B in the default build (the three
// feature-only fields' 56 nominal bytes net ~32 after alignment); the lossy
// feature adds two doubles (~312 B). The bounds carry a small margin so struct
// layout/alignment differences across platforms (Windows/macOS CI) don't trip
// the guard, while still catching a regression back toward the un-shrunk sizes
// (328 default / 344 lossy).
//
// Skip the guard under debug-STL configs (_GLIBCXX_DEBUG, libc++ debug,
// MSVC _ITERATOR_DEBUG_LEVEL>0): those inflate std::vector/iterator layout, so
// InfoNode legitimately exceeds these release-layout bounds without any real
// regression. The guard's job is to catch re-bloat in normal builds, not to
// constrain debug-iterator layouts.
#if !defined(_GLIBCXX_DEBUG) && !defined(_GLIBCXX_DEBUG_PEDANTIC)                            \
    && !defined(_LIBCPP_DEBUG) && !(defined(_LIBCPP_DEBUG_LEVEL) && _LIBCPP_DEBUG_LEVEL > 0) \
    && !(defined(_ITERATOR_DEBUG_LEVEL) && _ITERATOR_DEBUG_LEVEL > 0)
#if INFOMAP_FEATURE_LOSSY_MAP_EQUATION
static_assert(sizeof(InfoNode) <= 320, "InfoNode grew beyond the shrink target (lossy build)");
#else
static_assert(sizeof(InfoNode) <= 304, "InfoNode grew beyond the shrink target");
#endif
#endif

} // namespace infomap

#endif // INFONODE_H_
