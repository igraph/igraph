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
#include "iterators/infomapIterators.h"
#include "iterators/IterWrapper.h"
#include "../utils/MetaCollection.h"

#include <stdexcept>
#include <memory>
#include <iostream>
#include <vector>
#include <limits>

namespace infomap {

class InfomapBase;

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
  std::vector<int> metaData; // Categorical value for each meta data dimension

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
  MetaCollection metaCollection; // For modules
  std::vector<unsigned int> stateNodes; // For physically aggregated nodes

private:
  unsigned int m_childDegree = 0;
  bool m_childrenChanged = false;
  unsigned int m_numLeafMembers = 0;

  std::vector<InfoEdge*> m_outEdges;
  std::vector<InfoEdge*> m_inEdges;

  InfomapBase* m_infomap = nullptr;

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
        metaData(other.metaData),
        parent(other.parent),
        previous(other.previous),
        next(other.next),
        firstChild(other.firstChild),
        lastChild(other.lastChild),
        collapsedFirstChild(other.collapsedFirstChild),
        collapsedLastChild(other.collapsedLastChild),
        codelength(other.codelength),
        dirty(other.dirty),
        metaCollection(other.metaCollection),
        m_childDegree(other.m_childDegree),
        m_childrenChanged(other.m_childrenChanged),
        m_numLeafMembers(other.m_numLeafMembers) { }

  ~InfoNode() noexcept;

  InfoNode& operator=(const InfoNode& other)
  {
    data = other.data;
    index = other.index;
    stateId = other.stateId;
    physicalId = other.physicalId;
    layerId = other.layerId;
    metaData = other.metaData;
    parent = other.parent;
    previous = other.previous;
    next = other.next;
    firstChild = other.firstChild;
    lastChild = other.lastChild;
    collapsedFirstChild = other.collapsedFirstChild;
    collapsedLastChild = other.collapsedLastChild;
    codelength = other.codelength;
    dirty = other.dirty;
    metaCollection = other.metaCollection;
    m_childDegree = other.m_childDegree;
    m_childrenChanged = other.m_childrenChanged;
    m_numLeafMembers = other.m_numLeafMembers;
    return *this;
  }

  // ---------------------------- Getters ----------------------------

  unsigned int getMetaData(unsigned int dimension = 0) noexcept
  {
    if (dimension >= metaData.size()) {
      return 0;
    }
    auto meta = metaData[dimension];
    return meta < 0 ? 0 : static_cast<unsigned int>(meta);
  }

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
   * Release the children and store the child pointers for later expansion
   * @return the number of children collapsed
   */
  unsigned int collapseChildren() noexcept;

  /**
   * Expand collapsed children
   * @return the number of collapsed children expanded
   */
  unsigned int expandChildren();

  // ------ OLD -----

  // After change, set the child degree if known instead of lazily computing it by traversing the linked list
  void setChildDegree(unsigned int value) noexcept;

  void setNumLeafNodes(unsigned int value) noexcept { m_numLeafMembers = value; }

  void addChild(InfoNode* child) noexcept;

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

  void remove(bool removeChildren) noexcept;

  void deleteChildren() noexcept;

  void addOutEdge(InfoNode& target, double weight, double flow = 0.0) noexcept
  {
    auto* edge = new InfoEdge(*this, target, weight, flow);
    m_outEdges.push_back(edge);
    target.m_inEdges.push_back(edge);
  }

private:
  void calcChildDegree() noexcept;
};

} // namespace infomap

#endif // INFONODE_H_
