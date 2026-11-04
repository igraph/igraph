/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef OUTPUT_VIEW_H_
#define OUTPUT_VIEW_H_

#include <cstdint>
#include <vector>
#include <functional>
#include <map>
#include <string>
#include <utility>

namespace infomap {

class InfomapBase;
class InfoNode;
class StateNetwork;

enum class OutputLeafPolicy : std::uint8_t {
  //! Honour --hide-bipartite-nodes.
  HideBipartite,
  //! Keep the feature nodes whatever the option says, because the flow tree's link
  //! section refers to them. Decided per output rather than from the global
  //! printFlowTree flag: one writer serves both .tree and .ftree, so asking the flag
  //! meant requesting --ftree also un-hid the nodes in the plain .tree (#908).
  KeepBipartite
};

struct OutputLeafRow {
  const InfoNode& node;
  std::vector<unsigned int> path;
  unsigned int moduleId = 0;
  double modularCentrality = 0.0;
  unsigned int stateId = 0;
  unsigned int physicalId = 0;
  unsigned int layerId = 0;
  double flow = 0.0;
  std::string name;
};

struct OutputTreeRow {
  const InfoNode& node;
  unsigned int depth = 0;
  unsigned int stateId = 0;
  unsigned int physicalId = 0;
  double flow = 0.0;
};

struct OutputStateNodeTarget {
  InfoNode* parent = nullptr;
  unsigned int childIndex = 0;
};

using OutputModulePath = std::string;
using OutputModuleLink = std::pair<unsigned int, unsigned int>;
using OutputModuleLinkMap = std::map<OutputModuleLink, double>;
using OutputModuleLinks = std::map<OutputModulePath, OutputModuleLinkMap>;

struct OutputModuleRow {
  OutputModulePath jsonPath;
  OutputModulePath linkPathLabel;
  double enterFlow = 0.0;
  double exitFlow = 0.0;
  unsigned int numChildren = 0;
  double codelength = 0.0;
  const OutputModuleLinkMap& links;
};

class OutputView {
public:
  using LeafCallback = std::function<void(const OutputLeafRow&)>;
  using TreeCallback = std::function<void(const OutputTreeRow&)>;
  using ModuleCallback = std::function<void(const OutputModuleRow&)>;

  OutputView(InfomapBase& infomap, const StateNetwork& network, bool states);

  bool isStateLevel() const { return m_states; }
  bool isPhysicalLevel() const { return !m_states; }
  bool isHigherOrderPhysicalLevel() const;
  bool isMultilayer() const;
  bool hasMetaData() const;
  const char* nodeIdHeaderName() const;
  unsigned int leafId(const OutputLeafRow& row) const;
  unsigned int leafId(const OutputTreeRow& row) const;

  void forEachLeaf(int moduleIndexLevel, OutputLeafPolicy filter, const LeafCallback& callback);
  void forEachTreeNode(const TreeCallback& callback);
  // OutputModuleRow::links is valid only for the duration of the callback.
  void forEachModule(const ModuleCallback& callback);

  std::map<unsigned int, OutputStateNodeTarget> stateNodeTargets();
  OutputModuleLinks moduleLinks();

private:
  bool shouldIncludeLeaf(const InfoNode& node, OutputLeafPolicy filter) const;
  OutputLeafRow leafRow(const InfoNode& node,
                        const std::vector<unsigned int>& path,
                        unsigned int moduleId,
                        double modularCentrality) const;
  OutputTreeRow treeRow(const InfoNode& node, unsigned int depth) const;
  std::string nodeName(const InfoNode& node) const;

  InfomapBase& m_infomap;
  const StateNetwork& m_network;
  bool m_states = false;
};

} // namespace infomap

#endif // OUTPUT_VIEW_H_
