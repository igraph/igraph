/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "OutputView.h"
#include "../core/InfomapBase.h"
#include "../core/InfoNode.h"
#include "../core/iterators/InfomapIterator.h"
#include "../core/StateNetwork.h"
#include "../utils/convert.h"

namespace infomap {

OutputView::OutputView(InfomapBase& infomap, const StateNetwork& network, bool states)
    : m_infomap(infomap), m_network(network), m_states(states) {}

bool OutputView::isHigherOrderPhysicalLevel() const
{
  return m_infomap.haveMemory() && !m_states;
}

bool OutputView::isMultilayer() const
{
  return m_infomap.isMultilayerNetwork();
}

bool OutputView::hasMetaData() const
{
  return m_infomap.haveMetaData();
}

const char* OutputView::nodeIdHeaderName() const
{
  return m_states ? "state_id" : "node_id";
}

unsigned int OutputView::leafId(const OutputLeafRow& row) const
{
  return m_states ? row.stateId : row.physicalId;
}

unsigned int OutputView::leafId(const OutputTreeRow& row) const
{
  return m_states ? row.stateId : row.physicalId;
}

void OutputView::forEachLeaf(int moduleIndexLevel, OutputLeafPolicy filter, const LeafCallback& callback)
{
  if (isHigherOrderPhysicalLevel()) {
    for (auto it(m_infomap.iterTreePhysical(moduleIndexLevel)); !it.isEnd(); ++it) {
      InfoNode& node = *it;
      if (node.isLeaf() && shouldIncludeLeaf(node, filter)) {
        callback(leafRow(node, it.path(), it.moduleId(), it.modularCentrality()));
      }
    }
    return;
  }

  for (auto it(m_infomap.iterTree(moduleIndexLevel)); !it.isEnd(); ++it) {
    InfoNode& node = *it;
    if (node.isLeaf() && shouldIncludeLeaf(node, filter)) {
      callback(leafRow(node, it.path(), it.moduleId(), it.modularCentrality()));
    }
  }
}

void OutputView::forEachTreeNode(const TreeCallback& callback)
{
  if (isHigherOrderPhysicalLevel()) {
    for (auto it(m_infomap.iterTreePhysical()); !it.isEnd(); ++it) {
      callback(treeRow(*it, it.depth()));
    }
    return;
  }

  for (auto it(m_infomap.iterTree()); !it.isEnd(); ++it) {
    callback(treeRow(*it, it.depth()));
  }
}

void OutputView::forEachModule(const ModuleCallback& callback)
{
  auto linksByModulePath = moduleLinks();
  const OutputModuleLinkMap emptyLinks;
  for (auto it(m_infomap.iterModules()); !it.isEnd(); ++it) {
    const auto path = io::stringify(it.path(), ":");
    const auto& module = *it;
    const auto linksIt = linksByModulePath.find(path);
    const auto& links = linksIt == linksByModulePath.end() ? emptyLinks : linksIt->second;
    callback({
        io::stringify(it.path(), ","),
        path.empty() ? "root" : path,
        module.data.enterFlow,
        module.data.exitFlow,
        module.infomapChildDegree(),
        module.codelength,
        links,
    });
  }
}

std::map<unsigned int, OutputStateNodeTarget> OutputView::stateNodeTargets()
{
  std::map<unsigned int, OutputStateNodeTarget> targets;

  if (isHigherOrderPhysicalLevel()) {
    for (auto it(m_infomap.iterTreePhysical()); !it.isEnd(); ++it) {
      if (it->isLeaf()) {
        for (auto stateId : it->stateNodes) {
          targets[stateId] = { it->parent, it.childIndex() };
        }
      }
    }
    return targets;
  }

  for (auto it(m_infomap.iterTree()); !it.isEnd(); ++it) {
    if (it->isLeaf()) {
      targets[it->stateId] = { it->parent, it.childIndex() };
    }
  }

  return targets;
}

OutputModuleLinks OutputView::moduleLinks()
{
  auto stateTargets = stateNodeTargets();
  OutputModuleLinks moduleLinks;

  for (auto& leaf : m_infomap.leafNodes()) {
    for (auto& link : leaf->outEdges()) {
      const double flow = link->data.flow;
      auto& sourceTarget = stateTargets[link->source->stateId];
      auto& targetTarget = stateTargets[link->target->stateId];
      InfoNode* sourceParent = sourceTarget.parent;
      InfoNode* targetParent = targetTarget.parent;

      auto sourceDepth = sourceParent->calculatePath().size() + 1;
      auto targetDepth = targetParent->calculatePath().size() + 1;

      auto sourceChildIndex = sourceTarget.childIndex;
      auto targetChildIndex = targetTarget.childIndex;

      auto sourceParentIt = InfomapParentIterator(sourceParent);
      auto targetParentIt = InfomapParentIterator(targetParent);

      while (targetDepth > sourceDepth) {
        ++targetParentIt;
        --targetDepth;
      }

      while (sourceDepth > targetDepth) {
        ++sourceParentIt;
        --sourceDepth;
      }

      auto currentDepth = sourceDepth;
      while (currentDepth > 0) {
        if (sourceParentIt == targetParentIt && sourceChildIndex != targetChildIndex) {
          const auto parentId = io::stringify(sourceParentIt->calculatePath(), ":");
          auto& linkMap = moduleLinks[parentId];
          linkMap[std::make_pair(sourceChildIndex + 1, targetChildIndex + 1)] += flow;
        }

        sourceChildIndex = sourceParentIt->childIndex();
        targetChildIndex = targetParentIt->childIndex();

        ++sourceParentIt;
        ++targetParentIt;

        --currentDepth;
      }
    }
  }

  return moduleLinks;
}

bool OutputView::shouldIncludeLeaf(const InfoNode& node, OutputLeafPolicy filter) const
{
  const auto shouldHideBipartiteNodes = filter == OutputLeafPolicy::HideBipartite
      && m_infomap.isBipartite() && m_infomap.hideBipartiteNodes;

  return !shouldHideBipartiteNodes || node.physicalId < m_network.bipartiteStartId();
}

OutputLeafRow OutputView::leafRow(const InfoNode& node,
                                  const std::vector<unsigned int>& path,
                                  unsigned int moduleId,
                                  double modularCentrality) const
{
  return { node, path, moduleId, modularCentrality, node.stateId, node.physicalId, node.layerId, node.data.flow, nodeName(node) };
}

OutputTreeRow OutputView::treeRow(const InfoNode& node, unsigned int depth) const
{
  return { node, depth, node.stateId, node.physicalId, node.data.flow };
}

std::string OutputView::nodeName(const InfoNode& node) const
{
  const auto nameIt = m_network.names().find(node.physicalId);
  if (nameIt != m_network.names().end()) {
    return nameIt->second;
  }
  return io::stringify(node.physicalId);
}

} // namespace infomap
