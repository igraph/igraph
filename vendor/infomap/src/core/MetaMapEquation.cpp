/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "MetaMapEquation.h"
#include "FlowData.h"
#include "InfoNode.h"
#include "../utils/Log.h"

#include <memory>
#include <vector>
#include <utility>

namespace infomap {

double MetaMapEquation::getModuleCodelength() const
{
  return moduleCodelength + metaCodelength * metaDataRate;
}

double MetaMapEquation::getCodelength() const
{
  return codelength + metaCodelength * metaDataRate;
}

// ===================================================
// IO
// ===================================================

std::ostream& MetaMapEquation::print(std::ostream& out) const
{
  return out << indexCodelength << " + " << moduleCodelength << " + " << metaCodelength << " = " << io::toPrecision(getCodelength());
}

std::ostream& operator<<(std::ostream& out, const MetaMapEquation& mapEq)
{
  return mapEq.print(out);
}

// ===================================================
// Init
// ===================================================

void MetaMapEquation::init(const Config& config)
{
  Log(3) << "MetaMapEquation::init()...\n";
  Base::init(config);
  numMetaDataDimensions = config.numMetaDataDimensions;
  metaDataRate = config.metaDataRate;
  weightByFlow = !config.unweightedMetaData;
}

void MetaMapEquation::initTree(InfoNode& root)
{
  Log(3) << "MetaMapEquation::initTree()...\n";
  initMetaNodes(root);
}

void MetaMapEquation::initNetwork(InfoNode& root)
{
  Log(3) << "MetaMapEquation::initNetwork()...\n";
  Base::initNetwork(root);
  m_unweightedNodeFlow = 1.0 / root.childDegree();
}

void MetaMapEquation::initPartition(std::vector<InfoNode*>& nodes)
{
  initPartitionOfMetaNodes(nodes);
  calculateCodelength(nodes);
}

void MetaMapEquation::initMetaNodes(InfoNode& root)
{
  bool notInitiated = !root.firstChild->hasMetaCollection() || root.firstChild->metaCollection().empty();
  if (notInitiated) {
    Log(3) << "MetaMapEquation::initMetaNodes()...\n";

    for (auto it = root.begin_post_depth_first(); !it.isEnd(); ++it) {
      auto& node = *it;
      if (node.isRoot()) {
        break;
      }
      if (node.isLeaf()) {
        if (!node.metaData().empty()) {
          // TODO: Use flow here and move weightByFlow choice to metaCollection, using flowCount?
          double flow = weightByFlow ? node.data.flow : m_unweightedNodeFlow;
          node.parent->ensureMetaCollection().add(node.metaData()[0], flow);
        } else {
          throw std::length_error("A node is missing meta data using MetaMapEquation");
        }
      } else if (node.hasMetaCollection()) {
        node.parent->ensureMetaCollection().add(node.metaCollection());
      }
    }
  }
}

void MetaMapEquation::initPartitionOfMetaNodes(std::vector<InfoNode*>& nodes)
{
  Log(4) << "MetaMapEquation::initPartitionOfMetaNodes()...\n";
  m_moduleToMetaCollection.clear();

  for (auto& n : nodes) {
    InfoNode& node = *n;
    unsigned int moduleIndex = node.index; // Assume unique module index for all nodes in this initiation phase
    if (!node.hasMetaCollection() || node.metaCollection().empty()) {
      if (!node.metaData().empty()) {
        double flow = weightByFlow ? node.data.flow : m_unweightedNodeFlow;
        node.ensureMetaCollection().add(node.metaData()[0], flow);
      } else
        throw std::length_error("A node is missing meta data using MetaMapEquation");
    }
    m_moduleToMetaCollection[moduleIndex] = node.metaCollection();
  }
}

// ===================================================
// Codelength
// ===================================================

void MetaMapEquation::calculateCodelength(std::vector<InfoNode*>& nodes)
{
  calculateCodelengthTerms(nodes);

  calculateCodelengthFromCodelengthTerms();

  metaCodelength = 0.0;

  // Treat each node as a single module
  for (InfoNode* n : nodes) {
    InfoNode& node = *n;
    if (node.hasMetaCollection())
      metaCodelength += node.metaCollection().calculateEntropy();
  }
}

double MetaMapEquation::calcCodelength(const InfoNode& parent) const
{
  return parent.isLeafModule() ? calcCodelengthOnModuleOfLeafNodes(parent) : Base::calcCodelengthOnModuleOfModules(parent);
}

double MetaMapEquation::calcCodelengthOnModuleOfLeafNodes(const InfoNode& parent) const
{
  double indexLength = Base::calcCodelength(parent);

  // Meta addition
  MetaCollection metaCollection;
  for (const InfoNode& node : parent) {
    if (node.hasMetaCollection() && !node.metaCollection().empty())
      metaCollection.add(node.metaCollection());
    else
      metaCollection.add(node.metaData()[0], weightByFlow ? node.data.flow : m_unweightedNodeFlow); // TODO: Initiate to collection and use all dimensions
  }

  double _metaCodelength = metaCollection.calculateEntropy();

  return indexLength + metaDataRate * _metaCodelength;
}

INFOMAP_HOT double MetaMapEquation::getDeltaCodelengthOnMovingNode(InfoNode& current,
                                                                   DeltaFlow& oldModuleDelta,
                                                                   DeltaFlow& newModuleDelta,
                                                                   std::vector<FlowData>& moduleFlowData,
                                                                   std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  double deltaMetaL = 0.0;

  unsigned int oldModuleIndex = oldModuleDelta.module;
  unsigned int newModuleIndex = newModuleDelta.module;

  // Remove codelength of old and new module before changes
  deltaMetaL -= getCurrentModuleMetaCodelength(oldModuleIndex, current, 0);
  deltaMetaL -= getCurrentModuleMetaCodelength(newModuleIndex, current, 0);
  // Add codelength of old module with current node removed
  deltaMetaL += getCurrentModuleMetaCodelength(oldModuleIndex, current, -1);
  // Add codelength of old module with current node added
  deltaMetaL += getCurrentModuleMetaCodelength(newModuleIndex, current, 1);

  return deltaL + deltaMetaL * metaDataRate;
}

INFOMAP_HOT double MetaMapEquation::getDeltaCodelengthOnMovingNodeHoisted(InfoNode& current,
                                                                          DeltaFlow& oldModuleDelta,
                                                                          const OldSideTerms& oldSide,
                                                                          DeltaFlow& newModuleDelta,
                                                                          std::vector<FlowData>& moduleFlowData,
                                                                          std::vector<unsigned int>& moduleMembers)
{
  double deltaL = Base::getDeltaCodelengthOnMovingNodeHoisted(current, oldModuleDelta, oldSide, newModuleDelta, moduleFlowData, moduleMembers);

  double deltaMetaL = 0.0;

  unsigned int oldModuleIndex = oldModuleDelta.module;
  unsigned int newModuleIndex = newModuleDelta.module;

  // Remove codelength of old and new module before changes
  deltaMetaL -= getCurrentModuleMetaCodelength(oldModuleIndex, current, 0);
  deltaMetaL -= getCurrentModuleMetaCodelength(newModuleIndex, current, 0);
  // Add codelength of old module with current node removed
  deltaMetaL += getCurrentModuleMetaCodelength(oldModuleIndex, current, -1);
  // Add codelength of old module with current node added
  deltaMetaL += getCurrentModuleMetaCodelength(newModuleIndex, current, 1);

  return deltaL + deltaMetaL * metaDataRate;
}

double MetaMapEquation::getCurrentModuleMetaCodelength(unsigned int module, InfoNode& current, int addRemoveOrNothing)
{
  auto& currentMetaCollection = m_moduleToMetaCollection[module];
  const MetaCollection* nodeMeta = current.metaCollectionPtr();

  double moduleMetaCodelength = 0.0;

  if (addRemoveOrNothing == 0 || nodeMeta == nullptr) {
    moduleMetaCodelength = currentMetaCollection.calculateEntropy();
  }
  // If add or remove, do the change, calculate new codelength and then undo the change
  else if (addRemoveOrNothing == 1) {
    currentMetaCollection.add(*nodeMeta);
    moduleMetaCodelength = currentMetaCollection.calculateEntropy();
    currentMetaCollection.remove(*nodeMeta);
  } else {
    currentMetaCollection.remove(*nodeMeta);
    moduleMetaCodelength = currentMetaCollection.calculateEntropy();
    currentMetaCollection.add(*nodeMeta);
  }

  return moduleMetaCodelength;
}

// ===================================================
// Consolidation
// ===================================================

void MetaMapEquation::updateCodelengthOnMovingNode(InfoNode& current,
                                                   DeltaFlow& oldModuleDelta,
                                                   DeltaFlow& newModuleDelta,
                                                   std::vector<FlowData>& moduleFlowData,
                                                   std::vector<unsigned int>& moduleMembers)
{
  Base::updateCodelengthOnMovingNode(current, oldModuleDelta, newModuleDelta, moduleFlowData, moduleMembers);

  double deltaMetaL = 0.0;

  unsigned int oldModuleIndex = oldModuleDelta.module;
  unsigned int newModuleIndex = newModuleDelta.module;

  // Remove codelength of old and new module before changes
  deltaMetaL -= getCurrentModuleMetaCodelength(oldModuleIndex, current, 0);
  deltaMetaL -= getCurrentModuleMetaCodelength(newModuleIndex, current, 0);

  // Update meta data from moving node
  updateMetaData(current, oldModuleIndex, newModuleIndex);

  // Add codelength of old and new module after changes
  deltaMetaL += getCurrentModuleMetaCodelength(oldModuleIndex, current, 0);
  deltaMetaL += getCurrentModuleMetaCodelength(newModuleIndex, current, 0);

  metaCodelength += deltaMetaL;
}

void MetaMapEquation::updateMetaData(InfoNode& current, unsigned int oldModuleIndex, unsigned int bestModuleIndex)
{
  const MetaCollection* nodeMeta = current.metaCollectionPtr();
  if (nodeMeta == nullptr)
    return;

  // Remove meta id from old module (can be a set of meta ids when moving submodules in coarse tune)
  m_moduleToMetaCollection[oldModuleIndex].remove(*nodeMeta);

  // Add meta id to new module
  m_moduleToMetaCollection[bestModuleIndex].add(*nodeMeta);
}

void MetaMapEquation::consolidateModules(std::vector<InfoNode*>& modules)
{
  for (auto& module : modules) {
    if (module == nullptr)
      continue;
    module->ensureExtras().metaCollection = std::make_unique<MetaCollection>(m_moduleToMetaCollection[module->index]);
  }
}

// ===================================================
// Debug
// ===================================================

void MetaMapEquation::printDebug() const
{
  Log() << "MetaMapEquation\n";
  Base::printDebug();
}

} // namespace infomap
