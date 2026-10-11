/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "ClusterMap.h"
#include "SafeFile.h"
#include "../utils/Log.h"
#include "../utils/Console.h"
#include "../utils/FileURI.h"
#include "../utils/format.h"
#include <sstream>
#include "../utils/convert.h"

namespace infomap {

void ClusterMap::readClusterData(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId)
{
  m_clusterIds.clear();
  m_flowData.clear();
  m_treePaths.clear();
  m_extension.clear();
  m_isHigherOrder = false;
  m_hasTreeLeafIdType = false;
  m_treeLeafIdType = TreeLeafIdType::physical;

  FileURI file(filename);
  m_extension = file.getExtension();
  if (m_extension == "tree" || m_extension == "ftree") {
    return readTree(filename, includeFlow, layerNodeToStateId);
  }
  if (m_extension == "clu") {
    return readClu(filename, includeFlow, layerNodeToStateId);
  }
  throw std::runtime_error(fmt::format(FMT_STRING("Input cluster data from file '{}' is of unknown extension '{}'. Must be 'clu', 'tree' or 'ftree'."), filename, m_extension));
}

/**
 * Sample from .tree file (physical level)
# Codelength = 3.46227314 bits.
# path flow name node_id
1:1:1 0.0384615 "1" 1
1:1:2 0.025641 "2" 2
 *
 * Sample from .tree file (state level)
# path flow name state_id node_id
1:1 0.166667 "i" 1 1
2:1 0.166667 "i" 4 1
 */
void ClusterMap::readTree(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId)
{
  bool isMultilayer = layerNodeToStateId != nullptr;

  SafeInFile input(filename);
  std::string line;
  std::istringstream lineStream;
  std::istringstream pathStream;

  unsigned int lineNr = 0;

  while (!std::getline(input, line).fail()) {
    ++lineNr;
    if (line.empty())
      continue;
    if (line[0] == '#') {
      // Header lines like `# path flow name node_id [...]` are just
      // human-readable decoration — the row column count is authoritative.
      continue;
    }
    if (line[0] == '*') {
      break;
    }

    lineStream.clear();
    lineStream.str(line);

    std::string pathString;
    double flow;
    std::string name;
    unsigned int parsedId;
    unsigned int nodeId = 0;
    unsigned int layerId = 0;
    if (!(lineStream >> pathString))
      throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse tree path from line '{}'"), line));
    if (!(lineStream >> flow))
      throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse node flow from line '{}'"), line));
    // Everything between the first and the last quote, the same rule the network parser
    // uses for Pajek vertex names. Reading up to the first quote with getline truncated
    // a name that contains one and shifted the rest of the line, so Infomap could not
    // read back its own .tree: "gene "X", alias" failed as "Couldn't parse node id"
    // (#908). The written format is unchanged.
    const auto flowEnd = lineStream.tellg();
    const auto searchFrom = flowEnd < 0 ? std::size_t { 0 } : static_cast<std::size_t>(flowEnd);
    const auto nameStart = line.find_first_of('"', searchFrom);
    const auto nameEnd = line.find_last_of('"');
    if (nameStart == std::string::npos || nameEnd <= nameStart)
      throw std::runtime_error(fmt::format(FMT_STRING("Can't parse node name from line {} ('{}')."), lineNr, line));
    name = line.substr(nameStart + 1, nameEnd - nameStart - 1);
    lineStream.seekg(static_cast<std::streamoff>(nameEnd + 1));
    if (!(lineStream >> parsedId))
      throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse node id from line '{}'"), line));

    const auto hasExplicitNodeId = static_cast<bool>(lineStream >> nodeId);
    const auto inferredLeafIdType = hasExplicitNodeId ? TreeLeafIdType::state : TreeLeafIdType::physical;

    if (!m_hasTreeLeafIdType) {
      m_treeLeafIdType = inferredLeafIdType;
      m_hasTreeLeafIdType = true;
    } else if (inferredLeafIdType != m_treeLeafIdType) {
      // Earlier rows had a different column count. A file that mixes
      // physical-id (4 columns) and state-id (5+ columns) rows cannot be
      // parsed safely.
      throw std::runtime_error(fmt::format(FMT_STRING("Mixed state and physical tree ids are not supported in line '{}'."), line));
    }

    if (m_treeLeafIdType == TreeLeafIdType::state) {
      if (hasExplicitNodeId) {
        m_isHigherOrder = true;
      } else if (m_isHigherOrder) {
        throw std::runtime_error(fmt::format(FMT_STRING("Missing node id from line '{}'."), line));
      }
      if (isMultilayer) {
        if (!hasExplicitNodeId)
          throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse node key from line '{}'"), line));
        if (!(lineStream >> layerId))
          throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse layer id from line '{}'"), line));
      }
    }

    bool multilayerNodeFound = false;
    unsigned int stateId = parsedId;

    if (isMultilayer && m_treeLeafIdType == TreeLeafIdType::state) {
      // Re-map state id from layer/node ids in case the multilayer state ids differ
      auto it = layerNodeToStateId->find(layerId);
      if (it != layerNodeToStateId->end()) {
        auto nodeIdToStateId = it->second.find(nodeId);
        if (nodeIdToStateId != it->second.end()) {
          stateId = nodeIdToStateId->second;
          multilayerNodeFound = true;
        }
      }
      if (!multilayerNodeFound) {
        // Skip rows whose layer/node combination is not present in the network
        continue;
      }
    }

    pathStream.clear();
    pathStream.str(pathString);
    unsigned int childNumber;

    Path path;
    while (pathStream >> childNumber) {
      pathStream.get(); // Extract the delimiting character also
      if (childNumber == 0)
        throw std::runtime_error("There is a '0' in the tree path, lowest allowed integer is 1.");
      path.push_back(childNumber); // Keep 1-based indexing in path
    }

    m_treePaths.push_back({ stateId, path, m_treeLeafIdType });

    if (includeFlow)
      m_flowData[stateId] = flow;
  }
}

void ClusterMap::readClu(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId)
{
  auto isMultilayer = layerNodeToStateId != nullptr;

  Console::detail(1, "reading initial partition from '{}'", filename);
  SafeInFile input(filename);
  std::string line;
  std::istringstream lineStream;
  std::map<unsigned int, unsigned int> clusterData;

  while (!std::getline(input, line).fail()) {
    if (line.empty() || line[0] == '#' || line[0] == '*')
      continue;

    lineStream.clear();
    lineStream.str(line);
    // # state_id module flow node_id layer_id

    unsigned int stateId;
    unsigned int nodeId;
    unsigned int moduleId;
    unsigned int layerId;

    // Validated rather than read straight into the unsigneds: `>> unsigned`
    // accepts a leading '-' and stores the value modulo 2^32.
    std::string stateIdToken;
    std::string moduleIdToken;
    if (!(lineStream >> stateIdToken >> moduleIdToken))
      throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse node key and cluster id from line '{}'"), line));
    if (!io::parseNonNegativeInteger(stateIdToken, stateId) || !io::parseNonNegativeInteger(moduleIdToken, moduleId))
      throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse node key and cluster id from line '{}': expected two non-negative integers, got '{}' and '{}'"), line, stateIdToken, moduleIdToken));

    auto flow = 0.0;
    if (lineStream >> flow) {
      if (includeFlow)
        m_flowData[stateId] = flow;
    }

    auto multilayerNodeFound = false;
    if (isMultilayer) {
      if (!(lineStream >> nodeId))
        throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse node key from line '{}'"), line));

      if (!(lineStream >> layerId))
        throw std::runtime_error(fmt::format(FMT_STRING("Couldn't parse layer id from line '{}'"), line));

      // get new state id from map
      auto it = layerNodeToStateId->find(layerId);

      if (it != layerNodeToStateId->end()) {
        auto nodeIdToStateId = it->second.find(nodeId);
        if (nodeIdToStateId != it->second.end()) {
          stateId = nodeIdToStateId->second;
          multilayerNodeFound = true;
        }
      }
    }

    if (isMultilayer && !multilayerNodeFound) {
      continue;
    }

    m_clusterIds[stateId] = moduleId;
  }
}

} // namespace infomap
