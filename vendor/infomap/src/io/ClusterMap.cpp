/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#include "ClusterMap.h"
#include "SafeFile.h"
#include "InfomapError.h"
#include "../utils/Log.h"
#include "../utils/Console.h"
#include "../utils/FileURI.h"
#include "../utils/format.h"
#include <algorithm>
#include <limits>
#include <unordered_set>
#include <sstream>
#include "../utils/convert.h"

namespace infomap {

namespace {

  // A tree path is positive integers separated by ':' (1-based, kept as written).
  // Reading it with `>> unsigned` stopped at the first non-digit, so `x` and `:1` gave
  // an empty path that initTree then indexed, and `1:x` or `1.2` were cut or misread
  // without a word (#1208).
  bool parseTreePath(const std::string& pathString, Path& path)
  {
    std::size_t start = 0;
    while (true) {
      const auto end = std::min(pathString.find(':', start), pathString.size());
      if (end == start)
        return false;
      unsigned long long childNumber = 0;
      for (auto i = start; i < end; ++i) {
        const char c = pathString[i];
        if (c < '0' || c > '9')
          return false;
        childNumber = childNumber * 10 + static_cast<unsigned int>(c - '0');
        if (childNumber > std::numeric_limits<unsigned int>::max())
          return false;
      }
      if (childNumber == 0)
        return false;
      path.push_back(static_cast<unsigned int>(childNumber));
      if (end == pathString.size())
        return true;
      start = end + 1;
    }
  }

} // namespace

void ClusterMap::readClusterData(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId)
{
  m_clusterIds.clear();
  m_flowData.clear();
  m_treePaths.clear();
  m_extension.clear();
  m_filename = filename;
  m_duplicateClusterIds = {};
  m_isHigherOrder = false;
  m_hasTreeLeafIdType = false;
  m_treeLeafIdType = TreeLeafIdType::physical;

  // The readers throw a plain runtime_error on a missing or malformed file, and
  // FileURI an invalid_argument on a malformed name. Both are input errors (exit
  // 2), as for the network and the meta data, not internal ones (#1263). Other
  // logic errors pass through: they are bugs.
  try {
    FileURI file(filename);
    m_extension = file.getExtension();
    if (m_extension == "tree" || m_extension == "ftree") {
      return readTree(filename, includeFlow, layerNodeToStateId);
    }
    if (m_extension == "clu") {
      return readClu(filename, includeFlow, layerNodeToStateId);
    }
    throw std::runtime_error(fmt::format(FMT_STRING("Input cluster data from file '{}' is of unknown extension '{}'. Must be 'clu', 'tree' or 'ftree'."), filename, m_extension));
  } catch (const InfomapError&) {
    throw;
  } catch (const std::runtime_error& e) {
    throw InfomapError(ExitCode::InputError, e.what());
  } catch (const std::invalid_argument& e) {
    throw InfomapError(ExitCode::InputError, e.what());
  }
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

  unsigned int lineNr = 0;
  std::unordered_set<unsigned int> seenStateIds;

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
    Path path;
    if (!parseTreePath(pathString, path))
      throw std::runtime_error(fmt::format(FMT_STRING("Invalid tree path '{}' on line {} ('{}'). A path is positive integers separated by ':', such as 1:2:1."), pathString, lineNr, line));
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

    // A states tree has one row per state node. A repeated state id put the same
    // leaf into two modules, and a multi-level tree then crashed in initTree (#1335).
    // Physical rows repeat by design on higher-order networks, one per module of a
    // split physical node, and are paired with states in normalizeTreePaths.
    if (m_treeLeafIdType == TreeLeafIdType::state && !seenStateIds.insert(stateId).second) {
      const auto node = isMultilayer ? fmt::format(FMT_STRING("Node {} in layer {}"), nodeId, layerId) : fmt::format(FMT_STRING("State id {}"), stateId);
      throw std::runtime_error(fmt::format(FMT_STRING("{} appears more than once in the tree, again on line {} ('{}'). "
                                                      "A tree with state ids has one row per state node."),
                                           node,
                                           lineNr,
                                           line));
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
  // Repeat bookkeeping, kept only for ids that actually repeat: a well-formed clu
  // has one row per node (#1039), so keying this off the insert result below keeps
  // the common case at no extra storage and no extra lookup. Populated only from
  // rows that reach m_clusterIds, so multilayer rows skipped for a layer/node the
  // network does not have are never mistaken for duplicates. `conflicting` counts
  // the repeats that changed the module, which is a different thing from a row
  // repeated verbatim.
  struct RepeatCount {
    unsigned int rows = 0;
    unsigned int conflicting = 0;
  };
  std::map<unsigned int, RepeatCount> repeats;

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

    // insert, not find-then-assign: on a repeat its iterator already points at
    // what the id held, so the old value, the fact that it is a repeat, and the
    // overwrite all come out of a single traversal. A well-formed file pays exactly
    // one lookup per row, which is what a clu reader should cost.
    const auto insertion = m_clusterIds.insert({ stateId, moduleId });
    const auto entry = insertion.first;
    if (!insertion.second) {
      const bool changesModule = entry->second != moduleId;
      entry->second = moduleId; // the last row wins, which is what this reports on
      auto& repeat = repeats[stateId];
      if (repeat.rows == 0)
        repeat.rows = 1; // the row already in m_clusterIds
      ++repeat.rows;
      if (changesModule)
        ++repeat.conflicting;
    }
  }

  // The example id is there to illustrate the problem, so a conflicting id always
  // beats one that merely repeats a row verbatim, and within a tier the id with the
  // most rows wins. Without the first rule a file holding both kinds could name the
  // harmless one while the warning talks about changed modules.
  bool exampleIsConflicting = false;
  for (const auto& idRepeat : repeats) {
    const auto nodeId = idRepeat.first;
    const auto& repeat = idRepeat.second;
    m_duplicateClusterIds.rows += repeat.rows - 1;
    ++m_duplicateClusterIds.ids;

    const bool isConflicting = repeat.conflicting > 0;
    if (isConflicting) {
      m_duplicateClusterIds.conflictingRows += repeat.conflicting;
      ++m_duplicateClusterIds.conflictingIds;
    }

    const bool betterTier = isConflicting && !exampleIsConflicting;
    const bool sameTierMoreRows = isConflicting == exampleIsConflicting
        && repeat.rows > m_duplicateClusterIds.maxRowsForOneId;
    if (betterTier || sameTierMoreRows) {
      m_duplicateClusterIds.maxRowsForOneId = repeat.rows;
      m_duplicateClusterIds.exampleId = nodeId;
      exampleIsConflicting = isConflicting;
    }
  }
}

} // namespace infomap
