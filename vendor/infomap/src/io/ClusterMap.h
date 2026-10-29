/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef CLUSTER_MAP_H_
#define CLUSTER_MAP_H_

#include <cstdint>
#include <string>
#include <map>
#include <vector>

namespace infomap {

using Path = std::vector<unsigned int>; // 1-based indexing

using NodePath = std::pair<unsigned int, Path>;

using NodePaths = std::vector<NodePath>;

enum class TreeLeafIdType : std::uint8_t {
  physical,
  state,
};

struct TreePath {
  unsigned int nodeId = 0; // Either a physical id or a state id, see idType
  Path path;
  TreeLeafIdType idType = TreeLeafIdType::physical;
};

using TreePaths = std::vector<TreePath>;

// What a clu file's repeated node ids amounted to. A clu file has one row per
// node, so a repeat means the file is not a partition of the network it is being
// applied to: the reader keeps whichever row it saw last and the earlier module
// assignments are gone. The physical clu of a higher-order network is the case
// that matters -- it has one row per (physical node, module) pair, so overlapping
// modules make its ids repeat by construction (#1039).
struct DuplicateClusterIds {
  unsigned int rows = 0; // repeated rows, whether or not they changed anything
  unsigned int ids = 0; // distinct node ids seen more than once
  unsigned int exampleId = 0; // the worst one, named in the warning
  unsigned int maxRowsForOneId = 0; // how many rows that id had

  // Repeats that actually changed the id's module. A file can repeat a row
  // verbatim, and then nothing is replaced and the partition it describes is the
  // one the reader ends up with -- so only these justify saying the result is not
  // the file's partition.
  unsigned int conflictingRows = 0;
  unsigned int conflictingIds = 0;

  bool any() const noexcept { return rows > 0; }
  bool anyConflicting() const noexcept { return conflictingRows > 0; }
};

class ClusterMap {
public:
  void readClusterData(const std::string& filename, bool includeFlow = false, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);

  const std::map<unsigned int, unsigned int>& clusterIds() const noexcept
  {
    return m_clusterIds;
  }

  const TreePaths& treePaths() const noexcept { return m_treePaths; }

  const std::string& extension() const noexcept { return m_extension; }

  //! The file the last readClusterData() call parsed; empty before any read.
  const std::string& filename() const noexcept { return m_filename; }

  TreeLeafIdType treeLeafIdType() const noexcept { return m_treeLeafIdType; }

  // Only ever non-empty for a clu file. The tree reader has its own report for the
  // same situation, which it can make more precise because a tree row carries a
  // path it can pair with a state id -- see normalizeTreePaths.
  const DuplicateClusterIds& duplicateClusterIds() const noexcept { return m_duplicateClusterIds; }

private:
  void readTree(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);
  void readClu(const std::string& filename, bool includeFlow, const std::map<unsigned int, std::map<unsigned int, unsigned int>>* layerNodeToStateId = nullptr);

  std::map<unsigned int, unsigned int> m_clusterIds;
  std::map<unsigned int, double> m_flowData;
  TreePaths m_treePaths;
  std::string m_extension;
  std::string m_filename;
  DuplicateClusterIds m_duplicateClusterIds;
  bool m_isHigherOrder = false;
  bool m_hasTreeLeafIdType = false;
  TreeLeafIdType m_treeLeafIdType = TreeLeafIdType::physical;
};

} // namespace infomap

#endif // CLUSTER_MAP_H_
