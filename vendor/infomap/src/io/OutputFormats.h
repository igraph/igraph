/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef OUTPUT_FORMATS_H_
#define OUTPUT_FORMATS_H_

#include <cstdint>
#include <string>
#include <vector>

namespace infomap {

enum class OutputKind : std::uint8_t {
  Clu,
  Tree,
  FlowTree,
  Newick,
  Json,
  Csv,
  PajekNetwork,
  StateNetwork,
  FlowNetwork
};

struct OutputFileFormat {
  std::string resultKey;
  std::string displayName;
  bool states;
  std::string suffix;
  std::string extension;
  std::string mimeType;
  unsigned int resultOrder;
};

struct OutputFormat {
  std::string optionName;
  OutputKind kind;
  std::vector<OutputFileFormat> files;
};

const std::vector<OutputFormat>& outputFormats();
std::vector<std::string> outputFormatNames();
const OutputFormat* findOutputFormat(const std::string& optionName);
const OutputFileFormat* findOutputFileFormat(const std::string& resultKey);
std::string outputFilename(const std::string& basename, const OutputFileFormat& file);
std::string outputFilenameForResultKey(const std::string& basename, const std::string& resultKey);
std::string outputFormatsManifestJson();

} // namespace infomap

#endif // OUTPUT_FORMATS_H_
