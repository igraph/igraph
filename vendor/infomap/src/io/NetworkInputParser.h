/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef NETWORK_INPUT_PARSER_H_
#define NETWORK_INPUT_PARSER_H_

#include "SafeFile.h"
#include "../core/StateNetwork.h"
#include "../utils/Log.h"
#include "../utils/Console.h"
#include "../utils/convert.h"
#include "../utils/format.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace infomap {
namespace input {

  struct NetworkInputOptions {
    bool undirectedFlow = true;
    unsigned int matchableMultilayerIds = 0;
  };

  struct ParsedVertex {
    unsigned int id = 0;
    std::string name;
    bool hasWeight = false;
    double weight = 1.0;
  };

  struct ParsedLink {
    unsigned int source = 0;
    unsigned int target = 0;
    double weight = 1.0;
  };

  struct ParsedStateNode {
    StateNetwork::StateNode node;
    bool hasWeight = false;
  };

  struct ParsedMultilayerLink {
    unsigned int sourceLayer = 0;
    unsigned int sourceNode = 0;
    unsigned int targetLayer = 0;
    unsigned int targetNode = 0;
    double weight = 1.0;
  };

  struct ParsedMultilayerIntraLink {
    unsigned int layer = 0;
    unsigned int sourceNode = 0;
    unsigned int targetNode = 0;
    double weight = 1.0;
  };

  struct ParsedMultilayerInterLink {
    unsigned int sourceLayer = 0;
    unsigned int node = 0;
    unsigned int targetLayer = 0;
    double weight = 1.0;
  };

  namespace detail {

    using InsensitiveStringSet = std::set<std::string, io::InsensitiveCompare>;

    struct ParseProgress {
      unsigned int numPhysicalNodes = 0;
      unsigned int numLinks = 0;
      unsigned int numIntraLayerLinks = 0;
      unsigned int numInterLayerLinks = 0;
      std::set<unsigned int> layers;
    };

    inline bool isInputWhitespace(char c)
    {
      return c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v';
    }

    inline bool isDigit(char c)
    {
      return c >= '0' && c <= '9';
    }

    inline void skipWhitespace(const char*& p)
    {
      while (isInputWhitespace(*p)) {
        ++p;
      }
    }

    inline bool parseUnsigned(const char*& p, unsigned int& value)
    {
      skipWhitespace(p);
      if (*p == '+') {
        ++p;
      }
      if (!isDigit(*p)) {
        return false;
      }

      unsigned long long result = 0;
      const unsigned long long maxValue = std::numeric_limits<unsigned int>::max();
      do {
        result = result * 10 + static_cast<unsigned long long>(*p - '0');
        if (result > maxValue) {
          return false;
        }
        ++p;
      } while (isDigit(*p));

      value = static_cast<unsigned int>(result);
      return true;
    }

    inline bool parseOptionalDouble(const char*& p, double& value)
    {
      skipWhitespace(p);
      if (*p == '\0' || *p == '#') {
        return false;
      }

      char* end = nullptr;
      const double parsed = std::strtod(p, &end);
      if (end == p) {
        return false;
      }

      value = parsed;
      p = end;
      return true;
    }

    inline ParsedLink parseLink(const std::string& line)
    {
      const char* p = line.c_str();
      ParsedLink link;
      if (!parseUnsigned(p, link.source) || !parseUnsigned(p, link.target))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse link data from line '{}'"), line));
      if (!parseOptionalDouble(p, link.weight)) {
        link.weight = 1.0;
      }
      return link;
    }

    inline ParsedStateNode parseStateNode(const std::string& line)
    {
      const char* p = line.c_str();
      ParsedStateNode parsed;
      if (!parseUnsigned(p, parsed.node.id) || !parseUnsigned(p, parsed.node.physicalId))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse any state node from line '{}'"), line));

      auto nameStart = line.find_first_of('\"', static_cast<std::size_t>(p - line.c_str()));
      auto nameEnd = line.find_last_of('\"');
      if (nameStart < nameEnd) {
        parsed.node.name = std::string(line.begin() + nameStart + 1, line.begin() + nameEnd);
        p = line.c_str() + nameEnd + 1;
      }
      if (parseOptionalDouble(p, parsed.node.weight)) {
        parsed.hasWeight = true;
        if (parsed.node.weight < 0)
          throw std::runtime_error(fmt::format(FMT_STRING("Negative state node weight ({}) from line '{}'"), parsed.node.weight, line));
      }
      return parsed;
    }

    inline ParsedMultilayerLink parseMultilayerLink(const std::string& line)
    {
      const char* p = line.c_str();
      ParsedMultilayerLink link;
      if (!parseUnsigned(p, link.sourceLayer) || !parseUnsigned(p, link.sourceNode) || !parseUnsigned(p, link.targetLayer) || !parseUnsigned(p, link.targetNode))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse multilayer link data from line '{}'"), line));
      if (!parseOptionalDouble(p, link.weight)) {
        link.weight = 1.0;
      }
      return link;
    }

    inline ParsedMultilayerIntraLink parseMultilayerIntraLink(const std::string& line)
    {
      const char* p = line.c_str();
      ParsedMultilayerIntraLink link;
      if (!parseUnsigned(p, link.layer) || !parseUnsigned(p, link.sourceNode) || !parseUnsigned(p, link.targetNode))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse intra-multilayer link data from line '{}'"), line));
      if (!parseOptionalDouble(p, link.weight)) {
        link.weight = 1.0;
      }
      return link;
    }

    inline ParsedMultilayerInterLink parseMultilayerInterLink(const std::string& line)
    {
      const char* p = line.c_str();
      ParsedMultilayerInterLink link;
      if (!parseUnsigned(p, link.sourceLayer) || !parseUnsigned(p, link.node) || !parseUnsigned(p, link.targetLayer))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse inter-multilayer link data from line '{}'"), line));
      if (!parseOptionalDouble(p, link.weight)) {
        link.weight = 1.0;
      }
      if (link.sourceLayer == link.targetLayer)
        throw std::runtime_error(fmt::format(FMT_STRING("Inter-layer link from line '{}' doesn't go between different layers."), line));
      return link;
    }

    inline InsensitiveStringSet generalValidHeadings()
    {
      return { "*vertices", "*states", "*multilayer", "*multiplex", "*intra", "*inter", "*paths", "*edges", "*arcs", "*links", "*contexts", "*bipartite" };
    }

    inline InsensitiveStringSet generalIgnoredHeadings()
    {
      return { "*contexts" };
    }

    template <typename Sink>
    std::string parseVertices(std::ifstream& file, const std::string&, Sink& sink, ParseProgress& progress)
    {
      Console::detail(1, "parsing vertices");
      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line.empty() || line[0] == '#')
          continue;

        if (line[0] == '*')
          break;

        std::istringstream extractor(line);
        ParsedVertex vertex;
        std::string idToken;
        if (!(extractor >> idToken))
          throw std::runtime_error(fmt::format(FMT_STRING("Can't parse node id from line '{}'"), line));
        // Validated rather than read straight into the unsigned: `>> unsigned`
        // accepts "-2" and stores 4294967294.
        if (!io::parseNonNegativeInteger(idToken, vertex.id))
          throw std::runtime_error(fmt::format(FMT_STRING("Can't parse node id '{}' from line '{}': expected a non-negative integer no greater than {}"), idToken, line, std::numeric_limits<unsigned int>::max()));

        auto nameStart = line.find_first_of('\"');
        auto nameEnd = line.find_last_of('\"');
        if (nameStart < nameEnd) {
          vertex.name = std::string(line.begin() + nameStart + 1, line.begin() + nameEnd);
          line = line.substr(nameEnd + 1);
          extractor.clear();
          extractor.str(line);
        } else {
          if (!(extractor >> vertex.name))
            throw std::runtime_error(fmt::format(FMT_STRING("Can't parse node name from line '{}'"), line));
        }
        if ((extractor >> vertex.weight)) {
          vertex.hasWeight = true;
          if (vertex.weight < 0)
            throw std::runtime_error(fmt::format(FMT_STRING("Negative node weight ({}) from line '{}'"), vertex.weight, line));
        }

        sink.onPhysicalNode(vertex);
        ++progress.numPhysicalNodes;
      }
      Console::detail(1, "{} physical nodes added", progress.numPhysicalNodes);
      return line;
    }

    template <typename Sink>
    std::string parseStateNodes(std::ifstream& file, const std::string&, Sink& sink)
    {
      Console::detail(1, "parsing state nodes");
      unsigned int numStateNodesFound = 0;
      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line.empty() || line[0] == '#')
          continue;

        if (line[0] == '*')
          break;

        auto stateNode = parseStateNode(line);
        sink.onStateNode(stateNode);
        ++numStateNodesFound;
      }
      Console::detail(1, "{} state nodes added", numStateNodesFound);
      return line;
    }

    template <typename Sink>
    std::string parseLinks(std::ifstream& file, Sink& sink, ParseProgress& progress)
    {
      bool parsingLinks = false;
      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line.empty() || line[0] == '#')
          continue;

        if (line[0] == '*')
          break;

        if (!parsingLinks) {
          parsingLinks = true;
          Console::detail(1, "parsing links");
        }

        sink.onLink(parseLink(line));
        ++progress.numLinks;
      }
      if (parsingLinks)
        Console::detail(1, "{} links", progress.numLinks);
      return line;
    }

    template <typename Sink>
    std::string parseMultilayerLinks(std::ifstream& file, Sink& sink, const NetworkInputOptions& options, ParseProgress& progress)
    {
      Console::detail(1, "parsing multilayer links");

      if (options.matchableMultilayerIds > 0) {
        Console::detail(1, "creating matchable state ids: nodeId << (log2({}) + 1) | layerId", options.matchableMultilayerIds);
      }

      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line.empty() || line[0] == '#')
          continue;

        if (line[0] == '*')
          break;

        const auto link = parseMultilayerLink(line);
        sink.onMultilayerLink(link);
        progress.layers.insert(link.sourceLayer);
        progress.layers.insert(link.targetLayer);
        if (link.sourceLayer == link.targetLayer) {
          ++progress.numIntraLayerLinks;
        } else {
          ++progress.numInterLayerLinks;
        }
      }
      Console::detail(1, "{} links in {} layers", progress.numIntraLayerLinks + progress.numInterLayerLinks, progress.layers.size());
      Console::detail(1, "{} intra-layer links", progress.numIntraLayerLinks);
      Console::detail(1, "{} inter-layer links", progress.numInterLayerLinks);
      return line;
    }

    template <typename Sink>
    std::string parseMultilayerIntraLinks(std::ifstream& file, Sink& sink, const NetworkInputOptions& options, ParseProgress& progress)
    {
      Console::detail(1, "parsing intra-layer links");

      if (options.matchableMultilayerIds > 0) {
        Console::detail(1, "creating matchable state ids: nodeId << (log2({}) + 1) | layerId", options.matchableMultilayerIds);
      }

      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line.empty() || line[0] == '#')
          continue;

        if (line[0] == '*')
          break;

        const auto link = parseMultilayerIntraLink(line);
        sink.onMultilayerIntraLink(link);
        progress.layers.insert(link.layer);
        ++progress.numIntraLayerLinks;
      }
      Console::detail(1, "{} intra-layer links", progress.numIntraLayerLinks);
      return line;
    }

    template <typename Sink>
    std::string parseMultilayerInterLinks(std::ifstream& file, Sink& sink, ParseProgress& progress)
    {
      Console::detail(1, "parsing inter-layer links");
      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line.empty() || line[0] == '#')
          continue;

        if (line[0] == '*')
          break;

        const auto link = parseMultilayerInterLink(line);
        sink.onMultilayerInterLink(link);
        progress.layers.insert(link.sourceLayer);
        progress.layers.insert(link.targetLayer);
        ++progress.numInterLayerLinks;
      }
      Console::detail(1, "{} inter-layer links", progress.numInterLayerLinks);
      return line;
    }

    template <typename Sink>
    std::string parseBipartiteLinks(std::ifstream& file, const std::string& heading, Sink& sink)
    {
      Console::detail(1, "parsing bipartite links");
      std::istringstream extractor(heading);
      std::string tmp;
      unsigned int bipartiteStartId = 0;
      if (!(extractor >> tmp >> bipartiteStartId))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse bipartite start id from line '{}'"), heading));

      Console::detail(1, "using bipartite start id {}", bipartiteStartId);
      sink.onBipartiteStart(bipartiteStartId);
      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line.empty() || line[0] == '#')
          continue;

        if (line[0] == '*')
          break;

        sink.onBipartiteLink(line, parseLink(line));
      }
      return line;
    }

    inline std::string ignoreSection(std::ifstream& file, const std::string& heading)
    {
      Console::note(0, "Ignoring unrecognized section '{}'.", heading);
      std::string line;
      while (!std::getline(file, line).fail()) {
        if (line[0] == '*')
          break;
      }
      return line;
    }

    template <typename Sink>
    void parseNetwork(const std::string& filename,
                      const InsensitiveStringSet& validHeadings,
                      const InsensitiveStringSet& ignoreHeadings,
                      Sink& sink,
                      const NetworkInputOptions& options,
                      const std::string& startHeading = "")
    {
      SafeInFile input(filename);
      ParseProgress progress;

      std::string heading = !startHeading.empty() ? startHeading : parseLinks(input, sink, progress);

      while (!heading.empty() && heading[0] == '*') {
        std::string headingLowerCase = io::tolower(io::firstWord(heading));
        if (validHeadings.count(headingLowerCase) == 0) {
          throw std::runtime_error(fmt::format(FMT_STRING("Unrecognized heading in network file: '{}'."), headingLowerCase));
        }
        bool shouldIgnoreHeading = ignoreHeadings.count(headingLowerCase) > 0;
        if (!shouldIgnoreHeading && headingLowerCase == "*vertices") {
          heading = parseVertices(input, heading, sink, progress);
        } else if (!shouldIgnoreHeading && headingLowerCase == "*states") {
          heading = parseStateNodes(input, heading, sink);
        } else if (!shouldIgnoreHeading && headingLowerCase == "*edges") {
          if (!options.undirectedFlow) {
            Log() << "\n";
            Console::note(0, "Links marked as undirected but parsed as directed.");
          }
          heading = parseLinks(input, sink, progress);
        } else if (!shouldIgnoreHeading && headingLowerCase == "*arcs") {
          if (options.undirectedFlow) {
            Log() << "\n";
            Console::note(0, "Links marked as directed but parsed as undirected.");
          }
          heading = parseLinks(input, sink, progress);
        } else if (!shouldIgnoreHeading && headingLowerCase == "*links") {
          heading = parseLinks(input, sink, progress);
        } else if (!shouldIgnoreHeading && (headingLowerCase == "*multilayer" || headingLowerCase == "*multiplex")) {
          heading = parseMultilayerLinks(input, sink, options, progress);
        } else if (!shouldIgnoreHeading && headingLowerCase == "*intra") {
          heading = parseMultilayerIntraLinks(input, sink, options, progress);
        } else if (!shouldIgnoreHeading && headingLowerCase == "*inter") {
          heading = parseMultilayerInterLinks(input, sink, progress);
        } else if (!shouldIgnoreHeading && headingLowerCase == "*bipartite") {
          heading = parseBipartiteLinks(input, heading, sink);
        } else {
          heading = ignoreSection(input, headingLowerCase);
        }
      }

      Console::detail(1, "done");
    }

  } // namespace detail

  template <typename Sink>
  void parseNetworkInput(const std::string& filename, Sink& sink, const NetworkInputOptions& options)
  {
    Console::detail(1, "parsing {} network from {}", options.undirectedFlow ? "undirected" : "directed", filename);
    detail::parseNetwork(filename, detail::generalValidHeadings(), detail::generalIgnoredHeadings(), sink, options);
  }

  template <typename Sink>
  void parseMetaDataInput(const std::string& filename, Sink& sink)
  {
    Console::detail(1, "parsing meta data from {}", filename);
    SafeInFile input(filename);
    std::string line;
    unsigned int numMetaDataColumns = 0;
    unsigned int numMetaDataRows = 0;
    while (!std::getline(input, line).fail()) {
      if (line.empty() || line[0] == '#')
        continue;

      if (line[0] == '*')
        break;

      std::istringstream extractor(line);
      unsigned int nodeId = 0;
      std::string nodeIdToken;
      if (!(extractor >> nodeIdToken))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse node id from line '{}'"), line));
      if (!io::parseNonNegativeInteger(nodeIdToken, nodeId))
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse node id '{}' from line '{}': expected a non-negative integer no greater than {}"), nodeIdToken, line, std::numeric_limits<unsigned int>::max()));

      // Same treatment as the id above, and for the same reason: `>> unsigned`
      // takes "-1" as 4294967295, which then becomes -1 again in this vector<int>,
      // so two different files produced the same category with no complaint. The
      // upper bound is int, not unsigned, because that is what the vector holds.
      std::vector<int> metaData;
      std::string metaToken;
      while (extractor >> metaToken) {
        // An inline '#' ends the line, the same convention the link rows use
        // (parseLink stops at '#'). Reading the categories as unsigned used to give
        // this for free -- the stream failed on '#' and quietly ended the loop -- so
        // validating the tokens has to keep it explicitly or "1 2 # note" stops
        // parsing, which it did.
        const auto commentStart = metaToken.find('#');
        const bool commentFollows = commentStart != std::string::npos;
        if (commentFollows)
          metaToken.erase(commentStart);

        if (!metaToken.empty()) {
          unsigned int metaId = 0;
          if (!io::parseNonNegativeInteger(metaToken, metaId) || metaId > static_cast<unsigned int>(std::numeric_limits<int>::max()))
            throw std::runtime_error(fmt::format(FMT_STRING("Can't parse meta category '{}' from line '{}': expected a non-negative integer no greater than {}"), metaToken, line, std::numeric_limits<int>::max()));
          metaData.push_back(static_cast<int>(metaId));
        }

        if (commentFollows)
          break;
      }
      if (metaData.empty())
        throw std::runtime_error(fmt::format(FMT_STRING("Can't parse any meta data from line '{}'"), line));

      sink.onMetaData(nodeId, metaData);
      numMetaDataColumns = std::max<unsigned int>(numMetaDataColumns, metaData.size());
      ++numMetaDataRows;
    }
    Console::detail(1, "parsed {} columns of meta data for {} nodes", numMetaDataColumns, numMetaDataRows);
  }

} // namespace input
} // namespace infomap

#endif // NETWORK_INPUT_PARSER_H_
