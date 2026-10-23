/*******************************************************************************
 Infomap software package for multi-level network clustering
 Copyright (c) 2013, 2014 Daniel Edler, Anton Holmgren, Martin Rosvall

 This file is part of the Infomap software package.
 See file LICENSE_GPLv3.txt for full license details.
 For more information, see <http://www.mapequation.org>
 ******************************************************************************/

#ifndef JSON_NETWORK_INPUT_PARSER_H_
#define JSON_NETWORK_INPUT_PARSER_H_

// JSON network input parser for the infomap-network v1.0 format (RFC #645).
//
// Parsing is order- and emitter-independent via two SAX passes (no DOM):
//   pass 1 captures the root scalar discriminators (format/version/type/...),
//   pass 2 streams `nodes`/`edges` into the same core network sink the text
//   parser uses, so all network semantics are inherited from the core builder.
//
// Phase 1 wires the `standard` network type (and the omitted-type default).
// bipartite/multilayer/state, embedded `meta`, and embedded `path` are parsed
// into the pending structures but not yet emitted; later phases connect them.

#include "NetworkInputParser.h"
#include "../utils/Console.h"
#include "../utils/format.h"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace infomap {
namespace input {

  namespace json_detail {

    using Json = nlohmann::json;

    enum class NetworkType : std::uint8_t { Standard,
                                            Bipartite,
                                            Multilayer,
                                            State };

    enum class MultilayerMode : std::uint8_t { Full,
                                               Intra,
                                               IntraInter };

    struct JsonHeader {
      std::string format;
      std::string version;
      std::string type = "standard";
      std::string multilayer = "full";
      bool hasDirected = false;
      bool directed = false;
      bool hasBipartiteStartId = false;
      unsigned int bipartiteStartId = 0;
      bool hasNodes = false;
      bool hasStates = false;
    };

    // Integer coercion (RFC): integral-valued doubles are accepted, non-integral
    // values are rejected, ids are bounded by unsigned int (uint32) max.
    inline bool coerceUnsigned(double value, unsigned int& out)
    {
      if (std::floor(value) != value)
        return false;
      if (value < 0.0)
        return false;
      if (value > static_cast<double>(std::numeric_limits<unsigned int>::max()))
        return false;
      out = static_cast<unsigned int>(value);
      return true;
    }

    inline bool coerceUnsigned(std::int64_t value, unsigned int& out)
    {
      if (value < 0)
        return false;
      if (static_cast<std::uint64_t>(value) > std::numeric_limits<unsigned int>::max())
        return false;
      out = static_cast<unsigned int>(value);
      return true;
    }

    inline bool coerceUnsigned(std::uint64_t value, unsigned int& out)
    {
      if (value > std::numeric_limits<unsigned int>::max())
        return false;
      out = static_cast<unsigned int>(value);
      return true;
    }

    // A numeric value flowing through the SAX callbacks, normalized to double
    // plus an exactly-integral marker so id coercion can reject e.g. 1.5.
    struct Number {
      double value = 0.0;
      bool integral = false;
    };

    // Pass 1: capture root scalars, skip the big arrays without materializing them.
    struct HeaderSaxHandler {
      JsonHeader& header;
      std::size_t objectDepth = 0;
      std::size_t arrayDepth = 0;
      std::string currentKey;

      explicit HeaderSaxHandler(JsonHeader& h) : header(h) {}

      bool atRoot() const { return objectDepth == 1 && arrayDepth == 0; }

      bool null() { return true; }
      bool boolean(bool val)
      {
        if (atRoot() && currentKey == "directed") {
          header.hasDirected = true;
          header.directed = val;
        }
        return true;
      }
      bool number_integer(std::int64_t val) { return number_unsigned(static_cast<std::uint64_t>(val < 0 ? 0 : val)); }
      bool number_unsigned(std::uint64_t val)
      {
        if (atRoot() && currentKey == "bipartiteStartId") {
          unsigned int id = 0;
          if (coerceUnsigned(val, id)) {
            header.hasBipartiteStartId = true;
            header.bipartiteStartId = id;
          }
        }
        return true;
      }
      bool number_float(double, const std::string&) { return true; }
      bool string(std::string& val)
      {
        if (atRoot()) {
          if (currentKey == "format")
            header.format = val;
          else if (currentKey == "version")
            header.version = val;
          else if (currentKey == "type")
            header.type = val;
          else if (currentKey == "multilayer")
            header.multilayer = val;
        }
        return true;
      }
      bool binary(Json::binary_t&) { return true; }
      bool start_object(std::size_t)
      {
        ++objectDepth;
        return true;
      }
      bool key(std::string& val)
      {
        currentKey = val;
        return true;
      }
      bool end_object()
      {
        --objectDepth;
        return true;
      }
      bool start_array(std::size_t)
      {
        if (objectDepth == 1 && arrayDepth == 0) {
          if (currentKey == "nodes")
            header.hasNodes = true;
          else if (currentKey == "states")
            header.hasStates = true;
        }
        ++arrayDepth;
        return true;
      }
      bool end_array()
      {
        --arrayDepth;
        return true;
      }
      bool parse_error(std::size_t position, const std::string& last_token, const Json::exception& ex)
      {
        throw std::runtime_error(fmt::format(FMT_STRING("JSON parse error at byte {} (near '{}'): {}"), position, last_token, ex.what()));
      }
    };

    // Accumulated fields of the current `nodes[]` / `states[]` / `edges[]` element.
    struct PendingNode {
      bool hasId = false;
      unsigned int id = 0;
      std::string name;
      bool hasWeight = false;
      double weight = 1.0;
      bool hasMeta = false;
      int meta = 0;
      std::vector<unsigned int> path; // collected; emitted in a later phase
    };

    struct PendingEdge {
      bool hasSource = false;
      unsigned int source = 0;
      bool hasTarget = false;
      unsigned int target = 0;
      bool hasWeight = false;
      double weight = 1.0;
      std::vector<unsigned int> layers; // multilayer only
    };

    struct PendingState {
      bool hasId = false;
      unsigned int id = 0;
      bool hasNode = false;
      unsigned int node = 0;
      std::string name;
      bool hasWeight = false;
      double weight = 1.0;
      std::vector<unsigned int> path; // collected; emitted in a later phase
    };

    // Pass 2: stream nodes/edges into the core network sink.
    template <typename Sink>
    struct ContentSaxHandler {
      const JsonHeader& header;
      Sink& sink;
      const NetworkInputOptions& options;

      enum class Section : std::uint8_t { None,
                                          Nodes,
                                          Edges,
                                          States };
      Section section = Section::None;
      std::size_t objectDepth = 0;
      std::size_t arrayDepth = 0;
      std::string currentKey;
      std::set<std::string> warnedUnknown;

      NetworkType networkType = NetworkType::Standard;
      MultilayerMode multilayerMode = MultilayerMode::Full;

      PendingNode node;
      PendingEdge edge;
      PendingState state;
      std::set<unsigned int> emittedImplicitStates;

      ContentSaxHandler(const JsonHeader& h, Sink& s, const NetworkInputOptions& o)
          : header(h), sink(s), options(o)
      {
        if (header.type == "bipartite")
          networkType = NetworkType::Bipartite;
        else if (header.type == "multilayer")
          networkType = NetworkType::Multilayer;
        else if (header.type == "state")
          networkType = NetworkType::State;

        if (header.multilayer == "intra")
          multilayerMode = MultilayerMode::Intra;
        else if (header.multilayer == "intra-inter")
          multilayerMode = MultilayerMode::IntraInter;
      }

      bool inRootScalar() const { return objectDepth == 1 && arrayDepth == 0; }
      bool inElement() const { return objectDepth == 2; }
      bool inFieldArray() const { return objectDepth == 2 && arrayDepth == 2; }

      void warnUnknown(const std::string& context, const std::string& field)
      {
        const std::string label = context + "." + field;
        if (warnedUnknown.insert(label).second)
          Console::note(0, "Ignoring unrecognized JSON field '{}'.", label);
      }

      void onElementNumber(const Number& number)
      {
        if (section == Section::Nodes) {
          if (currentKey == "id") {
            if (!coerceUnsigned(number.value, node.id) || !number.integral)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid node id '{:g}' (must be a non-negative integer <= {})"), number.value, std::numeric_limits<unsigned int>::max()));
            node.hasId = true;
          } else if (currentKey == "weight") {
            // Same rule the text parser applies to a *Vertices weight, and the
            // one the normative schema states (weight.minimum: 0). Without it a
            // negative weight reached addPhysicalNode verbatim -- unlike addLink,
            // which is the documented funnel that rejects these -- and the run
            // only failed later, in the flow guard, with a message about
            // bipartite teleportation that has nothing to do with the cause.
            if (!(number.value >= 0) || !std::isfinite(number.value))
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid node weight '{:g}' (must be a finite non-negative number)"), number.value));
            node.weight = number.value;
            node.hasWeight = true;
          } else if (currentKey == "meta") {
            unsigned int meta = 0;
            // Bounded by int, not unsigned: the value is stored in a vector<int>, so
            // anything above INT_MAX came back out negative -- the same wrap the text
            // path had, and fixing one without the other would just move the
            // disagreement between the two parsers.
            if (!coerceUnsigned(number.value, meta) || !number.integral || meta > static_cast<unsigned int>(std::numeric_limits<int>::max()))
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid node meta '{:g}' (must be a non-negative integer no greater than {})"), number.value, std::numeric_limits<int>::max()));
            node.meta = static_cast<int>(meta);
            node.hasMeta = true;
          } else if (currentKey == "path" && inFieldArray()) {
            unsigned int module = 0;
            if (!coerceUnsigned(number.value, module) || !number.integral || module == 0)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid path module '{:g}' (must be a positive integer)"), number.value));
            node.path.push_back(module);
          }
        } else if (section == Section::States) {
          if (currentKey == "id") {
            if (!coerceUnsigned(number.value, state.id) || !number.integral)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid state id '{:g}' (must be a non-negative integer <= {})"), number.value, std::numeric_limits<unsigned int>::max()));
            state.hasId = true;
          } else if (currentKey == "node") {
            if (!coerceUnsigned(number.value, state.node) || !number.integral)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid state physical node '{:g}' (must be a non-negative integer)"), number.value));
            state.hasNode = true;
          } else if (currentKey == "weight") {
            if (!(number.value >= 0) || !std::isfinite(number.value))
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid state weight '{:g}' (must be a finite non-negative number)"), number.value));
            state.weight = number.value;
            state.hasWeight = true;
          } else if (currentKey == "path" && inFieldArray()) {
            unsigned int module = 0;
            if (!coerceUnsigned(number.value, module) || !number.integral || module == 0)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid path module '{:g}' (must be a positive integer)"), number.value));
            state.path.push_back(module);
          }
        } else if (section == Section::Edges) {
          if (currentKey == "source") {
            if (!coerceUnsigned(number.value, edge.source) || !number.integral)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid edge source '{:g}' (must be a non-negative integer <= {})"), number.value, std::numeric_limits<unsigned int>::max()));
            edge.hasSource = true;
          } else if (currentKey == "target") {
            if (!coerceUnsigned(number.value, edge.target) || !number.integral)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid edge target '{:g}' (must be a non-negative integer <= {})"), number.value, std::numeric_limits<unsigned int>::max()));
            edge.hasTarget = true;
          } else if (currentKey == "weight") {
            edge.weight = number.value;
            edge.hasWeight = true;
          } else if (currentKey == "layers" && inFieldArray()) {
            unsigned int layer = 0;
            if (!coerceUnsigned(number.value, layer) || !number.integral)
              throw std::runtime_error(fmt::format(FMT_STRING("Invalid layer id '{:g}' (must be a non-negative integer)"), number.value));
            edge.layers.push_back(layer);
          }
        }
      }

      void flushNode()
      {
        if (!node.hasId)
          throw std::runtime_error("JSON node is missing the required 'id' field");
        if (networkType == NetworkType::Multilayer && !node.path.empty())
          throw std::runtime_error("'path' on a node is invalid for a multilayer network; use type 'state' with states[].path");
        ParsedVertex vertex;
        vertex.id = node.id;
        vertex.name = node.name;
        vertex.hasWeight = node.hasWeight;
        vertex.weight = node.weight;
        sink.onPhysicalNode(vertex);
        // Embedded one-dimensional meta enables meta coding the same way the
        // in-memory set_meta_data / --meta-data inputs do (presence enables).
        if (node.hasMeta)
          sink.onMetaData(node.id, std::vector<int> { node.meta });
        if (!node.path.empty())
          sink.onInitialPartitionPath(node.id, node.path, /*stateLevel=*/false);
        node = PendingNode {};
      }

      void flushState()
      {
        if (!state.hasId || !state.hasNode)
          throw std::runtime_error("JSON state is missing the required 'id'/'node' field");
        ParsedStateNode parsed;
        parsed.node.id = state.id;
        parsed.node.physicalId = state.node;
        parsed.node.name = state.name;
        parsed.hasWeight = state.hasWeight;
        if (state.hasWeight)
          parsed.node.weight = state.weight;
        sink.onStateNode(parsed);
        emittedImplicitStates.insert(state.id);
        if (!state.path.empty())
          sink.onInitialPartitionPath(state.id, state.path, /*stateLevel=*/true);
        state = PendingState {};
      }

      // For state networks without an explicit states[] array, infer one state
      // per edge endpoint with an identity physical mapping (physicalId == id).
      void ensureImplicitState(unsigned int stateId)
      {
        if (networkType != NetworkType::State || header.hasStates)
          return;
        if (!emittedImplicitStates.insert(stateId).second)
          return;
        ParsedStateNode parsed;
        parsed.node.id = stateId;
        parsed.node.physicalId = stateId;
        sink.onStateNode(parsed);
      }

      ParsedLink toLink() const
      {
        ParsedLink link;
        link.source = edge.source;
        link.target = edge.target;
        link.weight = edge.hasWeight ? edge.weight : 1.0;
        return link;
      }

      void flushMultilayerEdge()
      {
        const double weight = edge.hasWeight ? edge.weight : 1.0;
        const std::size_t numLayers = edge.layers.size();

        if (multilayerMode == MultilayerMode::Intra) {
          if (numLayers != 1)
            throw std::runtime_error("multilayer 'intra' edge requires exactly one layer id");
          sink.onMultilayerIntraLink(ParsedMultilayerIntraLink { edge.layers[0], edge.source, edge.target, weight });
        } else if (multilayerMode == MultilayerMode::IntraInter) {
          if (numLayers == 1) {
            sink.onMultilayerIntraLink(ParsedMultilayerIntraLink { edge.layers[0], edge.source, edge.target, weight });
          } else if (numLayers == 2) {
            if (edge.source != edge.target)
              throw std::runtime_error("multilayer 'intra-inter' inter-layer edge requires source == target (use multilayer 'full' for non-diagonal inter-layer links)");
            sink.onMultilayerInterLink(ParsedMultilayerInterLink { edge.layers[0], edge.source, edge.layers[1], weight });
          } else {
            throw std::runtime_error("multilayer 'intra-inter' edge requires one or two layer ids");
          }
        } else { // Full
          if (numLayers != 2)
            throw std::runtime_error("multilayer 'full' edge requires exactly two layer ids");
          sink.onMultilayerLink(ParsedMultilayerLink { edge.layers[0], edge.source, edge.layers[1], edge.target, weight });
        }
      }

      void flushEdge()
      {
        if (networkType == NetworkType::Multilayer) {
          if (!edge.hasSource || !edge.hasTarget)
            throw std::runtime_error("JSON multilayer edge is missing the required 'source'/'target' field");
          if (edge.layers.empty())
            throw std::runtime_error("JSON multilayer edge is missing the required 'layers' field");
          flushMultilayerEdge();
          edge = PendingEdge {};
          return;
        }

        if (!edge.hasSource || !edge.hasTarget)
          throw std::runtime_error("JSON edge is missing the required 'source'/'target' field");

        if (networkType == NetworkType::Bipartite) {
          const ParsedLink link = toLink();
          sink.onBipartiteLink(fmt::format(FMT_STRING("{} {} {:g}"), link.source, link.target, link.weight), link);
        } else if (networkType == NetworkType::State) {
          ensureImplicitState(edge.source);
          ensureImplicitState(edge.target);
          sink.onLink(toLink());
        } else {
          sink.onLink(toLink());
        }
        edge = PendingEdge {};
      }

      bool routeNumber(const Number& number)
      {
        if (inRootScalar())
          return true; // header scalars already captured in pass 1
        if (objectDepth >= 2)
          onElementNumber(number);
        return true;
      }

      bool null() { return true; }
      bool boolean(bool) { return true; }
      bool number_integer(std::int64_t val)
      {
        return routeNumber(Number { static_cast<double>(val), true });
      }
      bool number_unsigned(std::uint64_t val)
      {
        return routeNumber(Number { static_cast<double>(val), true });
      }
      bool number_float(double val, const std::string&)
      {
        return routeNumber(Number { val, std::floor(val) == val });
      }
      bool string(std::string& val)
      {
        if (inElement() && currentKey == "name") {
          if (section == Section::Nodes)
            node.name = val;
          else if (section == Section::States)
            state.name = val;
        }
        return true;
      }
      bool binary(Json::binary_t&) { return true; }

      bool start_object(std::size_t)
      {
        ++objectDepth;
        return true;
      }
      bool key(std::string& val)
      {
        currentKey = val;
        if (objectDepth == 1 && arrayDepth == 0) {
          const bool knownRoot = val == "format" || val == "version" || val == "type"
              || val == "directed" || val == "multilayer" || val == "bipartiteStartId"
              || val == "nodes" || val == "states" || val == "edges";
          if (!knownRoot)
            warnUnknown("(root)", val);
        } else if (inElement()) {
          if (section == Section::Nodes && !(val == "id" || val == "name" || val == "weight" || val == "meta" || val == "path"))
            warnUnknown("node", val);
          else if (section == Section::States && !(val == "id" || val == "node" || val == "name" || val == "weight" || val == "path"))
            warnUnknown("state", val);
          else if (section == Section::Edges) {
            const bool knownEdgeField = val == "source" || val == "target" || val == "weight"
                || (networkType == NetworkType::Multilayer && val == "layers");
            if (!knownEdgeField)
              warnUnknown("edge", val);
          }
        }
        return true;
      }
      bool end_object()
      {
        if (objectDepth == 2) {
          if (section == Section::Nodes)
            flushNode();
          else if (section == Section::States)
            flushState();
          else if (section == Section::Edges)
            flushEdge();
        }
        --objectDepth;
        return true;
      }
      bool start_array(std::size_t)
      {
        ++arrayDepth;
        if (objectDepth == 1 && arrayDepth == 1) {
          if (currentKey == "nodes")
            section = Section::Nodes;
          else if (currentKey == "edges")
            section = Section::Edges;
          else if (currentKey == "states")
            section = Section::States;
        }
        return true;
      }
      bool end_array()
      {
        --arrayDepth;
        if (objectDepth == 1 && arrayDepth == 0)
          section = Section::None;
        return true;
      }
      bool parse_error(std::size_t position, const std::string& last_token, const Json::exception& ex)
      {
        throw std::runtime_error(fmt::format(FMT_STRING("JSON parse error at byte {} (near '{}'): {}"), position, last_token, ex.what()));
      }
    };

    inline void validateHeader(const JsonHeader& header, const NetworkInputOptions& options)
    {
      if (header.format != "infomap-network")
        throw std::runtime_error(fmt::format(FMT_STRING("Not an infomap-network document (format = '{}')"), header.format));
      if (header.version != "1.0")
        throw std::runtime_error(fmt::format(FMT_STRING("Unsupported infomap-network version '{}' (this build supports 1.0)"), header.version));

      const bool validType = header.type == "standard" || header.type == "bipartite"
          || header.type == "multilayer" || header.type == "state";
      if (!validType)
        throw std::runtime_error(fmt::format(FMT_STRING("Unknown infomap-network type '{}'"), header.type));

      // 'directed' is advisory (like *Edges/*Arcs in the text format): the run
      // flag governs the flow model. Multilayer flow handles direction itself.
      if (header.hasDirected && header.type != "multilayer") {
        const bool fileDirected = header.directed;
        const bool runDirected = !options.undirectedFlow;
        if (fileDirected != runDirected) {
          Log() << "\n";
          Console::note(0, "JSON 'directed': {} but flow parsed as {}.", fileDirected ? "true" : "false", runDirected ? "directed" : "undirected");
        }
      }

      if (header.type == "bipartite" && !header.hasBipartiteStartId)
        throw std::runtime_error("infomap-network type 'bipartite' requires a numeric 'bipartiteStartId'");

      if (header.type == "multilayer") {
        const bool validMode = header.multilayer == "full" || header.multilayer == "intra" || header.multilayer == "intra-inter";
        if (!validMode)
          throw std::runtime_error(fmt::format(FMT_STRING("Unknown multilayer mode '{}'"), header.multilayer));
      }
    }

  } // namespace json_detail

  // Detect infomap-network by content: first non-whitespace byte is '{'.
  inline bool looksLikeJsonNetwork(const std::string& filename)
  {
    std::ifstream file(filename);
    if (!file)
      return false;
    char c = 0;
    while (file.get(c)) {
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == '\f' || c == '\v')
        continue;
      return c == '{';
    }
    return false;
  }

  template <typename Sink>
  void parseJsonNetworkInput(const std::string& filename, Sink& sink, const NetworkInputOptions& options)
  {
    using namespace json_detail;
    Console::detail(1, "parsing infomap-network from {}", filename);

    JsonHeader header;
    {
      std::ifstream file(filename);
      if (!file)
        throw std::runtime_error(fmt::format(FMT_STRING("Can't open file '{}'"), filename));
      HeaderSaxHandler handler(header);
      Json::sax_parse(file, &handler);
    }

    validateHeader(header, options);
    Console::detail(1, "infomap-network type '{}'", header.type);

    // The bipartite boundary is a discriminator known from the header pass, so
    // emit it before any links regardless of nodes/edges key order.
    if (header.type == "bipartite")
      sink.onBipartiteStart(header.bipartiteStartId);

    {
      std::ifstream file(filename);
      if (!file)
        throw std::runtime_error(fmt::format(FMT_STRING("Can't open file '{}'"), filename));
      ContentSaxHandler<Sink> handler(header, sink, options);
      Json::sax_parse(file, &handler);
    }

    Console::detail(1, "done");
  }

} // namespace input
} // namespace infomap

#endif // JSON_NETWORK_INPUT_PARSER_H_
