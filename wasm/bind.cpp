#include <memory>
#include <iostream>
#include <functional>

#include "igraph.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>

using namespace emscripten;

#define IGRAPH_VECTOR_EMCC(name) \
    val get_##name(const igraph_t& graph) { \
        const auto& vector{graph.name}; \
        const auto size{(vector.stor_end - vector.stor_begin)}; \
        return val{typed_memory_view(size, vector.stor_begin)}; \
    }

bool get_directed(const igraph_t& graph) {
    return graph.directed;
}

igraph_integer_t get_n(const igraph_t& graph) {
    return graph.n;
}

val get_vector(const igraph_vector_int_t& v) {
    const auto size{(v.stor_end - v.stor_begin)};
    return val{typed_memory_view(size, v.stor_begin)};
}

igraph_vector_int_t create_vector(val arr) {
    igraph_vector_int_t v{};
    const auto data{convertJSArrayToNumberVector<igraph_integer_t>(arr)};

    igraph_vector_int_init(&v, data.size());
    for (size_t i{0}; i < data.size(); ++i) VECTOR(v)[i] = data[i];

    return v;
}

igraph_error_t famous(igraph_t* graph, std::string name) {
    return igraph_famous(graph, name.c_str());
}

IGRAPH_VECTOR_EMCC(from)
IGRAPH_VECTOR_EMCC(to)
IGRAPH_VECTOR_EMCC(oi)
IGRAPH_VECTOR_EMCC(ii)
IGRAPH_VECTOR_EMCC(is)
IGRAPH_VECTOR_EMCC(os)

EMSCRIPTEN_BINDINGS(IGraph)
{
  class_<igraph_t>("graph").
    constructor().
    property("n", &get_n).
    property("directed", &get_directed).
    property("from", &get_from).
    property("to", &get_to).
    property("oi", &get_oi).
    property("ii", &get_ii).
    property("is", &get_is).
    property("os", &get_os);

  class_<igraph_vector_int_t>("intVector").
    constructor().
    constructor(&create_vector, allow_raw_pointers()).
    property("begin", &get_vector);

  enum_<igraph_error_type_t>("igraph_error");

  function("empty", &igraph_empty, allow_raw_pointers());
  function("destroy", &igraph_destroy, allow_raw_pointers());
  function("create", &igraph_create, allow_raw_pointers());
  function("ring", &igraph_ring, allow_raw_pointers());
  function("edgelist", &igraph_get_edgelist, allow_raw_pointers());
  function("famous", &famous, allow_raw_pointers());
}
#endif
