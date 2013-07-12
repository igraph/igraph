/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_MATCHING_H
#define IGRAPHPP_MATCHING_H

namespace igraph {

class Graph;
class Vector;
class VectorBool;
class VectorLong;

void maximum_bipartite_matching(const Graph& graph, const VectorBool& types,
        integer_t* matching_size, real_t* matching_weight, VectorLong* matching,
        const Vector* weights, real_t eps);

}         // end of namespace

#endif    // IGRAPHPP_MATCHING_H

