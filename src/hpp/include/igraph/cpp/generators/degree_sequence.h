/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GENERATORS_DEGREE_SEQUENCE_H
#define IGRAPHPP_GENERATORS_DEGREE_SEQUENCE_H

#include <memory>

namespace igraph {

class Graph;

/// Generates an undirected random graph with a given degree sequence
std::auto_ptr<Graph> degree_sequence_game(const Vector& degrees,
        DegreeSequenceMethod method = IGRAPH_DEGSEQ_SIMPLE);

/// Generates a directed random graph with a given degree sequence
std::auto_ptr<Graph> degree_sequence_game(const Vector& outdegrees, const Vector& indegrees,
        DegreeSequenceMethod method = IGRAPH_DEGSEQ_SIMPLE);

}         // end of namespace

#endif    // IGRAPHPP_GENERATORS_DETERMINISTIC_H
