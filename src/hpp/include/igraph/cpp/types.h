/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_TYPES_H
#define IGRAPHPP_TYPES_H

#include <stdexcept>
#include <igraph/igraph_constants.h>

namespace igraph {

typedef igraph_integer_t integer_t;
typedef igraph_real_t real_t;
typedef igraph_bool_t bool_t;

typedef igraph_add_weights_t AddWeights;
typedef igraph_connectedness_t Connectedness;
typedef igraph_degseq_t DegreeSequenceMethod;
typedef igraph_edgeorder_type_t EdgeOrderType;
typedef igraph_neimode_t NeighborMode;

}       // end of namespaces

#endif  // IGRAPHPP_VECTOR_H


