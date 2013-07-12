/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_ERROR_H
#define IGRAPHPP_ERROR_H

#include <igraph/igraph_error.h>

#define IGRAPH_TRY(call) {\
    int __result = call;    \
    if (__result != 0)      \
        throw std::runtime_error("igraph error"); \
}

#define IGRAPHPP_TRY_NEW(variable, type)   \
    try {                                \
        variable = new type;             \
    } catch (const std::bad_alloc&) { \
        IGRAPH_ERROR("std::bad_alloc thrown in C++ code", IGRAPH_ENOMEM); \
    }

#endif    // IGRAPHPP_ERROR_H

